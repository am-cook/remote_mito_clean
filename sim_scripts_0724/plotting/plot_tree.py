#!/usr/bin/env python

# saves a static cas tree image for one urid / timepoint, colored by one or more
# clade attributes. no config file -- pure CLI args, called once per
# urid/timepoint/color-by combination you want an image for.

import argparse
from pathlib import Path

# perform_em lives in the sibling em/ directory
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'em'))
import perform_em as em


def parse_args():
    parser = argparse.ArgumentParser(
        description = 'Save a static cas tree image.',
        formatter_class = argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('--urid', required = True, type = int,
                        help = 'urid of the simulation run')
    parser.add_argument('--timepoint', required = True, type = int,
                        help = 'Timepoint to plot')
    parser.add_argument('--color-by', required = True, nargs = '+',
                        dest = 'color_by',
                        help = 'Clade attribute name(s) to color leaves by '
                               '(e.g. clade_celltype 20clade_topological_treedist_gt)')
    parser.add_argument('--output-dir', required = True,
                        dest = 'output_dir',
                        help = 'Directory to write the image into')
    parser.add_argument('--param-combo', default = None,
                        dest = 'param_combo',
                        help = 'scoremat short_name for reconstructed tree '
                               '(omit to use the GT tree)')
    parser.add_argument('--clade-colors', default = None,
                        dest = 'clade_colors',
                        help = 'String key for internal-node coloring '
                               '(e.g. induced_differentiation) or omit for none')
    parser.add_argument('--include-alive-ring', action = 'store_true',
                        dest = 'include_alive_ring',
                        help = 'Add an alive-status ring to each leaf')
    parser.add_argument('--image-format', default = 'png',
                        dest = 'image_format',
                        choices = ['png', 'pdf', 'svg'],
                        help = 'Output image format (default: png)')
    parser.add_argument('--cluster-vmin', type = float, default = -0.5,
                        dest = 'cluster_vmin')
    parser.add_argument('--cluster-vmax', type = float, default = 1.0,
                        dest = 'cluster_vmax')
    return parser.parse_args()


def main():
    args = parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents = True, exist_ok = True)

    # process_urid returns (timepoint_res_dict, model_dict) -- only the first is
    # needed here
    print(f'Processing urid={args.urid}, timepoint={args.timepoint} ...')
    timepoint_res_dict = em.process_urid(
        urid = args.urid,
        param_dict_list = em.param_dict_list_rand10 + [em.param_dict_gt5],
        timepoints = [args.timepoint],
        cluster_vmin = args.cluster_vmin,
        cluster_vmax = args.cluster_vmax,
    )[0]

    tpr = timepoint_res_dict[args.timepoint]
    em.make_clade_df(timepoint_res = tpr, urid = args.urid, timepoint = args.timepoint)

    print('Building tree figure ...')
    img_bytes = tpr.new_plot_cass_tree(
        color_by = args.color_by,
        return_html_string = False,
        clade_colors = args.clade_colors,
        include_alive_status_ring = args.include_alive_ring,
        interactive = False,
        param_combo = args.param_combo,
        image_format = args.image_format,
    )

    # filename records the tree source (gt or the scoremat combo) and what the
    # leaves were colored by, so repeated calls don't overwrite each other
    color_tag = '_'.join(args.color_by)
    tree_tag = args.param_combo if args.param_combo else 'gt'
    out_path = out_dir / f'{args.urid}_tp{args.timepoint}_{tree_tag}_{color_tag}.{args.image_format}'

    out_path.write_bytes(img_bytes)
    print(f'Saved: {out_path}')


if __name__ == '__main__':
    main()
