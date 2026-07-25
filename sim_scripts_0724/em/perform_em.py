#!/usr/bin/env python3
# longmito_36.py - auto-generated from mito_clonal_sub.ipynb cells 0, 1, 2, 29

import sys
import argparse
sys.setrecursionlimit(100000)

# Imports (cell 0)
import base64
try:
    import kaleido
    kaleido.get_chrome_sync()
except Exception:
    pass
from sklearn.metrics import confusion_matrix, precision_score, recall_score, f1_score, classification_report
import joblib
import networkx as nx
import re
import os
import json
import pandas as pd
from ete3 import Tree
from pathlib import Path
import matplotlib.pyplot as plt
import glob
import numpy as np
import copy
import plotly.graph_objects as go
from matplotlib.colors import rgb2hex
from matplotlib.colors import ListedColormap
import itertools
from scipy.spatial.distance import pdist, squareform, jensenshannon
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage
import tqdm
from cassiopeia.data import CassiopeiaTree
from cassiopeia.plotting import plot_plotly
import seaborn as sns
from collections import Counter, defaultdict
from sklearn.metrics import adjusted_rand_score, precision_recall_curve, auc, roc_curve, roc_auc_score
from typing import List, Mapping, Optional, Union
import plotly.graph_objs as go
from plotly.colors import hex_to_rgb
from plotly.express import imshow
from ipywidgets import widgets, interactive, VBox, HBox
from functools import lru_cache
import pickle
import ipywidgets as w
from IPython.display import display, HTML
import plotly.io as pio
from scipy.linalg import expm
from collections import defaultdict
import logging
import dendropy as dp
from scipy.special import logsumexp
import random 
from sklearn.mixture import GaussianMixture
from scipy.special import logit
import gc
from typing import Dict, Any, Optional
import matplotlib.lines as mlines
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from numpy.linalg import matrix_power
from scipy.special import rel_entr
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import accuracy_score
from concurrent.futures import ProcessPoolExecutor, as_completed

path_to_output_dir = Path('./output')
path_to_em_induction_prob = path_to_output_dir / 'em_induction_prob'

induced_color = 'darkred'

uninduced_color = 'steelblue'

# Path/param utilities

def extract_ints(path, joint = False):
    num_ints_pattern = r'.*_list_(\d+)_ints_.*'
    all_matches = re.findall(num_ints_pattern, path)

    if not joint:
        return int(all_matches[0])
    return [int(match) for match in all_matches]

def extract_frac(path, joint = False):
    frac_recover_pattern = r'.*_RP_([\d\.]+)_samp.*'
    all_matches = re.findall(frac_recover_pattern, path)

    if not joint:
        return float(all_matches[0])
    return [float(match) for match in all_matches]

def extract_bin(path, joint = False):
    bin_pattern = r'.*B_([TF]).*'
    all_matches = re.findall(bin_pattern, path)

    if not joint:
        return all_matches[0]
    return all_matches

def extract_af(path, joint = False):
    af_pattern = r'.*_AF_([\d\.]+)_B.*'
    all_matches = re.findall(af_pattern, path)

    if not joint:
        return float(all_matches[0])
    return [float(match) for match in all_matches]

def extract_cell_rec_fracs(path):
    samp_to_time = re.search(r'samp_(.*?)_time', path).group(1)
    sample_fracs = re.findall(r'[\w]+-(\d+\.?\d*)', samp_to_time)
    return sample_fracs

def extract_mt_or_bc(path):
    '''
    also used to determine whether a path is joint or not, depending on length of all_matches
    '''
    
    mt_or_bc_pattern = r'.*roc_(\w+)_list'
    all_matches = re.findall(mt_or_bc_pattern, path)

    return all_matches

def read_in_tree(path_to_tree):
    newick = ''.join(open(path_to_tree).read().split())
    # clean internal nodes
    clean_newick = re.sub(r'\)([^:);]+):', r'):', newick)
    tr = Tree(clean_newick, format = 1)

    return tr

def extract_timepoint(path):
    
    timepoint_pattern = r'.*time_([\d\.]+).*'
    all_matches = re.findall(timepoint_pattern, str(path))

    timept = float(all_matches[0])

    if timept.is_integer():
        return int(timept)
    
    return float(timept)

def get_sim_param_dict(urid, return_param_file_name = False):
    '''
    read the json simulation param file in as a dictionary using urid
    '''
    # path_to_run_specs_dir = Path(f'/dartfs/rc/lab/M/McKennaLab/projects/Aidan/simulations/r_sim_clean/output/run_specs/{urid}/')
    path_to_run_specs_dir = path_to_output_dir / f'run_specs/{urid}/'
    path_to_params_json = [str(path) for path in path_to_run_specs_dir.iterdir() if path.is_file()][0]
    with open(path_to_params_json, 'r') as f:
        params_dict = json.load(f)

    if return_param_file_name:
        return params_dict, path_to_params_json.name
    return params_dict

def get_run_name(urid):
    '''
    extract the run name from the json param file using urid 
    '''

    params_dict = get_sim_param_dict(urid)
    
    run_name = params_dict.get('savename', None)
    
    return run_name

def combo_name_to_short_name(combo_name):

    matches = re.findall(string = combo_name,
                         pattern = r'(?:bc|mt).*')

    # if this matches an actual mt/bc clade assignment, return it
    if len(matches) > 0:
        return matches[0]

    # else if this is some indicator like alive/not or induced/not, return None as an indicator to later skip
    return None

def get_tms(urid, path_to_output_dir = path_to_output_dir):

    # run_spec_dir_path = Path(f'/dartfs/rc/lab/M/McKennaLab/projects/Aidan/simulations/r_sim_clean/output/run_specs/{urid}/')
    run_spec_dir_path = path_to_output_dir / f'run_specs/{urid}'
    run_spec_file_names = [f for f in os.listdir(run_spec_dir_path) if f.endswith('.json')]
    with open(run_spec_dir_path / run_spec_file_names[0], 'r') as f:
        run_spec_dict = json.load(f)
    
    uninduced_tm = run_spec_dict['cell_type_dict']['uninduced_transition_matrix']
    induced_tm = run_spec_dict['cell_type_dict']['induced_transition_matrix']

    return uninduced_tm, induced_tm

# Data classes


class TimepointRes():
    '''
    stores ScorematRes objects with a unified cell_pop
    '''

    def __init__(self, urid, run_name, timepoint, scoremat_res_list = []):
        # path_to_cellpop = Path(f'/dartfs/rc/lab/M/McKennaLab/projects/Aidan/simulations/r_sim_clean/output/cell_populations/{urid}/cell_population_{run_name}_time_{timepoint}.json')
        path_to_cellpop = path_to_output_dir / f'cell_populations/{urid}/cell_population_{run_name}_time_{timepoint}.json'
        # filter cellpop to only include terminal alive cells at this timepoint:
        with open(path_to_cellpop, 'r') as f:
            self.cell_pop = json.load(f)

        # upon creation of timepoint res, assign differentiation_induced and editing_induced to their own clades in clade_assignments
        # this makes the structure compatible with new_plot_cass_tree()
        for cell, celldict in self.cell_pop.items():
            if 'clade_assignments' not in celldict.keys():
                celldict['clade_assignments'] = {}
            celldict['clade_assignments']['clade_differentiation_induced'] = celldict['induced_differentiation']
            celldict['clade_assignments']['clade_editing_induced'] = celldict['induced_editing']
            celldict['clade_assignments']['clade_alive'] = celldict['alive']
            celldict['clade_assignments']['clade_terminal'] = celldict['terminal']
            celldict['clade_assignments']['clade_celltype'] = celldict['celltype']

        alive_terminal_cellpop = {linstring: lindict for linstring, lindict in self.cell_pop.items() if lindict['terminal'] and lindict['alive']}
        self.alive_terminal_pop = alive_terminal_cellpop

        # path_to_gt_tree = Path(f'/dartfs/rc/lab/M/McKennaLab/projects/Aidan/simulations/r_sim_clean/output/processed_newicks/{urid}/ground_truth_tree_{run_name}_time_{timepoint}.newick')
        path_to_gt_tree = path_to_output_dir / f'processed_newicks/{urid}/ground_truth_tree_{run_name}_time_{timepoint}.newick'
        self.gttree_path = path_to_gt_tree
        self.gt_tree = read_in_tree(path_to_gt_tree)

        self.scoremat_res_list = []

        
    def new_plot_cass_tree(self,
                            color_by,
                           return_html_string = True,
                          include_alive_status_ring = False,
                          clade_colors = None,
                          interactive = True,
                          param_combo = None,
                          image_format = 'png'):
        
        if isinstance(color_by, str):
            color_by = [color_by]

        base_keep_cols = ['celltype', 'terminal', 'alive']
        keep_cols = base_keep_cols + color_by

        # if we're including alive_status on the tree plot, we DON'T want to filter out dead cells here
        if include_alive_status_ring:
            # find this data for all terminal cells regardless of alive status
            filt_dict = {
                linstring: [inner_dict.get(col, None) for col in base_keep_cols] + [inner_dict['clade_assignments'].get(col, f'NAclade_{col}') for col in color_by]
                for linstring, inner_dict in self.cell_pop.items() if inner_dict['terminal']
                # for linstring, inner_dict in self.alive_terminal_pop.items()
            }
        # if we DON'T want to include alive status as a feature on the plots, we filter out dead cells here in addition to the normal non-terminal filtering
        elif not include_alive_status_ring:
            filt_dict = {
                linstring: [inner_dict.get(col, None) for col in base_keep_cols] + [inner_dict['clade_assignments'].get(col, f'NAclade_{col}') for col in color_by]
                for linstring, inner_dict in self.cell_pop.items() if inner_dict['terminal'] and inner_dict['alive']
                # for linstring, inner_dict in self.alive_terminal_pop.items()
            }

        last_linstring = next(reversed(self.alive_terminal_pop))
        last_linstring_dict = self.alive_terminal_pop[last_linstring]
        
        filt_df = pd.DataFrame.from_dict(
            filt_dict, orient = 'index', columns = keep_cols
        )

        # get all unique clades across all cells
        groups = list(set([clade for cell, celldict in self.cell_pop.items() for clade in celldict['clade_assignments'].values()]))

        try:
            # sort groups alphabetically
            groups_sorted = sorted(groups)
        except TypeError as e: # arises if, e.g., we have boolean clade assignment labels
            groups_sorted = groups
            

        # per-column palettes: each color_by column gets its own maximally-spaced HLS palette
        group_to_color_dict = {}
        group_to_ind = {}
        idx = 0
        for col in color_by:
            col_vals = sorted(filt_df[col].unique(), key = str)
            col_palette = sns.color_palette('hls', len(col_vals))
            for v, c in zip(col_vals, col_palette):
                if v not in group_to_color_dict:
                    group_to_color_dict[v] = rgb2hex(c)
                    group_to_ind[v] = idx
                    idx += 1

        if include_alive_status_ring:
            alive_to_color_dict = {
                'alive': '#0000FF',
                'dead': '#000000'
            }

            # add alive binary column to filt_df:
            filt_df['inner_color'] = filt_df['alive'].map(lambda x: 'alive' if bool(x) else 'dead')
            
        else:
            alive_to_color_dict = None

        if param_combo is None:
            gt_tree_copy = read_in_tree(self.gttree_path)
        else:
            matching = [sr for sr in self.scoremat_res_list if sr.short_name == param_combo]
            if not matching:
                raise ValueError(
                    f'No ScorematRes found with short_name == {param_combo!r}. '
                    f'Available: {[sr.short_name for sr in self.scoremat_res_list]}'
                )
            gt_tree_copy = read_in_tree(matching[0].recontree_path)


     

        keep_linstrings = prune_dead_lineages(cellpop_dict = self.cell_pop)
        keep_nodes = [node for node in gt_tree_copy.traverse('preorder') if node.name in keep_linstrings]
        gt_tree_copy.prune(keep_nodes, preserve_branch_length = True)
        # print(f'pre prune, {filt_df.shape=}')
        filt_df = filt_df[filt_df.index.isin(keep_linstrings)]
        # print(f'post prune, {filt_df.shape=}')
        
      

        # if color internal nodes and descendants by some feature:
        if clade_colors is not None:
            if isinstance(clade_colors, dict):
                # pre-built linstring -> hex color dict: pass through directly
                pass
            else:
                # string key: try celldict first, then clade_assignments as fallback
                def _lookup_clade_val(celldict, key):
                    if key in celldict:
                        return celldict[key]
                    return celldict.get('clade_assignments', {}).get(key)

                linstring_to_val = {
                    linstring: _lookup_clade_val(celldict, clade_colors)
                    for linstring, celldict in self.cell_pop.items()
                    if linstring in keep_linstrings
                }
                unique_vals = sorted(
                    set(v for v in linstring_to_val.values() if v is not None),
                    key = str
                )
                if set(unique_vals) <= {True, False}:
                    # boolean key: existing red/black behavior
                    clade_colors = {ls: '#ff0000' if v else '#000000'
                                   for ls, v in linstring_to_val.items()}
                else:
                    # multi-value key: one unique hls color per distinct clade value
                    pal = sns.color_palette('hls', len(unique_vals))
                    val_to_color = {v: rgb2hex(c) for v, c in zip(unique_vals, pal)}
                    clade_colors = {ls: val_to_color.get(v, '#888888')
                                   for ls, v in linstring_to_val.items()}

        
        fig = build_cassiopeia_figure(
            tr = gt_tree_copy,
            cell_state_filt_df = filt_df,
            group_to_color_dict = group_to_color_dict,
            inner_color_dict = alive_to_color_dict,
            group_to_ind = group_to_ind,
            meta_columns = list(color_by),
            return_html_string = return_html_string,
            clade_colors = clade_colors,
            interactive = interactive,
            image_format = image_format,
        )

        # # Store a tiny bit of metadata on the FigureWidget for possible later use
        # fig._metadata = {
        #     '"group_to_color_dict': group_to_color_dict,
        #     'group_to_ind': group_to_ind,
        #     'cell_state_filt_df': filt_df,
        # }

        return fig


    def interactive_plot(self,
                         color_options,
                         default_selection = None,
                         include_alive_status_ring = False):
        
        default_set = set(default_selection or [])
        cb_width = '1000px' # change to bigger size if clade names are cut off
        checkboxes = [
            widgets.Checkbox(
                value = (opt in default_set),   # True  pre-checked
                description = opt,
                indent = False,
                layout = widgets.Layout(width = cb_width)
            )
            for opt in color_options
        ]
    
        # Put the checkboxes in a vertical box (you could also use HBox/GridBox)
        selector = widgets.VBox(checkboxes)
        out = widgets.Output()
    
        
        def currently_selected():
            return [cb.description for cb in checkboxes if cb.value]
    
        def refresh(change = None):
            selected = currently_selected()
            if not selected:
                with out:
                    out.clear_output()
                    print('At least one color must be given')
                return
    
            # Call the public method that builds the FigureWidget
            fig = self.new_plot_cass_tree(
                color_by = selected,
                include_alive_status_ring = include_alive_status_ring,
            )
            with out:
                out.clear_output()
                display(fig)

        # for each selected checkbox, refresh
        for cb in checkboxes:
            cb.observe(refresh, names = 'value')
    
        # show initial plot
        refresh()
    
        return widgets.VBox([selector, out])

class ScorematRes():
    '''
    stores results for a given score matrix
    '''

    def __init__(self, urid, timepoint, run_name, scoremat_path = None,  timepoint_res = None,
                tree = None, 
                 # cell_pop = None, alive_and_terminal_pop = None,
                recontree_path = None, gt_shortcut = False):
        '''
        initialze with run features and paths to results
        '''

        # abbreviated for gt results (the rest can be skipped):
        self.urid = urid
        self.timepoint = timepoint
        self.timepoint_res = timepoint_res

        # need dp version of gt tree even in shortcut
        # self.gttree_path = Path(f'/dartfs/rc/lab/M/McKennaLab/projects/Aidan/simulations/r_sim_clean/output/processed_newicks/{urid}/ground_truth_tree_{run_name}_time_{timepoint}.newick')
        self.gttree_path = path_to_output_dir / f'processed_newicks/{urid}/ground_truth_tree_{run_name}_time_{timepoint}.newick'
        self.gt_tree = read_in_tree(self.gttree_path)
        self.dp_gt_tree = dp.Tree.get(path = self.gttree_path,
                          schema = 'newick',
                         preserve_underscores = True)

        if not gt_shortcut:
        
            scoremat_path_stem = Path(scoremat_path).stem
    
            # use scoremat_path to extract run details
            mt_or_bc = extract_mt_or_bc(path = scoremat_path_stem)
            if len(mt_or_bc) == 1:
                self.mt_or_bc = mt_or_bc[0]
                self.joint = False
            else:
                self.mt_or_bc = 'joint'
                self.joint = True
            
            num_ints = extract_ints(path = scoremat_path_stem, joint = self.joint)
            rec_frac = extract_frac(path = scoremat_path_stem, joint = self.joint)

            cell_rec_fracs = extract_cell_rec_fracs(scoremat_path_stem)
    
    
            # only mt accepts af and bin as params
            if self.mt_or_bc != 'bc':
                af = extract_af(path = scoremat_path_stem, joint = self.joint)
                bin = extract_bin(path = scoremat_path_stem, joint = self.joint)
            else:
                af = 0
                bin = 'F'
    
            self.num_ints = num_ints
            self.rec_frac = rec_frac
            self.af = af
            self.bin = bin
    
            self.run_name = run_name
        
            self.scoremat = pd.read_csv(scoremat_path, index_col = 0)
    
            self.short_name = f"{self.mt_or_bc}_{self.num_ints}ints_{rec_frac}RP_{af}AF_{bin}B_{'-'.join(cell_rec_fracs)}"
            
            
            # get matching recon tree path:
            # path_to_recon_tree_dir = Path(f'/dartfs/rc/lab/M/McKennaLab/projects/Aidan/simulations/r_sim_clean/output/recon_trees/{urid}/')
            path_to_recon_tree_dir = path_to_output_dir / f'recon_trees/{urid}/'
            scoremat_path_stem = Path(scoremat_path).stem

            try:
                self.recontree_path = path_to_recon_tree_dir / Path(f'score_{str(scoremat_path_stem)}') / Path(f'score_{str(scoremat_path_stem)}.treefile')
                self.recon_tree = read_in_tree(self.recontree_path)
                
                # temp solution to acclereate score calculation: also read in dendropy trees for gt and recon
                self.dp_recon_tree = dp.Tree.get(path = self.recontree_path,
                                  schema = 'newick',
                                 preserve_underscores = True)
            except: # if there are no trees associated with this scoremat, that's ok
                self.recontree_path = None
                self.recon_tree = None
                self.dp_recon_tree = None

        
    def add_clade_dict_to_cellpop(self, cellpop_attr,
                                  cell_to_group_dict):

        # check to see if color_by is in the cell-level dicts yet (quick & dirty check):
        first_linstring, first_linstring_dict = next(iter(self.timepoint_res.cell_pop.items()))

        # in case we want to color by multiple metadata factors
        if isinstance(cellpop_attr, str):
            cellpop_attr = [cellpop_attr]
        for this_attr in cellpop_attr:
            if this_attr not in first_linstring_dict['clade_assignments'].keys():
                for cellname, inner_dict in self.timepoint_res.cell_pop.items():
                    inner_dict['clade_assignments'][this_attr] = cell_to_group_dict.get(cellname, None)
    
        
    
    def make_clade_clustermap(self, 
                              cell_clade_dict,
                             order_by_clade = True):

        '''
        Plot a clustermap of cell x unique mutations


        '''

        
        clade_counts = dict(Counter(cell_clade_dict.values()))
        clade_counts = {k:v for k, v in sorted(clade_counts.items(), key = lambda x: x[1], reverse = True)}
        # sorted_nbl_clades = sorted(list(set(med_nbl_clades.values())), key = lambda clade: int(clade.split('_')[1]))
        sorted_clades = sorted(list(set(cell_clade_dict.values())), key = lambda clade: clade_counts[clade], reverse = True)
        palette = sns.color_palette('muted', n_colors = len(sorted_clades))
        clade_to_color_dict = dict(zip(sorted_clades, palette))

        hex_dict = {k:rgb2hex(v) for k, v in clade_to_color_dict.items()}
        # print(f'in make_clade_clustermap, clade_to_color_dict = {hex_dict}')
        
        
        # map score mat indices (cell names) to colors:
        clade_series = pd.Series([cell_clade_dict[cell_name] for cell_name in self.scoremat.index], index = self.scoremat.index)
        row_cols_list = [clade_to_color_dict[clade] for clade in clade_series]
        row_cols_series = pd.Series(row_cols_list, index = self.scoremat.index)
            
        if order_by_clade:
            order = clade_series.loc[clade_series.map(clade_counts).sort_values(ascending = False).index].index
            scores_to_plot = self.scoremat.loc[order]
            row_colors_to_plot = row_cols_series.loc[order]
        else:
            scores_to_plot = self.scoremat
            row_colors_to_plot = row_cols_series
        
        g = sns.clustermap(scores_to_plot,
                           method = 'average',
                           metric = 'euclidean',
                           cmap = 'Blues',
                           row_colors = row_colors_to_plot,
                           col_cluster = False,
                          row_cluster = False)

        g.ax_heatmap.set_title(f'{self.__str__()}')

        
        return clade_to_color_dict


    def plot_clade_distribution_heatmap(self, cats, varying = 'dist', use_branch_length = True):

        '''
        if varying == 'dist', we are plotting clade distributions across dist thresholds
        elif varying == 'num_clades', we are plotting clade distributions across numbers of clades

        '''
        df_list = []
    
        for cat in cats:
            if varying == 'dist':
                clades = self.assign_topological_clades(num_clades = None,
                                              dist_cutoff = cat, 
                                              use_branch_length = use_branch_length)
            elif varying == 'num_clades':
                clades = self.assign_topological_clades(num_clades = cat)
                
        
            countdict = dict(Counter(clades.values()))
            
            clade_count_df = pd.DataFrame.from_dict(countdict, orient = 'index')
            clade_count_df.columns = [f'cutoff_{dist_cutoff}']
                
            df_list.append(clade_count_df)
        
        stacked_res = pd.concat(df_list, axis = 1).fillna(0)
        
        sorted_inds = sorted(list(stacked_res.index), key = lambda x: int(x.split('_')[1]))
        stacked_res = stacked_res.loc[sorted_inds, ].T
        
        norm_stacked_res = copy.deepcopy(stacked_res)
        
        rowsums = norm_stacked_res.sum(axis = 1)
        norm_stacked_res = norm_stacked_res.div(rowsums, axis = 0)
    
        plt.figure(figsize = (10, 10))
        sns.heatmap(norm_stacked_res, cmap = 'Blues')
        # plt.title(f'Norm distribution across clades; use_branch_length == {use_branch_length}')

    def plot_dists_from_root(self, use_branch_length = True, return_median = False, show_fig = True):

        root_node = self.recon_tree.get_tree_root()

        if use_branch_length:
            dists_from_root = [node.get_distance(root_node) for node in self.recon_tree.traverse('preorder')]
        else:
            dists_from_root = [node.get_distance(root_node, topology_only = True) for node in self.recon_tree.traverse('preorder')]

        med_dist_val = np.round(np.median(dists_from_root), 3)
        
        if show_fig:
            ax = sns.histplot(dists_from_root)
            ax.axvline(x = med_dist_val, color = 'red', linestyle = '--')
            ax.text(x = med_dist_val, y = ax.get_ylim()[1]*.96, s = f'  {med_dist_val}', color = 'red')
            plt.title(f'{self.__str__()} \ndists from root\nuse_branch_length = {use_branch_length}')

        if return_median:
            return med_dist_val

    def clustermap_tree_workflow(self, cellpop_clade_attr_name, method, 
                                 use_median = True, num_clades = None, distance_metric = None,
                                 use_branch_length = True, show_clustermap = False,
                                 order_by_clade = True):

        '''
        if num_clades is None, will use median distance from root metric for clade assignments!
        if method == 'topological':
            if num_clades is None:
                use median distance as threshold where all distances calculated according to use_branch_length
            else:
                num_clades dictates number of clades
        elif method == 'hierarchical':
            if distance_metric == "tree_dist":
                distance matrix calculated depending on use_branch_length         
            elif distance_metric == "hamming":
                distance matrix calculated from cell-level score mats
                
        '''

        if method not in ['topological', 'hierarchical']:
            raise Exception ("method must be in ['topological', 'hierarchical']")

        if method == 'topological':
            if num_clades is None:
                med_dist = self.plot_dists_from_root(use_branch_length = use_branch_length, return_median = True,
                                             show_fig = show_intermediate_figs)
                varying = 'dist'
                cats = [med_dist]
            else:
                varying = 'num_clades'
                cats = [num_clades]
                med_dist = None

            if show_clustermap:
                self.plot_clade_distribution_heatmap(cats = cats, varying = 'dist', 
                                                     use_branch_length = use_branch_length)

            cell_clade_dict = self.assign_topological_clades(num_clades = num_clades, dist_cutoff = med_dist, 
                                                               use_branch_length = use_branch_length)
        
            
        elif method == 'hierarchical':
            cell_clade_dict = self.assign_hierarchical_clades(distance_metric = distance_metric,
                                                          num_clusters = num_clades,
                                                          use_branch_length = use_branch_length)

                
        clade_to_color_dict = self.make_clade_clustermap(cell_clade_dict = cell_clade_dict,
                                                             order_by_clade = order_by_clade)
    
        # print(f'Counter(cell_clade_dict.values()) == {Counter(cell_clade_dict.values())}')
        hex_dict = {k: rgb2hex(v) for k, v in clade_to_color_dict.items()}
        # print(f'clade_to_color_dict = {hex_dict}')
        # self.plot_cass_tree(color_by = cellpop_clade_attr_name, 
        #                    # cell_to_group_dict = cell_clade_dict, 
        #                     group_to_color_dict = clade_to_color_dict)

        # update
        self.new_plot_cass_tree(color_by = cellpop_clade_attr_name)
        
        return cell_clade_dict
       

    
    def __str__(self):
        growing_str = ''
        if self.mt_or_bc == 'mt' or self.mt_or_bc == 'joint':
            growing_str += f'mt score mat with {self.num_ints} ints, AF = {self.af}, recovery fraction = {self.rec_frac}, bin = {self.bin}'
        if self.mt_or_bc == 'bc' or self.mt_or_bc == 'joint':
            growing_str += f'bc score mat with {self.num_ints} ints, AF = {self.af}, recovery fraction = {self.rec_frac}, bin = {self.bin}'
        return growing_str


    def assign_scoremat_clades(self,
                               attr_name,
                               return_clade_dict = True,
                              min_clade_size = 10,
                              return_defining_sig = False,
                              num_clades = None):


        for linstring, cell_dict in self.timepoint_res.cell_pop.items():
            if 'clade_assignments' not in cell_dict.keys():
                cell_dict['clade_assignments'] = {}

        
        scoremat = copy.deepcopy(self.scoremat)
        shortname = self.short_name

    
        row_signature = scoremat.apply(lambda r: tuple(r.astype(int)), axis = 1)
        
        clone_groups = row_signature.groupby(row_signature).groups
        signature_counts = {sig: len(cells) for sig, cells in clone_groups.items()}

        sorted_sigs = sorted(clone_groups.keys(), key = lambda clone: signature_counts[clone], reverse = True)
        
        cell_to_clone = {}
        # clone_id_to_cells = {}
        clone_id_to_cells = defaultdict(list)
        clade_to_sig = {}
        clade_freqs = {}

        if num_clades is not None:
            center_sigs = sorted_sigs[:num_clades]
            for clade_id, sig in enumerate(center_sigs, start = 1):
                clade_name = f'clade_{clade_id}'
                cells = clone_groups[sig]
                clade_to_sig[clade_name] = sig
                clade_freqs[clade_name] = signature_counts[sig]
                for cell in cells:
                    cell_to_clone[cell] = clade_name
                clone_id_to_cells[clade_name].extend(cells)
    
            for sig in sorted_sigs[num_clades:]:
                cells = clone_groups[sig]
                best_clade, best_dist, best_freq = None, None, -1
                for clade_name, center_sig in clade_to_sig.items():
                    dist = sum(a != b for a, b in zip(sig, center_sig))
                    if (best_dist is None or dist < best_dist or
                            (dist == best_dist and clade_freqs[clade_name] > best_freq)):
                        best_dist = dist
                        best_clade = clade_name
                        best_freq = clade_freqs[clade_name]
                for cell in cells:
                    cell_to_clone[cell] = best_clade
                clone_id_to_cells[best_clade].extend(cells)
                clade_freqs[best_clade] += signature_counts[sig]


        else:
            next_clade_id = 1
            for sig in sorted_sigs:
                cells = clone_groups[sig]
                n_cells = signature_counts[sig]
                if n_cells >= min_clade_size:
                    clade_name = f'clade_{next_clade_id}'
                    next_clade_id += 1
    
                    for cell in cells:
                        cell_to_clone[cell] = clade_name
    
                    # clone_id_to_cells[clade_name] = clone_id_to_cells.get(clade_name, []) + list(cells)
                    clone_id_to_cells[clade_name].extend(cells)
                    # print('post get')
                    clade_to_sig[clade_name] = sig
                    clade_freqs[clade_name] = n_cells
    
                    continue
    
                best_clade = None
                best_dist = None
                best_freq = -1
    
    
              
                if not clade_to_sig:
                    clade_name = f'clade_{next_clade_id}'
                    next_clade_id += 1
    
                    for cell in cells:
                        cell_to_clone[cell] = clade_name
    
                    # clone_id_to_cells[clade_name] = clone_id_to_cells.get(clade_name, []) + list(cells)
                    clone_id_to_cells[clade_name].extend(cells)
                    # print('post get')
                    clade_to_sig[clade_name] = sig
                    clade_freqs[clade_name] = n_cells
    
                    continue
    
                    
    
                for old_clade, old_sig in clade_to_sig.items():
    
                    
                    # dist given by hamming dist:
                    dist = sum(a != b for a, b in zip(sig, old_sig))
                    
                    if (best_dist is None) or (dist < best_dist) or (dist == best_dist and clade_freqs[old_clade] > best_freq):
                        best_dist = dist
                        best_clade = old_clade
                        best_freq = clade_freqs[old_clade]
    
                for cell in cells:
                    cell_to_clone[cell] = best_clade
                clone_id_to_cells[best_clade].extend(cells)
                clade_freqs[best_clade] += n_cells

        
        for linstring, clade in cell_to_clone.items():
            self.timepoint_res.cell_pop[linstring]['clade_assignments'][attr_name] = clade
            # print(f'{linstring}: {clade}')

        if return_clade_dict and return_defining_sig:
            return cell_to_clone, clade_to_sig
        elif return_clade_dict:
            return cell_to_clone
        elif return_defining_sig:
            return clade_to_sig
        else:
            return None

    def assign_random_clades(self,
                             attr_name,
                             num_clades,
                             tree = None,
                             return_clade_dict = True):

        '''
        map each leaf to one of num_clades clades randomly
        
        '''

        np.random.seed(42)
        
        if tree is None:
            tree = self.recon_tree

        # initialize an empty clade dict attr for each cell in cell_pop if not already present:
        for linstring, cell_dict in self.timepoint_res.cell_pop.items():
            if 'clade_assignments' not in cell_dict.keys():
                cell_dict['clade_assignments'] = {}
        
        all_leaf_linstrings = [leaf.name for leaf in tree.iter_leaves()]

        clade_names = [f'clade_{clade_num}' for clade_num in range(num_clades)]
        rand_clade_assignments = random.choices(clade_names, k = len(all_leaf_linstrings))

        terminal_leaf_to_clade = dict(zip(all_leaf_linstrings, rand_clade_assignments))

        for linstring, clade in terminal_leaf_to_clade.items():
            self.timepoint_res.cell_pop[linstring]['clade_assignments'][attr_name] = clade

        if return_clade_dict:
            return terminal_leaf_to_clade

   
    def assign_topological_clades(self, 
                                  attr_name,
                                  tree = None,
                                  return_clade_dict = True,
                                  num_clades = None,
                                  dist_cutoff = None, 
                                  use_branch_length = False):

        '''
        can define clades according to some distance metric or by specifying fixed number of clades
    
        '''

        if tree is None:
            tree = self.recon_tree

        # initialize an empty clade dict attr for each cell in cell_pop if not already present:
        for linstring, cell_dict in self.timepoint_res.cell_pop.items():
            if 'clade_assignments' not in cell_dict.keys():
                cell_dict['clade_assignments'] = {}
        
        # Avoid copy.deepcopy(tree): for trees with depth = thousands of divisions
        # deepcopy recurses to that depth and always hits Python's recursion limit.
        # Instead, use an id()-keyed dict so we never mutate the original tree at all.
        root_node = tree.get_tree_root()

        # if specifying epxlicit number of nodes
        if num_clades is not None:

            cut_nodes = [root_node]

            while len(cut_nodes) < num_clades:

                elig_split_nodes = [node for node in cut_nodes if not node.is_leaf()]

                if len(elig_split_nodes) == 0:
                    break

                # split the node with the most downstream leaves
                split_cand = max(elig_split_nodes, key = lambda node: len(node.get_leaves()))

                if split_cand is None:
                    break

                # replace split node with its daughters
                cut_nodes.remove(split_cand)
                cut_nodes.extend(split_cand.children)

            # ensure root node remains in cut_nodes
            if root_node not in cut_nodes:
                cut_nodes.append(root_node)

        else:
            cut_nodes = []

            for node in tree.traverse('preorder'):
                if node.is_root():
                    continue
                if use_branch_length:
                    root_dist = node.get_distance(root_node)
                    parent_dist = node.up.get_distance(root_node)
                else:
                    root_dist = node.get_distance(root_node, topology_only = True)
                    parent_dist = node.up.get_distance(root_node, topology_only = True)

                if parent_dist < dist_cutoff <= root_dist:
                    cut_nodes.append(node)

            if len(cut_nodes) == 0:
                cut_nodes = [root_node]

            if root_node not in cut_nodes:
                cut_nodes.append(root_node)

        # Map node object id -> clade label (avoids add_feature / tree mutation)
        cut_node_clade = {id(cut_node): f'clade_{cut_num}'
                          for cut_num, cut_node in enumerate(cut_nodes)}

        # Walk each leaf up to its nearest cut-node ancestor
        terminal_leaf_to_clade = {}
        for leaf in tree.iter_leaves():
            parent = leaf
            while parent is not None and id(parent) not in cut_node_clade:
                parent = parent.up
            cl = cut_node_clade.get(id(parent), 'unlabeled')
            terminal_leaf_to_clade[leaf.name] = cl

        for linstring, clade in terminal_leaf_to_clade.items():
            self.timepoint_res.cell_pop[linstring]['clade_assignments'][attr_name] = clade

        if return_clade_dict:
            return terminal_leaf_to_clade


    def get_hamming_dist_df(self):
        # Suppose `mut_mat` is a DataFrame (cells  mutation sites) of 0/1
        dist_arr = squareform(pdist(self.scoremat.values, metric = 'hamming'))
        dist_df = pd.DataFrame(dist_arr, index = self.scoremat.index, columns = self.scoremat.index)
    
        return dist_df

    
    # slow way:
    def get_tree_dist_df(self, tree = None, use_branch_length = True):
        if tree is None:
            tree = self.recon_tree
        leaves = [leaf.name for leaf in tree.iter_leaves()]
    
        dist_mat = np.zeros((len(leaves), len(leaves)))
    
        # for i, leaf_i in enumerate(leaves):
            # for j, leaf_j in enumerate(leaves):
        for i in range(1, len(leaves)):
            leaf_i = leaves[i]
            for j in range(i):
                leaf_j = leaves[j]
                # if j >= i:
                #     continue
                if use_branch_length:
                    dist = tree.get_distance(leaf_i, leaf_j)
                else:
                    dist = tree.get_distance(leaf_i, leaf_j, topology_only = True)
    
                dist_mat[i, j] = dist_mat[j, i] = dist
    
        dist_df = pd.DataFrame(dist_mat, index = leaves, columns = leaves)
    
        return dist_df

    def dp_get_tree_dist_df(self, tree_path = None, dp_tree = None, use_branch_length = True):

        # taxon_namespace = dp.TaxonNamespace()
    
        if tree_path is not None:
            dp_tree = dp.Tree.get(path = tree_path,
                                  schema = 'newick',
                                 preserve_underscores = True)
            
        pdm = dp_tree.phylogenetic_distance_matrix()
        
        # Gather leaf labels once
        leaf_tax = list(dp_tree.leaf_node_iter())
        leaf_labels = [leaf.taxon.label for leaf in dp_tree.leaf_node_iter()]
        n = len(leaf_labels)
        
        dist_mat = np.empty((n, n), dtype = float)
        
        for i, a in enumerate(leaf_tax):
            for j, b in enumerate(leaf_tax):
                # dist_mat[i, j] = pdm(a.taxon, b.taxon)
                dist_mat[i,j] = pdm.distance(a.taxon, b.taxon, is_weighted_edge_distances = use_branch_length)
        
        dist_df = pd.DataFrame(dist_mat, index = leaf_labels, columns = leaf_labels)
        
        return dist_df
    
    def cut_hierarchical(self, dist_df, num_clusters):
    
        condensed = squareform(dist_df.values)
        Z = linkage(condensed, method = 'average')
        cluster_labs = fcluster(Z, t = num_clusters, criterion = 'maxclust')
        cluster_lab_names = [f'clade_{int(lab)}' for lab in cluster_labs]
        clusters_dict = dict(zip(dist_df.columns, cluster_lab_names))
        return clusters_dict

    
    def assign_hierarchical_clades(self, 
                                  distance_metric,
                                  num_clusters,
                                   attr_name,
                                   tree = None,
                                   return_clade_dict = True,
                                  use_branch_length = True):
        '''
        assign clades by performing hierarchical clustering of distance matrix
        distance matrix is computed by distance_metric
        distance_metric can be hamming dist between cell-level scoremats or 
            topological or branch-length distances
        '''
        # print(f'beginning of assign_hierarch, id(self.timepoint_res) == {id(self.timepoint_res)}')

        # should not fire
        if tree is None:
            tree = self.recon_tree

        # initialize an empty clade dict attr for each cell in cell_pop if not already present:
        for linstring, cell_dict in self.timepoint_res.cell_pop.items():
            if 'clade_assignments' not in cell_dict.keys():
                cell_dict['clade_assignments'] = {}
        
        if distance_metric == 'tree_dist':
            dist_df = self.dp_get_tree_dist_df(use_branch_length = use_branch_length, dp_tree = tree)
        elif distance_metric == 'hamming':
            dist_df = self.get_hamming_dist_df()
        else:
            raise Exception('dist_df must be one of ["tree_dist" or "hamming"]')

        clade_dict = self.cut_hierarchical(dist_df = dist_df,
                                     num_clusters = num_clusters)

        # add elements of clade dict to TimepointRes under attr_name
        for linstring, clade in clade_dict.items():
            self.timepoint_res.cell_pop[linstring]['clade_assignments'][attr_name] = clade

        # print(ex_timepoint_res.cell_pop)

        if return_clade_dict:
            return clade_dict

# Clade and clustering helpers


def jaccard_dicts(dict_a, dict_b):
    
    same_clade = {linstring for linstring in list(dict_a.keys())+list(dict_b.keys()) \
                  if dict_a[linstring] == dict_b[linstring]}
    
    total_linstrings = dict_a.keys() | dict_b.keys()
    
    return len(same_clade) / len(total_linstrings) 

def pairwise_jaccards(dict_list1, dict_list2,
                     dict1_names = None, dict2_names = None):

    '''
    Find all pairwise jaccards between two lists of linstring:clade dicts
    Can always pass in the same list to both args
    '''
    j_arr = np.empty((len(dict_list1), len(dict_list2)))

    for i in range(len(dict_list1)):
        for j in range(len(dict_list2)):
            j_arr[i,j] = jaccard_dicts(dict_list1[i],
                                  dict_list2[j])

    rownames = dict1_names if dict1_names is not None else range(len(dict_list1))
    colnames = dict2_names if dict2_names is not None else range(len(dict_list2))
    return pd.DataFrame(j_arr, index = rownames, columns = colnames)

def fetch_clusters(dict1, dict2):
    sorted_keys = sorted(set(dict1) & set(dict2))
    dict1_clusters = np.array([dict1[k] for k in sorted_keys])
    dict2_clusters = np.array([dict2[k] for k in sorted_keys])
    return dict1_clusters, dict2_clusters

def pairwise_concordance(dict_list1, dict_list2,
                         dict1_names = None, dict2_names = None):

    rownames = dict1_names if dict1_names is not None else range(len(dict_list1))
    colnames = dict2_names if dict2_names is not None else range(len(dict_list2))


    j_arr = np.empty((len(dict_list1), len(dict_list2)))
    for i in range(len(dict_list1)):
        for j in range(len(dict_list2)):
            clusters1, clusters2 = fetch_clusters(dict_list1[i],
                                                  dict_list2[j])
            j_arr[i,j] = adjusted_rand_score(clusters1, clusters2)

    return pd.DataFrame(j_arr, index = rownames, columns = colnames)

def build_cassiopeia_figure(
        tr,
        cell_state_filt_df,
        group_to_color_dict,
        group_to_ind,
        meta_columns,
        inner_color_dict = None,
        return_html_string = True,
        clade_colors = None,
        interactive = True,
        image_format = 'png'):


    G = nx.DiGraph()
    def rename_internal_node(n):
        return n.name if n.name else f'internal_{id(n)}'
    def add_edges(p):
        for c in p.children:
            G.add_edge(rename_internal_node(p), rename_internal_node(c))
            add_edges(c)
    add_edges(tr)

    cass_tree = CassiopeiaTree(tree = G)


    if inner_color_dict is not None:
        if 'inner_color' not in cell_state_filt_df.columns:
            cell_state_filt_df['inner_color'] = cell_state_filt_df['alive'].map(
                lambda x: 'alive' if bool(x) else 'dead'
            )
    
    cass_tree.cell_meta = cell_state_filt_df
    
    custom_cmap = sns.color_palette('muted', n_colors = len(group_to_color_dict))

    # Drop constant-value columns - cassiopeia's continuous colorstrip crashes with min==max
    varying_meta = [col for col in meta_columns
                    if col not in cell_state_filt_df.columns
                    or cell_state_filt_df[col].nunique() > 1]
    fig = go.Figure(
        plot_plotly(cass_tree,
                    meta_data = varying_meta,
                    categorical_cmap = ListedColormap(group_to_color_dict.values()),
                    value_mapping = group_to_ind,
                   clade_colors = clade_colors)
    )
    fig.update_layout(template = 'plotly_white')


    if inner_color_dict is None:
        if interactive:
            if return_html_string:
                html_str = pio.to_html(fig, include_plotlyjs = 'cdn')
                return html_str
            return fig
    
    
        static_bytes = pio.to_image(fig, format = image_format)

        if return_html_string:
            # Embed the raw bytes in a base-64 data-uri.
            b64 = base64.b64encode(static_bytes).decode()

            html_img = f'<img src="data:application/pdf;base64,{b64}" alt="Cassiopeia tree (static)"/>'
    
            return html_img
        else:
            # Return bytes directly - the caller can write them to disk or whatever.
            return static_bytes

    inner_lab_to_ind = {lab: i for i, lab in enumerate(inner_color_dict.keys())}
    leaf_names = fig.data[0].text

    # strip formatting away from leaf names, if present
    leaf_names = [re.findall(r'.*<br>([\d_]+)', leaf_name)[0] for leaf_name in leaf_names if '<br>' in leaf_name]

    inner_inds = [inner_lab_to_ind[cell_state_filt_df.loc[leaf, 'inner_color']] for leaf in leaf_names]

    inner_cmap = ListedColormap(list(inner_color_dict.values()))

    x = fig.data[0].x
    y = fig.data[0].y

    fig.add_trace(
        go.Scatter(x = x, y = y, mode = 'markers',
                   marker = dict(color = inner_inds,
                                 colorscale = inner_cmap.colors,
                                 cmin = 0,
                                 cmax = len(inner_color_dict) - 1,
                                 size = 8,
                                 line = dict(width = 0)))
    )

    if interactive:
        if return_html_string:
            html_str = pio.to_html(fig, include_plotlyjs = 'cdn')
            return html_str
        return fig


    static_bytes = pio.to_image(fig, format = 'png')

    if return_html_string:
        # Embed the raw bytes in a base-64 data-uri.
        b64 = base64.b64encode(static_bytes).decode()

        html_img = f'<img src="data:application/pdf;base64,{b64}" alt="Cassiopeia tree (static)"/>'

        return html_img
    else:
        # Return bytes directly - the caller can write them to disk or whatever.
        return static_bytes

# Pipeline helpers


def make_clade_df(timepoint_res, urid, timepoint):
    # path_toclade_dfs = Path(f'/dartfs/rc/lab/M/McKennaLab/projects/Aidan/simulations/r_sim_clean/output/clade_assignment_dfs/{urid}/{timepoint}')
    path_toclade_dfs = path_to_output_dir / f'clade_assignment_dfs/{urid}/{timepoint}'
    os.makedirs(path_toclade_dfs, exist_ok = True)

    cell_pop = timepoint_res.cell_pop
    dict_of_dicts = {cell:celldict['clade_assignments'] for cell, celldict in cell_pop.items()}
    # removal of cells not assigned to clades:
    dict_of_dicts = {k:v for k, v in dict_of_dicts.items() if len(v) > 0}
    clade_df = pd.DataFrame.from_dict(dict_of_dicts, orient = 'index')

    clade_df.to_csv(os.path.join(path_toclade_dfs, f'clade_df_{urid}_t{timepoint}.csv'))

    return clade_df

def scoremat_assign(scoremat, min_clade_size):
    """
    Pure function: assign scoremat clades from a binary allele DataFrame.

    Mirrors ScorematRes.assign_scoremat_clades() but requires no cell_pop,
    no trees, and no ScorematRes wrapper - takes the matrix directly.

    Cells are grouped by their allele signature (row tuple). Groups with
    n >= min_clade_size become named clades; smaller groups are assigned
    to the nearest existing clade by Hamming distance, with ties broken
    by clade frequency.

    Parameters
    ----------
    scoremat       : pd.DataFrame  rows=linstrings, cols=alleles, values in {0,1}
    min_clade_size : int

    Returns
    -------
    dict  {linstring: clade_name}
    """
    row_sig = scoremat.apply(lambda r: tuple(r.astype(int)), axis = 1)
    clone_groups = row_sig.groupby(row_sig).groups
    sig_counts = {sig: len(cells) for sig, cells in clone_groups.items()}
    sorted_sigs = sorted(clone_groups, key = lambda s: sig_counts[s], reverse = True)

    cell_to_clade = {}
    clade_to_sig = {}
    clade_freqs = {}
    next_id = 1

    for sig in sorted_sigs:
        cells = clone_groups[sig]
        n = sig_counts[sig]

        if n >= min_clade_size:
            name = f'clade_{next_id}'; next_id += 1
            clade_to_sig[name] = sig
            clade_freqs[name]  = n
            for cell in cells:
                cell_to_clade[cell] = name
            continue

        # below threshold: force-start from the largest group if nothing formed yet
        if not clade_to_sig:
            name = f'clade_{next_id}'; next_id += 1
            clade_to_sig[name] = sig
            clade_freqs[name]  = n
            for cell in cells:
                cell_to_clade[cell] = name
            continue

        # assign to nearest existing clade by Hamming; tie-break by clade size
        best, best_d, best_f = None, None, -1
        for cn, cs in clade_to_sig.items():
            d = sum(a != b for a, b in zip(sig, cs))
            if best_d is None or d < best_d or (d == best_d and clade_freqs[cn] > best_f):
                best, best_d, best_f = cn, d, clade_freqs[cn]

        for cell in cells:
            cell_to_clade[cell] = best
        clade_freqs[best] += n

    return cell_to_clade

def scoremat_clade_sweep(
    urid_list,
    min_clade_sizes,
    timepoints = None,
    base_dir = path_to_output_dir,
):
    """
    Assign scoremat clades for every combination of:
        urid  timepoint  scoremat CSV (AF/RP)  min_clade_size

    Does not load any trees, does not run EM, and does not require
    process_urid() or TimepointRes. Works directly from the binary
    allele CSVs in output/{urid}/score_mats/{urid}/matrices/csvs/.

    Parameters
    ----------
    urid_list       : list of int or str
    min_clade_sizes : list of int  e.g. [5, 10, 20]
    timepoints      : list of int, optional  - None = all available
    base_dir        : Path  - root of per-urid output dirs (default path_to_output_dir)

    Returns
    -------
    pd.DataFrame with columns:
        urid, timepoint, param_combo, linstring, clade
    where param_combo encodes both the scoremat filter (RP/AF) and min_clade_size,
    e.g. "scoremat_min5_mt_250ints_1.0RP_0.0AF_TB_1-1-1-1-1-1"
    """
    records = []

    for urid in urid_list:
        urid_str = str(urid)
        csv_dir = Path(base_dir) / urid_str / 'score_mats' / urid_str / 'matrices' / 'csvs'

        if not csv_dir.exists():
            print(f'[scoremat_clade_sweep] No CSV dir found for {urid_str}, skipping')
            continue

        csv_paths = sorted(csv_dir.glob('*.csv'))
        if timepoints is not None:
            csv_paths = [p for p in csv_paths
                         if extract_timepoint(p.stem) in timepoints]

        for csv_path in csv_paths:
            stem = csv_path.stem
            tp = extract_timepoint(stem)

            try:
                rp = extract_frac(path = stem, joint = False)
                af = extract_af(path = stem, joint = False)
                ints = extract_ints(path = stem, joint = False)
                bin_ = extract_bin(path = stem, joint = False)
                fracs = extract_cell_rec_fracs(stem)
            except Exception as e:
                print(f'[scoremat_clade_sweep] Could not parse {stem}: {e}')
                continue

            tb_str = '-'.join(fracs)
            scoremat = pd.read_csv(csv_path, index_col = 0)

            for min_size in min_clade_sizes:
                param_combo = (
                    f'scoremat_min{min_size}_mt_{ints}ints_'
                    f'{rp}RP_{af}AF_{bin_}B_{tb_str}'
                )
                cell_to_clade = scoremat_assign(scoremat, min_clade_size = min_size)

                for linstring, clade in cell_to_clade.items():
                    records.append((urid_str, tp, param_combo, linstring, clade))

    return pd.DataFrame(
        records,
        columns = ['urid', 'timepoint', 'param_combo', 'linstring', 'clade']
    )

def subtree_has_alive(linstring, cellpop_dict, seen_dict = None):

    '''
    include seen_dict to prevent redundant recursions
    '''
    if seen_dict is None:
        seen_dict = {}

    if linstring in seen_dict:
        return seen_dict[linstring]

    celldict = cellpop_dict[linstring]

    # if cell is alive, keep it
    if celldict['alive']:
        seen_dict[linstring] = True
        return True

    # then, if cell has no descendants, remove it
    if not celldict['descendants']:
        seen_dict[linstring] = False
        return False

    # repeat process for each descendant of the node
    for child in celldict['descendants']:
        if subtree_has_alive(child, cellpop_dict, seen_dict):
            seen_dict[linstring] = True
            return True

    # if no descendants of this linstring are alive
    seen_dict[linstring] = False
    return False

def prune_dead_lineages(cellpop_dict):
    alive_or_descendants = {linstring: subtree_has_alive(linstring, cellpop_dict) for linstring in cellpop_dict.keys()}
    keep_nodes = [linstring for linstring, keep in alive_or_descendants.items() if keep]
    return keep_nodes

def make_clade_celltype_df(urid, timepoint, urid_timepoint_res_dict, path_to_output_dir = path_to_output_dir):
    
    # path_toclade_dfs = Path(f'/dartfs/rc/lab/M/McKennaLab/projects/Aidan/simulations/r_sim_clean/output/clade_assignment_dfs/{urid}/{timepoint}')
    path_toclade_dfs = path_to_output_dir / f'clade_assignment_dfs/{urid}/{timepoint}'
    
    ex_at_cell_pop = urid_timepoint_res_dict[urid][timepoint].alive_terminal_pop
    
    exclade_df = pd.read_csv(path_toclade_dfs / f'clade_df_{urid}_t{timepoint}.csv', index_col = 0).reset_index(names = 'cell')
    # exclade_df.rename(columns = {'index': 'cell'}, inplace = True)
    
    meltedclade_df = pd.melt(exclade_df, var_name = 'param_combo',
                             value_name = 'clade', ignore_index = False,
                             id_vars = 'cell')
    
    celltype_df = pd.DataFrame.from_dict({cell:celldict['celltype'] for cell, celldict in ex_at_cell_pop.items()},
                                          orient = 'index', columns = ['celltype']).reset_index(names = 'cell')
    
    full_df = pd.merge(meltedclade_df, celltype_df, on = 'cell', how = 'inner')

    return full_df

def process_urid(urid,
                 param_dict_list,
                 clade_color_by = None,
                # inner_tree = 'highest_acc',
                timepoints = None,
                cluster_vmin = None,
                cluster_vmax = None,
                include_ari_labels = True):

    '''

    if clade_color_by == 
    
    # if inner_tree == 'gt':
    #     plot ground truth tree before adding color_by clade labels
    # elif inner_tree == 'highest_acc':
    #     plot tree with highest reconstruction accuracy before adding color_by clade labels
    # elif inner_tree in 
    

    '''

    run_name = get_run_name(urid)

    # path_to_cluster_maps = Path(f'/dartfs/rc/lab/M/McKennaLab/projects/Aidan/simulations/r_sim_clean/output/cluster_maps/{urid}')
    path_to_cluster_maps = path_to_output_dir / f'cluster_maps/{urid}'
    os.makedirs(path_to_cluster_maps, exist_ok = True)


    # path_to_clade_trees = Path(f'/dartfs/rc/lab/M/McKennaLab/projects/Aidan/simulations/r_sim_clean/output/clade_trees/{urid}')
    path_to_clade_trees = path_to_output_dir / f'clade_trees/{urid}'
    os.makedirs(path_to_clade_trees, exist_ok = True)

    # path_to_scoremats = Path(f'/dartfs/rc/lab/M/McKennaLab/projects/Aidan/simulations/r_sim_clean/output/score_mats/{urid}/matrices/csvs/')
    path_to_scoremats = path_to_output_dir / f'score_mats/{urid}/matrices/csvs/'
    all_scoremat_paths = [str(path) for path in path_to_scoremats.iterdir() if path.is_file()]
    
    
    # extract available timepoints for this urid 
    if timepoints is None:
        timepoints = sorted(set([extract_timepoint(sm_path) for sm_path in all_scoremat_paths]))

    growing_clade_dict_collection = {}

    timepoint_res_dict = {}

    for timepoint in timepoints:


        timepoint_res = TimepointRes(urid = urid,
                                     run_name = run_name, 
                                     timepoint = timepoint)

        timepoint_scoremat_paths = [path for path in all_scoremat_paths if f'time_{timepoint}' in path]

        # initialze timepoint-specific dict
        growing_clade_dict_collection[timepoint] = {}

        # will map full clade combo X scoremat names --> clade combo (avoids downstream recovery with regex)
        clade_name_to_param_combo_dict = {}


        # first generate clade assignments on the ground truth tree for each clade_param_dict in param_dict_list:
        # do this by creating temporary scoremat object associated with gt data so that assign_hierarchical_clades() can be called on GT data
        # not the cleanest way of doing this, but that's ok for now
        # note that gt res cannot undergo hamming dist logic

        gt_scoremat = ScorematRes(urid = urid,
                                  timepoint = timepoint,
                                  timepoint_res = timepoint_res,
                                 gt_shortcut = True,
                                 run_name = run_name)
        
        for clade_param_dict in param_dict_list:


            distance_metric = clade_param_dict.get('distance_metric', None)
            num_clades = clade_param_dict.get('num_clades', None)
            scheme = clade_param_dict.get('scheme', None)
            use_branch_length = clade_param_dict.get('use_branch_length', None)
            clade_combo_name = clade_param_dict.get('name', None)
            dist_cutoff = clade_param_dict.get('dist_cutoff', None)

            if distance_metric == 'hamming':
                continue

            name = f'{clade_combo_name}_gt'

            # assign clades to cellpop using res+param combo names
            if scheme == 'hierarchical':

                # can accelerate tree distance matrix calculations using the dendropy workaround to precompute distances
                # tree argument is not used when distance_metric == 'hamming' so we can provide dp_tree in either case
            
                clade_dict = gt_scoremat.assign_hierarchical_clades(distance_metric = distance_metric,
                                                              num_clusters = num_clades,
                                                              use_branch_length = use_branch_length,
                                                                tree = gt_scoremat.dp_gt_tree,
                                                              attr_name = name,
                                                              return_clade_dict = True)
                growing_clade_dict_collection[timepoint][name] = clade_dict
                    
            elif scheme == 'topological':

                # for topological, still want to use the ete3 tree
                clade_dict = gt_scoremat.assign_topological_clades(attr_name = name,
                                                          return_clade_dict = True,
                                                          num_clades = num_clades,
                                                                         tree = timepoint_res.gt_tree,
                                                          dist_cutoff = dist_cutoff, 
                                                          use_branch_length = use_branch_length)
                growing_clade_dict_collection[timepoint][name] = clade_dict

            elif scheme == 'random':
                clade_dict = gt_scoremat.assign_random_clades(attr_name = name,
                                                             num_clades = num_clades,
                                                             tree = timepoint_res.gt_tree,
                                                             return_clade_dict = True)

            # can't use gt data for a scoremat approach
            elif scheme == 'scoremat':
                pass
            

        # then generate clade assignments on each reconstructed ground truth tree for each clade_param_dict in param_dict_list
        for scoremat_path in timepoint_scoremat_paths:

            try:
                this_scoremat_res = ScorematRes(scoremat_path = scoremat_path,
                                                run_name = run_name,
                                                urid = urid,
                                                timepoint_res = timepoint_res,
                                               timepoint = timepoint)
            except FileNotFoundError as e:
                print(f'Unable to create ScorematRes: {e}')
                continue
                

            timepoint_res.scoremat_res_list.append(this_scoremat_res)
            # print(f'this_scoremat_res.mt_or_bc == {this_scoremat_res.mt_or_bc}')


            # subject each scoremat res at this timepoint to each clade-building approach
            for clade_param_dict in param_dict_list:
    
                distance_metric = clade_param_dict.get('distance_metric', None)
                num_clades = clade_param_dict.get('num_clades', None)
                scheme = clade_param_dict.get('scheme', None)
                use_branch_length = clade_param_dict.get('use_branch_length', None)
                clade_combo_name = clade_param_dict.get('name', None)
                dist_cutoff = clade_param_dict.get('dist_cutoff', None)
                min_bc_clade_size = clade_param_dict.get('min_bc_clade_size', 1)
                min_mt_clade_size = clade_param_dict.get('min_mt_clade_size', 1)

                
                # assign clades according to each clade generation combo, for each modality-reconstructed tree
                # for res in timepoint_res.scoremat_res_list:

                # create identifier from scoremat params and clade combo
                # name = f"{this_scoremat_res.short_name}_{clade_combo_name}"
                name = f"{clade_combo_name}_{this_scoremat_res.short_name}"

                # the mapping of clade combo names --> names is constant across timepoints
                # so doing this for each timepoint is a bit redundant but ok
                clade_name_to_param_combo_dict[name] = clade_combo_name
                
                # assign clades to cellpop using res+param combo names
                if scheme == 'hierarchical':

                    # can accelerate tree distance matrix calculations using the dendropy workaround to precompute distances
                    # tree argument is not used when distance_metric == 'hamming' so we can provide dp_tree in either case
                    try:
                        clade_dict = this_scoremat_res.assign_hierarchical_clades(distance_metric = distance_metric,
                                                                      num_clusters = num_clades,
                                                                      use_branch_length = use_branch_length,
                                                                      attr_name = name,
                                                                      return_clade_dict = True,
                                                                        tree = this_scoremat_res.dp_recon_tree)
                        growing_clade_dict_collection[timepoint][name] = clade_dict
                    except Exception as e:
                        print(f'failed hierarchical recon scheme: {e}')
                        
                elif scheme == 'topological':
                    try:
                        clade_dict = this_scoremat_res.assign_topological_clades(attr_name = name,
                                                                  return_clade_dict = True,
                                                                  num_clades = num_clades,
                                                                  dist_cutoff = dist_cutoff, 
                                                                  use_branch_length = use_branch_length,
                                                                                tree = this_scoremat_res.recon_tree)
                        growing_clade_dict_collection[timepoint][name] = clade_dict
                    except Exception as e:
                        print(f'failed topological recon scheme: {e}')

                elif scheme == 'random':
                    try:
                        clade_dict = this_scoremat_res.assign_random_clades(attr_name = name,
                                                                             num_clades = num_clades,
                                                                             return_clade_dict = True,
                                                                           tree = this_scoremat_res.recon_tree)
                    except Exception as e:
                        print(f'failed random recon scheme: {e}')

                elif scheme == 'scoremat':
                    # try:

                    # prioritizing num_clades over min_clade_size
                    clade_dict = this_scoremat_res.assign_scoremat_clades(attr_name = name,
                                                                          return_clade_dict = True,
                                                                          num_clades = num_clades,
                                                                         min_clade_size = min_bc_clade_size if this_scoremat_res.mt_or_bc == 'bc' else min_mt_clade_size,
                                                                         return_defining_sig = False)
                    growing_clade_dict_collection[timepoint][name] = clade_dict
                    # except Exception as e:
                    #     print(f'failed scoremat recon scheme: {e}')
    
                    
                # print(f'intermediate growing_clade_dict_collection[timepoint].keys() == {growing_clade_dict_collection[timepoint].keys()}')


        # now for each clade combo, plot a cluster map comparing mt clusters to bc clusters under that combo
        all_clade_combo_names = [clade_param_dict.get('name') for clade_param_dict in param_dict_list]

        # isolate only the dicts for this clade combo and divide into mt and bc lists
        for clade_combo_name in all_clade_combo_names:

            # print(f'clade_combo_name = {clade_combo_name}')
            # print(f'growing_clade_dict_collection[timepoint].keys() = {growing_clade_dict_collection[timepoint].keys()}')
            
            
            mt_name_clade_dicts = [(clade_name, clade_dict)
                              for clade_name, clade_dict in growing_clade_dict_collection[timepoint].items() 
                                if '_mt_' in clade_name and clade_combo_name in clade_name]
            bc_name_clade_dicts = [(clade_name, clade_dict)
                              for clade_name, clade_dict in growing_clade_dict_collection[timepoint].items() 
                                if '_bc' in clade_name and clade_combo_name in clade_name]

            # print(f'len(mt_name_clade_dicts) == {len(mt_name_clade_dicts)}')
            # print(f'mt_name_clade_dicts[0] == {mt_name_clade_dicts[0]}')
            # print(f'len(bc_name_clade_dicts) == {len(bc_name_clade_dicts)}')
            # print(f'bc_name_clade_dicts[0] == {bc_name_clade_dicts[0]}')

            try:
                mt_names, mt_clade_dicts = zip(*mt_name_clade_dicts)
                
            except ValueError as e:
                print(f'mt modality not found; will not plot this pairwise clustermap')
                continue
                
            try:
                bc_names, bc_clade_dicts = zip(*bc_name_clade_dicts)
            except ValueError as e:
                # print(f'bc modality not found; will not plot this pairwise clustermap')
                continue
            
            pairwise_ari = pairwise_concordance(dict_list1 = mt_clade_dicts, 
                                                dict_list2 = bc_clade_dicts,
                                                dict1_names = mt_names,
                                                dict2_names = bc_names)
            
    
            g = sns.clustermap(pairwise_ari,
                              row_cluster = False,
                              col_cluster = False,
                              cmap = 'Blues',
                              vmin = cluster_vmin,
                              vmax = cluster_vmax)

            if include_ari_labels:

                ax = g.ax_heatmap
                ari_vals = pairwise_ari.values

                if cluster_vmin is None and cluster_vmax is None:
                    midpoint = 0.5
                else:
                    midpoint = (cluster_vmin + cluster_vmax) / 2
    
                for (i,j), val in np.ndenumerate(ari_vals):
                    lab_color = 'white' if val > midpoint else 'black'
    
                    ax.text(j + 0.5, i + 0.5,
                            f'{np.round(val, 2)}',
                            ha = 'center',
                            va = 'center',
                            color = lab_color)

            g.cax.set_title('ARI')
            g.ax_heatmap.set_xlabel('')
            g.ax_heatmap.set_ylabel('')
            plt.tight_layout()
            
            
            g.fig.savefig(path_to_cluster_maps / f'{clade_combo_name}_time{timepoint}.pdf')
    
            plt.close(g.fig)

        timepoint_res_dict[timepoint] = timepoint_res

        this_clade_df = make_clade_df(timepoint_res = timepoint_res, 
                                      urid = urid, 
                                      timepoint = timepoint)

        
    # print(f'end of process, clade_name_to_param_combo_dict = {clade_name_to_param_combo_dict}')

    return timepoint_res_dict, clade_name_to_param_combo_dict

# EM utilities


def checkstochastic(M: np.ndarray, name: str) -> None:
    """Small helper - raise if rows do not sum to 1 or contain <0."""
    if not np.allclose(M.sum(axis = 1), 1.0, atol = 1e-9):
        raise ValueError(f"{name} rows must sum to 1.")
    if np.any(M < -1e-12):
        raise ValueError(f"{name} contains negative entries.")

def identify_rare_cells(leaf2type_int, Phi0, Phi1, eps = 1e-12):
        rarity = {}
        for leaf_name, cell_type in leaf2type_int.items():
            prob0 = Phi0[:, cell_type].mean()
            prob1 = Phi1[:, cell_type].mean()
            rarity[leaf_name] = prob1/(prob0+eps)
        return rarity

def clade_prior_from_celltype_enrichment(clade_cells, leaf2type, celltype_induction_score_map):
        '''
        for one clade, generated weighted induction prior based on cell type --> induction score map
        '''
        ct_counts = Counter([leaf2type[cell] for cell in clade_cells])
        total = sum(ct_counts.values())
        if total == 0:
            return 0.5
    
        running_weight = 0
        for ct, ct_count in ct_counts.items():
            running_weight += (ct_count * celltype_induction_score_map[ct])
        return running_weight / total

def beta_parameters(prior_mean, strength = 10):
    alpha = prior_mean * strength
    beta = (1.0 - prior_mean) * strength
    return alpha, beta

def initialize_p_frac_from_rarity(leaf2type_int, leaf2clade, Phi0, Phi1, clades, eps = 1e-12):
    """
    Initialize p_frac based on the rarity of cell types in each clade.
    
    High rarity -> high p_c (clade likely induced)
    Low rarity -> low p_c (clade likely un-induced)
    """
    
    # Compute rarity score for each cell
    rarity_scores = {}
    for leaf, cell_type in leaf2type_int.items():
        # Average probability of this cell type under each matrix (over all parent states)
        prob0 = np.mean(Phi0[:, cell_type])
        prob1 = np.mean(Phi1[:, cell_type])
        rarity_scores[leaf] = prob1 / (prob0 + eps)
    
    # Compute per-clade rarity
    p_frac_init = {}
    for c in clades:
        leaves_in_c = [leaf for leaf, clade in leaf2clade.items() if clade == c]
        
        if not leaves_in_c:
            p_frac_init[c] = 0.5  # default
            continue
        
        # Average rarity in this clade
        avg_rarity = np.mean([rarity_scores[leaf] for leaf in leaves_in_c])
        
        # Map rarity to p_c
        # Rarity of 1.0 -> p_c = 0.5 (no signal)
        # Rarity of 5.0 -> p_c = 0.9 (strong signal for induced)
        # Rarity of 0.2 -> p_c = 0.1 (strong signal for un-induced)
        
        # Using a logistic function: p_c = 1 / (1 + (1/avg_rarity)^scale)
        # Or simpler: p_c = avg_rarity / (1 + avg_rarity)
        
        p_frac_init[c] = avg_rarity / (1.0 + avg_rarity)  # maps [0, inf) -> [0, 1)
        
        # Clip to avoid extremes
        p_frac_init[c] = np.clip(p_frac_init[c], 0.1, 0.9)
    
    return p_frac_init

# EM classes


class cp_TreeMixtureNodeEM2():
    """
    EM for a node-only mixture of two K-state stochastic matrices on a rooted
    phylogeny.  The *likelihood-ratio* weighting scheme is implemented as one
    of the possible ``weight_strategy`` options.
    """
    #  Constructor
    def __init__(self,
             Phi0: np.ndarray,
             Phi1: np.ndarray,
             root_type,
             tree,
             leaf2clade: Dict[str, str],
             leaf2type: Dict[str, Any],
             type_encoder: Optional[Dict[Any, int]] = None,
             weight_strategy: str = 'leaf_prop'):
        """
        Parameters
        ----------
        weight_strategy : str
            One of
                'leaf_prop'   - leaf-proportion weighting (original code)
                'uniform'     - uniform-per-node
                'inv_leaf'    - inverse-leaf weighting
                
                
                'lr+uniform'  - LR multiplied by uniform (default for LR)
                'lr+inv'      - LR multiplied by inverse-leaf
        """
        # 1  Store tree and attach leaf metadata
        self.tree = tree
        self.leaf2clade = leaf2clade
        self.leaf2type = leaf2type
        self.weight_strategy = weight_strategy

        leaf_names = {leaf.name for leaf in self.tree.iter_leaves()}
        missing = leaf_names - set(self.leaf2clade.keys())
        if missing:
            raise ValueError(f"The following leaves are missing from leaf2clade: {missing}")
        # 2  Encode cell types as integers 0 ... K-1
        if type_encoder is None:
            uniq = sorted(set(self.leaf2type.values()))
            self.type_encoder = {t: i for i, t in enumerate(uniq)}
        else:
            self.type_encoder = type_encoder
        self.leaf2type_int = {
            leaf: self.type_encoder[self.leaf2type[leaf]]
            for leaf in leaf_names
        }
        # 3  Store the two node-transition matrices (must be stochastic)
        self.Phi0 = np.asarray(Phi0, dtype = float)
        self.Phi1 = np.asarray(Phi1, dtype = float)
        self.K = self.Phi0.shape[0]
        assert self.Phi0.shape == (self.K, self.K)
        assert self.Phi1.shape == (self.K, self.K)
        checkstochastic(self.Phi0, 'Phi0')
        checkstochastic(self.Phi1, 'Phi1')
        # 4  Root prior (founder type is known)
        if isinstance(root_type, str):
            root_type = self.type_encoder[root_type]
        self.root_type = int(root_type)
        # 5  Initialise clade-specific mixture parameters
        self.clades = sorted(set(self.leaf2clade.values()))
        self.p_frac = {c: 0.5 for c in self.clades}

        # cache for fast node-to-clade look-up
        self._assign_node_clades()
    #  Helper methods
    def _check_tree_is_rooted(self) -> None:
        """Raise if the supplied tree has no root."""
        if not self.tree.get_tree_root():
            raise ValueError('The supplied tree must be rooted.')

    def _assign_node_clades(self) -> None:
        """
        Attach a ``clade_set`` attribute to every node:
        the set of clades that have at least one descendant leaf under that node.
        """
        for node in self.tree.traverse('postorder'):
            if node.is_leaf():
                node.clade_set = {self.leaf2clade[node.name]}
            else:
                node.clade_set = set()
                for child in node.children:
                    node.clade_set.update(child.clade_set)
    #  Node-weight factories
    def _node_weights_uniform(self, node) -> Dict[str, float]:
        """Uniform-per-node: weight = 1 for every clade covered by the node."""
        return {c: 1.0 for c in node.clade_set}

    def _node_weights_inverse_leaf(self, node) -> Dict[str, float]:
        """Inverse-leaf: weight = 1 / (# descendant leaves)."""
        n_leaves = len(list(node.iter_leaves()))
        if n_leaves == 0:
            return {}
        inv = 1.0 / n_leaves
        return {c: inv for c in node.clade_set}

    def _node_weights_leaf_prop(self, node):
        n_leaves = len(list(node.iter_leaves()))
        if n_leaves == 0:
            return {}
        inv = 1.0 / n_leaves
        return {c: n_leaves for c in node.clade_set}

    def _node_weights_lr(self, node) -> Dict[str, float]:
        """
        LR weighting: weight = LR_u for every clade covered by the node.
        ``node.lr_ratio`` is computed during the pruning pass (see below).
        """
        return {c: node.lr_ratio for c in node.clade_set}

    def _node_weights_lr_plus_uniform(self, node) -> Dict[str, float]:
        """LR  uniform (equivalent to pure LR because uniform = 1)."""
        return {c: node.lr_ratio for c in node.clade_set}

    def _node_weights_lr_plus_inv(self, node) -> Dict[str, float]:
        """LR  inverse-leaf."""
        n_leaves = len(list(node.iter_leaves()))
        if n_leaves == 0:
            return {}
        inv = 1.0 / n_leaves
        return {c: node.lr_ratio * inv for c in node.clade_set}

    def _node_weights_lr_plus_leaf_prop(self, node):
        n_leaves = len(list(node.iter_leaves()))
        if n_leaves == 0:
            return {}
        inv = 1.0 / n_leaves
        return {c: node.lr_ratio * n_leaves for c in node.clade_set}
        

    def pruninglikelihood(self) -> float:

        K = self.K
        eps = 1e-12
    
        total_log_factor = 0.0               # cumulative scaling term
    
        for node in self.tree.traverse('postorder'):
            if node.is_leaf():
                # leaf emission (one-hot)
                lk = np.zeros(K)
                lk[self.leaf2type_int[node.name]] = 1.0
                node.lik = lk
    
                # log vectors: log(1)=0 for the observed state, -inf otherwise
                node.logL0 = np.full(K, -np.inf)
                node.logL1 = np.full(K, -np.inf)
                obs = self.leaf2type_int[node.name]
                node.logL0[obs] = 0.0
                node.logL1[obs] = 0.0
    
                node.lr_ratio = 1.0               # dummy, never used for leaves
                node.log_factor = 0.0
                continue

            # 1 lift each child through Phi0 and Phi1 and sum the logs
            logL0 = np.zeros(K)
            logL1 = np.zeros(K)
            for child in node.children:

                prob0 = self.Phi0 @ child.lik
                prob1 = self.Phi1 @ child.lik
                # avoid log(0); eps is tiny enough not to affect the result
                prob0[prob0 == 0] = eps
                prob1[prob1 == 0] = eps
                logL0 += np.log(prob0)
                logL1 += np.log(prob1)
    
            # 2 remove a common offset (scaling)
            max0 = logL0.max()
            node.logL0 = logL0 - max0
            node.logL1 = logL1 - max0          # same offset for both components
            node.log_factor = max0
    
            # 3 LR scalar (for lr weighting)
            node.lr_ratio = np.exp(node.logL1.max()) / (np.exp(node.logL0.max()) + eps)
    
            # 4 clade-specific mixture probability p_u
            clade_weights = defaultdict(float)
            for leaf in node.iter_leaves():
                clade_weights[self.leaf2clade[leaf.name]] += 1.0
            total = sum(clade_weights.values())
            if total == 0:
                p_u = 0.0
            else:
                p_u = sum((w / total) * self.p_frac[c]
                          for c, w in clade_weights.items())
    
            # 5 mixture upward vector (still a probability vector)
            L0 = np.exp(node.logL0)            # already divided by exp(max0)
            L1 = np.exp(node.logL1)
            mixed = (1.0 - p_u) * L0 + p_u * L1
    
            mix_max = mixed.max()
            if mix_max == 0:
                mixed = np.full(K, eps)
                mix_max = eps
            node.lik = mixed / mix_max          # keep it normalised
            node.log_factor += np.log(mix_max)  # add the mixture scaling
    
            total_log_factor += node.log_factor
    
        # 6 root contribution
        root = self.tree.get_tree_root()
        root_prior = np.zeros(K)
        root_prior[self.root_type] = 1.0
    
        logL = np.log(np.dot(root_prior, root.lik) + eps) + total_log_factor
    
        # 7 diagnostics (only once)
        if not hasattr(self, '_diag_printed'):
            internal = sum(1 for n in self.tree.traverse() if not n.is_leaf())
            mixed_nodes = sum(1 for n in self.tree.traverse()
                              if not n.is_leaf() and getattr(n, 'has_mixture', False))
            self._diag_printed = True
    
        return logL

    
    def node_responsibility(self, node, parent_prior_vec = None):
        """
        Public version of the responsibility calculation.
        Returns the posterior probability _u that the division represented by ``node``
        used the induced transition matrix.
    
        Parameters
        ----------
        node : ete3.TreeNode
            An internal node (must not be a leaf or the root).
        parent_prior_vec : np.ndarray, optional
            The likelihood vector of the parent (parent.up.lik).
            If omitted, falls back to the root prior.
    
        Returns
        -------
        float
            Posterior responsibility _u  [0, 1].
        """
        eps = 1e-12
        # 1  Get the parent prior vector
        if parent_prior_vec is None:
            if node.up is None:
                parent_prior_vec = np.zeros(self.K)
                parent_prior_vec[self.root_type] = 1.0
            else:
                parent_prior_vec = node.up.lik
        # 2  Compute the clade-weighted prior mixture
        clade_weights = defaultdict(float)
        for leaf in node.iter_leaves():
            clade_weights[self.leaf2clade[leaf.name]] += 1.0
        total = sum(clade_weights.values())
        if total == 0:
            prior_mix = 0.0
        else:
            prior_mix = sum((w / total) * self.p_frac[c]
                            for c, w in clade_weights.items()) 
        # 3  Log-likelihood pieces that were stored by pruninglikelihood
        log_factor = node.log_factor
        
        #  CRITICAL FIX: weight by parent_prior_vec BEFORE taking logsumexp
        # We want log(parent_prior  L_vec), not just log(sum(L_vec))
        parent_prior_log = np.log(parent_prior_vec + eps)  # shape (K,)
        
        # log(parent_prior  exp(logL0 + log_factor))
        #  = log(parent_prior) + logL0 + log_factor
        log_term0 = logsumexp(parent_prior_log + node.logL0 + log_factor)
        
        # log(parent_prior  exp(logL1 + log_factor))
        log_term1 = logsumexp(parent_prior_log + node.logL1 + log_factor)
        # 4  Return the posterior _u using Bayes rule
        if prior_mix <= 0.0:
            return 0.0
        if prior_mix >= 1.0:
            return 1.0
        
        # _u = P(Z_u=1 | data)
        #     = prior_mix  P(data|Z_u=1) / [prior_mix  P(data|Z_u=1) + (1-prior_mix)  P(data|Z_u=0)]
        #     = prior_mix  exp(log_term1) / [prior_mix  exp(log_term1) + (1-prior_mix)  exp(log_term0)]
        
        # In log space:
        log_num = np.log(prior_mix + eps) + log_term1
        
        # logsumexp([log(prior_mix) + log_term1, log(1-prior_mix) + log_term0])
        log_den = logsumexp(np.array([
            np.log(prior_mix + eps) + log_term1,
            np.log(1.0 - prior_mix + eps) + log_term0
        ]))
        
        gamma = np.exp(log_num - log_den)
        return float(np.clip(gamma, 0.0, 1.0))

    def leaf_predictions(self,
                     threshold: float = 0.5,
                     use_path_responsibilities: bool = False) -> pd.DataFrame:
        """
        Convert the clade-level mixture weights (self.p_frac) into leaf-level
        predictions.
    
        Parameters
        ----------
        threshold : float, default 0.5
            Decision threshold for the hard label.  Leaves with
            ``prob_induced >= threshold`` are labelled induced.
        use_path_responsibilities : bool, default False
            If True, incorporate the posterior responsibilities (_u) of every
            internal node on the path from the leaf to the root.  The leaf
            probability is then
                P(induced | leaf) = _{u on path} [ p_parent_u + (1-p_parent)(1-_u) ].
            The product is evaluated in log-space for numerical stability.
            If False (default) the leaf probability is simply the clade mixture
            weight p_c.
    
        Returns
        -------
        pandas.DataFrame with columns:
            - leaf_name      : identifier of the leaf / cell
            - clade          : clade label taken from self.leaf2clade
            - p_clade        : posterior mixture weight for that clade (p_c)
            - prob_induced   : soft probability that the leaf is induced
            - pred_label     : hard label (induced or un-induced)
        """
        records = []
    
        for leaf in self.tree.iter_leaves():
            leaf_name = leaf.name
            clade = self.leaf2clade[leaf_name]
            p_c = float(self.p_frac[clade])      # clade-level posterior weight
            # A) Simple clade-only probability (fast)
            prob = p_c
            # B) Path-aware probability (if requested)
            if use_path_responsibilities:
                # Walk up from leaf to root, accumulating log-probability.
                node = leaf
                log_prob = np.log(p_c + 1e-12)      # start with the clade weight
                while not node.is_root():
                    parent = node.up
                    # _u for the division that produced `node`
                    gamma_u = self.node_responsibility(node)
    
                    # mixture weight for the parent clade (if the parent has a clade)
                    parent_clade = self.leaf2clade.get(parent.name, None)
                    if parent_clade is None:
                        p_parent = 0.5               # neutral prior when unknown
                    else:
                        p_parent = float(self.p_frac[parent_clade])
    
                    # Edge-level probability that this division was induced:
                    edge_prob = p_parent * gamma_u + (1.0 - p_parent) * (1.0 - gamma_u)
    
                    log_prob += np.log(edge_prob + 1e-12)
                    node = parent
                prob = np.exp(log_prob)          # back to probability space
            # Hard label
            label = 'induced' if prob >= threshold else 'un-induced'
    
            records.append({
                'leaf_name': leaf_name,
                'clade': clade,
                'p_clade': p_c,
                'prob_induced': prob,
                'pred_label': label
            })
    
        return pd.DataFrame.from_records(records)
    #  EM algorithm
    def fit(self,
            max_iter: int = 100,
            tol: float = 1e-2,
            verbose: bool = True,
            prior_alpha: float = 1.0,
            prior_beta: float = 1.0,
            prior_alpha_dict = None,
            prior_beta_dict = None,
            # reset_to_baseline = True
           ) -> list:
        """
        Run the EM algorithm.
        The ``weight_strategy`` chosen at construction determines how the
        per-node contributions are weighted in the M-step.
        """
        logL_history = []
        #  Choose the weight-lookup function once (so we do not test strings
        #  inside the inner loop).
        if self.weight_strategy == 'leaf_prop':
            weight_func = self._node_weights_leaf_prop  
        elif self.weight_strategy == 'uniform':
            weight_func = self._node_weights_uniform
        elif self.weight_strategy == 'inv_leaf':
            weight_func = self._node_weights_inverse_leaf
        elif self.weight_strategy == 'lr':
            weight_func = self._node_weights_lr
        elif self.weight_strategy == 'lr+uniform':
            weight_func = self._node_weights_lr_plus_uniform
        elif self.weight_strategy == 'lr+inv':
            weight_func = self._node_weights_lr_plus_inv
        elif self.weight_strategy == 'lr+leaf_prop':
            weight_func = self._node_weights_lr_plus_leaf_prop
        else:
            raise ValueError(f"Unknown weight_strategy: {self.weight_strategy}")

        # if reset_to_baseline:
        #     baseline_p_frac = {c: 0.5 for c in self.clades}
        #     self.p_frac = baseline_p_frac
        #     baseline_logL = self.pruninglikelihood()
        #     print(f"Baseline logL (all p_c=0.5): {baseline_logL:.8f}")
        # else:
        #     baseline_logL = None
    
        for it in tqdm.trange(max_iter, disable = not verbose):
            # E-step
            # **IMPORTANT** - call the *correct* pruning routine
            logL = self.pruninglikelihood()          # <-- fixed line
            logL_history.append(logL)
    
            # (optional) sanity check that vectors exist
            # for n in self.tree.traverse():
            #     if not n.is_leaf() and not n.is_root():
            #         assert n.L0_vec is not None and n.L1_vec is not None
    
            if verbose:
                logging.info(f"EM iteration {it+1:02d} - logL = {logL:.6f}")
    
            # Gather sufficient statistics
            clade_stats = {c: [0.0, 0.0] for c in self.clades}   # [numer, denom]
    
            # Root prior vector (used as parent prior for the root's children)
            root_prior = np.zeros(self.K)
            root_prior[self.root_type] = 1.0


            for node in self.tree.traverse('postorder'):
                if node.is_root() or node.is_leaf(): # added is_leaf() check
                    continue                       # root has no Z_u

                clade_weights = defaultdict(float)
                for leaf in node.iter_leaves():
                    clade_weights[self.leaf2clade[leaf.name]] += 1.0
                # if len(clade_weights) > 1:
                #     print(f"  Node {node.name}: split across {len(clade_weights)} clades: {dict(clade_weights)}")
                # elif len(clade_weights) == 1:
                #     print(f"  Node {node.name}: has clade_lengths == 1")
    
                # parent prior comes from the already-computed node.up.lik
                parent_vec = node.up.lik if node.up else root_prior
    
                # posterior responsibility
                gamma_u = self.node_responsibility(node, parent_vec)   # uses L0_vec/L1_vec
    
                # node weight according to strategy
                w_dict = weight_func(node)        # e.g. {'cladeA': 5.0, ...}
                # print(f'{w_dict = }')
    
                # accumulate numerator / denominator for each clade
                for c, w in w_dict.items():
                    clade_stats[c][0] += w * gamma_u    #  w
                    clade_stats[c][1] += w             #  w
    
            # M-step (MAP update)
            for c in self.clades:
                num, den = clade_stats[c]
                if den > 0:
                    # new_p = (num + prior_alpha - 1.0) / (den + prior_alpha + prior_beta - 2.0)
                    # self.p_frac[c] = np.clip(new_p, 1e-5, 1.0 - 1e-5)

                    raw_estimate = num/den

                    alpha = prior_alpha_dict.get(c, prior_alpha) if prior_alpha_dict else prior_alpha
                    beta = prior_beta_dict.get(c, prior_beta) if prior_beta_dict else prior_beta

                    new_p = (num + alpha - 1.0) / (den + alpha + beta - 2.0)

                    
                    self.p_frac[c] = np.clip(new_p, 1e-5, 1.0 - 1e-5) if np.isfinite(new_p) else 0.5

            # self.update_mixture()
    
            # convergence check
            logL_change = abs(logL_history[-1] - logL_history[-2]) if it > 0 else np.inf

            if verbose:
                if it > 0 and logL_change < tol:
                    logging.info(f"EM converged (log-likelihood change {logL_change:.2e} < {tol})")
                    break
            if it > 0 and abs(logL_history[-1] - logL_history[-2]) < tol:
                if verbose:
                    logging.info('EM converged (log-likelihood change < tol).')
                break


        # final_logL = logL_history[-1]
        # improvement = final_logL - baseline_logL
        # print(f"\nFinal logL:     {final_logL:.8f}")
        # print(f"Improvement:    {improvement:.8f}")
        # if improvement > 0:
        #     print(" Model improved over baseline")
        # else:
        #     print("  Model did NOT improve over baseline")
        self.logL_history_ = logL_history
        return logL_history
    #  Post-processing utilities
    def clade_posteriors(self, use_gmm = False) -> pd.DataFrame:
        """Return a DataFrame with one row per clade (p_induced & hard label)."""

        if use_gmm and len(self.clades) < 2:
            use_gmm = False  # GMM needs >= 2 samples; fall back to threshold

        if use_gmm:
            gmm_eps = 1e-6
            x = logit(np.clip(np.array([self.p_frac[c] for c in self.clades]), gmm_eps, 1-gmm_eps)).reshape(-1,1)
            
            gmm = GaussianMixture(n_components = 2, random_state = 42)
            gmm.fit(x)


            probs = gmm.predict_proba(x)
            labels = gmm.predict(x)

            # determine which component is induced and which is uninduced. induced assumed to have greater mean
            means = gmm.means_.flatten()
            induced_component = np.argmax(means)
        
            str_labels = ['induced' if label == induced_component else 'un-induced' for label in labels]
            
            # str_labels = ['induced' if label == 1 else 'uninduced' for label in labels]

            rows = [{'clade': self.clades[i], 
                     'p_induced': probs[i][induced_component],
                     'hard_label': str_labels[i]} for i in range(len(self.clades))]
            
        
        else:
            rows = []
            for c in self.clades:
            
                p = self.p_frac[c]
                rows.append({
                    'clade': c,
                    'p_induced': p,
                    'hard_label': 'induced' if p > 0.5 else 'un-induced'
                })
        return pd.DataFrame(rows) 

    
    def summary(self) -> None:
        """Print a quick textual summary."""
        df = self.clade_posteriors()
        print(df.to_string(index = False))
        if hasattr(self, 'logL_history_'):
            pass

    def diagnose_gamma_by_clade(self):
        """
        Inspect per-node _u values grouped by clade.
        Should be called AFTER fit() to see which clades prefer which matrix.
        """
        
        # First, recompute the pruning pass to populate L0_vec, L1_vec
        self.pruninglikelihood()
        
        root_prior = np.zeros(self.K)
        # root_prior[self.type_encoder[self.root_type]] = 1.0
        root_prior[self.root_type] = 1.0
        
        # Collect  values by clade
        gamma_by_clade = defaultdict(list)
        
        for node in self.tree.traverse('postorder'):
            if node.is_leaf() or node.is_root():
                continue
            
            parent_vec = node.up.lik if node.up else root_prior
            gamma = self.node_responsibility(node, parent_vec)
            
            # Determine which clade(s) this node belongs to
            clade_weights = defaultdict(float)
            for leaf in node.iter_leaves():
                clade_weights[self.leaf2clade[leaf.name]] += 1.0
            
            if clade_weights:
                # For nodes that span multiple clades, use the dominant one
                clade = max(clade_weights, key = clade_weights.get)
                gamma_by_clade[clade].append((node.name, gamma))
        
        # Print summary per clade
        for clade in sorted(self.clades):
            gammas = gamma_by_clade[clade]
            if gammas:
                gamma_array = np.array([g for _, g in gammas])
                mean_gamma = np.mean(gamma_array)
                std_gamma = np.std(gamma_array)
                min_gamma = np.min(gamma_array)
                max_gamma = np.max(gamma_array)
                
                
                # Show individual nodes if they vary significantly
                if std_gamma > 0.1:
                    for node_name, gamma in sorted(gammas, key = lambda x: x[1]):
                        pass
        
        # Overall statistics
        all_gammas = [g for gammas in gamma_by_clade.values() for _, g in gammas]
        
        return gamma_by_clade


# EM runners


def run_em_with_restarts(Phi0, Phi1, root_type, tree, leaf2clade, leaf2type,
                         type_encoder, weight_strategy, clade_alpha, clade_beta,
                         n_restarts = 5, max_iter = 200, tol = 1e-3, seed_start = 42,
                         verbose = False, p_frac_init = None,
                         alpha_values = None, Phi1_true = None):
    """
    Run EM multiple times with different random initializations.

    If alpha_values and Phi1_true are provided, sweeps over each alpha,
    blending Phi1_input(alpha) = (1-alpha)*Phi0 + alpha*Phi1_true, and
    returns a dict {alpha: best_result}.  Otherwise runs once with the
    supplied Phi1 and returns a single best_result dict.

    Parameters
    ----------
    alpha_values : list of float, optional
        Alpha values to sweep (0 = uninformative prior, 1 = exact prior).
    Phi1_true : np.ndarray, optional
        True induced transition matrix; required when alpha_values is given.
    """
    if alpha_values is not None:
        if Phi1_true is None:
            raise ValueError('Phi1_true must be supplied when alpha_values is given')
        alpha_results = {}
        for alpha in alpha_values:
            Phi1_input = (1.0 - alpha) * Phi0 + alpha * np.asarray(Phi1_true, dtype = float)
            _, best_result, _ = run_em_with_restarts(
                Phi0 = Phi0, Phi1 = Phi1_input,
                root_type = root_type, tree = tree,
                leaf2clade = leaf2clade, leaf2type = leaf2type,
                type_encoder = type_encoder, weight_strategy = weight_strategy,
                clade_alpha = clade_alpha, clade_beta = clade_beta,
                n_restarts = n_restarts, max_iter = max_iter, tol = tol,
                seed_start = seed_start, verbose = verbose,
                p_frac_init = p_frac_init,
            )
            alpha_results[alpha] = best_result
        return alpha_results

    results = []

    for restart_idx in range(n_restarts):
        seed = seed_start + restart_idx
        np.random.seed(seed)


        em = cp_TreeMixtureNodeEM2(
            Phi0 = Phi0,
            Phi1 = Phi1,
            root_type = root_type,
            tree = tree,
            leaf2clade = leaf2clade,
            leaf2type = leaf2type,
            type_encoder = type_encoder,
            weight_strategy = weight_strategy
        )


        # Compute base restart_init from priors (once per call, not once per restart)
        if restart_idx == 0:
            prior_mean_frac = {c: clade_alpha[c] / (clade_alpha[c] + clade_beta[c]) for c in em.clades}
            prior_mean_frac = {c: np.clip(p, 0.1, 0.9) for c, p in prior_mean_frac.items()}

        if p_frac_init is not None:
            # Caller supplied explicit init - use as-is for all restarts
            restart_init = p_frac_init
        elif restart_idx == 0:
            # First restart: use prior mean (deterministic anchor)
            restart_init = prior_mean_frac
        else:
            # Subsequent restarts: random perturbation around prior mean
            restart_init = {c: float(np.clip(
                                prior_mean_frac[c] + np.random.uniform(-0.4, 0.4), 0.05, 0.95))
                            for c in em.clades}

        for c, init_val in restart_init.items():
            em.p_frac[c] = init_val

         
        iteration_lr_val_dict = {}
        
        logL_hist = em.fit(
            max_iter = max_iter,
            tol = tol,
            verbose = verbose,  # suppress iteration details to reduce output
            prior_alpha_dict = clade_alpha,
            prior_beta_dict = clade_beta,
        )


        
        final_logL = logL_hist[-1]

        n_iterations = len(logL_hist)

        
        # Store results
        results.append({
            'restart': restart_idx,
            'seed': seed,
            'final_logL': final_logL,
            'n_iterations': n_iterations,
            'p_frac': em.p_frac.copy(),
            'clade_posteriors': em.clade_posteriors().copy(),
            'model': em,
            'clade_priors': p_frac_init,
            'iteration_lr_val_dict': iteration_lr_val_dict
        })
    
    # Find best solution
    best_idx = np.argmax([r['final_logL'] for r in results])
    best_result = results[best_idx]

    
    return results, best_result, p_frac_init

def evaluate_prediction_probs(em, leaf2clade, leaf2induced_gt, res_df, plot = True):
    cell_names = list(leaf2clade.keys())
    ground_truth_str = np.array([leaf2induced_gt[cell] for cell in cell_names])
    ground_truth = np.array([1 if label == 'induced' else 0 for label in ground_truth_str])


    

    clade_to_prob_dict = dict(zip(res_df['clade'], res_df['p_induced']))
    clade_to_hard_label = dict(zip(res_df['clade'], res_df['hard_label']))

    predicted_probs = np.array([clade_to_prob_dict[leaf2clade[cell]] for cell in cell_names])
    predicted_label = np.array([1 if clade_to_hard_label[leaf2clade[cell]] == 'induced' else 0 for cell in cell_names])

    precision, recall, pr_thresholds = precision_recall_curve(ground_truth, predicted_probs, pos_label = 1)
    pr_auc = auc(recall, precision)

    fpr, tpr, roc_thresholds = roc_curve(ground_truth, predicted_probs, pos_label = 1)
    if len(np.unique(ground_truth)) < 2 or np.any(np.isnan(predicted_probs)):
        roc_auc = float('nan')
    else:
        roc_auc = roc_auc_score(ground_truth, predicted_probs)

    
    # predicted_probs = np.array([em.p_frac[leaf2clade[cell]] for cell in cell_names])
    # predicted_label = np.array([1 if prob > 0.5 else 0 for prob in predicted_probs])

    
    # print(f'{predicted_label = }')
    # print(f'{ground_truth = }')
    accuracy = np.mean(predicted_label == ground_truth)

    # tp = np.sum((predicted_label == 'induced') & (ground_truth == 'induced'))
    # fp = np.sum((predicted_label == 'induced') & (ground_truth == 'un-induced'))
    # tn = np.sum((predicted_label == 'un-induced') & (ground_truth == 'un-induced'))
    # fn = np.sum((predicted_label == 'un-induced') & (ground_truth == 'induced'))

    tp = np.sum((predicted_label == 1) & (ground_truth == 1))
    fp = np.sum((predicted_label == 1) & (ground_truth == 0))
    tn = np.sum((predicted_label == 0) & (ground_truth == 0))
    fn = np.sum((predicted_label == 0) & (ground_truth == 1))

    # print(f'{tp = }')
    # print(f'{fp = }')
    # print(f'{tn = }')
    # print(f'{fn = }')
    
    precision_hard = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall_hard = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1_hard = 2 * (precision_hard * recall_hard) / (precision_hard + recall_hard) if (precision_hard + recall_hard) > 0 else 0
    
    metrics = {
        'pr_auc': pr_auc,
        'roc_auc': roc_auc,
        'accuracy': accuracy,
        'precision_hard': precision_hard,
        'recall_hard': recall_hard,
        'f1_hard': f1_hard,
        'tp': tp,
        'fp': fp,
        'tn': tn,
        'fn': fn,
    }
    
    
    if plot:
        fig, axes = plt.subplots(1, 2, figsize = (14, 5))
        
        # PR curve
        axes[0].plot(recall, precision, label = f'PR curve (AUC={pr_auc:.4f})')
        axes[0].fill_between(recall, precision, alpha = 0.2)
        axes[0].set_xlabel('Recall')
        axes[0].set_ylabel('Precision')
        axes[0].set_title('Precision-Recall Curve')
        axes[0].legend()
        axes[0].grid(alpha = 0.3)
        axes[0].set_xlim([0, 1])
        axes[0].set_ylim([0, 1])
        
        # ROC curve
        axes[1].plot(fpr, tpr, label = f'ROC curve (AUC={roc_auc:.4f})')
        axes[1].fill_between(fpr, tpr, alpha = 0.2)
        axes[1].plot([0, 1], [0, 1], 'k--', label = 'Random classifier')
        axes[1].set_xlabel('False Positive Rate')
        axes[1].set_ylabel('True Positive Rate')
        axes[1].set_title('ROC Curve')
        axes[1].legend()
        axes[1].grid(alpha = 0.3)
        axes[1].set_xlim([0, 1])
        axes[1].set_ylim([0, 1])
        
        plt.tight_layout()
        plt.savefig('evaluation_curves.png', dpi = 150, bbox_inches = 'tight')
        plt.close('all')

        # Add this to your code:
    
    # Histogram
    n_bins = min(50, len(np.unique(predicted_probs)))
    if n_bins > 1:
        plt.hist(predicted_probs, bins = n_bins)
    else:
        plt.text(0.5, 0.5, f"All probs = {predicted_probs[0]:.4f}", ha = 'center', transform = plt.gca().transAxes)
    plt.xlabel('Predicted Probability')
    plt.ylabel('Count')
    plt.close('all')

    return metrics, (precision, recall), (fpr, tpr), (predicted_probs, predicted_label, ground_truth)

def urid_to_induction_inference(urid_list, timepoints, param_dict_list, weight_strategy = 'leaf_prop',
                               clade_priors = None, upper_rarity_bound = 2,
                               include_param_combos = None, scoremat_filters = None,
                               verbose = False,
                               root_type = 'ct1', num_restarts = 1,
                               return_opt_model_dict = False, type_encoder = None,
                               fit_tol = 1e-6, enrichment_dict = None, strength = None,
                               use_gmm = False, p_frac_init = None,
                               path_to_output_dir = path_to_output_dir):

    '''
    if clade_priors == 'rarity':
        clade-specific priors determined from enrichment of rarest cell type
    elif clade_priors == 

    '''

    urid_timepoint_res_dict = {}
    
    for urid in tqdm.tqdm(urid_list):
        timepoint_res_dict, clade_name_to_param_combo_dict = process_urid(urid = urid, 
                                                                         param_dict_list = param_dict_list,
                                                                        cluster_vmin = -0.5,
                                                                        cluster_vmax = 1,
                                                                         timepoints = timepoints)
        urid_timepoint_res_dict[urid] = timepoint_res_dict
    
    for urid, timepoint_res_dict in urid_timepoint_res_dict.items():
        for tp, tpr in timepoint_res_dict.items():
            # if timepoints is not None:
            #     if tp not in timepoints:
            #         continue
            df = make_clade_df(timepoint_res = tpr,
                               urid = urid,
                               timepoint = tp)

    param_combo_induction_dict = {}
    full_em_model_dict = {}

    # new for loop ... 
    for urid in urid_timepoint_res_dict.keys():

        os.makedirs(path_to_em_induction_prob / str(urid), exist_ok = True)

        ex_uninduced_tm, ex_induced_tm = get_tms(urid, path_to_output_dir = path_to_output_dir)

        for timepoint in timepoints:
    
            full_df = make_clade_celltype_df(urid = urid, timepoint = timepoint, urid_timepoint_res_dict = urid_timepoint_res_dict, path_to_output_dir = path_to_output_dir)


            # full_subset_gt = full_df[full_df['param_combo'] == 'clade_differentiation_induced']
            # gt_linstring_to_induced_dict = dict(zip(full_subset_gt['cell'], ['induced' if clade else 'un-induced' for clade in full_subset_gt['clade']]))
            
            
            tpt_param_combo_induction_dict = {}
            tpt_em_model_dict = {}
            
            
            for param_combo in full_df['param_combo'].unique():
    
                
                full_df_copy = copy.deepcopy(full_df)
    
                
                # ex_combo_name = '5clade_hierarch_hamming_mt_250ints_1.0RP_0.2AF_TB'
                ex_short_name = combo_name_to_short_name(param_combo)

    
                # if this isn't a recon-tree-defined clade, assign the tree to be the gt tree
                if ex_short_name is None:

                    # additional pruning steps needed if working with gt tree
                    keep_linstrings = prune_dead_lineages(cellpop_dict = urid_timepoint_res_dict[urid][timepoint].cell_pop)
                    # Re-read from the newick file on disk: ete3's file reader is
                    # iterative, while deepcopy and tree.write() both recurse to
                    # tree depth and overflow the stack for late-timepoint trees.
                    tree = read_in_tree(urid_timepoint_res_dict[urid][timepoint].gttree_path)
                    keep_nodes = [node for node in tree.traverse('preorder') if node.name in keep_linstrings]
                    tree.prune(keep_nodes, preserve_branch_length = True)
                    filt_df = full_df[(full_df_copy['cell'].isin(keep_linstrings)) & (full_df_copy['param_combo'] == param_combo)]
                    
                
                else:


                    tree = [res.recon_tree for res in urid_timepoint_res_dict[urid][timepoint].scoremat_res_list 
                                                 if res.short_name == ex_short_name][0]
                
                
                    filt_df = full_df_copy[full_df_copy['param_combo'] == param_combo]

                if type_encoder is None:
                    type_encoder = dict(zip(sorted(filt_df['celltype'].unique()), range(filt_df['celltype'].nunique())))
                
                leaf2clade = dict(zip(filt_df['cell'], filt_df['clade']))
                leaf2type = dict(zip(filt_df['cell'], filt_df['celltype']))

                if include_param_combos is not None:
                    if param_combo not in include_param_combos:
                        continue
                if scoremat_filters is not None:
                    # gt combos (ending _gt) always pass; scoremat combos must
                    # contain at least one of the filter strings (e.g. '1.0RP_0.0AF')
                    if not param_combo.endswith('_gt'):
                        if not any(f in param_combo for f in scoremat_filters):
                            continue

                em = cp_TreeMixtureNodeEM2(
                    Phi0 = ex_uninduced_tm,
                    Phi1 = ex_induced_tm,
                    root_type = root_type,
                    tree = tree,
                    leaf2clade = leaf2clade,
                    leaf2type = leaf2type,
                    type_encoder = type_encoder,
                    weight_strategy = weight_strategy
                )


                if clade_priors is not None:
                    if clade_priors == 'rarity':
    
                        rarity_scores = identify_rare_cells(leaf2type_int = em.leaf2type_int, 
                                            Phi0 = em.Phi0, 
                                            Phi1 = em.Phi1, 
                                            eps = 1e-12)
                        clade_rarity = {}
                        clade_alpha = {}
                        clade_beta = {}
                        for c in em.clades:
                            leaves_in_c = [leaf for leaf in em.leaf2type.keys() if leaf2clade[leaf] == c]
                            mean_rarity = np.mean([rarity_scores[leaf] for leaf in leaves_in_c])
                            clade_rarity[c] = mean_rarity
        
                            # if mean_rarity > 1.5:
                            # if mean_rarity > 2:
                            if mean_rarity > upper_rarity_bound:
                                alpha, beta = 5.0, 1.0
                            elif mean_rarity > 1.2:
                                alpha, beta = 2.0, 1.0
                            else:
                                alpha, beta = 1.0, 1.0 # flat prior
    
                            clade_alpha[c] = alpha
                            clade_beta[c] = beta

    
                    elif clade_priors == 'enrichment':

                        clade_alpha = {}
                        clade_beta = {}
                        
                        clade2cells = {}
                        for leaf, clade in leaf2clade.items():
                            clade2cells[clade] = clade2cells.get(clade, []) + [leaf]
                        
                        
                        # clade_prior_param_dict = {}
                        for clade, cells in clade2cells.items():
                            
                            
                            # prior_mean = clade_prior_from_enrichment(clade_cells = cells, 
                            #                                             leaf2type = leaf2type, 
                            #                                             type_enrichment = type_prior)
                        
                            prior_mean = clade_prior_from_celltype_enrichment(clade_cells = cells, 
                                                                              leaf2type = leaf2type, 
                                                                              celltype_induction_score_map = enrichment_dict)
                        
                            beta_params = beta_parameters(prior_mean, strength = strength)

                            clade_alpha[clade] = beta_params[0]
                            clade_beta[clade] = beta_params[1]

                    elif isinstance(clade_priors, dict):
                        # assumes {clade_name: [alpha, beta]}
                        clade_alpha = {clade_name: vals[0] for clade_name, vals in clade_priors.items()}
                        clade_beta = {clade_name: vals[1] for clade_name, vals in clade_priors.items()}


                else:
                    clade_alpha = dict(zip(em.clades, [1 for _ in range(len(em.clades))]))
                    clade_beta = dict(zip(em.clades, [1 for _ in range(len(em.clades))]))
    
                            


    
                    


                assert np.allclose(np.array(ex_uninduced_tm).sum(axis = 1), 1.0, atol = 1e-3), 'uninduced should be transposed'
                assert np.allclose(np.array(ex_induced_tm).sum(axis = 1), 1.0, atol = 1e-3), 'induced should be transposed'

                results, best, init_p = run_em_with_restarts(
                    Phi0 = ex_uninduced_tm,
                    Phi1 = ex_induced_tm,
                    root_type = root_type,
                    tree = tree,
                    leaf2clade = leaf2clade,
                    leaf2type = leaf2type,
                    type_encoder = type_encoder,
                    weight_strategy = weight_strategy,
                    clade_alpha = clade_alpha,
                    clade_beta = clade_beta,
                    n_restarts = num_restarts,
                    max_iter = 200,
                    tol = 1e-3,
                    seed_start = 42,
                    verbose = verbose,
                    p_frac_init = p_frac_init,
                )
            
                # Access the best model
                best_model = best['model']
                
                
                em = best_model
                em.results = results
                
                
                df_post = em.clade_posteriors(use_gmm = use_gmm)
                
                
                
                tpt_param_combo_induction_dict[param_combo] = df_post
                tpt_em_model_dict[param_combo] = em


                # a few added steps post-addition to tpt_param_combo_induction_dict for evaluating auc, roc, etc.
                df_post['clade_name'] = param_combo
    
            param_combo_induction_dict[timepoint] = tpt_param_combo_induction_dict
            full_em_model_dict[timepoint] = tpt_em_model_dict

    # return param_combo_induction_dict

    growing_df_list = []
    for timepoint, timepoint_dict in param_combo_induction_dict.items():
        for clade_name, clade_df in timepoint_dict.items():
            clade_df['clade_name'] = clade_name
            clade_df['timepoint'] = timepoint
            growing_df_list.append(clade_df)
    param_combo_induction_df = pd.concat(growing_df_list, axis = 0)    

    final_df = param_combo_induction_df.reset_index(drop = True)
    final_df_copy = copy.deepcopy(final_df)

    eval_dict_per_combo = {}

    for timepoint in final_df['timepoint'].unique():
        timepoint_df = final_df[final_df['timepoint'] == timepoint]

        full_df = make_clade_celltype_df(urid = urid, timepoint = timepoint, urid_timepoint_res_dict = urid_timepoint_res_dict, path_to_output_dir = path_to_output_dir)
        full_subset_gt = full_df[full_df['param_combo'] == 'clade_differentiation_induced']
        gt_linstring_to_induced_dict = dict(zip(full_subset_gt['cell'], ['induced' if clade else 'un-induced' for clade in full_subset_gt['clade']]))

        for clade in timepoint_df['clade_name'].unique():
            temp_df = timepoint_df[timepoint_df['clade_name'] == clade]

            this_em = full_em_model_dict[timepoint][clade]


            display(temp_df)

            metrics, (precision, recall), (fpr, tpr), (predicted_probs, predicted_label, ground_truth) = evaluate_prediction_probs(em = this_em,
                                                                                 leaf2clade = this_em.leaf2clade,
                                                                                 leaf2induced_gt = gt_linstring_to_induced_dict,
                                                                                 plot = True,
                                                                                                                                  res_df = temp_df)
            eval_dict_per_combo[clade] = {'metrics': metrics,
                         'precision': precision,
                         'recall': recall,
                         'fpr': fpr,
                         'tpr': tpr,
                         'predicted_probs': predicted_probs,
                         'predicted_label': predicted_label,
                         'ground_truth': ground_truth}

    if return_opt_model_dict:
        return final_df_copy, full_em_model_dict, eval_dict_per_combo

    return final_df_copy

# Analysis and plotting


def matrix_separation(Phi0, Phi1, prior = None, max_steps = 60):
    """
    Compare two transition matrices across step counts.

    Metrics per step n:
      marginal_tv  - TV between marginals prior @ Phi0^n vs prior @ Phi1^n.
                     Collapses source type; the original metric.
      expected_tv  - E_src[ TV(Phi0^n[src,:], Phi1^n[src,:]) ] under prior.
                     Keeps co-occurrence with source: how distinguishable are
                     the two matrices starting from a typical cell?
      max_tv       - max over source types of TV(Phi0^n[src,:], Phi1^n[src,:]).
                     Which starting type maximises separation at step n?
      path2_tv     - TV between 2-step path distributions
                     P(X_0=i, X_1=j, X_n=k) under Phi0 vs Phi1.
                     Captures co-occurrence across intermediate steps.
    """
    K = Phi0.shape[0]
    if prior is None:
        prior = np.ones(K) / K

    records = []
    for n in range(1, max_steps + 1):
        Pn0 = matrix_power(Phi0, n)
        Pn1 = matrix_power(Phi1, n)

        # original: marginal TV
        marginal_tv = 0.5 * np.abs(prior @ Pn0 - prior @ Pn1).sum()

        # per-source TV, weighted by prior
        tv_per_src = 0.5 * np.abs(Pn0 - Pn1).sum(axis = 1)   # shape (K,)
        expected_tv = float(prior @ tv_per_src)
        max_tv = float(tv_per_src.max())

        # 2-step path joint: P(X_0=i, X_1=j, X_n=k) = prior[i]*Phi[i,j]*Phi^(n-1)[j,k]
        if n > 1:
            Pm0 = matrix_power(Phi0, n - 1)
            Pm1 = matrix_power(Phi1, n - 1)
            # broadcasting: (K,1,1)*(K,K,1)*(1,K,K) -> (K,K,K)
            path0 = prior[:, None, None] * Phi0[:, :, None] * Pm0[None, :, :]
            path1 = prior[:, None, None] * Phi1[:, :, None] * Pm1[None, :, :]
        else:
            path0 = prior[:, None] * Phi0
            path1 = prior[:, None] * Phi1
        path2_tv = float(0.5 * np.abs(path0 - path1).sum())

        records.append((n, marginal_tv, expected_tv, max_tv, path2_tv))

    return pd.DataFrame(records,
                        columns = ['n_steps', 'marginal_tv', 'expected_tv', 'max_tv', 'path2_tv'])

def plot_matrix_sep(urid, path_to_output_dir = path_to_output_dir, cell_type_names = None,
                    save_path = None):
    phi0, phi1 = map(np.array, get_tms(urid, path_to_output_dir = path_to_output_dir))
    K = phi0.shape[0]
    if cell_type_names is None:
        cell_type_names = [f'ct{i+1}' for i in range(K)]

    sep_df = matrix_separation(phi0, phi1)

    # Row 0: TV metrics over step count
    fig, axes = plt.subplots(1, 3, figsize = (14, 4))

    cols = ['marginal_tv', 'expected_tv', 'max_tv', 'path2_tv']
    colors = ['#95a5a6',     '#2980b9',     '#c0392b', '#27ae60']
    labels = ['marginal TV (old)', 'expected TV', 'max TV', '2-step path TV']

    ax = axes[0]
    for col, color, label in zip(cols, colors, labels):
        ax.plot(sep_df['n_steps'], sep_df[col], color = color, linewidth = 1.6, label = label)
    peak_n = int(sep_df.loc[sep_df['expected_tv'].idxmax(), 'n_steps'])
    ax.axvline(peak_n, color = 'black', linestyle = ':', linewidth = 0.9, label = f"peak n={peak_n}")
    ax.set_xlabel('Steps (n)'); ax.set_ylabel('TV distance')
    ax.set_title('Separation metrics vs step count')
    ax.legend(fontsize = 7, frameon = False)

    # Row 1 & 2: source-level TV at peak n
    from numpy.linalg import matrix_power as _mp
    Pn0 = _mp(phi0, peak_n)
    Pn1 = _mp(phi1, peak_n)
    tv_per_src = 0.5 * np.abs(Pn0 - Pn1).sum(axis = 1)

    ax = axes[1]
    bar_colors = plt.cm.RdYlGn(tv_per_src / tv_per_src.max())
    ax.bar(cell_type_names, tv_per_src, color = bar_colors, edgecolor = 'white')
    ax.set_xlabel('Source cell type'); ax.set_ylabel(f"TV at n={peak_n}")
    ax.set_title(f"Per-source separation at peak n={peak_n}")

    # Row 3: diff heatmap at peak n
    ax = axes[2]
    diff = phi0 - phi1   # 1-step difference; most interpretable
    vmax = np.abs(diff).max()
    im = ax.imshow(diff, cmap = 'RdBu_r', vmin = -vmax, vmax = vmax, aspect = 'auto')
    ax.set_xticks(range(K)); ax.set_xticklabels(cell_type_names, rotation = 45, ha = 'right')
    ax.set_yticks(range(K)); ax.set_yticklabels(cell_type_names)
    ax.set_xlabel('Target type'); ax.set_ylabel('Source type')
    ax.set_title('Phi0  Phi1  (1-step, red = Phi0 higher)')
    plt.colorbar(im, ax = ax, fraction = 0.046, pad = 0.04)

    fig.suptitle(f"Matrix separation - URID {urid}", fontweight = 'bold')
    plt.tight_layout()
    
    if save_path:
        fig.savefig(save_path, bbox_inches = 'tight')
    plt.close(fig)

    return sep_df, peak_n

def build_annotated_dfs(tpt_df_dict, urid, timepoints, urid_timepoint_res_dict,
                        param_combo = '20clade_topological_treedist_gt',
                        path_to_output_dir = path_to_output_dir):
    """
    Per-timepoint clade df augmented with majority-vote GT label.
    Returns {tp: DataFrame(clade, p_induced, gt_label)}.
    """
    annotated = {}
    for tp in timepoints:
        full_df = make_clade_celltype_df(urid = urid, timepoint = tp,
                                         urid_timepoint_res_dict = urid_timepoint_res_dict,
                                         path_to_output_dir = path_to_output_dir)
        target = full_df[full_df['param_combo'] == param_combo][['cell', 'clade']]
        gt_raw = (full_df[full_df['param_combo'] == 'clade_differentiation_induced']
                  [['cell', 'clade']].rename(columns = {'clade': 'induced'}))
        merged = target.merge(gt_raw, on = 'cell', how = 'inner')
        gt_map = (merged.groupby('clade')['induced']
                  .apply(lambda x: 'induced' if np.mean(x.astype(float)) > 0.5 else 'uninduced')
                  .to_dict())
        df = tpt_df_dict[tp][tpt_df_dict[tp]['clade_name'] == param_combo].copy()
        df['gt_label'] = df['clade'].map(gt_map)
        annotated[tp] = df
    return annotated

def build_cell_level_dfs(tpt_df_dict, urid, timepoints, urid_timepoint_res_dict,
                         param_combo = '20clade_topological_treedist_gt',
                         path_to_output_dir = path_to_output_dir):
    """
    Expand clade-level EM predictions to individual cells.
    Returns {tp: DataFrame(clade, p_induced, y_true)} where y_true in {0, 1}.
    """
    cell_dfs = {}
    for tp in timepoints:
        full_df = make_clade_celltype_df(urid = urid, timepoint = tp,
                                          urid_timepoint_res_dict = urid_timepoint_res_dict,
                                          path_to_output_dir = path_to_output_dir)
        clade_to_p = (tpt_df_dict[tp][tpt_df_dict[tp]['clade_name'] == param_combo]
                      .set_index('clade')['p_induced'].to_dict())
        cell_to_clade = (full_df[full_df['param_combo'] == param_combo]
                         .set_index('cell')['clade'].to_dict())
        cell_to_gt = (full_df[full_df['param_combo'] == 'clade_differentiation_induced']
                      .set_index('cell')['clade'].to_dict())
        cells = list(set(cell_to_clade) & set(cell_to_gt))
        cdf = pd.DataFrame({
            'clade':     [cell_to_clade[c] for c in cells],
            'p_induced': [clade_to_p.get(cell_to_clade[c], np.nan) for c in cells],
            'y_true':    [int(float(cell_to_gt[c]) > 0.5) for c in cells],
        }).dropna(subset = ['p_induced']).reset_index(drop = True)
        cell_dfs[tp] = cdf
    return cell_dfs

def plot_p_induced_per_clade(annotated_dfs, timepoints,
                              induced_color = induced_color, uninduced_color = uninduced_color,
                              title = 'EM p(induced) per clade', save_path = None):
    """Sorted scatter of clade p_induced coloured by GT label."""
    fig, axes = plt.subplots(1, len(timepoints), figsize = (10, 3.2), sharey = True)
    if len(timepoints) == 1:
        axes = [axes]
    for ax, tp in zip(axes, timepoints):
        df = annotated_dfs[tp].sort_values('p_induced')
        colors = [induced_color if g == 'induced' else uninduced_color for g in df['gt_label']]
        ax.scatter(range(len(df)), df['p_induced'], c = colors, s = 40, zorder = 3, edgecolors = 'none')
        ax.axhline(0.5, color = '#7f8c8d', linewidth = 0.8, linestyle = '--', zorder = 2)
        ax.set_title(f'Timepoint {tp}')
        ax.set_xlabel('Clade (sorted by p)')
        ax.set_ylim(0, 1); ax.set_xticks([])
    axes[0].set_ylabel('p(induced)')
    legend_elements = [
        Line2D([0],[0], marker = 'o', color = 'w', markerfacecolor = induced_color,  markersize = 7, label = 'Induced (GT)'),
        Line2D([0],[0], marker = 'o', color = 'w', markerfacecolor = uninduced_color, markersize = 7, label = 'Uninduced (GT)'),
    ]
    axes[-1].legend(handles = legend_elements, loc = 'upper left', frameon = False)
    if title:
        fig.suptitle(title, y = 1.02)
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, bbox_inches = 'tight')
    return fig

def plot_roc_curves(eval_dict_bytp, timepoints, title = 'ROC - EM classifier', save_path = None):
    """ROC curve per timepoint on a single axes."""
    validtps = [tp for tp in timepoints if eval_dict_bytp.get(tp, {}).get('fpr') is not None]
    tp_colors = plt.cm.viridis(np.linspace(0.1, 0.9, max(len(validtps), 1)))
    fig, ax = plt.subplots(figsize = (4, 4))
    for tp, col in zip(validtps, tp_colors):
        ed = eval_dict_bytp[tp]
        ax.plot(ed['fpr'], ed['tpr'], color = col, linewidth = 1.5,
                label = f"tp {tp}  (AUC={auc(ed['fpr'], ed['tpr']):.2f})")
    ax.plot([0,1],[0,1], color = '#bdc3c7', linewidth = 0.8, linestyle = '--')
    ax.set_xlabel('False positive rate'); ax.set_ylabel('True positive rate')
    ax.set_title(title); ax.legend(frameon = False)
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, bbox_inches = 'tight')
    return fig

def plot_pr_curves(eval_dict_bytp, timepoints, baselinetp = None,
                   title = 'Precision-Recall - EM classifier', save_path = None):
    """Precision-Recall curve per timepoint with no-skill baseline."""
    validtps = [tp for tp in timepoints if eval_dict_bytp.get(tp, {}).get('recall') is not None]
    tp_colors = plt.cm.viridis(np.linspace(0.1, 0.9, max(len(validtps), 1)))
    fig, ax = plt.subplots(figsize = (4, 4))
    for tp, col in zip(validtps, tp_colors):
        ed = eval_dict_bytp[tp]
        ax.plot(ed['recall'], ed['precision'], color = col, linewidth = 1.5,
                label = "tp {tp}  (AUC={auc(ed['recall'], ed['precision']):.2f})")
    reftp = baselinetp if baselinetp is not None else (validtps[0] if validtps else None)
    gt = eval_dict_bytp[reftp]['ground_truth'] if reftp is not None else []
    baseline = sum(1 for g in gt if g == 1) / len(gt) if len(gt) > 0 else 0.5
    ax.axhline(baseline, color = '#bdc3c7', linewidth = 0.8, linestyle = '--',
               label = f'Baseline ({baseline:.2f})')
    ax.set_xlabel('Recall'); ax.set_ylabel('Precision')
    ax.set_title(title); ax.set_xlim(0,1); ax.set_ylim(0,1.05)
    ax.legend(frameon = False)
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, bbox_inches = 'tight')
    return fig

def plot_p_induced_distribution(annotated_dfs, timepoints,
                                 induced_color = induced_color, uninduced_color = uninduced_color,
                                 title = 'EM p(induced) distribution by ground-truth label',
                                 save_path = None):
    """Violin + jitter of clade p_induced grouped by GT label."""
    fig, axes = plt.subplots(1, len(timepoints), figsize = (10, 3.2), sharey = True)
    if len(timepoints) == 1:
        axes = [axes]
    for ax, tp in zip(axes, timepoints):
        df = annotated_dfs[tp]
        for label, col, xpos in [('induced', induced_color, 1), ('uninduced', uninduced_color, 0)]:
            vals = df[df['gt_label'] == label]['p_induced'].values
            if len(vals) < 2 or np.std(vals) == 0:
                continue
            parts = ax.violinplot(vals, positions = [xpos], widths = 0.6, showmedians = True, showextrema = False)
            for pc in parts['bodies']:
                pc.set_facecolor(col); pc.set_alpha(0.6); pc.set_edgecolor('none')
            parts['cmedians'].set_color(col); parts['cmedians'].set_linewidth(2)
            ax.scatter(np.random.normal(xpos, 0.05, len(vals)), vals,
                       color = col, s = 18, alpha = 0.7, zorder = 3, edgecolors = 'none')
        ax.axhline(0.5, color = '#7f8c8d', linewidth = 0.8, linestyle = '--')
        ax.set_title(f'Timepoint {tp}')
        ax.set_xticks([0,1]); ax.set_xticklabels(['Uninduced','Induced'], rotation = 20, ha = 'right')
        ax.set_ylim(0, 1)
    axes[0].set_ylabel('p(induced)')
    if title:
        fig.suptitle(title, y = 1.02)
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, bbox_inches = 'tight')
    return fig

def plot_cell_concordance(cell_level_dfs, timepoints, threshold = 0.5, n_bins = 5,
                          induced_color = induced_color, uninduced_color = uninduced_color,
                          title = 'Cell-level concordance: EM p(induced) vs ground truth',
                          save_path = None):
    """
    3-row x len(timepoints) figure:
      Row 0 - Waterfall bars sorted by p_induced, coloured by GT, wrong predictions outlined
      Row 1 - Calibration reliability diagram
      Row 2 - Confusion matrix at threshold with raw counts and fraction-of-total
    """
    fig, axes = plt.subplots(3, len(timepoints), figsize = (13, 8.5),
                             gridspec_kw = {'height_ratios':[2.5,1.5,1.5], 'hspace':0.55, 'wspace':0.3})
    if title:
        fig.suptitle(title, fontsize = 10, fontweight = 'bold', y = 1.01)

    for col, tp in enumerate(timepoints):
        cdf = cell_level_dfs[tp].sort_values('p_induced').reset_index(drop = True)
        y_pred = (cdf['p_induced'] >= threshold).astype(int)
        n = len(cdf)
        bar_colors = [induced_color if y else uninduced_color for y in cdf['y_true']]

        # row 0: waterfall
        ax = axes[0, col]
        ax.bar(np.arange(n), cdf['p_induced'], color = bar_colors, width = 1.0, linewidth = 0)
        wrong = cdf.index[y_pred != cdf['y_true']].tolist()
        if wrong:
            ax.bar(wrong, cdf.loc[wrong, 'p_induced'],
                   color = [bar_colors[i] for i in wrong],
                   width = 1.0, linewidth = 0.6, edgecolor = '#111111')
        ax.axhline(threshold, color = '#7f8c8d', linewidth = 0.8, linestyle = '--')
        ax.set_xlim(-1, n); ax.set_ylim(0, 1)
        ax.set_xticks([]); ax.set_title(f'Timepoint {tp}', fontsize = 9)
        if col == 0:
            ax.set_ylabel('p(induced)')
        acc = (y_pred == cdf['y_true']).mean()
        ax.text(0.97, 0.95, f'acc={acc:.0%}', transform = ax.transAxes,
                ha = 'right', va = 'top', fontsize = 7.5, color = '#333333')

        # row 1: calibration
        ax = axes[1, col]
        bins = np.linspace(0, 1, n_bins + 1)
        bin_centers, frac_pos, bin_counts = [], [], []
        for lo, hi in zip(bins[:-1], bins[1:]):
            mask = (cdf['p_induced'] >= lo) & (cdf['p_induced'] < hi)
            if mask.sum() > 0:
                frac_pos.append(cdf.loc[mask, 'y_true'].mean())
                bin_centers.append((lo + hi) / 2)
                bin_counts.append(mask.sum())
        ax.plot([0,1],[0,1], color = '#bdc3c7', linewidth = 0.8, linestyle = '--', zorder = 1)
        ax.scatter(bin_centers, frac_pos, s = [10 + 3*c for c in bin_counts], color = '#27ae60', zorder = 3)
        ax.plot(bin_centers, frac_pos, color = '#27ae60', linewidth = 1.3, zorder = 2)
        ax.set_xlim(0,1); ax.set_ylim(0,1)
        ax.set_xticks([0,0.5,1]); ax.set_yticks([0,0.5,1])
        ax.tick_params(labelsize = 7)
        if col == 0:
            ax.set_ylabel('Frac. truly induced', fontsize = 8)
        ax.set_xlabel('Mean p(induced)', fontsize = 8)

        # row 2: confusion matrix
        ax = axes[2, col]
        cm = confusion_matrix(cdf['y_true'], y_pred, labels = [1, 0])
        cm_norm = cm.astype(float) / cm.sum()
        ax.imshow(cm_norm, cmap = 'Blues', vmin = 0, vmax = cm_norm.max()*1.2, aspect = 'auto')
        for i, j in itertools.product(range(2), range(2)):
            ax.text(j, i, f'{cm[i,j]}\n({cm_norm[i,j]:.0%})',
                    ha = 'center', va = 'center', fontsize = 7.5,
                    color = 'white' if cm_norm[i,j] > 0.45 else 'black')
        ax.set_xticks([0,1]); ax.set_xticklabels(['Pred +','Pred '], fontsize = 7.5)
        ax.set_yticks([0,1]); ax.set_yticklabels(['GT +','GT '], fontsize = 7.5)
        if col == 0:
            ax.set_ylabel('Ground truth', fontsize = 8)

    legend_handles = [
        Patch(facecolor = induced_color,   label = 'Induced (GT)'),
        Patch(facecolor = uninduced_color, label = 'Uninduced (GT)'),
        Patch(facecolor = '#aaaaaa', edgecolor = '#111111', linewidth = 0.8, label = 'Wrong prediction'),
    ]
    axes[0,-1].legend(handles = legend_handles, loc = 'upper left', frameon = False, fontsize = 7)
    if save_path:
        fig.savefig(save_path, bbox_inches = 'tight')
    return fig

def make_all_em_plots(param_combo_list,
                           full_modelres_dict,
                           urid_timepoint_res_dict,
                           return_plotting_dfs = False):
    
    full_plotting_df_dict = {}
    
    for urid in list(full_modelres_dict.keys()):
        
        full_plotting_df_dict[urid] = {}
    
        
        path_to_output_dir = Path('./output')
        path_to_em_fig_dir = path_to_output_dir / 'em_figs' / str(urid)

        os.makedirs(path_to_em_fig_dir, exist_ok = True)

        plot_matrix_sep(urid, path_to_output_dir = path_to_output_dir,
                        save_path = path_to_em_fig_dir / 'matrix_sep')


        for param_combo in param_combo_list:

            path_to_param_combo_dir = path_to_em_fig_dir / param_combo

            os.makedirs(path_to_param_combo_dir, exist_ok = True)

            timepoints = sorted(list(full_modelres_dict[urid].keys()))
            
            tpt_df_dict = {}
            eval_dict_bytp = {}
            
            for tp in timepoints:
                tpt_df_dict[tp] = full_modelres_dict[urid][tp]['tpt_induction_dict']
                ed = full_modelres_dict[urid][tp]['eval_dict'].get(param_combo)
                if ed is None:
                    print(f'WARNING: no eval_dict for urid={urid} tp={tp} combo={param_combo}', flush = True)
                    ed = {}
                eval_dict_bytp[tp] = ed


            annotated_dfs = build_annotated_dfs(
                tpt_df_dict, urid = urid, timepoints = timepoints,
                urid_timepoint_res_dict = urid_timepoint_res_dict,
                param_combo = param_combo,
                path_to_output_dir = path_to_output_dir)

            fig1 = plot_p_induced_per_clade(annotated_dfs, timepoints,
                                            save_path = path_to_param_combo_dir / f'em_p_induced_per_clade_{urid}_{param_combo}.pdf')
            plt.close(fig1)

            fig2 = plot_roc_curves(eval_dict_bytp, timepoints, save_path = path_to_param_combo_dir / f'em_roc_{urid}_{param_combo}.pdf')
            plt.close(fig2)

            fig3 = plot_pr_curves(eval_dict_bytp, timepoints, save_path = path_to_param_combo_dir / f'em_pr_{urid}_{param_combo}.pdf')
            plt.close(fig3)

            fig4 = plot_p_induced_distribution(annotated_dfs, timepoints,
                                                save_path = path_to_param_combo_dir / f'em_p_induced_distribution_{urid}_{param_combo}.pdf')
            plt.close(fig4)

            cell_level_dfs = build_cell_level_dfs(
                tpt_df_dict, urid = urid, timepoints = timepoints,
                urid_timepoint_res_dict = urid_timepoint_res_dict,
                param_combo = param_combo,
                path_to_output_dir = path_to_output_dir)

            fig5 = plot_cell_concordance(cell_level_dfs, timepoints,
                                        save_path = path_to_param_combo_dir / f'em_cell_concordance_{urid}_{param_combo}.pdf')
            plt.close(fig5)
            
            full_plotting_df_dict[urid][param_combo] = annotated_dfs
            
    if return_plotting_dfs:
        return full_plotting_df_dict

# CLI workers


def worker_process_urid(urid):
    """Phase-1 worker: load & cluster one urid."""
    sys.setrecursionlimit(100000)
    global path_to_output_dir, path_to_em_induction_prob
    path_to_output_dir = Path('./output')
    path_to_em_induction_prob = path_to_output_dir / 'em_induction_prob'
    timepoint_res_dict, _ = process_urid(
        urid = urid,
        param_dict_list = longmito_36_first_pass_list,
        cluster_vmin = -0.5,
        cluster_vmax = 1,
        timepoints = [8, 9, 10, 11, 12],
    )
    return urid, timepoint_res_dict

def worker_inference(args):
    """Phase-3 worker: run EM inference for one (urid, timepoint) pair."""
    sys.setrecursionlimit(100000)
    urid, timepoint = args
    global path_to_output_dir, path_to_em_induction_prob
    path_to_output_dir = Path('./output')
    path_to_em_induction_prob = path_to_output_dir / 'em_induction_prob'

    print(f'[EM start] urid={urid} tp={timepoint}', flush = True)

    hsc_enc = {'LT': 0, 'ST': 1, 'MyMPP': 2, 'LyMPP': 3, 'Mye': 4, 'Lym': 5}

    tpt_induction_dict, opt_mod_res, eval_dict = urid_to_induction_inference(
        [urid], [timepoint],
        param_dict_list = longmito_36_first_pass_list,
        weight_strategy = 'uniform',
        clade_priors = None,
        upper_rarity_bound = 1.5,
        include_param_combos = None,
        scoremat_filters = scoremat_filters,
        verbose = False,
        root_type = 'LT',
        num_restarts = 1,
        return_opt_model_dict = True,
        type_encoder = hsc_enc,
        enrichment_dict = None,
        strength = 10,
        use_gmm = True,
        path_to_output_dir = path_to_output_dir,
    )
    open(f'./progress_txt_files/progress_{urid}_{timepoint}.txt', 'w').close()
    # Strip tree references from EM models before pickling across the queue.
    # make_all_em_plots only uses tpt_induction_dict and eval_dict;
    # opt_mod_res is saved for inspection but the tree is not accessed again.
    for tp_dict in opt_mod_res.values():
        for em in tp_dict.values():
            if hasattr(em, 'tree'):
                em.tree = None
    return urid, timepoint, tpt_induction_dict, opt_mod_res, eval_dict

def alpha_get_tms(urid, output_dir = None):
    """Return (Phi0, Phi1_true, founder_type, type_order) from run_specs JSON."""
    od = Path(output_dir) if output_dir else Path('./output')
    candidates = [
        od / str(urid) / 'run_specs' / str(urid),
        od / str(urid) / 'run_specs',
        od / 'run_specs' / str(urid),
    ]
    jsons = []
    for p in candidates:
        jsons = list(p.glob('*.json'))
        if jsons:
            break
    if not jsons:
        raise FileNotFoundError(f"No run_specs JSON for urid={urid}; tried {candidates}")
    d = json.loads(jsons[0].read_text())
    ct = d['cell_type_dict']
    founder = ct.get('founder_cell_type', 'ct1')
    type_order = list(ct['cell_type_params'].keys())
    return (np.array(ct['uninduced_transition_matrix'], dtype = float),
            np.array(ct['induced_transition_matrix'],   dtype = float),
            founder, type_order)

def alpha_load_clade_df(urid, timepoint, output_dir = None):
    od = Path(output_dir) if output_dir else Path('./output')
    path = od / f"clade_assignment_dfs/{urid}/{timepoint}/clade_df_{urid}_t{timepoint}.csv"
    df = pd.read_csv(path, index_col = 0)
    df.index.name = 'cell'
    return df.reset_index()

def alpha_get_gt_tree_path(urid, timepoint, output_dir = None):
    od = Path(output_dir) if output_dir else Path('./output')
    urid_dir = od / f"processed_newicks/{urid}"
    matches = list(urid_dir.glob(f"ground_truth_tree_*_time_{timepoint}.newick"))
    if not matches:
        raise FileNotFoundError(f"GT tree not found for urid={urid} tp={timepoint}")
    return str(matches[0])

def alpha_get_recon_tree_path(urid, timepoint, output_dir = None):
    od = Path(output_dir) if output_dir else Path('./output')
    candidates = [
        od / 'recon_trees' / str(urid),                 # cluster: output/recon_trees/{urid}/
        od / str(urid) / 'recon_trees' / str(urid),     # local: output/{urid}/recon_trees/{urid}/
    ]
    glob_pat = f"*RP_1_samp*_time_{timepoint}_CD_T_AF_0_B_T/*.treefile"
    for recon_dir in candidates:
        matches = list(recon_dir.glob(glob_pat))
        if matches:
            return str(matches[0])
    return None

def alpha_build_leaf_maps(clade_df, param_combo):
    if param_combo not in clade_df.columns:
        return None, None, None
    sub = clade_df[
        ['cell', param_combo, 'clade_celltype', 'clade_differentiation_induced']
    ].dropna(subset = [param_combo, 'clade_celltype'])
    if len(sub) < 2:
        return None, None, None
    leaf2clade = dict(zip(sub['cell'], sub[param_combo].astype(str)))
    leaf2type = dict(zip(sub['cell'], sub['clade_celltype'].astype(str)))
    leaf2induced_gt = {
        cell: ('induced' if bool(val) else 'un-induced')
        for cell, val in zip(sub['cell'], sub['clade_differentiation_induced'])
    }
    return leaf2clade, leaf2type, leaf2induced_gt

def worker_alpha_sweep(args):
    """ProcessPoolExecutor worker: alpha sweep for one (urid, timepoint, param_combo)."""
    (urid, timepoint, param_combo, tree_path,
     leaf2clade, leaf2type, leaf2induced_gt,
     Phi0, Phi1_true, alpha_values, type_encoder,
     n_restarts, root_type, weight_strategy) = args

    # load and prune tree (done once; all alphas reuse it)
    raw_newick = ''.join(open(tree_path).read().split())
    clean_nwk = re.sub(r'\)([^:);]+):', r'):', raw_newick)
    tree = Tree(clean_nwk, format = 1)

    keep = set(leaf2clade.keys()) & {n.name for n in tree.iter_leaves()}
    if len(keep) < 4:
        print(f"[alpha_sweep] {urid} tp={timepoint} {param_combo}: {len(keep)} leaves, skip", flush = True)
        return []
    tree.prune(list(keep), preserve_branch_length = True)

    leaf2clade = {c: v for c, v in leaf2clade.items()      if c in keep}
    leaf2type = {c: v for c, v in leaf2type.items()       if c in keep}
    leaf2induced_gt = {c: v for c, v in leaf2induced_gt.items() if c in keep}

    clades = sorted(set(leaf2clade.values()))
    clade_alpha = {c: 1.0 for c in clades}
    clade_beta = {c: 1.0 for c in clades}

    cell_names = list(leaf2clade.keys())
    gt_arr = np.array([1 if leaf2induced_gt[c] == 'induced' else 0 for c in cell_names])

    print(f"[alpha_sweep] {urid} tp={timepoint} {param_combo}: "
          f"{len(cell_names)} cells, {len(alpha_values)} alpha values", flush = True)

    try:
        alpha_results = run_em_with_restarts(
            Phi0 = Phi0, Phi1 = Phi0,           # Phi1 unused when alpha_values provided
            root_type = root_type, tree = tree,
            leaf2clade = leaf2clade, leaf2type = leaf2type,
            type_encoder = type_encoder, weight_strategy = weight_strategy,
            clade_alpha = clade_alpha, clade_beta = clade_beta,
            n_restarts = n_restarts, max_iter = 200, tol = 1e-3,
            seed_start = 42, verbose = False,
            alpha_values = alpha_values, Phi1_true = Phi1_true,
        )
    except Exception as exc:
        print(f"[alpha_sweep] EM error {urid} tp={timepoint} {param_combo}: {exc}", flush = True)
        return []

    rows = []
    for alpha, best in alpha_results.items():
        em = best['model']
        df_post = em.clade_posteriors(use_gmm = False)

        clade_to_p = dict(zip(df_post['clade'].astype(str), df_post['p_induced']))
        clade_to_hl = dict(zip(df_post['clade'].astype(str), df_post['hard_label']))

        pred_prob = np.array([clade_to_p.get(str(leaf2clade[c]), 0.5) for c in cell_names])
        pred_hard = np.array([
            1 if clade_to_hl.get(str(leaf2clade[c]), 'un-induced') == 'induced' else 0
            for c in cell_names
        ])

        if len(np.unique(gt_arr)) < 2:
            roc_auc_val = float('nan')
            pr_auc_val = float('nan')
        else:
            roc_auc_val = roc_auc_score(gt_arr, pred_prob)
            if roc_auc_val < 0.5:
                pred_prob = 1.0 - pred_prob
                pred_hard = 1   - pred_hard
                roc_auc_val = 1.0 - roc_auc_val
            prec_c, rec_c, _ = precision_recall_curve(gt_arr, pred_prob, pos_label = 1)
            pr_auc_val = auc(rec_c, prec_c)

        tp_n = int(np.sum((pred_hard == 1) & (gt_arr == 1)))
        fp_n = int(np.sum((pred_hard == 1) & (gt_arr == 0)))
        tn_n = int(np.sum((pred_hard == 0) & (gt_arr == 0)))
        fn_n = int(np.sum((pred_hard == 0) & (gt_arr == 1)))

        prec_h = tp_n / (tp_n + fp_n) if (tp_n + fp_n) > 0 else 0.0
        rec_h = tp_n / (tp_n + fn_n) if (tp_n + fn_n) > 0 else 0.0
        f1_h = 2 * prec_h * rec_h / (prec_h + rec_h) if (prec_h + rec_h) > 0 else 0.0

        best_logL = best.get('final_logL', best.get('logL', float('nan')))

        rows.append(dict(
            urid = urid, timepoint = timepoint, param_combo = param_combo, alpha = alpha,
            roc_auc = round(roc_auc_val, 6),
            pr_auc = round(pr_auc_val, 6),
            f1_hard = round(f1_h, 6),
            precision = round(prec_h, 6),
            recall = round(rec_h, 6),
            tp = tp_n, fp = fp_n, tn = tn_n, fn = fn_n,
            pred_pos_frac = round((tp_n + fp_n) / max(tp_n + fp_n + tn_n + fn_n, 1), 4),
            n_cells = len(cell_names),
            n_clades = len(df_post),
            best_logL = (round(best_logL, 4) if not np.isnan(best_logL) else float('nan')),
        ))

    print(f"[alpha_sweep done] {urid} tp={timepoint} {param_combo} -> {len(rows)} rows", flush = True)
    return rows

# Module-level parameters and globals


Phi0 = np.array([[0.9, 0.05, 0.03, 0.02],
                 [0.04, 0.9, 0.04, 0.02],
                 [0.04, 0.04, 0.9, 0.02],
                 [0.02, 0.02, 0.02, 0.94]])

Phi1 = np.array([[0.9, 0.08, 0.0, 0.02],
                 [0.0, 0.9, 0.08, 0.02],
                 [0.08, 0.0, 0.9, 0.02],
                 [0.02, 0.02, 0.02, 0.94]])




scoremat_min1 = {'name': 'scoremat_min1',
                 'scheme': 'scoremat',
                 'use_branch_length': None,
                 'num_clades': None,
                 'distance_metric': None,
                 'dist_cutoff': None,
                 'include_in_plot': True,
           'min_bc_clade_size': 1,
           'min_mt_clade_size': 1}

scoremat_min5 = {'name': 'scoremat_min5',
                 'scheme': 'scoremat',
                 'use_branch_length': None,
                 'num_clades': None,
                 'distance_metric': None,
                 'dist_cutoff': None,
                 'include_in_plot': True,
           'min_bc_clade_size': 1,
           'min_mt_clade_size': 5}

scoremat_min10 = {'name': 'scoremat_min10',
                 'scheme': 'scoremat',
                 'use_branch_length': None,
                 'num_clades': None,
                 'distance_metric': None,
                 'dist_cutoff': None,
                 'include_in_plot': True,
           'min_bc_clade_size': 1,
           'min_mt_clade_size': 10}

scoremat_min20 = {'name': 'scoremat_min20',
                 'scheme': 'scoremat',
                 'use_branch_length': None,
                 'num_clades': None,
                 'distance_metric': None,
                 'dist_cutoff': None,
                 'include_in_plot': True,
           'min_bc_clade_size': 1,
           'min_mt_clade_size': 20}

scoremat_2clades = {'name': 'scoremat_2clades',
                 'scheme': 'scoremat',
                 'use_branch_length': None,
                 'num_clades': 2,
                 'distance_metric': None,
                 'dist_cutoff': None,
                 'include_in_plot': True,
           'min_bc_clade_size': None,
           'min_mt_clade_size': None}

scoremat_10clades = {'name': 'scoremat_10clades',
                 'scheme': 'scoremat',
                 'use_branch_length': None,
                 'num_clades': 10,
                 'distance_metric': None,
                 'dist_cutoff': None,
                 'include_in_plot': True,
           'min_bc_clade_size': None,
           'min_mt_clade_size': None}

scoremat_20clades = {'name': 'scoremat_20clades',
                 'scheme': 'scoremat',
                 'use_branch_length': None,
                 'num_clades': 20,
                 'distance_metric': None,
                 'dist_cutoff': None,
                 'include_in_plot': True,
           'min_bc_clade_size': None,
           'min_mt_clade_size': None}

param_dict_gt5 = {'name': '5clade_topological_treedist',
                 'scheme': 'topological',
                 'use_branch_length': True,
                 'num_clades': 5,
                 'distance_metric': None,
                 'dist_cutoff': None}

param_dict_gt10 = {'name': '10clade_topological_treedist',
                 'scheme': 'topological',
                 'use_branch_length': True,
                 'num_clades': 10,
                 'distance_metric': None,
                 'dist_cutoff': None}

param_dict_gt20 = {'name': '20clade_topological_treedist',
                 'scheme': 'topological',
                 'use_branch_length': True,
                 'num_clades': 20,
                 'distance_metric': None,
                 'dist_cutoff': None}

param_dict_gt100 = {'name': '100clade_topological_treedist',
                 'scheme': 'topological',
                 'use_branch_length': True,
                 'num_clades': 100,
                 'distance_metric': None,
                 'dist_cutoff': None}

param_dict_gt50 = {'name': '50clade_topological_treedist',
                 'scheme': 'topological',
                 'use_branch_length': True,
                 'num_clades': 50,
                 'distance_metric': None,
                 'dist_cutoff': None}

param_dict_gt2 = {'name': '2clade_topological_treedist',
                 'scheme': 'topological',
                 'use_branch_length': True,
                 'num_clades': 2,
                 'distance_metric': None,
                 'dist_cutoff': None}

random_dict_10clades = {'name': '10clade_random',
                        'scheme': 'random',
                        'use_branch_length': None,
                        'num_clades': 10,
                        'distance_metric': None,
                        'dust_cutoff': None}

scoremat_param_dict_list = [scoremat_min5, scoremat_min10, scoremat_min20]

only_topological_clades = [param_dict_gt5, param_dict_gt10, param_dict_gt20]

longmito_36_first_pass_list = [scoremat_min5, scoremat_min10, scoremat_min20, param_dict_gt5, param_dict_gt20]

param_dict_list_100clades = [param_dict_gt100]

param_dict_list_50clades = [param_dict_gt50]

param_dict_list_2clades = [param_dict_gt2]

param_dict_list_hsc = [param_dict_gt2, param_dict_gt20, param_dict_gt50]

param_dict_list_rand10 = [random_dict_10clades, param_dict_gt2, param_dict_gt20, param_dict_gt50, scoremat_min10, scoremat_min20, scoremat_min5, scoremat_min1]

param_dict_score_only = [scoremat_min10, scoremat_min20, scoremat_min5, scoremat_min1, scoremat_2clades, scoremat_10clades, scoremat_20clades]

alpha_gt_combo = {
    'gt5':  '5clade_topological_treedist_gt',
    'gt20': '20clade_topological_treedist_gt',
}

alpha_scoremat_combo = {
    'scoremat_min5':  'scoremat_min5_mt_250ints_1.0RP_0.0AF_TB_1-1-1-1-1-1',
    'scoremat_min10': 'scoremat_min10_mt_250ints_1.0RP_0.0AF_TB_1-1-1-1-1-1',
    'scoremat_min20': 'scoremat_min20_mt_250ints_1.0RP_0.0AF_TB_1-1-1-1-1-1',
}


n_workers = 45   # set to e.g. 8 to cap usage

scoremat_filters = ['clade_topological_treedist_mt_250ints_1.0RP_0.0AF']

gt_combos = ['20clade_topological_treedist_gt', '5clade_topological_treedist_gt']

mt_scoremat_clades = ['5clade_topological_treedist_mt_250ints_1.0RP_0.0AF_TB_1-1-1-1-1-1',
             '20clade_topological_treedist_mt_250ints_1.0RP_0.0AF_TB_1-1-1-1-1-1']

default_timepoints = [8, 9, 10, 11, 12]


# Entry point

if __name__ == '__main__':

    # CLI
    parser = argparse.ArgumentParser(
        prog = 'longmito_36.py',
        description = 'Run clade assignment and/or EM induction inference.',
    )
    parser.add_argument(
        '--mode',
        choices = ['full', 'scoremat_sweep', 'no_em', 'alpha_sweep'],
        default = 'full',
        help = (
            'full          - run all phases (process_urid -> EM -> plots) [default]\n'
            'scoremat_sweep - run scoremat_clade_sweep only, no trees, no EM\n'
            'no_em         - run phases 1-2 only (process_urid + make_clade_df); stop before EM\n'
            'alpha_sweep   - sweep Phi1 prior strength via alpha; runs EM in parallel with ProcessPoolExecutor'
        ),
    )
    parser.add_argument(
        '--urids',
        nargs = '+',
        type = int,
        default = None,
        metavar = 'URID',
        help = 'Space-separated list of urids to process (default: long_mito_36_urids)',
    )
    parser.add_argument(
        '--timepoints',
        nargs = '+',
        type = int,
        default = None,
        metavar = 'T',
        help = 'Timepoints to include (default: default_timepoints constant = 8 9 10 11 12)',
    )
    # scoremat_sweep-specific
    parser.add_argument(
        '--min-clade-sizes',
        nargs = '+',
        type = int,
        default = [5, 10, 20],
        metavar = 'N',
        help = '(scoremat_sweep only) min_clade_size values to sweep (default: 5 10 20)',
    )
    parser.add_argument(
        '--output',
        default = 'scoremat_clade_sweep_df.csv',
        metavar = 'PATH',
        help = '(scoremat_sweep only) path for output CSV (default: scoremat_clade_sweep_df.csv)',
    )
    parser.add_argument(
        '--clade-output',
        default = None,
        metavar = 'PATH',
        help = 'Path for combined clade-assignments CSV (written after Phase 2, before EM, '
             'in any mode). If omitted, no CSV is saved unless --mode no_em also sets --output.',
    )
    parser.add_argument(
        '--joblib-dir',
        default = '.',
        metavar = 'DIR',
        help = 'Directory to write joblib checkpoint and final result files (default: current dir)',
    )
    parser.add_argument(
        '--joblib-prefix',
        default = 'longmito_results',
        metavar = 'PREFIX',
        help = 'Filename stem for joblib outputs (default: longmito_results). '
             'Files are written as <DIR>/<PREFIX>_model_dict_ckpt_N.joblib, '
             '<DIR>/<PREFIX>_model_dict_FINAL.joblib, <DIR>/<PREFIX>_roc_dict_FINAL.joblib',
    )
    # alpha_sweep-specific arguments
    parser.add_argument(
        '--alpha-values',
        nargs = '+',
        type = float,
        default = [0.0, 0.1, 0.25, 0.5, 0.75, 1.0],
        metavar = 'A',
        help = '(alpha_sweep only) alpha values to sweep (default: 0 0.1 0.25 0.5 0.75 1.0). '
             'alpha=0 means uninformative prior, alpha=1 means exact prior.',
    )
    parser.add_argument(
        '--urid-file',
        default = None,
        metavar = 'PATH',
        help = '(alpha_sweep only) text file with one urid per line; '
             'lines starting with # are ignored (alternative to --urids).',
    )
    parser.add_argument(
        '--param-combos',
        nargs = '+',
        default = ['gt5', 'scoremat_min5'],
        choices = ['gt5', 'gt20', 'scoremat_min5', 'scoremat_min10', 'scoremat_min20'],
        metavar = 'COMBO',
        help = '(alpha_sweep only) clade-definition schemes to run '
             '(default: gt scoremat_min5).',
    )
    parser.add_argument(
        '--n-restarts',
        type = int,
        default = 3,
        metavar = 'N',
        help = '(alpha_sweep only) EM random restarts per alpha value (default: 3).',
    )
    parser.add_argument(
        '--alpha-out',
        default = 'alpha_sweep_results.csv',
        metavar = 'PATH',
        help = '(alpha_sweep only) output CSV path (default: alpha_sweep_results.csv).',
    )
    args = parser.parse_args()

    # resolve urids - alpha_sweep also accepts --urid-file
    if args.urid_file:
        with open(args.urid_file) as fh:
            file_urids = [int(l.strip()) for l in fh
                           if l.strip() and not l.strip().startswith('#')]
        urids = file_urids
    else:
        urids = args.urids if args.urids is not None else long_mito_36_urids

    timepoints = args.timepoints if args.timepoints is not None else default_timepoints

    joblib_dir = args.joblib_dir
    os.makedirs(joblib_dir, exist_ok = True)
    joblib_prefix = args.joblib_prefix

    # scoremat_sweep mode: no trees, no EM
    if args.mode == 'scoremat_sweep':
        print(
            f'=== scoremat_sweep  urids={len(urids)}  '
            f'timepoints={timepoints}  min_clade_sizes={args.min_clade_sizes} ==='
        )
        sweep_df = scoremat_clade_sweep(
            urid_list = urids,
            min_clade_sizes = args.min_clade_sizes,
            timepoints = timepoints,
        )
        sweep_df.to_csv(args.output, index = False)
        print(f'Saved {len(sweep_df):,} rows to {args.output}')
        sys.exit(0)

    # alpha_sweep mode: parallel EM over alpha values
    if args.mode == 'alpha_sweep':

        alpha_out = Path(args.alpha_out)
        alpha_values = args.alpha_values
        param_combos = args.param_combos
        n_restarts = args.n_restarts

        print(
            f'=== alpha_sweep  urids={len(urids)}  tp={timepoints}  '
            f'alpha={alpha_values}  combos={param_combos}  '
            f'n_restarts={n_restarts}  workers={n_workers} ==='
        )

        # build task list by loading TMs and leaf maps in the main process
        # (tree loading happens inside each worker to avoid pickling ete3 trees)
        tasks = []
        skip_count = 0

        for urid in urids:
            try:
                Phi0, Phi1_true, founder, type_order = alpha_get_tms(str(urid))
            except Exception as _e:
                print(f'[alpha_sweep] Could not load TMs for {urid}: {_e}')
                skip_count += len(timepoints) * len(param_combos)
                continue

            type_encoder = {t: i for i, t in enumerate(type_order)}

            for tp in timepoints:
                try:
                    clade_df = alpha_load_clade_df(str(urid), tp)
                except FileNotFoundError as _e:
                    print(f'[alpha_sweep] clade_df not found {urid} tp={tp}: {_e}')
                    skip_count += len(param_combos)
                    continue

                for pc_key in param_combos:
                    if pc_key in ('gt5', 'gt20'):
                        param_combo_full = alpha_gt_combo[pc_key]
                        try:
                            tree_path = alpha_get_gt_tree_path(str(urid), tp)
                        except FileNotFoundError as _e:
                            print(f'[alpha_sweep] GT tree not found {urid} tp={tp}: {_e}')
                            skip_count += 1
                            continue
                    else:
                        param_combo_full = alpha_scoremat_combo.get(pc_key)
                        tree_path = alpha_get_recon_tree_path(str(urid), tp)
                        if tree_path is None:
                            print(f'[alpha_sweep] recon tree not found {urid} tp={tp} {pc_key}')
                            skip_count += 1
                            continue

                    l2c, l2t, l2igt = alpha_build_leaf_maps(clade_df, param_combo_full)
                    if l2c is None:
                        print(f'[alpha_sweep] column absent {urid} tp={tp} {pc_key}')
                        skip_count += 1
                        continue
                    if len(l2c) < 4:
                        print(f'[alpha_sweep] too few cells {urid} tp={tp} {pc_key}')
                        skip_count += 1
                        continue

                    tasks.append((
                        str(urid), tp, param_combo_full, tree_path,
                        l2c, l2t, l2igt,
                        Phi0, Phi1_true, alpha_values, type_encoder,
                        n_restarts, founder, 'leaf_prop',
                    ))

                del clade_df
                gc.collect()

        print(f'[alpha_sweep] {len(tasks)} tasks queued, {skip_count} skipped')

        n_rows = 0
        with ProcessPoolExecutor(max_workers = n_workers) as pool:
            futures = {pool.submit(worker_alpha_sweep, t): t for t in tasks}
            for fut in tqdm.tqdm(as_completed(futures), total = len(tasks), desc = 'alpha_sweep'):
                try:
                    rows = fut.result()
                except Exception as exc:
                    task = futures[fut]
                    print(f'ERROR for {task[:3]}: {exc}')
                    continue
                if rows:
                    chunk = pd.DataFrame(rows)
                    chunk.to_csv(
                        alpha_out, mode = 'a', index = False,
                        header = not alpha_out.exists(),
                    )
                    n_rows += len(rows)

        print(f'[alpha_sweep] done - {n_rows} rows written to {alpha_out}')
        sys.exit(0)

    # Phase 1: load & cluster all urids (serial)
    # Parallelism is not viable here: TimepointRes objects contain deep ete3 Tree
    # objects whose pickle recursion depth exceeds any practical sys limit, causing
    # _ForkingPickler to crash when the worker tries to send results back over the
    # process queue.  Phase 3 (EM inference) is where the real compute lives.
    print('=== Phase 1: process_urid (serial) ===')
    longmito_36_resdict = {}
    for urid in tqdm.tqdm(urids):
        path_to_output_dir = Path('./output')
        path_to_em_induction_prob = path_to_output_dir / 'em_induction_prob'
        timepoint_res_dict, _ = process_urid(
            urid = urid,
            param_dict_list = longmito_36_first_pass_list,
            cluster_vmin = -0.5,
            cluster_vmax = 1,
            timepoints = timepoints,
        )
        longmito_36_resdict[urid] = timepoint_res_dict

    # Phase 2: make_clade_df (fast; keep serial)
    print('=== Phase 2: make_clade_df ===')
    for urid, timepoint_res_dict in longmito_36_resdict.items():
        for tp, tpr in timepoint_res_dict.items():
            make_clade_df(timepoint_res = tpr, urid = urid, timepoint = tp)

    # Save clade assignments CSV if requested (works in any mode)
    clade_csv_path = args.clade_output or (args.output if args.mode == 'no_em' else None)
    if clade_csv_path:
        print('=== Building combined clade DataFrame ===')
        records = []
        for urid, tpr_dict in longmito_36_resdict.items():
            for tp, tpr in tpr_dict.items():
                for linstring, celldict in tpr.cell_pop.items():
                    for param_combo, clade in celldict.get('clade_assignments', {}).items():
                        records.append((urid, tp, param_combo, linstring, clade))
        combined_df = pd.DataFrame(
            records,
            columns = ['urid', 'timepoint', 'param_combo', 'linstring', 'clade']
        )
        combined_df.to_csv(clade_csv_path, index = False)
        print(f'Saved {len(combined_df):,} rows to {clade_csv_path}')

    if args.mode == 'no_em':
        sys.exit(0)

    # Phase 3: EM inference over all (urid, timepoint) pairs in parallel
    print('=== Phase 3: urid_to_induction_inference (parallel) ===')
    tasks = [
        (urid, tp)
        for urid in urids
        for tp in timepoints
    ]

    model_longmito_36_dict = {urid: {} for urid in urids}
    roc_longmito_36_dict = {urid: {} for urid in urids}
    completed = 0

    with ProcessPoolExecutor(max_workers = n_workers) as pool:
        futures = {pool.submit(worker_inference, task): task for task in tasks}
        for future in tqdm.tqdm(as_completed(futures), total = len(tasks)):
            try:
                urid, timepoint, tpt_induction_dict, opt_mod_res, eval_dict = future.result()
            except Exception as exc:
                task = futures[future]
                print(f'ERROR for {task}: {exc}')
                continue

            model_longmito_36_dict[urid][timepoint] = {
                'tpt_induction_dict': tpt_induction_dict,
                'opt_mod_res':        opt_mod_res,
                'eval_dict':          eval_dict,
            }
            roc_longmito_36_dict[urid][timepoint] = {
                combo: ed['metrics']['roc_auc'] for combo, ed in eval_dict.items()
            }
            roc_summary = {c: f'{r:.3f}' for c, r in roc_longmito_36_dict[urid][timepoint].items()}
            print(f'[EM done]  urid={urid} tp={timepoint} roc={roc_summary}', flush = True)

            completed += 1
            # checkpoint every 50 completed tasks
            if completed % 50 == 0:
                joblib.dump(model_longmito_36_dict, f'{joblib_dir}/{joblib_prefix}_model_dict_ckpt_{completed}.joblib')
                # with open(f'./no_gt_combos_model_longmito_36_dict_ckpt_{completed}.pkl', 'wb') as f:
                #     pickle.dump(model_longmito_36_dict, f)

    # Save final results
    joblib.dump(model_longmito_36_dict, f'{joblib_dir}/{joblib_prefix}_model_dict_FINAL.joblib')
    joblib.dump(roc_longmito_36_dict,   f'{joblib_dir}/{joblib_prefix}_roc_dict_FINAL.joblib')

    # Phase 4: plots (serial; needs all results)
    print('=== Phase 4: make_all_em_plots ===')
    plotting_dfs = make_all_em_plots(
        param_combo_list = gt_combos + mt_scoremat_clades,
        full_modelres_dict = model_longmito_36_dict,
        urid_timepoint_res_dict = longmito_36_resdict,
        return_plotting_dfs = True,
    )
