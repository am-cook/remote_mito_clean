#!/usr/bin/env python3
# longmito_36.py — auto-generated from mito_clonal_sub.ipynb cells 0, 1, 2, 29

import sys
import argparse
sys.setrecursionlimit(100000)

# ── Imports (cell 0) ──────────────────────────────────────────────────────────
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

# ── Definitions (cell 1) ─────────────────────────────────────────────────────


path_to_output_dir = Path('./output')
path_to_em_induction_prob = path_to_output_dir / 'em_induction_prob'


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
        self.gt_tree_path = path_to_gt_tree
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
            filt_dict, orient="index", columns=keep_cols
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
            col_vals = sorted(filt_df[col].unique(), key=str)
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
            gt_tree_copy = read_in_tree(self.gt_tree_path)
        else:
            matching = [sr for sr in self.scoremat_res_list if sr.short_name == param_combo]
            if not matching:
                raise ValueError(
                    f'No ScorematRes found with short_name == {param_combo!r}. '
                    f'Available: {[sr.short_name for sr in self.scoremat_res_list]}'
                )
            gt_tree_copy = read_in_tree(matching[0].recon_tree_path)


        # # filter out dead&terminal cells so that internal cells are retained but we don't get info-less leaves
        # all_nodes_to_keep = [cell for cell in gt_tree_copy.traverse('preorder') if not (self.cell_pop[cell.name]['terminal'] and not self.cell_pop[cell.name]['alive'])]
        # gt_tree_copy.prune(all_nodes_to_keep, preserve_branch_length = True)
        
        ################### why did i do this
        # only plot alive and terminal cells in the gt tree:
        # gt_tree_copy = copy.deepcopy(self.gt_tree)
        # alive_terminal_leaves = list(self.alive_terminal_pop.keys())
        # leaves_to_keep = [leaf for leaf in gt_tree_copy.iter_leaves() if leaf.name in alive_terminal_leaves]

        # non_leaves = [node for node in gt_tree_copy.traverse('preorder') if not node.is_leaf()]
        # keep_nodes = leaves_to_keep + non_leaves
        
        # # prune gt tree to only include terminal/alive leaves
        # # gt_tree_copy.prune(leaves_to_keep, preserve_branch_length = True)


        # ########
        keep_linstrings = prune_dead_lineages(cellpop_dict = self.cell_pop)
        keep_nodes = [node for node in gt_tree_copy.traverse('preorder') if node.name in keep_linstrings]
        gt_tree_copy.prune(keep_nodes, preserve_branch_length = True)
        # print(f'pre prune, {filt_df.shape=}')
        filt_df = filt_df[filt_df.index.isin(keep_linstrings)]
        # print(f'post prune, {filt_df.shape=}')
        
        # ########
        ################### why did i do this

        # if color internal nodes and descendants by some feature:
        if clade_colors is not None:
            if isinstance(clade_colors, dict):
                # pre-built linstring → hex color dict: pass through directly
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
                    key=str
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
            cell_state_filt_df=filt_df,
            group_to_color_dict=group_to_color_dict,
            inner_color_dict = alive_to_color_dict,
            group_to_ind=group_to_ind,
            meta_columns=list(color_by),
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
                value=(opt in default_set),   # True ⇢ pre‑checked
                description=opt,
                indent=False,
                layout = widgets.Layout(width = cb_width)
            )
            for opt in color_options
        ]
    
        # Put the checkboxes in a vertical box (you could also use HBox/GridBox)
        selector = widgets.VBox(checkboxes)
        out = widgets.Output()
    
        
        def currently_selected():
            return [cb.description for cb in checkboxes if cb.value]
    
        def refresh(change=None):
            selected = currently_selected()
            if not selected:
                with out:
                    out.clear_output()
                    print('At least one color must be given')
                return
    
            # Call the public method that builds the FigureWidget
            fig = self.new_plot_cass_tree(
                color_by=selected,
                include_alive_status_ring=include_alive_status_ring,
            )
            with out:
                out.clear_output()
                display(fig)

        # for each selected checkbox, refresh
        for cb in checkboxes:
            cb.observe(refresh, names= 'value')
    
        # show initial plot
        refresh()
    
        return widgets.VBox([selector, out])


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

    cass_tree = CassiopeiaTree(tree=G)



    if inner_color_dict is not None:
        if 'inner_color' not in cell_state_filt_df.columns:
            cell_state_filt_df['inner_color'] = cell_state_filt_df['alive'].map(
                lambda x: 'alive' if bool(x) else 'dead'
            )
    
    cass_tree.cell_meta = cell_state_filt_df
    
    custom_cmap = sns.color_palette('muted', n_colors=len(group_to_color_dict))

    # Drop constant-value columns — cassiopeia's continuous colorstrip crashes with min==max
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
    fig.update_layout(template='plotly_white')


    
        
    if inner_color_dict is None:
        if interactive:
            if return_html_string:
                html_str = pio.to_html(fig, include_plotlyjs='cdn')
                return html_str
            return fig
    
    
        static_bytes = pio.to_image(fig, format=image_format)

        if return_html_string:
            # Embed the raw bytes in a base‑64 data‑uri.
            b64 = base64.b64encode(static_bytes).decode()

            html_img = f'<img src="data:application/pdf;base64,{b64}" alt="Cassiopeia tree (static)"/>'
    
            return html_img
        else:
            # Return bytes directly – the caller can write them to disk or whatever.
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
            html_str = pio.to_html(fig, include_plotlyjs='cdn')
            return html_str
        return fig


    static_bytes = pio.to_image(fig, format='png')

    if return_html_string:
        # Embed the raw bytes in a base‑64 data‑uri.
        b64 = base64.b64encode(static_bytes).decode()

        html_img = f'<img src="data:application/pdf;base64,{b64}" alt="Cassiopeia tree (static)"/>'

        return html_img
    else:
        # Return bytes directly – the caller can write them to disk or whatever.
        return static_bytes
                                
    
        

class ScorematRes():
    '''
    stores results for a given score matrix
    '''

    def __init__(self, urid, timepoint, run_name, scoremat_path = None,  timepoint_res = None,
                tree = None, 
                 # cell_pop = None, alive_and_terminal_pop = None,
                recon_tree_path = None, gt_shortcut = False):
        '''
        initialze with run features and paths to results
        '''

        # abbreviated for gt results (the rest can be skipped):
        self.urid = urid
        self.timepoint = timepoint
        self.timepoint_res = timepoint_res

        # need dp version of gt tree even in shortcut
        # self.gt_tree_path = Path(f'/dartfs/rc/lab/M/McKennaLab/projects/Aidan/simulations/r_sim_clean/output/processed_newicks/{urid}/ground_truth_tree_{run_name}_time_{timepoint}.newick')
        self.gt_tree_path = path_to_output_dir / f'processed_newicks/{urid}/ground_truth_tree_{run_name}_time_{timepoint}.newick'
        self.gt_tree = read_in_tree(self.gt_tree_path)
        self.dp_gt_tree = dp.Tree.get(path = self.gt_tree_path,
                          schema = 'newick',
                         preserve_underscores=True)

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
    
            self.short_name = f'{self.mt_or_bc}_{self.num_ints}ints_{rec_frac}RP_{af}AF_{bin}B_{"-".join(cell_rec_fracs)}'
            
            
            # get matching recon tree path:
            # path_to_recon_tree_dir = Path(f'/dartfs/rc/lab/M/McKennaLab/projects/Aidan/simulations/r_sim_clean/output/recon_trees/{urid}/')
            path_to_recon_tree_dir = path_to_output_dir / f'recon_trees/{urid}/'
            scoremat_path_stem = Path(scoremat_path).stem

            try:
                self.recon_tree_path = path_to_recon_tree_dir / Path(f'score_{str(scoremat_path_stem)}') / Path(f'score_{str(scoremat_path_stem)}.treefile')
                self.recon_tree = read_in_tree(self.recon_tree_path)
                
                # temp solution to acclereate score calculation: also read in dendropy trees for gt and recon
                self.dp_recon_tree = dp.Tree.get(path = self.recon_tree_path,
                                  schema = 'newick',
                                 preserve_underscores=True)
            except: # if there are no trees associated with this scoremat, that's ok
                self.recon_tree_path = None
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
    
        
    
    # def plot_cass_tree(self, color_by, 
    #                cell_to_group_dict = None,
    #                group_to_color_dict = None,
    #                inner_clade_dict_label = None):
    #     '''
    #     plot the reconstructed tree, colored by some metadata feature

        
    #     color_by is a string. check to see if that attribute already exists in cellpop; if it does, plot by it
    #         if it does not, cell_to_group_dict will be added to cellpop under attribute name color_by
    #     group_to_color_dict is a mapping of groups that are in cell_to_group_dict.values() to plot colors
        
    #     '''

    #     tr = copy.deepcopy(self.recon_tree)
        
    #     # # check to see if color_by is in the cell-level dicts yet (quick & dirty check):
    #     # first_linstring, first_linstring_dict = next(iter(self.timepoint_res.cell_pop.items()))

    #     # # in case we want to color by multiple metadata factors
    #     # if isinstance(color_by, str):
    #     #     color_by = [color_by]
            
    #     # for this_color_by in color_by:
    #     #     if this_color_by not in first_linstring_dict.keys():
    #     #         for cellname, inner_dict in self.timepoint_res.cell_pop.items():
    #     #             inner_dict[this_color_by] = cell_to_group_dict.get(cellname, None)

    #     # if the attrs still have to be added to the cell pop
    #     # print('made it here 1')
    #     if cell_to_group_dict is not None:
    #         self.add_clade_dict_to_cellpop(cellpop_attr = color_by,
    #                                        cell_to_group_dict = cell_to_group_dict)
    #     else:
    #         # print('defining cell to group dict')
    #         cell_to_group_dict = {linstring:celldict.get(inner_clade_dict_label, None) for linstring, 
    #                               celldict in self.timepoint_res.cell_pop.items()}
    #     # print('made it here 2')

    #     # self.add_clade_dict_to_cellpop(cell_to_group_dict = )
    
    #     base_keep_cols = ['celltype', 'terminal', 'alive']


    #     if isinstance(color_by, str):
    #         color_by = [color_by]
    #     keep_cols = base_keep_cols + color_by

    #     # print(f'keep_cols = {keep_cols}')
    
    #     # this is the full cell pop data dataframe
    #     # cell_state_full_dict = {linstring: [inner_dict[col] for col in keep_cols] for linstring, inner_dict in self.timepoint_res.cell_pop.items()}
    #     cell_state_full_dict = {linstring: [inner_dict.get(col, None) for col in keep_cols] for linstring, inner_dict in self.timepoint_res.cell_pop.items()}
    #     cell_state_full_df = pd.DataFrame.from_dict(cell_state_full_dict, orient = 'index', columns = keep_cols)
        
    #     # this is the dataframe corresponding to only cells that are alive, which will be used to annotate the tree
    #     # cell_state_filt_dict = {linstring: [inner_dict[col] for col in keep_cols] for linstring, inner_dict in self.timepoint_res.cell_pop.items() if inner_dict['alive'] and inner_dict['terminal']}
    #     cell_state_filt_dict = {linstring: [inner_dict.get(col, None) for col in keep_cols] for linstring, inner_dict in self.timepoint_res.cell_pop.items() if inner_dict['alive'] and inner_dict['terminal']}
    #     cell_state_filt_df = pd.DataFrame.from_dict(cell_state_filt_dict, orient = 'index', columns = keep_cols)
        
    #     if group_to_color_dict is None:
    #         # sort groups by their frequency descending
    #         if cell_to_group_dict is None:
    #             # group_counts = cell_state_filt_df[color_by[0]].value_counts(ascending = False).to_dict()
    #             groups = cell_state_filt_df[color_by[0]]
    #             group_counts = dict(Counter(groups))
    #             group_counts = {k:v for k, v in sorted(group_counts.items(),
    #                                                    key = lambda x: x[1], reverse = True)}
    #             # print(f'in cell to group dict none, group_counts = {group_counts}')
                
    #         else:
    #             # print('in the else for cell to group dict')
                
    #             group_counts = {k:v for k, v in sorted(dict(Counter(cell_to_group_dict.values())).items(),
    #                                                    key = lambda x: x[1],
    #                                                    reverse = True)}
    #             # print(f'in the else, group_counts = {group_counts}')
                
    #         group_palette = sns.color_palette('muted', n_colors = len(group_counts))
    #         group_palette_hex = [rgb2hex(col) for col in group_palette]
    #         group_to_color_dict = dict(zip(group_counts.keys(), group_palette_hex))

    #     else:
    #         ex_clade, ex_col = next(iter(group_to_color_dict.items()))

    #         # convert rgb tuples to hex codes if necessary
    #         if isinstance(ex_col, tuple):
    #             group_to_color_dict = {k:rgb2hex(v) for k, v in group_to_color_dict.items()}
                
            
            

    #     # add color data to meta
    #     cell_state_filt_df['node_color'] = cell_state_filt_df.index.map(
    #         lambda linstring: group_to_color_dict.get(cell_to_group_dict.get(linstring, None), 'black')
    #     )
    #     # linstring_to_color_dict = {linstring: rgb2hex(group_to_color_dict.get(cell_to_group_dict.get(linstring, None), 'black')) for linstring in cell_state_filt_df.index}


    #     # print(f'group_to_color_dict = {group_to_color_dict}')

        
    #     custom_cmap = ListedColormap(group_to_color_dict.values())
    #     group_to_ind = {group: idx for idx, group in enumerate(group_to_color_dict.keys())}
        

        

    #     # print(f"cell_state_filt_df['node_color'].value_counts() = {cell_state_filt_df['node_color'].value_counts()}")
        
    #     G = nx.DiGraph()
    #     def rename_internal_node(n): 
    #         if n.name:
    #             return n.name
    #         else: # should not fire since linstring info is stored in .treefile
    #             return f'internal_{id(n)}'
        
    #     def add_edges(p):
    #         for c in p.children:
    #             G.add_edge(rename_internal_node(p), rename_internal_node(c))
    #             add_edges(c)
    #     add_edges(tr) 

       
    #     # create cassio tree and assign metadat
    #     cass_tree = CassiopeiaTree(tree = G)
    #     cass_tree.cell_meta = cell_state_filt_df
        
        
        
    #     # node_colors = []
    #     # for n in G.nodes():
    #     #     node_colors.append(group_to_color_dict.get(cell_to_group_dict.get(n, None), 'black'))
        
    #     fig = plot_plotly(cass_tree, meta_data = color_by, 
    #                       categorical_cmap = custom_cmap,
    #                      value_mapping = group_to_ind)
    #     # fig.data[0].marker.color = node_colors

    #     fig.update_layout(template = 'plotly_white')
    #     fig.show()
    
    
    
            
    
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
        # if (num_clades is None and not use_median) or (num_clades is not None and use_median):
        #     raise Exception('exactly one of num_clades and use_median must be used')
        
        # if use_median:
        #     med_dist = self.plot_dists_from_root(use_branch_length = use_branch_length, return_median = True,
        #                                      show_fig = show_intermediate_figs)
        #     varying = 'dist'
        #     cats = [med_dist]
            
        # else:
        #     varying = 'num_clades'
        #     cats = [num_clades]
        #     med_dist = None
            
        # if show_intermediate_figs:
        #     self.plot_clade_distribution_heatmap(cats = [med_dist], varying = 'dist', use_branch_length = use_branch_length)
            
        # cell_clade_dict = self.assign_clades(num_clades = num_clades, dist_cutoff = med_dist, 
        #                                                        use_branch_length = use_branch_length)
        
        # clade_to_color_dict = self.make_clade_clustermap(cell_clade_dict = cell_clade_dict,
        #                                                  order_by_clade = order_by_clade)

        # print(f'Counter(cell_clade_dict.values()) == {Counter(cell_clade_dict.values())}')
        # hex_dict = {k: rgb2hex(v) for k, v in clade_to_color_dict.items()}
        # print(f'clade_to_color_dict = {hex_dict}')
        # self.plot_cass_tree(color_by = cellpop_clade_attr_name, 
        #            cell_to_group_dict = cell_clade_dict, 
        #             group_to_color_dict = clade_to_color_dict)
        
        # return cell_clade_dict

    
        
    

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
            for clade_id, sig in enumerate(center_sigs, start=1):
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
    
    
                ########################### new part
                # if there are no existing clades yet (first sig, if below min_clade_size),
                # manually create the first clade (which will be the biggest based on how the sorting was performed) 
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
    
                ########################### new part
                    
    
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

        # Map node object id → clade label (avoids add_feature / tree mutation)
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
        # Suppose `mut_mat` is a DataFrame (cells × mutation sites) of 0/1
        dist_arr = squareform(pdist(self.scoremat.values, metric='hamming'))
        dist_df = pd.DataFrame(dist_arr, index=self.scoremat.index, columns= self.scoremat.index)
    
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
                                 preserve_underscores=True)
            
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
        
        dist_df = pd.DataFrame(dist_mat, index=leaf_labels, columns=leaf_labels)
        
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
            #### if using tree_dist metric, use 
            # dist_df = self.get_tree_dist_df(use_branch_length = use_branch_length, tree = tree)
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




# def assign_gt_hierarchical_clades(self, urid, timepoint, timepoint_res):
#     # create temporary scoremat res object so that assign_hierarchical_clades() can be called on GT data
#     gt_scoremat = ScorematRes(urid = 

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

def make_clade_df(timepoint_res, urid, timepoint):
    # path_to_clade_dfs = Path(f'/dartfs/rc/lab/M/McKennaLab/projects/Aidan/simulations/r_sim_clean/output/clade_assignment_dfs/{urid}/{timepoint}')
    path_to_clade_dfs = path_to_output_dir / f'clade_assignment_dfs/{urid}/{timepoint}'
    os.makedirs(path_to_clade_dfs, exist_ok = True)

    cell_pop = timepoint_res.cell_pop
    dict_of_dicts = {cell:celldict['clade_assignments'] for cell, celldict in cell_pop.items()}
    # removal of cells not assigned to clades:
    dict_of_dicts = {k:v for k, v in dict_of_dicts.items() if len(v) > 0}
    clade_df = pd.DataFrame.from_dict(dict_of_dicts, orient = 'index')

    clade_df.to_csv(os.path.join(path_to_clade_dfs, f'clade_df_{urid}_t{timepoint}.csv'))

    return clade_df


def _scoremat_assign(scoremat, min_clade_size):
    """
    Pure function: assign scoremat clades from a binary allele DataFrame.

    Mirrors ScorematRes.assign_scoremat_clades() but requires no cell_pop,
    no trees, and no ScorematRes wrapper — takes the matrix directly.

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
    row_sig     = scoremat.apply(lambda r: tuple(r.astype(int)), axis=1)
    clone_groups = row_sig.groupby(row_sig).groups
    sig_counts  = {sig: len(cells) for sig, cells in clone_groups.items()}
    sorted_sigs = sorted(clone_groups, key=lambda s: sig_counts[s], reverse=True)

    cell_to_clade = {}
    clade_to_sig  = {}
    clade_freqs   = {}
    next_id       = 1

    for sig in sorted_sigs:
        cells = clone_groups[sig]
        n     = sig_counts[sig]

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
    timepoints=None,
    base_dir=path_to_output_dir,
):
    """
    Assign scoremat clades for every combination of:
        urid × timepoint × scoremat CSV (AF/RP) × min_clade_size

    Does not load any trees, does not run EM, and does not require
    process_urid() or TimepointRes. Works directly from the binary
    allele CSVs in output/{urid}/score_mats/{urid}/matrices/csvs/.

    Parameters
    ----------
    urid_list       : list of int or str
    min_clade_sizes : list of int  e.g. [5, 10, 20]
    timepoints      : list of int, optional  — None = all available
    base_dir        : Path  — root of per-urid output dirs (default path_to_output_dir)

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
        csv_dir  = Path(base_dir) / urid_str / 'score_mats' / urid_str / 'matrices' / 'csvs'

        if not csv_dir.exists():
            print(f'[scoremat_clade_sweep] No CSV dir found for {urid_str}, skipping')
            continue

        csv_paths = sorted(csv_dir.glob('*.csv'))
        if timepoints is not None:
            csv_paths = [p for p in csv_paths
                         if extract_timepoint(p.stem) in timepoints]

        for csv_path in csv_paths:
            stem = csv_path.stem
            tp   = extract_timepoint(stem)

            try:
                rp    = extract_frac(path=stem, joint=False)
                af    = extract_af(path=stem, joint=False)
                ints  = extract_ints(path=stem, joint=False)
                bin_  = extract_bin(path=stem, joint=False)
                fracs = extract_cell_rec_fracs(stem)
            except Exception as e:
                print(f'[scoremat_clade_sweep] Could not parse {stem}: {e}')
                continue

            tb_str   = '-'.join(fracs)
            scoremat = pd.read_csv(csv_path, index_col=0)

            for min_size in min_clade_sizes:
                param_combo = (
                    f'scoremat_min{min_size}_mt_{ints}ints_'
                    f'{rp}RP_{af}AF_{bin_}B_{tb_str}'
                )
                cell_to_clade = _scoremat_assign(scoremat, min_clade_size=min_size)

                for linstring, clade in cell_to_clade.items():
                    records.append((urid_str, tp, param_combo, linstring, clade))

    return pd.DataFrame(
        records,
        columns=['urid', 'timepoint', 'param_combo', 'linstring', 'clade']
    )


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


def get_candidate_cols(urid, param_dict_list, timepoint):
    tpt_res_dict, clade_name_to_param_combo_dict = process_urid(urid,
                                                                param_dict_list = param_dict_list,
                                                                cluster_vmin = -0.5,
                                                                cluster_vmax = 1)
    # print(f'clade_name_to_param_combo_dict = {clade_name_to_param_combo_dict}')
    timepoint_res = tpt_res_dict[timepoint]
    # candidate_cols = list(set([clade for cell, celldict in timepoint_res.cell_pop.items() for clade in celldict['clade_assignments'].keys()]))

    candidate_cols = sorted(list(clade_name_to_param_combo_dict.keys()))

    # print(f'in python, candidate_cols = {candidate_cols}')
    return timepoint_res, candidate_cols, clade_name_to_param_combo_dict


def urid_to_interactive(urid, param_dict_list, timepoint, selected_cols, timepoint_res):

    
    fig = timepoint_res.new_plot_cass_tree(
                color_by=selected_cols,
                include_alive_status_ring=False,
    )

    return fig

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
    
    
def get_ground_truth_inductions(urid, induction_type = 'differentiation'):

    if induction_type not in ['differentiation', 'editing']:
        raise ValueError('induction_type must be one of "differentiation" or "editing"')
    
    # induction_details_path = Path(f'/dartfs/rc/lab/M/McKennaLab/projects/Aidan/simulations/r_sim_clean/output/induction_details/{urid}/{induction_type}_induction_df.csv')
    induction_details_path = path_to_output_dir / f'induction_details/{urid}/{induction_type}_induction_df.csv'

    induction_df = pd.read_csv(induction_details_path, index_col = 0)

    return induction_df
    

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
    
    # path_to_clade_dfs = Path(f'/dartfs/rc/lab/M/McKennaLab/projects/Aidan/simulations/r_sim_clean/output/clade_assignment_dfs/{urid}/{timepoint}')
    path_to_clade_dfs = path_to_output_dir / f'clade_assignment_dfs/{urid}/{timepoint}'
    
    ex_at_cell_pop = urid_timepoint_res_dict[urid][timepoint].alive_terminal_pop
    
    ex_clade_df = pd.read_csv(path_to_clade_dfs / f'clade_df_{urid}_t{timepoint}.csv', index_col = 0).reset_index(names = 'cell')
    # ex_clade_df.rename(columns = {'index': 'cell'}, inplace = True)
    
    melted_clade_df = pd.melt(ex_clade_df, var_name = 'param_combo',
                             value_name = 'clade', ignore_index = False,
                             id_vars = 'cell')
    
    celltype_df = pd.DataFrame.from_dict({cell:celldict['celltype'] for cell, celldict in ex_at_cell_pop.items()},
                                          orient = 'index', columns = ['celltype']).reset_index(names = 'cell')
    
    full_df = pd.merge(melted_clade_df, celltype_df, on = 'cell', how = 'inner')

    return full_df


def plot_counts_for_param_combo(cell_counts_df, param_combo, urid, savename = ''):
    
    # path_to_celltype_bars = Path(f'/dartfs/rc/lab/M/McKennaLab/projects/Aidan/simulations/r_sim_clean/output/celltype_bars/{urid}')
    path_to_celltype_bars = path_to_output_dir / f'/celltype_bars/{urid}'
    os.makedirs(path_to_celltype_bars, exist_ok = True)
    
    combo_df = cell_counts_df[cell_counts_df['param_combo'] == param_combo]
    plt.figure(figsize = (10, 5))
    ax = sns.barplot(data = combo_df,
                     x = 'clade',
                     y = 'cell_count', 
                     hue = 'celltype', 
                     palette = 'muted')
    ax.set_xlabel('Clade')
    ax.set_ylabel('Number of cells')
    ax.set_title(f'Clade param combo = {param_combo}')

    ax.legend(title = 'Cell Type', loc = 'upper left', bbox_to_anchor = (1,1))
    plt.tight_layout()
    plt.savefig(path_to_celltype_bars / f'celltype_bars{savename}.pdf')
    plt.close('all')


def plot_celltypes_across_clades(urid, timepoint, urid_timepoint_res_dict, path_to_output_dir = path_to_output_dir):

    

    ct_clade_df = make_clade_celltype_df(urid = urid, 
                                        timepoint = timepoint,
                                        urid_timepoint_res_dict = urid_timepoint_res_dict,
                                        path_to_output_dir = path_to_output_dir)

    cellcount_df = ct_clade_df.groupby(['param_combo', 'clade', 'celltype'])['cell'].nunique().reset_index(name = 'cell_count')
    
    for param_combo in cellcount_df['param_combo'].unique():
        plot_counts_for_param_combo(cell_counts_df = cellcount_df, 
                                    param_combo = param_combo,
                                    urid = urid)


def plot_induced_prob_by_clade(df, urid, timepoint, include_cols = None, savename = ''):

    '''
    only retain param combos that contain substring(s) given in include_cols
    '''

    
    # path_to_induction_probs = Path(f'/dartfs/rc/lab/M/McKennaLab/projects/Aidan/simulations/r_sim_clean/output/em_induction_prob/{urid}')
    path_to_induction_probs = path_to_output_dir / f'em_induction_prob/{urid}'
    os.makedirs(path_to_induction_probs, exist_ok = True)
    
    tp_df = df[df['timepoint'] == timepoint]

    if include_cols is not None:
        retain_cols = list(set([col for col in tp_df['clade_name'].unique() for allowed_combo in include_cols if 
                                allowed_combo in col]))
        tp_df = tp_df[tp_df['clade_name'].isin(retain_cols)]

    plt.figure(figsize = (30, 10))
    # plt.figure()
    # ax = sns.scatterplot(data = tp_df,
    #                      x = 'clade_name',
    #                      y = 'p_induced',
    #                      hue = 'clade')

    ax = sns.barplot(data = tp_df,
                         x = 'clade_name',
                         y = 'p_induced',
                         hue = 'clade')
    ax.axhline(y = 0.5, color = 'red', linestyle = '--')
    ax.set_ylabel('Probability induced')
    ax.set_xlabel('Clade name')
    ax.set_title(f'{urid}: cladewise induction probability at time {timepoint}')
    ax.set_ylim(0, 1)
    

    plt.xticks(rotation = 45, ha = 'right')
    # ax.legend(title = 'Clade', bbox_to_anchor = (1.02, 1))
    if ax.get_legend():
        ax.get_legend().remove()
    plt.tight_layout()
    
    plt.savefig(path_to_induction_probs / f'induction_probs{savename}.pdf')
    plt.close('all')


# class TreeMixtureEM():
#     """
#     Maximum‑likelihood mixture EM on a rooted, time‑scaled tree.

#     Parameters
#     ----------
#     newick_path : str
#         Path to the Newick file that contains the lineage tree.
#     df_map : pd.DataFrame
#         Must contain columns ['cell_name', 'clade', 'cell_type'].
#     Q0, Q1 : np.ndarray, shape (K, K)
#         Un‑induced and induced CTMC rate matrices.
#     root_type : int or str
#         Known type of the founder cell (must match coding used for cell_type).
#     type_encoder : optional, default None
#         If cell_type are strings, supply a dict mapping type → integer (0‑based).
#         If None, the class will infer a mapping from the observed types.
#     """

#     def __init__(self,
#                  Q0,
#                  Q1,
#                  root_type,
#                  tree,
#                  leaf2clade,
#                  leaf2type,
#                  type_encoder=None):

#         # -------------------------------------------------
#         # 1️⃣  Load tree and attach leaf metadata
#         # -------------------------------------------------
#         # self.tree = Tree(newick_path, format=1)   # format=1 = with branch lengths
#         self.tree = tree
#         self._check_tree_is_rooted()

#         # # create a dictionary: leaf_name → (clade, observed_type)
#         # self.leaf2clade = {}
#         # self.leaf2type  = {}
#         self.leaf2clade = leaf2clade
#         self.leaf2type = leaf2type
#         # for _, row in df_map.iterrows():
#         #     name = str(row['cell_name'])
#         #     self.leaf2clade[name] = str(row['clade'])
#         #     self.leaf2type[name]  = row['cell_type']

#         # make sure every leaf in the tree has a mapping
#         leaf_names = {leaf.name for leaf in self.tree.iter_leaves()}
#         missing = leaf_names - set(self.leaf2clade.keys())
#         if missing:
#             raise ValueError(f"The following leaves are missing from df_map: {missing}")

#         # -------------------------------------------------
#         # 2️⃣  Encode cell types as integers 0 … K‑1
#         # -------------------------------------------------
#         if type_encoder is None:
#             uniq = sorted(set(self.leaf2type.values()))
#             self.type_encoder = {t: i for i, t in enumerate(uniq)}
#         else:
#             self.type_encoder = type_encoder

#         # convert to integers
#         self.leaf2type_int = {
#             leaf: self.type_encoder[self.leaf2type[leaf]]
#             for leaf in leaf_names
#         }

#         # store K
#         Q0 = np.asarray(Q0)
#         Q1 = np.asarray(Q1)
#         self.K = Q0.shape[0]
#         assert Q0.shape == (self.K, self.K) and Q1.shape == (self.K, self.K)

#         # -------------------------------------------------
#         # 3️⃣  Store rate matrices and pre‑compute branch‑length matrices
#         # -------------------------------------------------
#         self.Q0 = Q0.astype(float)
#         self.Q1 = Q1.astype(float)

#         # pre‑compute transition matrices for each *unique* branch length
#         self._precompute_branch_matrices()

#         # -------------------------------------------------
#         # 4️⃣  Root prior (founder type is known)
#         # -------------------------------------------------
#         if isinstance(root_type, str):
#             root_type = self.type_encoder[root_type]
#         self.root_type = int(root_type)

#         # -------------------------------------------------
#         # 5️⃣  Initialise clade‑specific mixture parameters
#         # -------------------------------------------------
#         self.clades = sorted(set(self.leaf2clade.values()))
#         # π_c = prior probability that clade c is induced (binary indicator)
#         self.pi = {c: 0.5 for c in self.clades}       # start at 0.5 for all clades
#         # for the *fractional* mixture version you could initialise p_c similarly
#         self.p_frac = {c: 0.5 for c in self.clades}   # optional, see later

#         # cache for fast node‑to‑clade lookup
#         self._assign_node_clades()

#     # ------------------------------------------------------------------
#     #  Helper methods
#     # ------------------------------------------------------------------
#     def _check_tree_is_rooted(self):
#         # if self.tree.is_rooted:
#         #     return
#         if self.tree.get_tree_root():
#             return
            
#         raise ValueError("The supplied tree must be rooted (Newick file should contain a root).")

#     def _precompute_branch_matrices(self):
#         """
#         Build dictionaries:
#             self.P0[length] = expm(Q0 * length)
#             self.P1[length] = expm(Q1 * length)
#         We store them for *exact* floating point values of branch length.
#         """
#         unique_lengths = {node.dist for node in self.tree.traverse()}
#         self.P0 = {}
#         self.P1 = {}
#         for L in unique_lengths:
#             self.P0[L] = expm(self.Q0 * L)
#             self.P1[L] = expm(self.Q1 * L)

#     def _assign_node_clades(self):
#         """
#         Attach a `.clade_set` attribute to every node: the *set* of clades that
#         have at least one descendant leaf under that node.
#         """
#         for node in self.tree.traverse("postorder"):
#             if node.is_leaf():
#                 node.clade_set = {self.leaf2clade[node.name]}
#             else:
#                 node.clade_set = set()
#                 for child in node.children:
#                     node.clade_set.update(child.clade_set)

#     # ------------------------------------------------------------------
#     #  Likelihood calculation (core of the EM)
#     # ------------------------------------------------------------------
#     def _branch_transition_matrix(self, node, clade):
#         """
#         Return the appropriate K×K transition matrix for the edge leading
#         *into* `node`, given that the edge belongs to `clade`.

#         The mixture model we use here is the **fractional** version:
#             P = (1 - p_c) * P0 + p_c * P1
#         where p_c = self.p_frac[clade].

#         If you prefer the binary indicator version (all-or‑none), replace
#         the mixture by either P0 or P1 depending on the current estimate of
#         the binary state.
#         """
#         length = node.dist
#         # fetch pre‑computed matrices
#         P0 = self.P0[length]
#         P1 = self.P1[length]

#         # mixture weight for this clade (fraction of cells that are induced)
#         p = self.p_frac[clade]          # 0 ≤ p ≤ 1
#         return (1.0 - p) * P0 + p * P1

#     def _pruning_likelihood(self):
#         """
#         Perform a single Felsenstein pruning pass on the *whole* tree,
#         using the current `p_frac` values for each clade.
#         Returns the total log‑likelihood (float).
#         Also stores the per‑node conditional likelihood vectors in `node.lik`
#         for possible downstream use (e.g., posterior decoding).
#         """
#         K = self.K
#         # bottom‑up traversal
#         for node in self.tree.traverse("postorder"):
#             if node.is_leaf():
#                 lk = np.zeros(K)
#                 lk[self.leaf2type_int[node.name]] = 1.0
#                 node.lik = lk
#             else:
#                 # start with a vector of ones (neutral element for element‑wise product)
#                 lk = np.ones(K)
#                 for child in node.children:
#                     # which clade(s) does this child belong to?
#                     # In the simplest formulation we assume *one* clade per leaf,
#                     # and that the whole branch inherits the clade of its descendant leaves.
#                     # If the child has multiple clades under it, we take a *weighted*
#                     # average of the mixture matrices (the weights are proportional
#                     # to the number of descendant leaves per clade).
#                     clade_weights = defaultdict(float)
#                     for leaf in child.iter_leaves():
#                         clade_weights[self.leaf2clade[leaf.name]] += 1.0
#                     total = sum(clade_weights.values())
#                     # build the effective transition matrix for this edge
#                     P_eff = np.zeros((K, K))
#                     for c, w in clade_weights.items():
#                         P_eff += (w/total) * self._branch_transition_matrix(child, c)

#                     # multiply conditional likelihoods
#                     lk = lk * (P_eff @ child.lik)
#                 node.lik = lk

#         # root prior: a degenerate distribution that forces the founder type
#         root = self.tree.get_tree_root()
#         root_prior = np.zeros(K)
#         root_prior[self.root_type] = 1.0
#         logL = np.log(np.dot(root_prior, root.lik))
#         return logL

#     # ------------------------------------------------------------------
#     #  EM algorithm
#     # ------------------------------------------------------------------
#     def fit(self,
#             max_iter=100,
#             tol=1e-6,
#             verbose=True):
#         """
#         Run the EM algorithm to estimate clade‑specific mixture fractions
#         (p_frac) and the total log‑likelihood.

#         Returns
#         -------
#         logL_history : list of log‑likelihood values (one per iteration)
#         """
#         logL_history = []
#         for it in tqdm.trange(max_iter, disable=not verbose):
#             # ---------- E‑step ----------
#             # Here the “complete data” consists of the hidden binary indicator
#             # Z_c for each clade.  Because we treat the mixture fraction p_c as a
#             # *continuous* latent variable, the posterior expectation of Z_c is
#             # simply p_c itself (the current estimate).
#             # If you wanted a true binary indicator, you would compute:
#             #   γ_c = Pr(Z_c = 1 | data) = p_c * L1_c / [p_c * L1_c + (1-p_c) * L0_c]
#             # where L1_c and L0_c are the clade‑specific likelihoods under Q1
#             # and Q0, respectively.  For simplicity we keep the fractional version.
#             # ------------------------------------
#             # Compute the log‑likelihood under the *current* parameters
#             logL = self._pruning_likelihood()
#             logL_history.append(logL)

#             if verbose:
#                 logging.info(f"EM iteration {it+1:02d} – logL = {logL:.6f}")

#             # ---------- M‑step ----------
#             # Update p_c for each clade by maximizing the expected complete‑data
#             # log‑likelihood.  In the fractional model the update has a closed‑form:
#             #
#             #        p_c_new = ( Σ_{branches in clade c}  w_branch * (partial likelihood of Q1) )
#             #                  ---------------------------------------------------------------
#             #                     Σ_{branches in clade c}  w_branch
#             #
#             # where w_branch is the posterior probability that a cell on that
#             # branch used the induced matrix.  A cheap (and often sufficient)
#             # approximation is to set p_c_new equal to the *average* of the
#             # posterior responsibilities that we can obtain from the E‑step.
#             #
#             # Below we compute a *soft* responsibility for each branch by
#             # evaluating the likelihood of that branch under Q0 and Q1 separately.
#             # -----------------------------------------------------------------
#             # First, collect per‑clade sufficient statistics.
#             #
#             clade_num = defaultdict(float)   # denominator Σ w_branch
#             clade_num_ind = defaultdict(float)  # numerator Σ w_branch * responsibility

#             for node in self.tree.traverse("postorder"):
#                 if node.is_root():
#                     continue
#                 # identify the clade(s) this edge belongs to (same logic as in pruning)
#                 clade_weights = defaultdict(float)
#                 for leaf in node.iter_leaves():
#                     clade_weights[self.leaf2clade[leaf.name]] += 1.0
#                 total_leafs = sum(clade_weights.values())

#                 # compute branch‑specific likelihoods under both matrices
#                 #   L0 = parent_lik * P0 * child.lik
#                 #   L1 = parent_lik * P1 * child.lik
#                 parent = node.up
#                 # get conditional vectors from previous pruning pass
#                 parent_lik = parent.lik
#                 child_lik  = node.lik

#                 # retrieve pre‑computed P0 and P1 for this edge length
#                 P0 = self.P0[node.dist]
#                 P1 = self.P1[node.dist]

#                 # compute the *joint* likelihood contributions (scalar)
#                 #   L = π_parent^T * P * child_lik
#                 #   where π_parent = parent.lik (already conditional on the subtree)
#                 L0 = np.dot(parent_lik, P0 @ child_lik)
#                 L1 = np.dot(parent_lik, P1 @ child_lik)

#                 # responsibility = probability that this *branch* used Q1
#                 #    = L1 / (L0 + L1)   (Bayes rule, assuming equal prior)
#                 # To avoid division‑by‑zero when both are tiny, add a tiny epsilon.
#                 eps = 1e-12
#                 resp = L1 / (L0 + L1 + eps)

#                 # distribute this responsibility to the clades proportionally
#                 for c, w in clade_weights.items():
#                     weight = w / total_leafs
#                     clade_num[c] += weight                     # denominator
#                     clade_num_ind[c] += weight * resp           # numerator

#             # finally update p_frac
#             for c in self.clades:
#                 if clade_num[c] > 0:
#                     new_p = clade_num_ind[c] / clade_num[c]
#                     # keep values inside (0,1) to avoid numerical issues
#                     self.p_frac[c] = np.clip(new_p, 1e-5, 1 - 1e-5)

            

#         self.logL_history_ = logL_history
#         return logL_history

#     # ------------------------------------------------------------------
#     #  Post‑processing utilities
#     # ------------------------------------------------------------------
#     def clade_posteriors(self):
#         """
#         Return a DataFrame with one row per clade:
#             clade | p_induced | hard_label
#         `p_induced` is the EM‑estimated mixture fraction (probability that
#         the clade uses the induced matrix).  `hard_label` is “induced” if
#         p_induced > 0.5, otherwise “un‑induced”.
#         """
#         rows = []
#         for c in self.clades:
#             p = self.p_frac[c]
#             rows.append({'clade': c,
#                          'p_induced': p,
#                          'hard_label': 'induced' if p > 0.5 else 'un-induced'})
#         return pd.DataFrame(rows)

#     def summary(self):
#         """Print a quick textual summary."""
#         df = self.clade_posteriors()
#         print(df.to_string(index=False))
#         print(f"\nFinal log‑likelihood = {self.logL_history_[-1]:.6f}")



import logging
from collections import defaultdict
from typing import Dict, Any, Optional

import numpy as np
import pandas as pd
import tqdm

# ----------------------------------------------------------------------
# Helper: sanity‑check that a matrix is stochastic (rows sum to 1, all ≥0)
# ----------------------------------------------------------------------
def _check_stochastic(P: np.ndarray, name: str) -> None:
    if not np.allclose(P.sum(axis=1), 1.0, atol=1e-9):
        raise ValueError(f"{name} rows must sum to 1 (found {P.sum(axis=1)})")
    if np.any(P < -1e-12):
        raise ValueError(f"{name} contains negative entries")


# # ----------------------------------------------------------------------
# # Main class
# # ----------------------------------------------------------------------
# class TreeMixtureNodeEM():
#     """
#     EM for a **node‑only** mixture of two K‑state stochastic matrices
#     on a rooted phylogeny.

#     Parameters
#     ----------
#     Phi0, Phi1 : np.ndarray, shape (K, K)
#         Two *node‑transition* probability matrices.
#         Phi0 – “un‑induced” process,
#         Phi1 – “induced” process.
#         Rows must sum to 1 and contain only non‑negative entries.
#     root_type : int or str
#         Known founder cell type (must match the integer encoding of cell types).
#     tree : ete3.Tree (or any object that behaves like it)
#         The rooted lineage tree (branch lengths are ignored by the model).
#     leaf2clade : dict[str, str]
#         Mapping leaf name → clade label.
#     leaf2type : dict[str, Any]
#         Mapping leaf name → observed cell‑type label (string or int).
#     type_encoder : dict[Any, int] , optional
#         If cell types are strings, supply a mapping to integers 0…K‑1.
#         If None the encoder is inferred from ``leaf2type``.
#     """

#     # ------------------------------------------------------------------
#     #  Constructor
#     # ------------------------------------------------------------------
#     def __init__(self,
#                  Phi0: np.ndarray,
#                  Phi1: np.ndarray,
#                  root_type,
#                  tree,
#                  leaf2clade: Dict[str, str],
#                  leaf2type: Dict[str, Any],
#                  type_encoder: Optional[Dict[Any, int]] = None,
#                 weight_strategy = 'leaf_prop'):
#         # --------------------------------------------------------------
#         # 1️⃣  Store tree and attach leaf metadata
#         # --------------------------------------------------------------
#         self.tree = tree
#         self._check_tree_is_rooted()

#         self.leaf2clade = leaf2clade
#         self.leaf2type  = leaf2type

#         self.weight_strategy = weight_strategy

#         leaf_names = {leaf.name for leaf in self.tree.iter_leaves()}
#         missing = leaf_names - set(self.leaf2clade.keys())
#         if missing:
#             raise ValueError(f"The following leaves are missing from leaf2clade: {missing}")

#         # --------------------------------------------------------------
#         # 2️⃣  Encode cell types as integers 0 … K‑1
#         # --------------------------------------------------------------
#         if type_encoder is None:
#             uniq = sorted(set(self.leaf2type.values()))
#             self.type_encoder = {t: i for i, t in enumerate(uniq)}
#         else:
#             self.type_encoder = type_encoder

#         self.leaf2type_int = {
#             leaf: self.type_encoder[self.leaf2type[leaf]]
#             for leaf in leaf_names
#         }

#         # --------------------------------------------------------------
#         # 3️⃣  Store the two node‑transition matrices (must be stochastic)
#         # --------------------------------------------------------------
#         self.Phi0 = np.asarray(Phi0, dtype=float)
#         self.Phi1 = np.asarray(Phi1, dtype=float)

#         self.K = self.Phi0.shape[0]
#         assert self.Phi0.shape == (self.K, self.K)
#         assert self.Phi1.shape == (self.K, self.K)

#         _check_stochastic(self.Phi0, "Phi0")
#         _check_stochastic(self.Phi1, "Phi1")

#         # --------------------------------------------------------------
#         # 4️⃣  Root prior (founder type is known)
#         # --------------------------------------------------------------
#         if isinstance(root_type, str):
#             root_type = self.type_encoder[root_type]
#         self.root_type = int(root_type)

#         # --------------------------------------------------------------
#         # 5️⃣  Initialise clade‑specific mixture parameters
#         # --------------------------------------------------------------
#         self.clades = sorted(set(self.leaf2clade.values()))
#         self.p_frac = {c: 0.5 for c in self.clades}   # start at 0.5 for every clade

#         # cache for fast node‑to‑clade look‑up (same as in the original code)
#         self._assign_node_clades()

#     # ------------------------------------------------------------------
#     #  Helper methods
#     # ------------------------------------------------------------------
#     def _check_tree_is_rooted(self) -> None:
#         """Raise if the supplied tree has no root."""
#         if not self.tree.get_tree_root():
#             raise ValueError("The supplied tree must be rooted.")

#     def _assign_node_clades(self) -> None:
#         """
#         Attach a ``clade_set`` attribute to every node:
#         the set of clades that have at least one descendant leaf under that node.
#         """
#         for node in self.tree.traverse("postorder"):
#             if node.is_leaf():
#                 node.clade_set = {self.leaf2clade[node.name]}
#             else:
#                 node.clade_set = set()
#                 for child in node.children:
#                     node.clade_set.update(child.clade_set)

#     def _effective_node_matrix(self, node) -> np.ndarray:
#         """
#         Return the convex combination
#                 Φ_node = (1‑p_c) Φ0 + p_c Φ1
#         where the weight p_c is a leaf‑proportion weighted average of the
#         clades that appear under ``node`` (exactly the same weighting that
#         the original edge‑wise code used).
#         """
#         # 1) count how many descendant leaves belong to each clade
#         clade_weights = defaultdict(float)
#         for leaf in node.iter_leaves():
#             clade_weights[self.leaf2clade[leaf.name]] += 1.0
#         total = sum(clade_weights.values())

#         # 2) build the weighted mixture
#         Phi_eff = np.zeros((self.K, self.K))
#         for c, w in clade_weights.items():
#             p = self.p_frac[c]                # mixture weight for clade c
#             Phi_eff += (w / total) * ((1.0 - p) * self.Phi0 + p * self.Phi1)
#         return Phi_eff

#     # ------------------------------------------------------------------
#     #  Likelihood calculation (core of the EM)
#     # ------------------------------------------------------------------
#     def _pruning_likelihood(self) -> float:
#         """
#         Single Felsenstein pruning pass under the *node‑only* model.
#         Stores a K‑dimensional conditional likelihood vector ``node.lik``
#         at every node and returns the total log‑likelihood.
#         """
#         K = self.K

#         # Bottom‑up traversal
#         for node in self.tree.traverse("postorder"):
#             if node.is_leaf():
#                 # leaf emission = one‑hot vector
#                 lk = np.zeros(K)
#                 lk[self.leaf2type_int[node.name]] = 1.0
#                 node.lik = lk
#             else:
#                 # children (binary tree assumed; if more children, the code still works)
#                 child_vec = np.ones(K)               # product of all children
#                 for child in node.children:
#                     child_vec *= child.lik            # element‑wise product

#                 # effective node‑transition matrix (Eq. 4)
#                 Phi_eff = self._effective_node_matrix(node)

#                 # propagate upward (Eq. 5)
#                 node.lik = Phi_eff.T @ child_vec      # note the transpose

#         # Root prior is degenerate (founder type known)
#         root = self.tree.get_tree_root()
#         root_prior = np.zeros(K)
#         root_prior[self.root_type] = 1.0
#         logL = np.log(np.dot(root_prior, root.lik) + 1e-12)   # tiny epsilon for safety
#         return logL

#     # def pruninglikelihood(self) -> float:
#     #             # leaf emission = one‑hot vector
#     #             lk = np.zeros(K)
#     #             lk[self.leaf2type_int[node.name]] = 1.0
#     #             node.lik = lk
#     #             # leaf nodes have no L0/L1 (they are not division points)
#     #             node.L0 = node.L1 = None
#     #             node.lr_ratio = 1.0                # dummy; never used for leaves
#     #         else:
#     #             # product of children likelihood vectors (element‑wise)
#     #             child_prod = np.ones(K)
#     #             for child in node.children:
#     #                 child_prod *= child.lik

#     #             # ----- compute the *pure* likelihood contributions -----
#     #             # These are the scalar L0 and L1 that we will need for the LR.
#     #             # The parent prior vector that reaches this node is NOT known
#     #             # yet (it will be the node.up.lik after the upward pass), so
#     #             # for the LR we simply use the *local* part:
#     #             #     L0_u = (Phi0.T @ child_prod)   (a K‑vector)
#     #             #     L1_u = (Phi1.T @ child_prod)
#     #             # Later, when we compute the responsibility, we will dot these
#     #             # with the *incoming* parent prior.
#     #             node.L0_vec = self.Phi0.T @ child_prod   # shape (K,)
#     #             node.L1_vec = self.Phi1.T @ child_prod

#     #             # store a scalar version of the LR that will be used as a
#     #             # weight irrespective of the parent prior:
#     #             #    LR_u = max(L1_vec) / (max(L0_vec) + eps)
#     #             # Using the max (or sum) works because all entries are
#     #             # non‑negative and we only need a *relative* magnitude.
#     #             node.lr_ratio = node.L1_vec.max() / (node.L0_vec.max() + eps)

#     #             # ----- build the *mixture* matrix for the upward pass -----
#     #             # The mixture weight for this node is a leaf‑proportion weighted
#     #             # average of the clade‑specific p_c (exactly the original code).
#     #             clade_weights = defaultdict(float)
#     #             for leaf in node.iter_leaves():
#     #                 clade_weights[self.leaf2clade[leaf.name]] += 1.0
#     #             total = sum(clade_weights.values())
#     #             Phi_eff = np.zeros((K, K))
#     #             for c, w in clade_weights.items():
#     #                 p = self.p_frac[c]
#     #                 Phi_eff += (w / total) * ((1.0 - p) * self.Phi0 + p * self.Phi1)

#     #             # propagate upward
#     #             node.lik = Phi_eff.T @ child_prod

#     #     # Root prior (founder type known)
#     #     root = self.tree.get_tree_root()
#     #     root_prior = np.zeros(K)
#     #     root_prior[self.root_type] = 1.0
#     #     logL = np.log(np.dot(root_prior, root.lik) + eps)
#     #     return logL

#     def _node_weights_uniform(self, node) -> Dict[str, float]:
#         """
#         Return a dict   weight[c] = 1  if clade c appears under ``node``,
#         otherwise 0.  This implements the “inverse‑leaf‑count” (uniform‑per‑node)
#         weighting described in the answer.
#         """
#         weights = {}
#         for c in node.clade_set:            # clade_set already contains every clade
#             weights[c] = 1.0                # every present clade gets the same weight
#         return weights

#     def _node_weights_inverse_leaf(self, node) -> Dict[str, float]:
#         """
#         weight[c] = 1 / (number of descendant leaves under node)   if clade c appears.
#         """
#         n_leaves = len(list(node.iter_leaves()))   # total leaves under node
#         if n_leaves == 0:
#             return {}
#         inv = 1.0 / n_leaves
#         return {c: inv for c in node.clade_set}

#     def _node_weights_lr(self, node):
#         return {c:node.lr_ratio for c in node.clade_set}

#     def _node_weights_lr_plus_inv(self, node):
#         # _node_weights_lr_plus_uniform would be identical to _node_weights_lr because of the uniform weights
#         n_leaves = len(list(node.iter_leaves()))   # total leaves under node
#         if n_leaves == 0:
#             return {}
#         inv = 1.0 / n_leaves
#         return {c: inv*node.lr_ratio for c in node.clade_set} # only difference between this adn inv leaf weights is lr producdt


#     # ------------------------------------------------------------------
#     #  EM algorithm
#     # ------------------------------------------------------------------
#     def fit(self,
#             max_iter: int = 100,
#             tol: float = 1e-6,
#             verbose: bool = True,
#             prior_alpha: float = 1.0,
#             prior_beta: float = 1.0) -> list:
#         """
#         Run the EM algorithm.

#         Parameters
#         ----------
#         max_iter : int
#             Maximum number of EM iterations.
#         tol : float
#             Convergence tolerance on the log‑likelihood.
#         verbose : bool
#             If True, tqdm progress bar + logging output.
#         prior_alpha, prior_beta : float
#             Hyper‑parameters of a **Beta(α,β)** prior on every p_c.
#             Setting both to 1 yields a flat prior (the version in the
#             original code).  Larger α/(α+β) pushes the estimate toward 1,
#             larger β pushes it toward 0.

#         Returns
#         -------
#         logL_history : list of float
#             Log‑likelihood after each EM iteration.
#         """
#         logL_history = []

#         # ------------------------------------------------------------------
#         #  Helper: compute the posterior responsibility γ_u for a given node
#         # ------------------------------------------------------------------
#         def _node_responsibility(node, parent_prior_vec):
#             """
#             Compute γ_u = Pr(Z_u = 1 | data) for an internal node ``node``.
#             ``parent_prior_vec`` is the conditional likelihood vector that
#             arrives from the parent (i.e. the parent’s *lik* before it is
#             multiplied by its own Φ).  The formula is Eq. 8′.
#             """
#             # effective node matrix for the two pure processes (no mixture)
#             # R0 = Φ0ᵀ (children_product), R1 = Φ1ᵀ (children_product)
#             child_product = np.ones(self.K)
#             for child in node.children:
#                 child_product *= child.lik

#             R0 = self.Phi0.T @ child_product    # shape (K,)
#             R1 = self.Phi1.T @ child_product

#             # scalar prior weight for the clade(s) of this node
#             # (same weighting that was used to build Φ_eff)
#             clade_weights = defaultdict(float)
#             for leaf in node.iter_leaves():
#                 clade_weights[self.leaf2clade[leaf.name]] += 1.0
#             total = sum(clade_weights.values())
#             prior_mix = 0.0
#             for c, w in clade_weights.items():
#                 prior_mix += (w / total) * self.p_frac[c]

#             # Now apply Bayes rule (Eq. 8′)
#             num = prior_mix * np.dot(parent_prior_vec, R1)
#             den = (prior_mix * np.dot(parent_prior_vec, R1) +
#                    (1.0 - prior_mix) * np.dot(parent_prior_vec, R0) + 1e-12)
#             return num / den

#         for it in tqdm.trange(max_iter, disable=not verbose):
#             # ------------------------------------------------------------------
#             #  E‑step : pruning pass + compute responsibilities γ_u
#             # ------------------------------------------------------------------
#             logL = self._pruning_likelihood()
#             logL_history.append(logL)

#             if verbose:
#                 logging.info(f"EM iteration {it+1:02d} – logL = {logL:.6f}")

#             if self.weight_strategy == 'leaf_prop':
#                 # ------------------------------------------------------------------
#                 #  Gather sufficient statistics (M‑step needs them)
#                 # ------------------------------------------------------------------
#                 # clade_stats[c] = (Σ w_u^{(c)} * γ_u , Σ w_u^{(c)})
#                 clade_stats = {c: [0.0, 0.0] for c in self.clades}
#                 root = self.tree.get_tree_root()
#                 # The root has no incoming Z, but its children do; we start the
#                 # recursion with the *root prior* as the incoming vector.
#                 parent_vec = np.zeros(self.K)
#                 parent_vec[self.root_type] = 1.0
    
#                 # Post‑order traversal again (children already have .lik from pruning)
#                 for node in self.tree.traverse("postorder"):
#                     if node.is_root():
#                         continue                      # nothing to update at the root
#                     # compute γ_u using the *parent* conditional vector
#                     # (parent is node.up)
#                     parent_node = node.up
#                     gamma_u = _node_responsibility(node, parent_node.lik)
    
#                     # weight each clade proportionally to the number of descendant leaves
#                     clade_weights = defaultdict(float)
#                     for leaf in node.iter_leaves():
#                         clade_weights[self.leaf2clade[leaf.name]] += 1.0
#                     total = sum(clade_weights.values())
    
#                     for c, w in clade_weights.items():
#                         weight = w / total
#                         clade_stats[c][0] += weight * gamma_u   # numerator Σ w·γ
#                         clade_stats[c][1] += weight             # denominator Σ w
    
#                 # ------------------------------------------------------------------
#                 #  M‑step : MAP update of the mixture fractions p_c
#                 # ------------------------------------------------------------------
#                 for c in self.clades:
#                     num, den = clade_stats[c]            # num = Σ w·γ , den = Σ w
#                     if den > 0:
#                         # Beta(α,β) prior → MAP (Eq. 12′)
#                         α = prior_alpha
#                         β = prior_beta
#                         new_p = (num + α - 1.0) / (den + α + β - 2.0)
#                         # keep away from exact 0/1 for numerical stability
#                         self.p_frac[c] = np.clip(new_p, 1e-5, 1 - 1e-5)
    
#                 # ------------------------------------------------------------------
#                 #  Convergence test
#                 # ------------------------------------------------------------------
#                 if it > 0 and abs(logL_history[-1] - logL_history[-2]) < tol:
#                     if verbose:
#                         logging.info("EM converged (log‑likelihood change < tol).")
#                     break

#             else:
#                 clade_stats = {c: [0.0, 0.0] for c in self.clades}   # [numer, denom]

#                 for node in self.tree.traverse("postorder"):
#                     if node.is_root():
#                         continue
#                     gamma_u = _node_responsibility(node, node.up.lik)   # already defined
        
#                     # uniform weight for each clade covered by the node
#                     if self.weight_strategy == 'uniform':
#                         w_dict = self._node_weights_uniform(node)           # {c: 1.0}
#                     elif self.weight_strategy == 'inv_leaf':
#                         w_dict = self._node_weights_inverse_leaf(node)
#                     for c, w in w_dict.items():
#                         clade_stats[c][0] += w * gamma_u   # numerator Σ γ
#                         clade_stats[c][1] += w            # denominator Σ 1
        
#                 # ---------- MAP update ----------
#                 for c in self.clades:
#                     num, den = clade_stats[c]
#                     if den > 0:
#                         new_p = (num + prior_alpha - 1.0) / (den + prior_alpha + prior_beta - 2.0)
#                         self.p_frac[c] = np.clip(new_p, 1e-5, 1.0 - 1e-5)
        
#                     # ---------- Convergence check ----------
#                     if it > 0 and abs(logL_history[-1] - logL_history[-2]) < tol:
#                         if verbose:
#                             logging.info("EM converged (log‑likelihood change < tol).")
#                         break

            

#         self.logL_history_ = logL_history
#         return logL_history

#     # ------------------------------------------------------------------
#     #  Post‑processing utilities
#     # ------------------------------------------------------------------
#     def clade_posteriors(self) -> pd.DataFrame:
#         """
#         Return a DataFrame with one row per clade:
#             clade | p_induced | hard_label
#         """
#         rows = []
#         for c in self.clades:
#             p = self.p_frac[c]
#             rows.append({
#                 'clade': c,
#                 'p_induced': p,
#                 'hard_label': 'induced' if p > 0.5 else 'un-induced'
#             })
#         return pd.DataFrame(rows)

#     def summary(self) -> None:
#         """Print a quick textual summary."""
#         df = self.clade_posteriors()
#         print(df.to_string(index=False))
#         if hasattr(self, "logL_history_"):
#             print(f"\nFinal log‑likelihood = {self.logL_history_[-1]:.6f}")



import logging
from collections import defaultdict
from typing import Dict, Any, Optional

import numpy as np
import pandas as pd
import tqdm


def checkstochastic(M: np.ndarray, name: str) -> None:
    """Small helper – raise if rows do not sum to 1 or contain <0."""
    if not np.allclose(M.sum(axis=1), 1.0, atol=1e-9):
        raise ValueError(f"{name} rows must sum to 1.")
    if np.any(M < -1e-12):
        raise ValueError(f"{name} contains negative entries.")


class TreeMixtureNodeEM2():
    """
    EM for a node‑only mixture of two K‑state stochastic matrices on a rooted
    phylogeny.  The *likelihood‑ratio* weighting scheme is implemented as one
    of the possible ``weight_strategy`` options.
    """

    # ------------------------------------------------------------------
    #  Constructor
    # ------------------------------------------------------------------
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
                'leaf_prop'   – leaf‑proportion weighting (original code)
                'uniform'     – uniform‑per‑node
                'inv_leaf'    – inverse‑leaf weighting
                
                
                'lr+uniform'  – LR multiplied by uniform (default for LR)
                'lr+inv'      – LR multiplied by inverse‑leaf
        """
        # --------------------------------------------------------------
        # 1️⃣  Store tree and attach leaf metadata
        # --------------------------------------------------------------
        self.tree = tree
        # self._check_tree_is_rooted()
        self.leaf2clade = leaf2clade
        self.leaf2type = leaf2type
        self.weight_strategy = weight_strategy

        leaf_names = {leaf.name for leaf in self.tree.iter_leaves()}
        missing = leaf_names - set(self.leaf2clade.keys())
        if missing:
            raise ValueError(f"The following leaves are missing from leaf2clade: {missing}")

        # --------------------------------------------------------------
        # 2️⃣  Encode cell types as integers 0 … K‑1
        # --------------------------------------------------------------
        if type_encoder is None:
            uniq = sorted(set(self.leaf2type.values()))
            self.type_encoder = {t: i for i, t in enumerate(uniq)}
        else:
            self.type_encoder = type_encoder
        self.leaf2type_int = {
            leaf: self.type_encoder[self.leaf2type[leaf]]
            for leaf in leaf_names
        }

        # --------------------------------------------------------------
        # 3️⃣  Store the two node‑transition matrices (must be stochastic)
        # --------------------------------------------------------------
        self.Phi0 = np.asarray(Phi0, dtype=float)
        self.Phi1 = np.asarray(Phi1, dtype=float)
        self.K = self.Phi0.shape[0]
        assert self.Phi0.shape == (self.K, self.K)
        assert self.Phi1.shape == (self.K, self.K)
        checkstochastic(self.Phi0, "Phi0")
        checkstochastic(self.Phi1, "Phi1")

        # --------------------------------------------------------------
        # 4️⃣  Root prior (founder type is known)
        # --------------------------------------------------------------
        if isinstance(root_type, str):
            root_type = self.type_encoder[root_type]
            # root_type = self.type_encoder.get(root_type, max(self.type_encoder.values())+1)
        self.root_type = int(root_type)

        # --------------------------------------------------------------
        # 5️⃣  Initialise clade‑specific mixture parameters
        # --------------------------------------------------------------
        self.clades = sorted(set(self.leaf2clade.values()))
        self.p_frac = {c: 0.5 for c in self.clades}

        # cache for fast node‑to‑clade look‑up
        self._assign_node_clades()

    # ------------------------------------------------------------------
    #  Helper methods
    # ------------------------------------------------------------------
    def _check_tree_is_rooted(self) -> None:
        """Raise if the supplied tree has no root."""
        if not self.tree.get_tree_root():
            raise ValueError("The supplied tree must be rooted.")

    def _assign_node_clades(self) -> None:
        """
        Attach a ``clade_set`` attribute to every node:
        the set of clades that have at least one descendant leaf under that node.
        """
        for node in self.tree.traverse("postorder"):
            if node.is_leaf():
                node.clade_set = {self.leaf2clade[node.name]}
            else:
                node.clade_set = set()
                for child in node.children:
                    node.clade_set.update(child.clade_set)

    # ------------------------------------------------------------------
    #  Node‑weight factories
    # ------------------------------------------------------------------
    def _node_weights_uniform(self, node) -> Dict[str, float]:
        """Uniform‑per‑node: weight = 1 for every clade covered by the node."""
        return {c: 1.0 for c in node.clade_set}

    def _node_weights_inverse_leaf(self, node) -> Dict[str, float]:
        """Inverse‑leaf: weight = 1 / (# descendant leaves)."""
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
        """LR × uniform (equivalent to pure LR because uniform = 1)."""
        return {c: node.lr_ratio for c in node.clade_set}

    def _node_weights_lr_plus_inv(self, node) -> Dict[str, float]:
        """LR × inverse‑leaf."""
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
        
    # 528
    # def pruninglikelihood(self) -> float:

    #     K = self.K
    #     eps = 1e-12
    
    #     total_log_factor = 0.0               # cumulative scaling term

    #     all_node_lr_ratios = {}
    
    #     for node in self.tree.traverse("postorder"):
    #         if node.is_leaf():
    #             # ----- leaf emission (one‑hot) -----
    #             lk = np.zeros(K)
    #             lk[self.leaf2type_int[node.name]] = 1.0
    #             node.lik = lk
    
    #             # log vectors: log(1)=0 for the observed state, -inf otherwise
    #             node.logL0 = np.full(K, -np.inf)
    #             node.logL1 = np.full(K, -np.inf)
    #             obs = self.leaf2type_int[node.name]
    #             node.logL0[obs] = 0.0
    #             node.logL1[obs] = 0.0
    
    #             node.lr_ratio = 1.0               # dummy, never used for leaves
    #             node.log_factor = 0.0
    #             continue

    #         # ----- 1️⃣ lift each child through Φ0 and Φ1 and sum the logs ----------
    #         logL0 = np.zeros(K)
    #         logL1 = np.zeros(K)
    #         for child in node.children:

    #             # 528
    #             prob0 = self.Phi0.T @ child.lik      # shape (K,)
    #             prob1 = self.Phi1.T @ child.lik
    #             # 528
                
    #             # # 528
    #             # prob0 = self.Phi0 @ child.lik      # shape (K,)
    #             # prob1 = self.Phi1 @ child.lik
    #             # # 528
                
    #             # avoid log(0); eps is tiny enough not to affect the result
    #             prob0[prob0 == 0] = eps
    #             prob1[prob1 == 0] = eps
    #             logL0 += np.log(prob0)
    #             logL1 += np.log(prob1)
    
    #         # ----- 2️⃣ remove a common offset (scaling) -----------------------------
    #         max0 = logL0.max()
    #         node.scale_factor = max0 # 528
    #         node.logL0 = logL0 - max0
    #         node.logL1 = logL1 - max0          # same offset for both components
    #         node.log_factor = max0
    
    #         # ----- 3️⃣ LR scalar (for “lr” weighting) -----------------------------
    #         node.lr_ratio = np.exp(node.logL1.max()) / (np.exp(node.logL0.max()) + eps) 
    #         all_node_lr_ratios[node] = node.lr_ratio
            
    
    #         # ----- 4️⃣ clade‑specific mixture probability p_u ----------------------
    #         clade_weights = defaultdict(float)
    #         for leaf in node.iter_leaves():
    #             clade_weights[self.leaf2clade[leaf.name]] += 1.0
    #         total = sum(clade_weights.values())
    #         if total == 0:
    #             p_u = 0.0
    #         else:
    #             p_u = sum((w / total) * self.p_frac[c]
    #                       for c, w in clade_weights.items())
    
    #         # ----- 5️⃣ mixture upward vector (still a probability vector) ----------
    #         L0 = np.exp(node.logL0)            # already divided by exp(max0)
    #         L1 = np.exp(node.logL1) 
    #         # node.lik0 = L0.copy() # 528
    #         # node.lik1 = L1.copy() # 528
    #         mixed = (1.0 - p_u) * L0 + p_u * L1
    
    #         mix_max = mixed.max()
    #         if mix_max == 0:
    #             mixed = np.full(K, eps)
    #             mix_max = eps
    #         node.lik = mixed / mix_max          # keep it normalised
    #         node.log_factor += np.log(mix_max)  # add the mixture scaling
    
    #         total_log_factor += node.log_factor
    
    #     # ----- 6️⃣ root contribution ------------------------------------------------
    #     root = self.tree.get_tree_root()
    #     root_prior = np.zeros(K)
    #     root_prior[self.root_type] = 1.0
    
    #     logL = np.log(np.dot(root_prior, root.lik) + eps) + total_log_factor
    
    #     # ----- 7️⃣ diagnostics (only once) -----------------------------------------
    #     if not hasattr(self, "_diag_printed"):
    #         internal = sum(1 for n in self.tree.traverse() if not n.is_leaf())
    #         mixed_nodes = sum(1 for n in self.tree.traverse()
    #                           if not n.is_leaf() and getattr(n, "has_mixture", False))
    #         print(f"[Mixture applied at {mixed_nodes}/{internal} internal nodes]")
    #         self._diag_printed = True
    
    #     return logL, all_node_lr_ratios
    # 528

    def pruninglikelihood(self) -> float:
        K   = self.K
        eps = 1e-12
        total_log_factor = 0.0
    
        for node in self.tree.traverse("postorder"):
            # ----------------------------------------------------------
            # 1️⃣  LEAVES – pure one‑hot probability vector
            # ----------------------------------------------------------
            if node.is_leaf():
                lk = np.zeros(K)
                lk[self.leaf2type_int[node.name]] = 1.0
                node.L0_vec = lk.copy()                # pure un‑induced
                node.L1_vec = lk.copy()                # pure induced
                node.lik    = lk.copy()                # mixture (= identical)
                node.logL0  = np.where(lk > 0, 0.0, -np.inf)
                node.logL1  = np.where(lk > 0, 0.0, -np.inf)
                node.scale_factor = 0.0
                node.log_factor   = 0.0
                node.lr_ratio     = 1.0
                continue
    
            # ----------------------------------------------------------
            # 2️⃣  INTERNAL NODES – propagate each *component* separately
            # ----------------------------------------------------------
            L0 = np.ones(K)
            L1 = np.ones(K)
            for child in node.children:
                # Pure un‑induced propagation – uses the child's *pure* L0 vector
                prob0 = self.Phi0.T @ child.L0_vec
                prob0[prob0 == 0] = eps
                L0 *= prob0
    
                # Pure induced propagation – uses the child's *pure* L1 vector
                prob1 = self.Phi1.T @ child.L1_vec
                prob1[prob1 == 0] = eps
                L1 *= prob1
    
            # ----------------------------------------------------------
            # 3️⃣  Scale for numerical stability (per‑node scalar)
            # ----------------------------------------------------------
            max0 = max(L0.max(), L1.max())
            L0 /= max0
            L1 /= max0
            node.scale_factor = np.log(max0)         # the *only* scaling that
                                                      # belongs to *this* node
    
            # Store the *pure* component log‑vectors (no eps added)
            with np.errstate(divide='ignore'):
                node.logL0 = np.log(L0)
                node.logL1 = np.log(L1)
    
            node.L0_vec = L0
            node.L1_vec = L1
    
            # ----------------------------------------------------------
            # 4️⃣  LR scalar (used only for ‘lr*’ weighting strategies)
            # ----------------------------------------------------------
            node.lr_ratio = L1.max() / (L0.max() + eps)
    
            # ----------------------------------------------------------
            # 5️⃣  Clade‑specific mixture probability p_u
            # ----------------------------------------------------------
            clade_weights = defaultdict(float)
            for leaf in node.iter_leaves():
                clade_weights[self.leaf2clade[leaf.name]] += 1.0
            total = sum(clade_weights.values())
            p_u   = 0.0 if total == 0 else sum(
                        (w / total) * self.p_frac[c] for c, w in clade_weights.items())
    
            # ----------------------------------------------------------
            # 6️⃣  Build the mixture upward vector (this is what the PARENT will see)
            # ----------------------------------------------------------
            mixed = (1.0 - p_u) * L0 + p_u * L1
            mix_max = mixed.max()
            if mix_max == 0:
                mixed = np.full(K, eps)
                mix_max = eps
            node.lik = mixed / mix_max               # normalised mixture vector
            node.log_factor = node.scale_factor + np.log(mix_max)   # cumulative
    
            total_log_factor += node.log_factor
    
        # --------------------------------------------------------------
        # Root contribution
        # --------------------------------------------------------------
        root = self.tree.get_tree_root()
        root_prior = np.zeros(K)
        root_prior[self.root_type] = 1.0
    
        logL = np.log(np.dot(root_prior, root.lik) + eps) + total_log_factor
    
        # Diagnostic (kept for compatibility)
        if not hasattr(self, "_diag_printed"):
            internal = sum(1 for n in self.tree.traverse() if not n.is_leaf())
            mixed_nodes = sum(1 for n in self.tree.traverse()
                              if not n.is_leaf() and getattr(n, "has_mixture", False))
            self._diag_printed = True
    
        return logL

    
    # def pruninglikelihood(self) -> float:
    #     """
    #     Bottom‑up (post‑order) Felsenstein pass.
    #     Stores for every internal node:
    #         - L0_vec, L1_vec          – pure‑component vectors
    #         - lr_ratio                – LR weight for the “lr” strategies
    #         - lik                     – mixture likelihood sent upward
    #     """
    #     K = self.K
    #     eps = 1e-12
    
    #     # diagnostics
    #     nodes_with_mixture = 0
    #     nodes_with_empty_clade_weights = 0
    
    #     for node in self.tree.traverse("postorder"):
    #         if node.is_leaf():
    #             # ----- leaf emission (one‑hot) -----
    #             lk = np.zeros(K)
    #             lk[self.leaf2type_int[node.name]] = 1.0
    #             node.lik = lk
    
    #             # dummy values – never used in the E‑step
    #             node.L0_vec = node.L1_vec = None
    #             node.lr_ratio = 1.0
    #             node.child_prod = None
    #             continue
    
    #         # ----- 1️⃣  accumulate lifted likelihoods for the two components -----
    #         L0_vec = np.ones(K)               # product_{children} (Φ0ᵀ @ child.lik)
    #         L1_vec = np.ones(K)               # product_{children} (Φ1ᵀ @ child.lik)
    
    #         for child in node.children:
    #             L0_vec *= self.Phi0.T @ child.lik
    #             L1_vec *= self.Phi1.T @ child.lik
    
    #         node.L0_vec = L0_vec
    #         node.L1_vec = L1_vec
    
    #         # ----- 2️⃣  LR scalar used by the “lr” weight strategies -----
    #         node.lr_ratio = L1_vec.max() / (L0_vec.max() + eps)
    
    #         # ----- 3️⃣  clade‑specific mixture weight p_u ----------------------
    #         clade_weights = defaultdict(float)
    #         for leaf in node.iter_leaves():
    #             clade_weights[self.leaf2clade[leaf.name]] += 1.0
    #         total = sum(clade_weights.values())
    
    #         if not clade_weights:                       # sanity check – should never happen
    #             print(f"⚠️  WARNING: Node {node.name} has EMPTY clade_weights!")
    #             nodes_with_empty_clade_weights += 1
    #             p_u = 0.0
    #         else:
    #             p_u = sum((w / total) * self.p_frac[c] for c, w in clade_weights.items())
    
    #         # ----- 4️⃣  mixture likelihood that is passed upward ---------------
    #         #   ℓ_u = (1 - p_u)·L0_vec + p_u·L1_vec
    #         node.lik = (1.0 - p_u) * L0_vec + p_u * L1_vec
    #         node.has_mixture = p_u > 0.0                # for diagnostics
    #         if node.has_mixture:
    #             nodes_with_mixture += 1
    
    #     # ----- root log‑likelihood -------------------------------------------------
    #     root = self.tree.get_tree_root()
    #     root_prior = np.zeros(K)
    #     root_prior[self.root_type] = 1.0
    #     logL = np.log(np.dot(root_prior, root.lik) + eps)
    
    #     # ----- diagnostic summary (run only once) ----------------------------------
    #     if not hasattr(self, "_diag_printed"):
    #         total_internal = sum(1 for n in self.tree.traverse() if not n.is_leaf())
    #         print(f"[Mixture applied at {nodes_with_mixture}/{total_internal} internal nodes]")
    #         if nodes_with_empty_clade_weights:
    #             print(f"[WARNING: {nodes_with_empty_clade_weights} nodes had empty clade_weights]")
    #         self._diag_printed = True
    
    #     return logL

    # 5/7
    # # ------------------------------------------------------------------
    # #  Likelihood calculation (core of the EM)
    # # ------------------------------------------------------------------
    # def pruninglikelihood(self) -> float:
    #     """
    #     Bottom‑up (post‑order) Felsenstein pass.
    #     For every **internal** node we store:
    #         - node.L0_vec  : Φ0ᵀ @ (product of child likelihoods)   (K‑vector)
    #         - node.L1_vec  : Φ1ᵀ @ (product of child likelihoods)   (K‑vector)
    #         - node.lr_ratio: scalar L1 / (L0 + eps)  (used as LR weight)
    #     Leaves get dummy values (they will never be used in the responsibility
    #     calculation).  The root also gets a dummy lr_ratio = 1.0 – it is never
    #     used because the root has no incoming Z variable.
    #     """
    #     K = self.K
    #     eps = 1e-12
        
    #     # Track diagnostics
    #     nodes_with_mixture = 0
    #     nodes_with_empty_clade_weights = 0
    
    #     for node in self.tree.traverse("postorder"):

    #         # new addition
    #         # node.L0_vec = node.L1_vec = None
            
    #         if node.is_leaf():
    #             # ----- leaf emission -----
    #             lk = np.zeros(K)
    #             lk[self.leaf2type_int[node.name]] = 1.0
    #             node.lik = lk
    
    #             # dummy values for leaves (they are never division points)
    #             node.L0_vec = node.L1_vec = None
    #             node.lr_ratio = 1.0
    #             node.child_prod = None
    #         else:
    #             # ----- product of children likelihood vectors -----
    #             child_prod = np.ones(K)
    #             for child in node.children:
    #                 child_prod *= child.lik
    #             node.child_prod = child_prod
    
    #             # ----- pure (un‑mixed) likelihood contributions -----
    #             node.L0_vec = self.Phi0.T @ child_prod      # shape (K,)
    #             node.L1_vec = self.Phi1.T @ child_prod      # shape (K,)
    
    #             # ----- LR scalar (used later as a weight) -----
    #             node.lr_ratio = node.L1_vec.max() / (node.L0_vec.max() + eps)
    
    #             # ===== BUILD MIXTURE MATRIX (with diagnostics) =====
    #             clade_weights = defaultdict(float)
    #             for leaf in node.iter_leaves():
    #                 clade_weights[self.leaf2clade[leaf.name]] += 1.0
    #             total = sum(clade_weights.values())
    
    #             # DIAGNOSTIC: Check if clade_weights is empty
    #             if not clade_weights:
    #                 print(f"⚠️  WARNING: Node {node.name} has EMPTY clade_weights!")
    #                 nodes_with_empty_clade_weights += 1
    
    #             # Build Phi_eff with diagnostics
    #             Phi_eff = np.zeros((K, K))
    #             Phi_eff_is_mixture = False  # track if Phi_eff differs from Phi0
    
    #             for c, w in clade_weights.items():
    #                 p = self.p_frac[c]
    #                 # print(f"   DEBUG: clade {c}, p_frac[{c}] = {p}, p > 1e-6 = {p > 1e-6}")
    #                 component = (w / total) * ((1.0 - p) * self.Phi0 + p * self.Phi1)
    #                 Phi_eff += component

    #                 baseline_component = (w/total)*self.Phi0
    #                 component_diff = np.linalg.norm(component - baseline_component)
    #                 # print(f"      component differs from baseline by: {component_diff:.10f}")
    
    #                 # Check if this component differs from Phi0
    #                 if p > 1e-6:  # if p is noticeably different from 0
    #                     Phi_eff_is_mixture = True
    
    #             # Track if mixture was applied
    #             if Phi_eff_is_mixture:
    #                 nodes_with_mixture += 1
                
    #             node.has_mixture = Phi_eff_is_mixture
    
    #             # ----- propagate upward -----
    #             node.lik = Phi_eff.T @ child_prod
    
    #     # Count total internal nodes
    #     total_internal = sum(1 for n in self.tree.traverse() if not n.is_leaf())
    
    #     # ===== DIAGNOSTIC: Summary of mixture usage =====
    #     if not hasattr(self, '_diag_printed'):
    #         print(f"   [Mixture applied at {nodes_with_mixture}/{total_internal} internal nodes]")
    #         if nodes_with_empty_clade_weights > 0:
    #             print(f"   [WARNING: {nodes_with_empty_clade_weights} nodes had empty clade_weights]")
    #         self._diag_printed = True
    
    #     # ----- root prior and total log‑likelihood -----
    #     root = self.tree.get_tree_root()
    #     root_prior = np.zeros(K)
    #     root_prior[self.root_type] = 1.0

    #     print(f'{root.lik = }')
    #     logL = np.log(np.dot(root_prior, root.lik) + eps)
        
    #     return logL 
    # 5/7

    
    # def pruninglikelihood(self) -> float:
    #     """
    #     Bottom‑up (post‑order) Felsenstein pass.
    #     For every **internal** node we store:
    #         - node.L0_vec  : Φ0ᵀ @ (product of child likelihoods)   (K‑vector)
    #         - node.L1_vec  : Φ1ᵀ @ (product of child likelihoods)   (K‑vector)
    #         - node.lr_ratio: scalar L1 / (L0 + eps)  (used as LR weight)
    #     Leaves get dummy values (they will never be used in the responsibility
    #     calculation).  The root also gets a dummy lr_ratio = 1.0 – it is never
    #     used because the root has no incoming Z variable.
    #     """
    #     K = self.K
    #     eps = 1e-12
    
    #     for node in self.tree.traverse("postorder"):
    #         if node.is_leaf():
    #             # ----- leaf emission -----
    #             lk = np.zeros(K)
    #             lk[self.leaf2type_int[node.name]] = 1.0
    #             node.lik = lk
    
    #             # dummy values for leaves (they are never division points)
    #             node.L0_vec = node.L1_vec = None
    #             node.lr_ratio = 1.0
    #             node.child_prod = None
    #         else:
    #             # ----- product of children likelihood vectors -----
    #             child_prod = np.ones(K)
    #             for child in node.children:
    #                 child_prod *= child.lik
    #             node.child_prod = child_prod
    
    #             # ----- pure (un‑mixed) likelihood contributions -----
    #             node.L0_vec = self.Phi0.T @ child_prod      # shape (K,)
    #             node.L1_vec = self.Phi1.T @ child_prod      # shape (K,)
    
    #             # ----- LR scalar (used later as a weight) -----
    #             node.lr_ratio = node.L1_vec.max() / (node.L0_vec.max() + eps)
    
    #             # ----- build the mixture matrix for the upward pass -----
    #             clade_weights = defaultdict(float)
    #             for leaf in node.iter_leaves():
    #                 clade_weights[self.leaf2clade[leaf.name]] += 1.0
    #             total = sum(clade_weights.values())
    
    #             Phi_eff = np.zeros((K, K))
    #             for c, w in clade_weights.items():
    #                 p = self.p_frac[c]
    #                 Phi_eff += (w / total) * ((1.0 - p) * self.Phi0 + p * self.Phi1)
    
    #             # ----- propagate upward -----
    #             node.lik = Phi_eff.T @ child_prod
    
    #     # ----- root prior and total log‑likelihood -----
    #     root = self.tree.get_tree_root()
    #     root_prior = np.zeros(K)
    #     root_prior[self.root_type] = 1.0
    #     logL = np.log(np.dot(root_prior, root.lik) + eps)


    #     return logL


    # def update_mixture(self):
    #     """Re‑compute node.lik after self.p_frac changes."""
    #     K = self.K
    #     for node in self.tree.traverse("postorder"):
    #         if node.is_leaf():
    #             continue
    
    #         # recompute lifted products exactly as in pruninglikelihood
    #         L0_vec = np.ones(K)
    #         L1_vec = np.ones(K)
    #         for child in node.children:
    #             L0_vec *= self.Phi0.T @ child.lik
    #             L1_vec *= self.Phi1.T @ child.lik
    #         node.L0_vec, node.L1_vec = L0_vec, L1_vec
    
    #         # clade‑specific mixture weight
    #         clade_weights = defaultdict(float)
    #         for leaf in node.iter_leaves():
    #             clade_weights[self.leaf2clade[leaf.name]] += 1.0
    #         total = sum(clade_weights.values())
    #         p_u = sum((w / total) * self.p_frac[c] for c, w in clade_weights.items()) \
    #               if total > 0 else 0.0
    
    #         node.lik = (1.0 - p_u) * L0_vec + p_u * L1_vec
    # 5/7
    # def update_mixture(self):
    #     """Re‑compute node likelihoods after self.p_frac changes."""
    #     K = self.K
    #     # Post-order traversal ensures children are updated before parents
    #     for node in self.tree.traverse("postorder"):
    #         if node.is_leaf():
    #             continue
            
    #         # RECOMPUTE child_prod from fresh child likelihoods
    #         child_prod = np.ones(K)
    #         for child in node.children:
    #             child_prod *= child.lik
    #         node.child_prod = child_prod  # Update the cache
            
    #         # Compute clade weights for this node
    #         clade_weights = defaultdict(float)
    #         for leaf in node.iter_leaves():
    #             clade_weights[self.leaf2clade[leaf.name]] += 1.0
    #         total = sum(clade_weights.values())
            
    #         if total == 0:
    #             node.lik = self.Phi0.T @ child_prod
    #             continue
            
    #         # Build the mixture matrix with updated p_frac
    #         Phi_eff = np.zeros((K, K))
    #         for c, w in clade_weights.items():
    #             p = self.p_frac[c]  # CURRENT (updated) value
    #             Phi_eff += (w / total) * ((1.0 - p) * self.Phi0 + p * self.Phi1)
            
    #         # Update the node likelihood
    #         node.lik = Phi_eff.T @ child_prod
    # 5/7

    # def update_mixture(self):
    #     """Re‑compute Φ_eff after self.p_frac changes."""
    #     K = self.K
    #     for node in self.tree.traverse("postorder"):
    #         if node.is_leaf():
    #             continue
    #         # reuse already‑computed child_prod
    #         child_prod = node.child_prod if node.child_prod is not None else np.ones(K)

    #         # recompute clade weights for this node
    #         clade_weights = defaultdict(float)
    #         for leaf in node.iter_leaves():
    #             clade_weights[self.leaf2clade[leaf.name]] += 1.0
    #         total = sum(clade_weights.values())
    #         if total == 0:
    #             node.lik = self.Phi0.T @ child_prod
    #             continue

    #         Phi_eff = np.zeros((K, K))
    #         for c, w in clade_weights.items():
    #             p = self.p_frac[c]                 # fresh value
    #             Phi_eff += (w / total) * ((1.0 - p) * self.Phi0 + p * self.Phi1)

    #         node.lik = Phi_eff.T @ child_prod

    # def recompute_mixture(self):
    #     """Re‑compute Φ_eff for all internal nodes after self.p_frac changes."""
    #     K = self.K
    #     eps = 1e-12
    #     for node in self.tree.traverse("postorder"):
    #         if node.is_leaf():
    #             continue
    
    #         # product of child likelihood vectors (already cached)
    #         child_prod = node.child_prod if node.child_prod is not None else np.ones(K)
    
    #         # compute clade weights for this node
    #         clade_weights = defaultdict(float)
    #         for leaf in node.iter_leaves():
    #             clade_weights[self.leaf2clade[leaf.name]] += 1.0
    #         total = sum(clade_weights.values())
    #         if total == 0:
    #             # No leaves under this node – treat it as pure Φ0
    #             node.lik = self.Phi0.T @ child_prod
    #             continue
    
    #         # build the mixture matrix
    #         Phi_eff = np.zeros((K, K))
    #         for c, w in clade_weights.items():
    #             p = self.p_frac[c]                # <-- uses the UPDATED p_frac
    #             Phi_eff += (w / total) * ((1.0 - p) * self.Phi0 + p * self.Phi1)
    
    #         # store the new upward likelihood
    #         node.lik = Phi_eff.T @ child_prod


    ######### comment 5/14
    
    # def node_responsibility(self, node, parent_prior_vec=None):
    #     """
    #     Public version of the responsibility calculation that was previously
    #     nested inside ``fit``.  It returns the posterior probability γ_u that the
    #     division represented by ``node`` used the induced transition matrix.
    
    #     Parameters
    #     ----------
    #     node : ete3.TreeNode
    #         An internal node (must not be a leaf or the root).
    #     parent_prior_vec : np.ndarray, optional
    #         The likelihood vector of the parent (parent.up.lik).  If omitted the
    #         method will retrieve it from the tree; if the tree has not been
    #         evaluated yet it falls back to the root prior.
    
    #     Returns
    #     -------
    #     float
    #         Posterior responsibility γ_u ∈ [0, 1].
    #     """
    #     # ------------------------------------------------------------------
    #     # 1️⃣  Get the parent prior vector
    #     # ------------------------------------------------------------------
    #     if parent_prior_vec is None:
    #         if node.up is None:
    #             # root prior (one‑hot on the known root type)
    #             parent_prior_vec = np.zeros(self.K)
    #             parent_prior_vec[self.root_type] = 1.0
    #         else:
    #             parent_prior_vec = node.up.lik
    
    #     # ------------------------------------------------------------------
    #     # 2️⃣  Compute the clade‑weighted prior mixture (same as in fit)
    #     # ------------------------------------------------------------------
    #     clade_weights = defaultdict(float)
    #     for leaf in node.iter_leaves():
    #         clade_weights[self.leaf2clade[leaf.name]] += 1.0
    #     total = sum(clade_weights.values())
    #     if total == 0:
    #         prior_mix = 0.0
    #     else:
    #         prior_mix = sum((w / total) * self.p_frac[c]
    #                         for c, w in clade_weights.items())
    
    #     # ------------------------------------------------------------------
    #     # 3️⃣  Log‑likelihood pieces that were stored by pruninglikelihood
    #     # ------------------------------------------------------------------
    #     log_factor = node.log_factor                # scaling factor from the log‑space pass
    #     log_pi_ell1 = logsumexp(node.logL1 + log_factor)
    #     log_pi_ell0 = logsumexp(node.logL0 + log_factor)

    
    #     # ------------------------------------------------------------------
    #     # 4️⃣  Return the posterior γ_u
    #     # ------------------------------------------------------------------
    #     if prior_mix <= 0.0:
    #         return 0.0
    #     if prior_mix >= 1.0:
    #         return 1.0
    
    #     log_num = np.log(prior_mix) + log_pi_ell1
    #     term1   = np.log(prior_mix)       + log_pi_ell1
    #     term2   = np.log(1.0 - prior_mix) + log_pi_ell0
    #     log_den = logsumexp(np.array([term1, term2]))
    #     gamma   = np.exp(log_num - log_den)
    #     return float(np.clip(gamma, 0.0, 1.0))

    ######### comment 5/14

    # 528
    # def node_responsibility(self, node, parent_prior_vec=None):
    #     """
    #     Public version of the responsibility calculation.
    #     Returns the posterior probability γ_u that the division represented by ``node``
    #     used the induced transition matrix.
    
    #     Parameters
    #     ----------
    #     node : ete3.TreeNode
    #         An internal node (must not be a leaf or the root).
    #     parent_prior_vec : np.ndarray, optional
    #         The likelihood vector of the parent (parent.up.lik).
    #         If omitted, falls back to the root prior.
    
    #     Returns
    #     -------
    #     float
    #         Posterior responsibility γ_u ∈ [0, 1].
    #     """
    #     eps = 1e-12
        
    #     # ------------------------------------------------------------------
    #     # 1️⃣  Get the parent prior vector
    #     # ------------------------------------------------------------------
    #     if parent_prior_vec is None:
    #         if node.up is None:
    #             parent_prior_vec = np.zeros(self.K)
    #             parent_prior_vec[self.root_type] = 1.0
    #         else:
    #             parent_prior_vec = node.up.lik # 528
    #             # parent_prior_vec = np.zeros(self.K)
    #             # parent_prior_vec[self.root_type] = 1.0
        
    #     # ------------------------------------------------------------------
    #     # 2️⃣  Compute the clade‑weighted prior mixture
    #     # ------------------------------------------------------------------
    #     clade_weights = defaultdict(float)
    #     for leaf in node.iter_leaves():
    #         clade_weights[self.leaf2clade[leaf.name]] += 1.0
    #     total = sum(clade_weights.values())
    #     if total == 0:
    #         prior_mix = 0.0
    #     else:
    #         prior_mix = sum((w / total) * self.p_frac[c]
    #                         for c, w in clade_weights.items()) 
        
    #     # ------------------------------------------------------------------
    #     # 3️⃣  Log‑likelihood pieces that were stored by pruninglikelihood
    #     # ------------------------------------------------------------------
    #     ############## 528
    #     log_factor = node.log_factor
        
    #     # 🔑 CRITICAL FIX: weight by parent_prior_vec BEFORE taking logsumexp
    #     # We want log(parent_prior · L_vec), not just log(sum(L_vec))
    #     parent_prior_log = np.log(parent_prior_vec + eps)  # shape (K,)
        
    #     # log(parent_prior · exp(logL0 + log_factor))
    #     #  = log(parent_prior) + logL0 + log_factor
    #     log_term0 = logsumexp(parent_prior_log + node.logL0 + log_factor)
        
    #     # log(parent_prior · exp(logL1 + log_factor))
    #     log_term1 = logsumexp(parent_prior_log + node.logL1 + log_factor)
    #     ############## 528


    #     # ############ new 528
    #     # L0 = np.exp(node.logL0 + node.log_factor)
    #     # L1 = np.exp(node.logL1 + node.log_factor)
        
    #     # term0 = np.dot(parent_prior_vec, L0)
    #     # term1 = np.dot(parent_prior_vec, L1)
        
    #     # log_term0 = np.log(term0 + eps)
    #     # log_term1 = np.log(term1 + eps)
    #     # ############ new 528
        
    #     # ------------------------------------------------------------------
    #     # 4️⃣  Return the posterior γ_u using Bayes rule
    #     # ------------------------------------------------------------------
    #     if prior_mix <= 0.0:
    #         return 0.0
    #     if prior_mix >= 1.0:
    #         return 1.0
        
    #     # γ_u = P(Z_u=1 | data)
    #     #     = prior_mix · P(data|Z_u=1) / [prior_mix · P(data|Z_u=1) + (1-prior_mix) · P(data|Z_u=0)]
    #     #     = prior_mix · exp(log_term1) / [prior_mix · exp(log_term1) + (1-prior_mix) · exp(log_term0)]
        
    #     # In log space:
    #     log_num = np.log(prior_mix + eps) + log_term1
        
    #     # logsumexp([log(prior_mix) + log_term1, log(1-prior_mix) + log_term0])
    #     log_den = logsumexp(np.array([
    #         np.log(prior_mix + eps) + log_term1,
    #         np.log(1.0 - prior_mix + eps) + log_term0
    #     ]))
        
    #     gamma = np.exp(log_num - log_den)
    #     return float(np.clip(gamma, 0.0, 1.0))
    # 528

    # 528
    # def node_responsibility(self, node, parent_prior_vec=None):
    #     """
    #     γ_u = P(Z_u = 1 | data) for an internal, non‑root node.
    #     Bug‑fixed version (see notes in the answer).
    #     """
    #     eps = 1e-12
    
    #     # ------------------------------------------------------------------
    #     # 1️⃣  Parent prior – use the *probability* vector directly, but turn
    #     #     zero entries into -inf in log space so they contribute nothing.
    #     # ------------------------------------------------------------------
    #     if parent_prior_vec is None:
    #         if node.up is None:
    #             parent_prior_vec = np.zeros(self.K)
    #             parent_prior_vec[self.root_type] = 1.0
    #         else:
    #             parent_prior_vec = node.up.lik
    
    #     # log of the probability vector (zero entries → -inf, no eps)
    #     with np.errstate(divide='ignore'):
    #         log_parent_prior = np.log(parent_prior_vec)
    #     log_parent_prior[parent_prior_vec == 0] = -np.inf
    
    #     # ------------------------------------------------------------------
    #     # 2️⃣  Per‑node scaling factor (NOT the cumulative one!)
    #     # ------------------------------------------------------------------
    #     # node.logL0 and node.logL1 have already been scaled by the *child*
    #     # offsets during pruning.  The scalar that *must* be added back is the
    #     # offset that was removed when this node’s own L0/L1 vectors were
    #     # generated – stored in node.scale_factor (we create it below).
    #     log_factor = getattr(node, "scale_factor", 0.0)
    
    #     # ------------------------------------------------------------------
    #     # 3️⃣  log( parent_prior · exp(logL_k + log_factor) )
    #     # ------------------------------------------------------------------
    #     log_term0 = logsumexp(log_parent_prior + node.logL0 + log_factor)
    #     log_term1 = logsumexp(log_parent_prior + node.logL1 + log_factor)
    
    #     # ------------------------------------------------------------------
    #     # 4️⃣  Prior mixture weight (the same formula you already use)
    #     # ------------------------------------------------------------------
    #     clade_weights = defaultdict(float)
    #     for leaf in node.iter_leaves():
    #         clade_weights[self.leaf2clade[leaf.name]] += 1.0
    #     total = sum(clade_weights.values())
    #     if total == 0:
    #         prior_mix = 0.0
    #     else:
    #         prior_mix = sum((w / total) * self.p_frac[c]
    #                         for c, w in clade_weights.items())
    
    #     # ------------------------------------------------------------------
    #     # 5️⃣  Bayes rule in log space
    #     # ------------------------------------------------------------------
    #     if prior_mix <= 0.0:
    #         return 0.0
    #     if prior_mix >= 1.0:
    #         return 1.0
    
    #     log_num = np.log(prior_mix) + log_term1
    #     log_den = logsumexp(np.array([
    #         np.log(prior_mix) + log_term1,
    #         np.log(1.0 - prior_mix) + log_term0
    #     ]))
    
    #     gamma = np.exp(log_num - log_den)
    #     return float(np.clip(gamma, 0.0, 1.0))
    # 528

    def node_responsibility(self, node, parent_prior_vec=None):
        eps = 1e-12
    
        # 1️⃣ parent prior (probability vector, zeros stay as zeros)
        if parent_prior_vec is None:
            if node.up is None:
                parent_prior_vec = np.zeros(self.K)
                parent_prior_vec[self.root_type] = 1.0
            else:
                parent_prior_vec = node.up.lik
    
        with np.errstate(divide='ignore'):
            log_parent_prior = np.log(parent_prior_vec)
        log_parent_prior[parent_prior_vec == 0] = -np.inf
    
        # 2️⃣ per‑node scalar
        log_factor = getattr(node, "scale_factor", 0.0)
    
        # 3️⃣ log( parent_prior · L_k ) for each component
        log_term0 = logsumexp(log_parent_prior + np.log(node.L0_vec + eps) + log_factor)
        log_term1 = logsumexp(log_parent_prior + np.log(node.L1_vec + eps) + log_factor)
    
        # 4️⃣ clade‑weighted prior mixture
        clade_weights = defaultdict(float)
        for leaf in node.iter_leaves():
            clade_weights[self.leaf2clade[leaf.name]] += 1.0
        total = sum(clade_weights.values())
        prior_mix = 0.0 if total == 0 else sum(
                        (w / total) * self.p_frac[c] for c, w in clade_weights.items())
    
        if prior_mix <= 0.0:
            return 0.0
        if prior_mix >= 1.0:
            return 1.0
    
        log_num = np.log(prior_mix) + log_term1
        log_den = logsumexp(np.array([log_num,
                                      np.log(1.0 - prior_mix) + log_term0]))
    
        gamma = np.exp(log_num - log_den)
        return float(np.clip(gamma, 0.0, 1.0))

    def leaf_predictions(self,
                     threshold: float = 0.5,
                     use_path_responsibilities: bool = False) -> pd.DataFrame:
        """
        Convert the clade‑level mixture weights (self.p_frac) into leaf‑level
        predictions.
    
        Parameters
        ----------
        threshold : float, default 0.5
            Decision threshold for the hard label.  Leaves with
            ``prob_induced >= threshold`` are labelled “induced”.
        use_path_responsibilities : bool, default False
            If True, incorporate the posterior responsibilities (γ_u) of every
            internal node on the path from the leaf to the root.  The leaf
            probability is then
                P(induced | leaf) = ∏_{u on path} [ p_parent·γ_u + (1-p_parent)·(1-γ_u) ].
            The product is evaluated in log‑space for numerical stability.
            If False (default) the leaf probability is simply the clade mixture
            weight p_c.
    
        Returns
        -------
        pandas.DataFrame with columns:
            - leaf_name      : identifier of the leaf / cell
            - clade          : clade label taken from self.leaf2clade
            - p_clade        : posterior mixture weight for that clade (p_c)
            - prob_induced   : soft probability that the leaf is induced
            - pred_label     : hard label (“induced” or “un‑induced”)
        """
        records = []
    
        for leaf in self.tree.iter_leaves():
            leaf_name = leaf.name
            clade = self.leaf2clade[leaf_name]
            p_c   = float(self.p_frac[clade])      # clade‑level posterior weight
    
            # --------------------------------------------------------------
            # A) Simple clade‑only probability (fast)
            # --------------------------------------------------------------
            prob = p_c
    
            # --------------------------------------------------------------
            # B) Path‑aware probability (if requested)
            # --------------------------------------------------------------
            if use_path_responsibilities:
                # Walk up from leaf to root, accumulating log‑probability.
                node = leaf
                log_prob = np.log(p_c + 1e-12)      # start with the clade weight
                while not node.is_root():
                    parent = node.up
                    # γ_u for the division that produced `node`
                    gamma_u = self.node_responsibility(node)
    
                    # mixture weight for the parent clade (if the parent has a clade)
                    parent_clade = self.leaf2clade.get(parent.name, None)
                    if parent_clade is None:
                        p_parent = 0.5               # neutral prior when unknown
                    else:
                        p_parent = float(self.p_frac[parent_clade])
    
                    # Edge‑level probability that this division was induced:
                    edge_prob = p_parent * gamma_u + (1.0 - p_parent) * (1.0 - gamma_u)
    
                    log_prob += np.log(edge_prob + 1e-12)
                    node = parent
                prob = np.exp(log_prob)          # back to probability space
    
            # --------------------------------------------------------------
            # Hard label
            # --------------------------------------------------------------
            label = "induced" if prob >= threshold else "un‑induced"
    
            records.append({
                "leaf_name": leaf_name,
                "clade": clade,
                "p_clade": p_c,
                "prob_induced": prob,
                "pred_label": label
            })
    
        return pd.DataFrame.from_records(records)
    # ------------------------------------------------------------------
    #  EM algorithm
    # ------------------------------------------------------------------
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
        per‑node contributions are weighted in the M‑step.
        """
        logL_history = []

        # ------------------------------------------------------------------
        #  Helper: posterior responsibility γ_u for a node
        # ------------------------------------------------------------------
        # def node_responsibility(node, parent_prior_vec):
        #     """
        #     Posterior probability that the division represented by ``node`` used the
        #     induced matrix (Z_u = 1).  ``parent_prior_vec`` must be the conditional
        #     likelihood vector coming from the node's parent (i.e. node.up.lik).
        #     The function is **never called for the root**, but we keep a guard
        #     anyway for safety.
        #     """
        #     # Safety check – the root has no Z variable
        #     if node.is_root():
        #         return 0.0
        
        #     # L0_vec / L1_vec are guaranteed to exist for internal nodes
        #     L0_vec = node.L0_vec          # shape (K,)
        #     L1_vec = node.L1_vec          # shape (K,)

        #     if L1_vec is None or L0_vec is None:
        #         print(f"\n❌ ERROR: Node {node.name} has L1_vec={L1_vec is None}, L0_vec={L0_vec is None}")
        #         print(f"   is_leaf={node.is_leaf()}, is_root={node.is_root()}")
        #         print(f"   children={len(node.children)}")
        #         # Print nearby nodes to help debugging
        #         for n in node.up.children if node.up else []:
        #             print(f"   sibling {n.name}: is_leaf={n.is_leaf()}, L1_vec={n.L1_vec is not None}")
        #         raise RuntimeError(f"Node {node.name} missing L1_vec or L0_vec after pruning!")
        
        #     # ----- mixture weight for the clade(s) covered by this node -----
        #     clade_weights = defaultdict(float)
        #     for leaf in node.iter_leaves():
        #         clade_weights[self.leaf2clade[leaf.name]] += 1.0
        #     total = sum(clade_weights.values())
        
        #     prior_mix = 0.0
        #     for c, w in clade_weights.items():
        #         prior_mix += (w / total) * self.p_frac[c]
        
        #     # ----- Bayes rule (Eq. 8′) -----
        #     # print(f'{prior_mix = }')
        #     # print(f'{L1_vec = }')
        #     num = prior_mix * np.dot(parent_prior_vec, L1_vec)
        #     den = (prior_mix * np.dot(parent_prior_vec, L1_vec) +
        #            (1.0 - prior_mix) * np.dot(parent_prior_vec, L0_vec) + 1e-12)
        #     return num / den


        # ####### most current 5/12 but making class method
        # def node_responsibility(node, parent_prior_vec):
        #     """
        #     Posterior probability that the division represented by ``node`` used the
        #     induced matrix (Z_u = 1).  Operates entirely in log‑space, matching the
        #     vectors stored by ``pruninglikelihood``.
        #     """
        #     # ------------------------------------------------------------------
        #     # 0️⃣  The root never has a Z variable
        #     # ------------------------------------------------------------------
        #     if node.is_root():
        #         return 0.0
        
        #     # ------------------------------------------------------------------
        #     # 1️⃣  Compute the *clade‑specific* mixture weight for this node
        #     # ------------------------------------------------------------------
        #     clade_weights = defaultdict(float)
        #     for leaf in node.iter_leaves():
        #         # leaf → clade mapping lives on the EM object (self.leaf2clade)
        #         clade_weights[self.leaf2clade[leaf.name]] += 1.0
        
        #     total = sum(clade_weights.values())
        #     if total == 0:                # should never happen for an internal node
        #         prior_mix = 0.0
        #     else:
        #         # weighted average of the current p_frac values
        #         prior_mix = sum((w / total) * self.p_frac[c]
        #                         for c, w in clade_weights.items())
        
        #     # ------------------------------------------------------------------
        #     # 2️⃣  Log‑space “dot‑products” π ⋅ ℓ₀ and π ⋅ ℓ₁
        #     # ------------------------------------------------------------------
        #     # `node.logL0` and `node.logL1` are the *scaled* log‑likelihood vectors.
        #     # The common scaling factor that was subtracted during pruning is stored
        #     # in `node.log_factor`.  Adding it back gives the true log‑vectors.
        #     log_factor = node.log_factor
        
        #     # log(π ⋅ ℓ₁)  =  logsumexp( logπ + logℓ₁ )  =  logsumexp( node.logL1 + log_factor )
        #     log_pi_ell1 = logsumexp(node.logL1 + log_factor)
        
        #     # log(π ⋅ ℓ₀)  =  logsumexp( node.logL0 + log_factor )
        #     log_pi_ell0 = logsumexp(node.logL0 + log_factor)
        
        #     # ------------------------------------------------------------------
        #     # 3️⃣  Build the log‑numerator and log‑denominator of γ_u
        #     # ------------------------------------------------------------------
        #     # Guard the extreme cases where the mixture weight is exactly 0 or 1.
        #     if prior_mix <= 0.0:
        #         return 0.0
        #     if prior_mix >= 1.0:
        #         return 1.0
        
        #     log_num = np.log(prior_mix) + log_pi_ell1
        
        #     # Two additive terms that must be summed in log‑space:
        #     term1 = np.log(prior_mix)      + log_pi_ell1
        #     term2 = np.log(1.0 - prior_mix) + log_pi_ell0
        #     log_den = logsumexp(np.array([term1, term2]))
        
        #     # ------------------------------------------------------------------
        #     # 4️⃣  Return the responsibility γ_u ∈ [0, 1]
        #     # ------------------------------------------------------------------
        #     gamma = np.exp(log_num - log_den)
        #     # Clip only to protect against tiny floating‑point round‑off errors.
        #     return float(np.clip(gamma, 0.0, 1.0))

        # ####### most current 5/12 but making class method

        # ------------------------------------------------------------------
        #  Choose the weight‑lookup function once (so we do not test strings
        #  inside the inner loop).
        # ------------------------------------------------------------------
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

        iteration_lr_val_dict = {}
    
        for it in tqdm.trange(max_iter, disable=not verbose):
            # ------------------- E‑step -------------------
            # **IMPORTANT** – call the *correct* pruning routine
            # logL, lr_vals = self.pruninglikelihood()          # <-- fixed line
            logL = self.pruninglikelihood()
            # iteration_lr_val_dict[it] = lr_vals
            logL_history.append(logL)
            
    
            # (optional) sanity check that vectors exist
            # for n in self.tree.traverse():
            #     if not n.is_leaf() and not n.is_root():
            #         assert n.L0_vec is not None and n.L1_vec is not None
    
            if verbose:
                logging.info(f"EM iteration {it+1:02d} – logL = {logL:.6f}")
    
            # ------------------- Gather sufficient statistics -------------------
            clade_stats = {c: [0.0, 0.0] for c in self.clades}   # [numer, denom]
    
            # Root prior vector (used as parent prior for the root's children)
            root_prior = np.zeros(self.K)
            root_prior[self.root_type] = 1.0


            for node in self.tree.traverse("postorder"):
                if node.is_root() or node.is_leaf(): # added is_leaf() check
                    continue                       # root has no Z_u

                clade_weights = defaultdict(float)
                for leaf in node.iter_leaves():
                    clade_weights[self.leaf2clade[leaf.name]] += 1.0
                # if len(clade_weights) > 1:
                #     print(f"  Node {node.name}: split across {len(clade_weights)} clades: {dict(clade_weights)}")
                # elif len(clade_weights) == 1:
                #     print(f"  Node {node.name}: has clade_lengths == 1")
    
                # parent prior comes from the already‑computed node.up.lik
                parent_vec = node.up.lik if node.up else root_prior # temp comment 528
                # parent_vec = root_prior
    
                # ------------------- posterior responsibility -------------------
                gamma_u = self.node_responsibility(node, parent_vec)   # uses L0_vec/L1_vec
    
                # ------------------- node weight according to strategy -------
                w_dict = weight_func(node)        # e.g. {'cladeA': 5.0, ...}
                # print(f'{w_dict = }')
    
                # accumulate numerator / denominator for each clade
                for c, w in w_dict.items():
                    clade_stats[c][0] += w * gamma_u    # Σ w·γ
                    clade_stats[c][1] += w             # Σ w
    
            # ------------------- M‑step (MAP update) -------------------
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
    
            # ------------------- convergence check -------------------
            logL_change = abs(logL_history[-1] - logL_history[-2]) if it > 0 else np.inf

            if verbose:
                if it > 0 and logL_change < tol:
                    logging.info(f"EM converged (log‑likelihood change {logL_change:.2e} < {tol})")
                    break
            if it > 0 and abs(logL_history[-1] - logL_history[-2]) < tol:
                if verbose:
                    logging.info("EM converged (log‑likelihood change < tol).")
                break


        # final_logL = logL_history[-1]
        # improvement = final_logL - baseline_logL
        # print(f"\nFinal logL:     {final_logL:.8f}")
        # print(f"Improvement:    {improvement:.8f}")
        # if improvement > 0:
        #     print("✓ Model improved over baseline")
        # else:
        #     print("⚠️  Model did NOT improve over baseline")
        self.logL_history_ = logL_history
        return logL_history, iteration_lr_val_dict
    # ------------------------------------------------------------------
    #  Post‑processing utilities
    # ------------------------------------------------------------------
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
        print(df.to_string(index=False))
        if hasattr(self, "logL_history_"):
            pass

    def diagnose_gamma_by_clade(self):
        """
        Inspect per-node γ_u values grouped by clade.
        Should be called AFTER fit() to see which clades prefer which matrix.
        """
        
        # First, recompute the pruning pass to populate L0_vec, L1_vec
        self.pruninglikelihood()
        
        root_prior = np.zeros(self.K)
        # root_prior[self.type_encoder[self.root_type]] = 1.0
        root_prior[self.root_type] = 1.0
        
        # Collect γ values by clade
        gamma_by_clade = defaultdict(list)
        
        for node in self.tree.traverse("postorder"):
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
                clade = max(clade_weights, key=clade_weights.get)
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
                    for node_name, gamma in sorted(gammas, key=lambda x: x[1]):
                        pass
        
        # Overall statistics
        all_gammas = [g for gammas in gamma_by_clade.values() for _, g in gammas]
        
        return gamma_by_clade



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
    
# ------------------------------------------------------------
# 3️⃣  Convert mean to Beta parameters (strength = pseudo‑samples)
# ------------------------------------------------------------
def beta_parameters(prior_mean, strength=10):
    alpha = prior_mean * strength
    beta  = (1.0 - prior_mean) * strength
    return alpha, beta

def urid_to_induction_inference(urid_list, timepoints, param_dict_list, weight_strategy = 'leaf_prop',
                               clade_priors = None, upper_rarity_bound = 2,
                               include_param_combos = None, scoremat_filters = None,
                               verbose = False,
                               root_type = 'ct1', num_restarts = 1,
                               return_opt_model_dict = False, type_encoder = None,
                               fit_tol = 1e-6, enrichment_dict = None, strength = None,
                               use_gmm = False, p_frac_init = None, use_cp_class = True,
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
                    tree = read_in_tree(urid_timepoint_res_dict[urid][timepoint].gt_tree_path)
                    keep_nodes = [node for node in tree.traverse('preorder') if node.name in keep_linstrings]
                    tree.prune(keep_nodes, preserve_branch_length = True)
                    # print(f'pre prune, {filt_df.shape=}')
                    filt_df = full_df[(full_df_copy['cell'].isin(keep_linstrings)) & (full_df_copy['param_combo'] == param_combo)]
                    
                
                else:


            
                    tree = [res.recon_tree for res in urid_timepoint_res_dict[urid][timepoint].scoremat_res_list 
                                                 if res.short_name == ex_short_name][0]
                
                
    
                    # print(f'{param_combo}_{ex_short_name}')
                
                    # only look at full_df clade/cell type rows that match this short_name/param_combo combination
                    filt_df = full_df_copy[full_df_copy['param_combo'] == param_combo]

                if type_encoder is None:
                    type_encoder = dict(zip(sorted(filt_df['celltype'].unique()), range(filt_df['celltype'].nunique())))
                
                leaf2clade = dict(zip(filt_df['cell'], filt_df['clade']))
                leaf2type = dict(zip(filt_df['cell'], filt_df['celltype']))

                # print(f'{filt_df.shape = }')
                # print(f'{len(leaf2clade) = }')
                # print(f'{len(leaf2type) = }')

                # num_leaves = len(list(tree.iter_leaves()))
                
                
                # em = TreeMixtureEM(
                #                  Q0 = ex_uninduced_tm,
                #                  Q1 = ex_induced_tm,
                #                  root_type = 'ct1',
                #                  tree = tree,
                #                  leaf2clade = leaf2clade,
                #                  leaf2type = leaf2type,
                #                  type_encoder=type_encoder)

                # em = TreeMixtureNodeEM(Phi0 = ex_uninduced_tm,
                #                      Phi1 = ex_induced_tm,
                #                      root_type = 'ct1',
                #                      tree = tree,
                #                      leaf2clade = leaf2clade,
                #                      leaf2type = leaf2type,
                #                      type_encoder=type_encoder,
                #                       weight_strategy = weight_strategy)

                # if param_combo not in ['20clade_topological_treedist_gt', 'clade_differentiation_induced']:
                # if param_combo not in ['20clade_topological_treedist_gt']:
                if include_param_combos is not None:
                    if param_combo not in include_param_combos:
                        continue
                if scoremat_filters is not None:
                    # gt combos (ending _gt) always pass; scoremat combos must
                    # contain at least one of the filter strings (e.g. '1.0RP_0.0AF')
                    if not param_combo.endswith('_gt'):
                        if not any(f in param_combo for f in scoremat_filters):
                            continue

                if not use_cp_class:
                    em = og_TreeMixtureNodeEM2(Phi0 = ex_uninduced_tm,
                                         Phi1 = ex_induced_tm,
                                         root_type = root_type,
                                         tree = tree,
                                         leaf2clade = leaf2clade,
                                         leaf2type = leaf2type,
                                         type_encoder = type_encoder,
                                          weight_strategy = weight_strategy)
                else:
                    em = cp_TreeMixtureNodeEM2(Phi0 = ex_uninduced_tm,
                                         Phi1 = ex_induced_tm,
                                         root_type = root_type,
                                         tree = tree,
                                         leaf2clade = leaf2clade,
                                         leaf2type = leaf2type,
                                         type_encoder = type_encoder,
                                          weight_strategy = weight_strategy)
                

                # em = new514_TreeMixtureEM(Phi0 = ex_uninduced_tm,
                #                         Phi1 = ex_induced_tm,
                #                         root_type = root_type,
                #                         tree = tree,
                #                         leaf2clade = leaf2clade,
                #                         leaf2type = leaf2type,
                #                         type_encoder= type_encoder)


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
    
                            
                            # clade_prior_param_dict[clade] = beta_params



                    ####################
                    # print(f'clade_rarity dict = {clade_rarity}')
                    ####################
    
                    
                    # # ============================================================
                    # # DIAGNOSTIC 1: Test if p_c affects likelihood
                    # # ============================================================
                    # print("\n" + "="*60)
                    # print("DIAGNOSTIC 1: Testing if p_c affects logL")
                    # print("="*60)
                    
                    # original_p_frac = {c: em.p_frac[c] for c in em.clades}
                    
                    # # Set all p_c = 0.0 (all un-induced)
                    # for c in em.clades:
                    #     em.p_frac[c] = 0.0
                    # # em.update_mixture()
                    # logL_all_uninduced = em.pruninglikelihood()
                    # print(f"logL (all p_c=0.0): {logL_all_uninduced:.10f}")
                    
                    # # Set all p_c = 1.0 (all induced)
                    # for c in em.clades:
                    #     em.p_frac[c] = 1.0
                    # # em.update_mixture()
                    # logL_all_induced = em.pruninglikelihood()
                    # print(f"logL (all p_c=1.0): {logL_all_induced:.10f}")
                    
                    # # Set all p_c = 0.5 (mixed)
                    # for c in em.clades:
                    #     em.p_frac[c] = 0.5
                    # # em.update_mixture()
                    # logL_mixed = em.pruninglikelihood()
                    # print(f"logL (all p_c=0.5): {logL_mixed:.10f}")
                    
                    # # Restore original
                    # for c in em.clades:
                    #     em.p_frac[c] = original_p_frac[c]
                    
                    # print(f"\nDifference (p_c=0 vs p_c=1): {abs(logL_all_induced - logL_all_uninduced):.10f}")
                    
                    # if logL_all_uninduced == logL_all_induced == logL_mixed:
                    #     print("❌ PROBLEM: logL does NOT change with p_c!")
                    #     print("   The mixture parameters are not affecting the likelihood.")
                    # else:
                    #     print("✓ logL varies with p_c (em is working)")
                    
                    # print("="*60)
                    # print()


                    ############5/7
                    # logL_hist = em.fit(max_iter=200, tol=1e-6, verbose=True,
                    #                   prior_alpha_dict = clade_alpha,
                    #                   prior_beta_dict = clade_beta)
                    ############5/7

                assert np.allclose(np.array(ex_uninduced_tm).sum(axis=1), 1.0, atol=1e-3), 'uninduced should be transposed'
                assert np.allclose(np.array(ex_induced_tm).sum(axis=1), 1.0, atol=1e-3), 'induced should be transposed'

                results, best, init_p = run_em_with_restarts(
                    Phi0=ex_uninduced_tm,
                    Phi1=ex_induced_tm,
                    root_type=root_type,
                    tree=tree,
                    leaf2clade=leaf2clade,
                    leaf2type=leaf2type,
                    type_encoder=type_encoder,
                    weight_strategy=weight_strategy,
                    clade_alpha=clade_alpha,
                    clade_beta=clade_beta,
                    n_restarts=num_restarts,
                    max_iter=200,
                    tol=1e-3,
                    seed_start=42,
                    verbose = verbose,
                    p_frac_init = p_frac_init,
                    use_cp_class = use_cp_class
                )
            
                # Access the best model
                best_model = best['model']
                
                ####################
                # # Access individual restart results
                # for i, result in enumerate(results):
                #     print(f"\nRestart {i+1}:")
                #     print(result['clade_posteriors'])
                ####################
                
                em = best_model
                em.results = results
                # else:
                #     logL_hist = em.fit(max_iter=200, tol=1e-3, verbose=verbose)
                
                
                # em.summary()          # prints clade table and final log‑likelihood
                df_post = em.clade_posteriors(use_gmm = use_gmm)
                
                
                ####################
                # print(f'initial_p estimates = {init_p}')

                # print("\nFinal clade estimates:")
                # print(df)
                
                # # Print raw p_frac values
                # print("\nRaw p_frac values:")
                # for c in em.clades:
                #     print(f"  {c}: {em.p_frac[c]:.6f}")
                ####################
                
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

            ####################
            # print(f'{this_em.leaf2clade = }')
            # print(f'{gt_linstring_to_induced_dict = }')
            ####################

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
             
def plot_overlap_fracs(full_df, param_combo1, param_combo2):
    sub_df = full_df[full_df['param_combo'].isin([param_combo1, param_combo2])]

    pivot_sub = sub_df.pivot(index = 'cell',
                             columns = 'param_combo',
                             values = 'clade')

    df_1_to_2 = pd.crosstab(
        pivot_sub[param_combo1],
        pivot_sub[param_combo2],
        normalize = 'index')

    plt.figure(figsize = (12, 12))
    sns.heatmap(df_1_to_2, cmap = 'Blues', annot_kws = {'color': 'orange', 'fontsize': 6}, annot = True)
    plt.title(f'Fraction of {param_combo1} in each {param_combo2}')
    plt.close('all')

    df_2_to_1 = pd.crosstab(
        pivot_sub[param_combo2],
        pivot_sub[param_combo1],
        normalize = 'index')

    plt.figure(figsize = (12, 12))
    sns.heatmap(df_2_to_1, cmap = 'Blues', annot_kws = {'color': 'orange', 'rotation': 90, 'fontsize': 6}, annot = True)
    plt.title(f'Fraction of {param_combo2} in each {param_combo1}')
    plt.close('all')

    return df_1_to_2, df_2_to_1


def initialize_p_frac_from_rarity(leaf2type_int, leaf2clade, Phi0, Phi1, clades, eps=1e-12):
    """
    Initialize p_frac based on the rarity of cell types in each clade.
    
    High rarity → high p_c (clade likely induced)
    Low rarity → low p_c (clade likely un-induced)
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
        # Rarity of 1.0 → p_c = 0.5 (no signal)
        # Rarity of 5.0 → p_c = 0.9 (strong signal for induced)
        # Rarity of 0.2 → p_c = 0.1 (strong signal for un-induced)
        
        # Using a logistic function: p_c = 1 / (1 + (1/avg_rarity)^scale)
        # Or simpler: p_c = avg_rarity / (1 + avg_rarity)
        
        p_frac_init[c] = avg_rarity / (1.0 + avg_rarity)  # maps [0, inf) → [0, 1)
        
        # Clip to avoid extremes
        p_frac_init[c] = np.clip(p_frac_init[c], 0.1, 0.9)
    
    return p_frac_init

def evaluate_prediction_probs(em, leaf2clade, leaf2induced_gt, res_df, plot = True):
    cell_names = list(leaf2clade.keys())
    ground_truth_str = np.array([leaf2induced_gt[cell] for cell in cell_names])
    ground_truth = np.array([1 if label == 'induced' else 0 for label in ground_truth_str])

    # print('in evaluate pred probs')
    # print(f'{predicted_probs = }')
    # print(f'{ground_truth = }')

    

    

    # predicted_label = np.array(['induced' if prob > 0.5 else 'un-induced' for prob in predicted_probs])

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
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        
        # PR curve
        axes[0].plot(recall, precision, label=f'PR curve (AUC={pr_auc:.4f})')
        axes[0].fill_between(recall, precision, alpha=0.2)
        axes[0].set_xlabel('Recall')
        axes[0].set_ylabel('Precision')
        axes[0].set_title('Precision-Recall Curve')
        axes[0].legend()
        axes[0].grid(alpha=0.3)
        axes[0].set_xlim([0, 1])
        axes[0].set_ylim([0, 1])
        
        # ROC curve
        axes[1].plot(fpr, tpr, label=f'ROC curve (AUC={roc_auc:.4f})')
        axes[1].fill_between(fpr, tpr, alpha=0.2)
        axes[1].plot([0, 1], [0, 1], 'k--', label='Random classifier')
        axes[1].set_xlabel('False Positive Rate')
        axes[1].set_ylabel('True Positive Rate')
        axes[1].set_title('ROC Curve')
        axes[1].legend()
        axes[1].grid(alpha=0.3)
        axes[1].set_xlim([0, 1])
        axes[1].set_ylim([0, 1])
        
        plt.tight_layout()
        plt.savefig('evaluation_curves.png', dpi=150, bbox_inches='tight')
        plt.close('all')

        # Add this to your code:
    
    # Histogram
    n_bins = min(50, len(np.unique(predicted_probs)))
    if n_bins > 1:
        plt.hist(predicted_probs, bins=n_bins)
    else:
        plt.text(0.5, 0.5, f"All probs = {predicted_probs[0]:.4f}", ha="center", transform=plt.gca().transAxes)
    plt.xlabel('Predicted Probability')
    plt.ylabel('Count')
    plt.close('all')

    return metrics, (precision, recall), (fpr, tpr), (predicted_probs, predicted_label, ground_truth)

def run_em_with_restarts(Phi0, Phi1, root_type, tree, leaf2clade, leaf2type,
                         type_encoder, weight_strategy, clade_alpha, clade_beta,
                         n_restarts=5, max_iter=200, tol=1e-3, seed_start=42,
                         verbose=False, p_frac_init=None, use_cp_class=True,
                         alpha_values=None, Phi1_true=None):
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
            raise ValueError("Phi1_true must be supplied when alpha_values is given")
        alpha_results = {}
        for alpha in alpha_values:
            Phi1_input = (1.0 - alpha) * Phi0 + alpha * np.asarray(Phi1_true, dtype=float)
            _, best_result, _ = run_em_with_restarts(
                Phi0=Phi0, Phi1=Phi1_input,
                root_type=root_type, tree=tree,
                leaf2clade=leaf2clade, leaf2type=leaf2type,
                type_encoder=type_encoder, weight_strategy=weight_strategy,
                clade_alpha=clade_alpha, clade_beta=clade_beta,
                n_restarts=n_restarts, max_iter=max_iter, tol=tol,
                seed_start=seed_start, verbose=verbose,
                p_frac_init=p_frac_init, use_cp_class=use_cp_class,
            )
            alpha_results[alpha] = best_result
        return alpha_results

    results = []

    for restart_idx in range(n_restarts):
        seed = seed_start + restart_idx
        np.random.seed(seed)


        ####################
        # print(f"\n{'='*70}")
        # print(f"RESTART {restart_idx + 1}/{n_restarts} (seed={seed})")
        # print(f"{'='*70}")
        ####################

        if not use_cp_class:
            # # Initialize with random p_frac values
            em = og_TreeMixtureNodeEM2(
                Phi0=Phi0,
                Phi1=Phi1,
                root_type=root_type,
                tree=tree,
                leaf2clade=leaf2clade,
                leaf2type=leaf2type,
                type_encoder=type_encoder,
                weight_strategy=weight_strategy
            )

        else:
            # # Initialize with random p_frac values
            em = cp_TreeMixtureNodeEM2(
                Phi0=Phi0,
                Phi1=Phi1,
                root_type=root_type,
                tree=tree,
                leaf2clade=leaf2clade,
                leaf2type=leaf2type,
                type_encoder=type_encoder,
                weight_strategy=weight_strategy
            )

        # em = new514_TreeMixtureEM(Phi0 = Phi0,
        #                                 Phi1 = Phi1,
        #                                 root_type = root_type,
        #                                 tree = tree,
        #                                 leaf2clade = leaf2clade,
        #                                 leaf2type = leaf2type,
        #                                 type_encoder= type_encoder)

        ############5/7
        # em = TreeMixtureNodeEM3(
        #     Phi0=Phi0,
        #     Phi1=Phi1,
        #     root_type=root_type,
        #     tree=tree,
        #     leaf2clade=leaf2clade,
        #     leaf2type=leaf2type,
        #     type_encoder=type_encoder,
        #     weight_strategy=weight_strategy
        # )
        ############5/7
        
        # # Random initialization of p_frac
        # for c in em.clades:
        #     em.p_frac[c] = np.random.uniform(0.1, 0.9)
        # print(f"Initial p_frac (random):")
        # for c in em.clades[:3]:  # show first 3
        #     print(f"  {c}: {em.p_frac[c]:.4f}")

        # p_frac_init = initialize_p_frac_from_rarity(leaf2type_int = em.leaf2type_int,
        #                                             leaf2clade = em.leaf2clade,
        #                                             Phi0 = em.Phi0,
        #                                             Phi1 = em.Phi1,
        #                                             clades = em.clades)

        # Compute base init from priors (once per call, not once per restart)
        if restart_idx == 0:
            _base_init = {c: clade_alpha[c] / (clade_alpha[c] + clade_beta[c]) for c in em.clades}
            _base_init = {c: np.clip(p, 0.1, 0.9) for c, p in _base_init.items()}

        if p_frac_init is not None:
            # Caller supplied explicit init — use as-is for all restarts
            _init = p_frac_init
        elif restart_idx == 0:
            # First restart: use prior mean (deterministic anchor)
            _init = _base_init
        else:
            # Subsequent restarts: random perturbation around prior mean
            _init = {c: float(np.clip(
                        _base_init[c] + np.random.uniform(-0.4, 0.4), 0.05, 0.95))
                     for c in em.clades}

        for c, init_val in _init.items():
            em.p_frac[c] = init_val

         


        iteration_lr_val_dict = {}
        
        # Run EM
        # logL_hist, iteration_lr_val_dict = em.fit(
        logL_hist = em.fit(
            max_iter=max_iter,
            tol=tol,
            verbose=verbose,  # suppress iteration details to reduce output
            prior_alpha_dict=clade_alpha,
            prior_beta_dict=clade_beta,
            # reset_to_baseline = False
        )

        for c in em.clades:
            pass

        
        final_logL = logL_hist[-1]

        n_iterations = len(logL_hist)

        ####################
        # print(f"\nResults:")
        # print(f"  Final logL: {final_logL:.8f}")
        # print(f"  Iterations: {n_iterations}")
        # print(f"  Converged: {'Yes' if n_iterations < max_iter else 'No'}")
        ####################
        
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

    ####################
    # print(f"\n{'='*70}")
    # print(f"SUMMARY: Best solution from restart {best_idx + 1} (seed={best_result['seed']})")
    # print(f"Final logL: {best_result['final_logL']:.8f}")
    # print(f"{'='*70}\n")
    
    # print("Best solution posteriors:")
    # print(best_result['clade_posteriors'].to_string())
    
    # # Compare all restarts
    # print(f"\n{'='*70}")
    # print("Comparison of all restarts:")
    # print(f"{'='*70}")
    # for r in results:
    #     n_induced = sum(1 for c in r['p_frac'] if r['p_frac'][c] > 0.5)
    #     print(f"Restart {r['restart']+1} (seed={r['seed']:3d}): "
    #           f"logL={r['final_logL']:10.6f}, "
    #           f"iterations={r['n_iterations']:3d}, "
    #           f"n_induced={n_induced}")
    ####################
    
    return results, best_result, p_frac_init


def js_between_tms(tm0, tm1):

    if isinstance(tm0, list):
        tm0 = np.array(tm0)
    if isinstance(tm1, list):
        tm1 = np.array(tm1)
    
    if not tm0.shape == tm1.shape:
        raise Exception('transition matrices have different dimensions')

    js_dists = []
    for i in range(tm0.shape[0]):
        js_dist = jensenshannon(tm0[i, :], tm1[i, :])
        js_dists.append(js_dist)

    return np.mean(js_dists)

def matrix_diff_rare_types(tm0, tm1, rare_types = None):

    if isinstance(tm0, list):
        tm0 = np.array(tm0)
    if isinstance(tm1, list):
        tm1 = np.array(tm1)
    
    # rare types are those with lower probability in tm0
    if rare_types is None:
        mean_prob_tm0 = tm0.mean(axis=0)
        rare_types = np.where(mean_prob_tm0 < np.median(mean_prob_tm0))[0]
    
    K = tm0.shape[0]
    
    # For each parent state, compute average boost in rare-type probability
    boosts = []
    for i in range(K):
        boost = np.mean([tm1[i, j] - tm0[i, j] for j in rare_types])
        boosts.append(boost)
    
    return {
        'per_state_boost': np.array(boosts),
        'mean_boost': np.mean(boosts),
        'max_boost': np.max(boosts),
        'rare_types': rare_types,
    }


import numpy as np
import pandas as pd
from scipy.special import logsumexp
from collections import defaultdict
import tqdm


class new514_TreeMixtureEM():
    """
    EM for a binary mixture of two transition matrices on a rooted tree.

    Each internal node u has a latent variable:
        Z_u ~ Bernoulli(p_clade)

    Z_u determines whether transitions below that node use:
        Phi0 or Phi1
    """

    def __init__(
        self,
        Phi0,
        Phi1,
        root_type,
        tree,
        leaf2clade,
        leaf2type,
        type_encoder=None,
    ):

        self.tree = tree
        self.leaf2clade = leaf2clade
        self.leaf2type = leaf2type

        self.Phi0 = np.asarray(Phi0, dtype=float)
        self.Phi1 = np.asarray(Phi1, dtype=float)

        self.K = self.Phi0.shape[0]

        if type_encoder is None:
            uniq = sorted(set(leaf2type.values()))
            self.type_encoder = {x: i for i, x in enumerate(uniq)}
        else:
            self.type_encoder = type_encoder

        self.leaf2type_int = {
            k: self.type_encoder[v]
            for k, v in leaf2type.items()
        }

        if isinstance(root_type, str):
            root_type = self.type_encoder[root_type]

        self.root_type = int(root_type)

        self.clades = sorted(set(leaf2clade.values()))

        # initialize clade mixture probs
        self.p_frac = {
            c: 0.5
            for c in self.clades
        }

        self._assign_node_clades()

    # ============================================================
    # TREE HELPERS
    # ============================================================

    def _assign_node_clades(self):

        for node in self.tree.traverse("postorder"):

            if node.is_leaf():
                node.clade_set = {
                    self.leaf2clade[node.name]
                }

            else:
                node.clade_set = set()

                for ch in node.children:
                    node.clade_set.update(ch.clade_set)

    def _node_prior(self, node):

        counts = defaultdict(float)

        for leaf in node.iter_leaves():
            counts[self.leaf2clade[leaf.name]] += 1

        total = sum(counts.values())

        if total == 0:
            return 0.5

        p = 0.0

        for c, w in counts.items():
            p += (w / total) * self.p_frac[c]

        return np.clip(p, 1e-8, 1 - 1e-8)

    # ============================================================
    # UPWARD PASS
    # ============================================================

    def upward_pass(self): 

        K = self.K

        for node in self.tree.traverse("postorder"):

            # ---------------------------------------------
            # LEAF
            # ---------------------------------------------
            if node.is_leaf():

                obs = self.leaf2type_int[node.name]

                logL = np.full(K, -np.inf)
                logL[obs] = 0.0

                node.log_up = logL

                continue

            # ---------------------------------------------
            # INTERNAL NODE
            # ---------------------------------------------

            logL0 = np.zeros(K)
            logL1 = np.zeros(K)

            for child in node.children:

                child_log = child.log_up

                # shape: (K_parent, K_child)
                term0 = np.log(self.Phi0 + 1e-300) + child_log[None, :]
                term1 = np.log(self.Phi1 + 1e-300) + child_log[None, :]

                contrib0 = logsumexp(term0, axis=1)
                contrib1 = logsumexp(term1, axis=1)

                logL0 += contrib0
                logL1 += contrib1

            p = self._node_prior(node)

            node.logL0 = logL0
            node.logL1 = logL1

            node.log_up = logsumexp(
                np.vstack([
                    np.log(1 - p) + logL0,
                    np.log(p) + logL1
                ]),
                axis=0
            )

    # ============================================================
    # TREE LOG LIKELIHOOD
    # ============================================================

    def loglikelihood(self):

        self.upward_pass()

        root = self.tree.get_tree_root()

        root_prior = np.full(self.K, -np.inf)
        root_prior[self.root_type] = 0.0

        return logsumexp(root_prior + root.log_up)

    # ============================================================
    # DOWNWARD PASS
    # ============================================================

    def downward_pass(self):

        root = self.tree.get_tree_root()

        root.log_down = np.full(self.K, -np.inf)
        root.log_down[self.root_type] = 0.0

        for node in root.traverse("preorder"):

            if node.is_leaf():
                continue

            for child in node.children:

                siblings = [
                    s for s in node.children
                    if s != child
                ]

                log_msg = np.copy(node.log_down)

                for sib in siblings:

                    term = (
                        np.log(self.Phi0 + 1e-300)
                        + sib.log_up[None, :]
                    )

                    log_msg += logsumexp(term, axis=1)

                child.log_down = log_msg

    # ============================================================
    # RESPONSIBILITIES
    # ============================================================

    def compute_responsibilities(self):

        self.downward_pass()

        responsibilities = {}

        for node in self.tree.traverse():

            if node.is_root() or node.is_leaf():
                continue

            p = self._node_prior(node)

            log_num1 = logsumexp(
                node.log_down
                + np.log(p)
                + node.logL1
            )

            log_num0 = logsumexp(
                node.log_down
                + np.log(1 - p)
                + node.logL0
            )

            gamma = np.exp(
                log_num1
                - logsumexp([log_num0, log_num1])
            )

            responsibilities[node] = gamma

        return responsibilities

    # ============================================================
    # EM
    # ============================================================

    def fit(
        self,
        max_iter=100,
        tol=1e-4,
        verbose=True,
        prior_alpha=1.0,
        prior_beta=1.0,
    ):

        logL_history = []

        for it in tqdm.trange(max_iter):

            # ---------------------------------------------
            # E STEP
            # ---------------------------------------------

            logL = self.loglikelihood()

            responsibilities = self.compute_responsibilities()

            logL_history.append(logL)

            # ---------------------------------------------
            # M STEP
            # ---------------------------------------------

            numer = defaultdict(float)
            denom = defaultdict(float)

            for node, gamma in responsibilities.items():

                for c in node.clade_set:

                    numer[c] += gamma
                    denom[c] += 1.0

            for c in self.clades:

                a = prior_alpha
                b = prior_beta

                self.p_frac[c] = (
                    numer[c] + a - 1
                ) / (
                    denom[c] + a + b - 2
                )

                self.p_frac[c] = np.clip(
                    self.p_frac[c],
                    1e-5,
                    1 - 1e-5
                )

            # ---------------------------------------------
            # convergence
            # ---------------------------------------------

            if verbose:


                for c in sorted(list(self.clades), key = lambda x: int(x.split('_')[1])):
                    pass

            if it > 0:

                delta = abs(
                    logL_history[-1]
                    - logL_history[-2]
                )

                if delta < tol:
                    break

        self.logL_history_ = logL_history

        return logL_history

    # ============================================================
    # OUTPUT
    # ============================================================

    def clade_posteriors(self):

        rows = []

        for c in self.clades:

            p = self.p_frac[c]

            rows.append({
                "clade": c,
                "p_induced": p,
            })

        return pd.DataFrame(rows)


def pc_df_to_label_dict(pc_df):

    p = pc_df['p_induced'].values
    
    eps = 1e-6
    x = logit(np.clip(p, eps, 1-eps)).reshape(-1,1)
    
    gmm = GaussianMixture(n_components=2, random_state=0)
    gmm.fit(x)
    
    labels = gmm.predict(x)

    # determine which component is induced and which is uninduced. induced assumed to have greater mean
    means = gmm.means_.flatten()
    induced_component = np.argmax(means)


    for i in range(2):
        vals = p[labels == i]

    str_labels = ['induced' if label == induced_component else 'un-induced' for label in labels]
    
    # probs = gmm.predict_proba(x)
    # str_labels = ['induced' if label == 1 else 'uninduced' for label in labels]

    cp_df = copy.deepcopy(pc_df)
    cp_df.insert(loc = pc_df.shape[1], column = 'lab', value = str_labels)
    cp_df.insert(loc = pc_df.shape[1], column = 'x', value = x)
    
    
    res_dict = dict(zip(cp_df['clade'], str_labels))

    return res_dict, cp_df
    


from collections import defaultdict
import numpy as np

def clade_scores_percentile(model, percentile=95):

    gamma_by_clade = defaultdict(list)

    root_prior = np.zeros(model.K)
    root_prior[model.root_type] = 1.0

    for node in model.tree.traverse("postorder"):

        if node.is_leaf() or node.is_root():
            continue

        parent_vec = node.up.lik if node.up else root_prior

        gamma = model.node_responsibility(node, parent_vec)

        clade = max(
            set(model.leaf2clade[l.name] for l in node.iter_leaves()),
            key=lambda c: sum(
                1 for l in node.iter_leaves()
                if model.leaf2clade[l.name] == c
            )
        )

        gamma_by_clade[clade].append(gamma)

    rows = []

    for clade, gammas in gamma_by_clade.items():

        score = np.percentile(gammas, percentile)

        rows.append({
            'clade': clade,
            'score': score
        })

    return (
        pd.DataFrame(rows)
        .sort_values('score', ascending=False)
    )



class cp_TreeMixtureNodeEM2():
    """
    EM for a node‑only mixture of two K‑state stochastic matrices on a rooted
    phylogeny.  The *likelihood‑ratio* weighting scheme is implemented as one
    of the possible ``weight_strategy`` options.
    """

    # ------------------------------------------------------------------
    #  Constructor
    # ------------------------------------------------------------------
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
                'leaf_prop'   – leaf‑proportion weighting (original code)
                'uniform'     – uniform‑per‑node
                'inv_leaf'    – inverse‑leaf weighting
                
                
                'lr+uniform'  – LR multiplied by uniform (default for LR)
                'lr+inv'      – LR multiplied by inverse‑leaf
        """
        # --------------------------------------------------------------
        # 1️⃣  Store tree and attach leaf metadata
        # --------------------------------------------------------------
        self.tree = tree
        # self._check_tree_is_rooted()
        self.leaf2clade = leaf2clade
        self.leaf2type = leaf2type
        self.weight_strategy = weight_strategy

        leaf_names = {leaf.name for leaf in self.tree.iter_leaves()}
        missing = leaf_names - set(self.leaf2clade.keys())
        if missing:
            raise ValueError(f"The following leaves are missing from leaf2clade: {missing}")

        # --------------------------------------------------------------
        # 2️⃣  Encode cell types as integers 0 … K‑1
        # --------------------------------------------------------------
        if type_encoder is None:
            uniq = sorted(set(self.leaf2type.values()))
            self.type_encoder = {t: i for i, t in enumerate(uniq)}
        else:
            self.type_encoder = type_encoder
        self.leaf2type_int = {
            leaf: self.type_encoder[self.leaf2type[leaf]]
            for leaf in leaf_names
        }

        # --------------------------------------------------------------
        # 3️⃣  Store the two node‑transition matrices (must be stochastic)
        # --------------------------------------------------------------
        self.Phi0 = np.asarray(Phi0, dtype=float)
        self.Phi1 = np.asarray(Phi1, dtype=float)
        self.K = self.Phi0.shape[0]
        assert self.Phi0.shape == (self.K, self.K)
        assert self.Phi1.shape == (self.K, self.K)
        checkstochastic(self.Phi0, "Phi0")
        checkstochastic(self.Phi1, "Phi1")

        # --------------------------------------------------------------
        # 4️⃣  Root prior (founder type is known)
        # --------------------------------------------------------------
        if isinstance(root_type, str):
            root_type = self.type_encoder[root_type]
            # root_type = self.type_encoder.get(root_type, max(self.type_encoder.values())+1)
        self.root_type = int(root_type)

        # --------------------------------------------------------------
        # 5️⃣  Initialise clade‑specific mixture parameters
        # --------------------------------------------------------------
        self.clades = sorted(set(self.leaf2clade.values()))
        self.p_frac = {c: 0.5 for c in self.clades}

        # cache for fast node‑to‑clade look‑up
        self._assign_node_clades()

    # ------------------------------------------------------------------
    #  Helper methods
    # ------------------------------------------------------------------
    def _check_tree_is_rooted(self) -> None:
        """Raise if the supplied tree has no root."""
        if not self.tree.get_tree_root():
            raise ValueError("The supplied tree must be rooted.")

    def _assign_node_clades(self) -> None:
        """
        Attach a ``clade_set`` attribute to every node:
        the set of clades that have at least one descendant leaf under that node.
        """
        for node in self.tree.traverse("postorder"):
            if node.is_leaf():
                node.clade_set = {self.leaf2clade[node.name]}
            else:
                node.clade_set = set()
                for child in node.children:
                    node.clade_set.update(child.clade_set)

    # ------------------------------------------------------------------
    #  Node‑weight factories
    # ------------------------------------------------------------------
    def _node_weights_uniform(self, node) -> Dict[str, float]:
        """Uniform‑per‑node: weight = 1 for every clade covered by the node."""
        return {c: 1.0 for c in node.clade_set}

    def _node_weights_inverse_leaf(self, node) -> Dict[str, float]:
        """Inverse‑leaf: weight = 1 / (# descendant leaves)."""
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
        """LR × uniform (equivalent to pure LR because uniform = 1)."""
        return {c: node.lr_ratio for c in node.clade_set}

    def _node_weights_lr_plus_inv(self, node) -> Dict[str, float]:
        """LR × inverse‑leaf."""
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
    
        for node in self.tree.traverse("postorder"):
            if node.is_leaf():
                # ----- leaf emission (one‑hot) -----
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

            # ----- 1️⃣ lift each child through Φ0 and Φ1 and sum the logs ----------
            logL0 = np.zeros(K)
            logL1 = np.zeros(K)
            for child in node.children:
                # prob0 = self.Phi0.T @ child.lik      # shape (K,)
                # prob1 = self.Phi1.T @ child.lik

                prob0 = self.Phi0 @ child.lik
                prob1 = self.Phi1 @ child.lik
                # avoid log(0); eps is tiny enough not to affect the result
                prob0[prob0 == 0] = eps
                prob1[prob1 == 0] = eps
                logL0 += np.log(prob0)
                logL1 += np.log(prob1)
    
            # ----- 2️⃣ remove a common offset (scaling) -----------------------------
            max0 = logL0.max()
            node.logL0 = logL0 - max0
            node.logL1 = logL1 - max0          # same offset for both components
            node.log_factor = max0
    
            # ----- 3️⃣ LR scalar (for “lr” weighting) -----------------------------
            node.lr_ratio = np.exp(node.logL1.max()) / (np.exp(node.logL0.max()) + eps)
    
            # ----- 4️⃣ clade‑specific mixture probability p_u ----------------------
            clade_weights = defaultdict(float)
            for leaf in node.iter_leaves():
                clade_weights[self.leaf2clade[leaf.name]] += 1.0
            total = sum(clade_weights.values())
            if total == 0:
                p_u = 0.0
            else:
                p_u = sum((w / total) * self.p_frac[c]
                          for c, w in clade_weights.items())
    
            # ----- 5️⃣ mixture upward vector (still a probability vector) ----------
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
    
        # ----- 6️⃣ root contribution ------------------------------------------------
        root = self.tree.get_tree_root()
        root_prior = np.zeros(K)
        root_prior[self.root_type] = 1.0
    
        logL = np.log(np.dot(root_prior, root.lik) + eps) + total_log_factor
    
        # ----- 7️⃣ diagnostics (only once) -----------------------------------------
        if not hasattr(self, "_diag_printed"):
            internal = sum(1 for n in self.tree.traverse() if not n.is_leaf())
            mixed_nodes = sum(1 for n in self.tree.traverse()
                              if not n.is_leaf() and getattr(n, "has_mixture", False))
            self._diag_printed = True
    
        return logL

    
   

    def node_responsibility(self, node, parent_prior_vec=None):
        """
        Public version of the responsibility calculation.
        Returns the posterior probability γ_u that the division represented by ``node``
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
            Posterior responsibility γ_u ∈ [0, 1].
        """
        eps = 1e-12
        
        # ------------------------------------------------------------------
        # 1️⃣  Get the parent prior vector
        # ------------------------------------------------------------------
        if parent_prior_vec is None:
            if node.up is None:
                parent_prior_vec = np.zeros(self.K)
                parent_prior_vec[self.root_type] = 1.0
            else:
                parent_prior_vec = node.up.lik
        
        # ------------------------------------------------------------------
        # 2️⃣  Compute the clade‑weighted prior mixture
        # ------------------------------------------------------------------
        clade_weights = defaultdict(float)
        for leaf in node.iter_leaves():
            clade_weights[self.leaf2clade[leaf.name]] += 1.0
        total = sum(clade_weights.values())
        if total == 0:
            prior_mix = 0.0
        else:
            prior_mix = sum((w / total) * self.p_frac[c]
                            for c, w in clade_weights.items()) 
        
        # ------------------------------------------------------------------
        # 3️⃣  Log‑likelihood pieces that were stored by pruninglikelihood
        # ------------------------------------------------------------------
        log_factor = node.log_factor
        
        # 🔑 CRITICAL FIX: weight by parent_prior_vec BEFORE taking logsumexp
        # We want log(parent_prior · L_vec), not just log(sum(L_vec))
        parent_prior_log = np.log(parent_prior_vec + eps)  # shape (K,)
        
        # log(parent_prior · exp(logL0 + log_factor))
        #  = log(parent_prior) + logL0 + log_factor
        log_term0 = logsumexp(parent_prior_log + node.logL0 + log_factor)
        
        # log(parent_prior · exp(logL1 + log_factor))
        log_term1 = logsumexp(parent_prior_log + node.logL1 + log_factor)
        
        # ------------------------------------------------------------------
        # 4️⃣  Return the posterior γ_u using Bayes rule
        # ------------------------------------------------------------------
        if prior_mix <= 0.0:
            return 0.0
        if prior_mix >= 1.0:
            return 1.0
        
        # γ_u = P(Z_u=1 | data)
        #     = prior_mix · P(data|Z_u=1) / [prior_mix · P(data|Z_u=1) + (1-prior_mix) · P(data|Z_u=0)]
        #     = prior_mix · exp(log_term1) / [prior_mix · exp(log_term1) + (1-prior_mix) · exp(log_term0)]
        
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
        Convert the clade‑level mixture weights (self.p_frac) into leaf‑level
        predictions.
    
        Parameters
        ----------
        threshold : float, default 0.5
            Decision threshold for the hard label.  Leaves with
            ``prob_induced >= threshold`` are labelled “induced”.
        use_path_responsibilities : bool, default False
            If True, incorporate the posterior responsibilities (γ_u) of every
            internal node on the path from the leaf to the root.  The leaf
            probability is then
                P(induced | leaf) = ∏_{u on path} [ p_parent·γ_u + (1-p_parent)·(1-γ_u) ].
            The product is evaluated in log‑space for numerical stability.
            If False (default) the leaf probability is simply the clade mixture
            weight p_c.
    
        Returns
        -------
        pandas.DataFrame with columns:
            - leaf_name      : identifier of the leaf / cell
            - clade          : clade label taken from self.leaf2clade
            - p_clade        : posterior mixture weight for that clade (p_c)
            - prob_induced   : soft probability that the leaf is induced
            - pred_label     : hard label (“induced” or “un‑induced”)
        """
        records = []
    
        for leaf in self.tree.iter_leaves():
            leaf_name = leaf.name
            clade = self.leaf2clade[leaf_name]
            p_c   = float(self.p_frac[clade])      # clade‑level posterior weight
    
            # --------------------------------------------------------------
            # A) Simple clade‑only probability (fast)
            # --------------------------------------------------------------
            prob = p_c
    
            # --------------------------------------------------------------
            # B) Path‑aware probability (if requested)
            # --------------------------------------------------------------
            if use_path_responsibilities:
                # Walk up from leaf to root, accumulating log‑probability.
                node = leaf
                log_prob = np.log(p_c + 1e-12)      # start with the clade weight
                while not node.is_root():
                    parent = node.up
                    # γ_u for the division that produced `node`
                    gamma_u = self.node_responsibility(node)
    
                    # mixture weight for the parent clade (if the parent has a clade)
                    parent_clade = self.leaf2clade.get(parent.name, None)
                    if parent_clade is None:
                        p_parent = 0.5               # neutral prior when unknown
                    else:
                        p_parent = float(self.p_frac[parent_clade])
    
                    # Edge‑level probability that this division was induced:
                    edge_prob = p_parent * gamma_u + (1.0 - p_parent) * (1.0 - gamma_u)
    
                    log_prob += np.log(edge_prob + 1e-12)
                    node = parent
                prob = np.exp(log_prob)          # back to probability space
    
            # --------------------------------------------------------------
            # Hard label
            # --------------------------------------------------------------
            label = "induced" if prob >= threshold else "un‑induced"
    
            records.append({
                "leaf_name": leaf_name,
                "clade": clade,
                "p_clade": p_c,
                "prob_induced": prob,
                "pred_label": label
            })
    
        return pd.DataFrame.from_records(records)
    # ------------------------------------------------------------------
    #  EM algorithm
    # ------------------------------------------------------------------
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
        per‑node contributions are weighted in the M‑step.
        """
        logL_history = []

     

        # ------------------------------------------------------------------
        #  Choose the weight‑lookup function once (so we do not test strings
        #  inside the inner loop).
        # ------------------------------------------------------------------
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
    
        for it in tqdm.trange(max_iter, disable=not verbose):
            # ------------------- E‑step -------------------
            # **IMPORTANT** – call the *correct* pruning routine
            logL = self.pruninglikelihood()          # <-- fixed line
            logL_history.append(logL)
    
            # (optional) sanity check that vectors exist
            # for n in self.tree.traverse():
            #     if not n.is_leaf() and not n.is_root():
            #         assert n.L0_vec is not None and n.L1_vec is not None
    
            if verbose:
                logging.info(f"EM iteration {it+1:02d} – logL = {logL:.6f}")
    
            # ------------------- Gather sufficient statistics -------------------
            clade_stats = {c: [0.0, 0.0] for c in self.clades}   # [numer, denom]
    
            # Root prior vector (used as parent prior for the root's children)
            root_prior = np.zeros(self.K)
            root_prior[self.root_type] = 1.0


            for node in self.tree.traverse("postorder"):
                if node.is_root() or node.is_leaf(): # added is_leaf() check
                    continue                       # root has no Z_u

                clade_weights = defaultdict(float)
                for leaf in node.iter_leaves():
                    clade_weights[self.leaf2clade[leaf.name]] += 1.0
                # if len(clade_weights) > 1:
                #     print(f"  Node {node.name}: split across {len(clade_weights)} clades: {dict(clade_weights)}")
                # elif len(clade_weights) == 1:
                #     print(f"  Node {node.name}: has clade_lengths == 1")
    
                # parent prior comes from the already‑computed node.up.lik
                parent_vec = node.up.lik if node.up else root_prior
    
                # ------------------- posterior responsibility -------------------
                gamma_u = self.node_responsibility(node, parent_vec)   # uses L0_vec/L1_vec
    
                # ------------------- node weight according to strategy -------
                w_dict = weight_func(node)        # e.g. {'cladeA': 5.0, ...}
                # print(f'{w_dict = }')
    
                # accumulate numerator / denominator for each clade
                for c, w in w_dict.items():
                    clade_stats[c][0] += w * gamma_u    # Σ w·γ
                    clade_stats[c][1] += w             # Σ w
    
            # ------------------- M‑step (MAP update) -------------------
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
    
            # ------------------- convergence check -------------------
            logL_change = abs(logL_history[-1] - logL_history[-2]) if it > 0 else np.inf

            if verbose:
                if it > 0 and logL_change < tol:
                    logging.info(f"EM converged (log‑likelihood change {logL_change:.2e} < {tol})")
                    break
            if it > 0 and abs(logL_history[-1] - logL_history[-2]) < tol:
                if verbose:
                    logging.info("EM converged (log‑likelihood change < tol).")
                break


        # final_logL = logL_history[-1]
        # improvement = final_logL - baseline_logL
        # print(f"\nFinal logL:     {final_logL:.8f}")
        # print(f"Improvement:    {improvement:.8f}")
        # if improvement > 0:
        #     print("✓ Model improved over baseline")
        # else:
        #     print("⚠️  Model did NOT improve over baseline")
        self.logL_history_ = logL_history
        return logL_history
    # ------------------------------------------------------------------
    #  Post‑processing utilities
    # ------------------------------------------------------------------
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
        print(df.to_string(index=False))
        if hasattr(self, "logL_history_"):
            pass

    def diagnose_gamma_by_clade(self):
        """
        Inspect per-node γ_u values grouped by clade.
        Should be called AFTER fit() to see which clades prefer which matrix.
        """
        
        # First, recompute the pruning pass to populate L0_vec, L1_vec
        self.pruninglikelihood()
        
        root_prior = np.zeros(self.K)
        # root_prior[self.type_encoder[self.root_type]] = 1.0
        root_prior[self.root_type] = 1.0
        
        # Collect γ values by clade
        gamma_by_clade = defaultdict(list)
        
        for node in self.tree.traverse("postorder"):
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
                clade = max(clade_weights, key=clade_weights.get)
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
                    for node_name, gamma in sorted(gammas, key=lambda x: x[1]):
                        pass
        
        # Overall statistics
        all_gammas = [g for gammas in gamma_by_clade.values() for _, g in gammas]
        
        return gamma_by_clade


class og_TreeMixtureNodeEM2():
    """
    EM for a node‑only mixture of two K‑state stochastic matrices on a rooted
    phylogeny.  The *likelihood‑ratio* weighting scheme is implemented as one
    of the possible ``weight_strategy`` options.
    """

    # ------------------------------------------------------------------
    #  Constructor
    # ------------------------------------------------------------------
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
                'leaf_prop'   – leaf‑proportion weighting (original code)
                'uniform'     – uniform‑per‑node
                'inv_leaf'    – inverse‑leaf weighting
                
                
                'lr+uniform'  – LR multiplied by uniform (default for LR)
                'lr+inv'      – LR multiplied by inverse‑leaf
        """
        # --------------------------------------------------------------
        # 1️⃣  Store tree and attach leaf metadata
        # --------------------------------------------------------------
        self.tree = tree
        # self._check_tree_is_rooted()
        self.leaf2clade = leaf2clade
        self.leaf2type = leaf2type
        self.weight_strategy = weight_strategy

        leaf_names = {leaf.name for leaf in self.tree.iter_leaves()}
        missing = leaf_names - set(self.leaf2clade.keys())
        if missing:
            raise ValueError(f"The following leaves are missing from leaf2clade: {missing}")

        # --------------------------------------------------------------
        # 2️⃣  Encode cell types as integers 0 … K‑1
        # --------------------------------------------------------------
        if type_encoder is None:
            uniq = sorted(set(self.leaf2type.values()))
            self.type_encoder = {t: i for i, t in enumerate(uniq)}
        else:
            self.type_encoder = type_encoder
        self.leaf2type_int = {
            leaf: self.type_encoder[self.leaf2type[leaf]]
            for leaf in leaf_names
        }

        # --------------------------------------------------------------
        # 3️⃣  Store the two node‑transition matrices (must be stochastic)
        # --------------------------------------------------------------
        self.Phi0 = np.asarray(Phi0, dtype=float)
        self.Phi1 = np.asarray(Phi1, dtype=float)
        self.K = self.Phi0.shape[0]
        assert self.Phi0.shape == (self.K, self.K)
        assert self.Phi1.shape == (self.K, self.K)
        checkstochastic(self.Phi0, "Phi0")
        checkstochastic(self.Phi1, "Phi1")

        # --------------------------------------------------------------
        # 4️⃣  Root prior (founder type is known)
        # --------------------------------------------------------------
        if isinstance(root_type, str):
            root_type = self.type_encoder[root_type]
            # root_type = self.type_encoder.get(root_type, max(self.type_encoder.values())+1)
        self.root_type = int(root_type)

        # --------------------------------------------------------------
        # 5️⃣  Initialise clade‑specific mixture parameters
        # --------------------------------------------------------------
        self.clades = sorted(set(self.leaf2clade.values()))
        self.p_frac = {c: 0.5 for c in self.clades}

        # cache for fast node‑to‑clade look‑up
        self._assign_node_clades()

    # ------------------------------------------------------------------
    #  Helper methods
    # ------------------------------------------------------------------
    def _check_tree_is_rooted(self) -> None:
        """Raise if the supplied tree has no root."""
        if not self.tree.get_tree_root():
            raise ValueError("The supplied tree must be rooted.")

    def _assign_node_clades(self) -> None:
        """
        Attach a ``clade_set`` attribute to every node:
        the set of clades that have at least one descendant leaf under that node.
        """
        for node in self.tree.traverse("postorder"):
            if node.is_leaf():
                node.clade_set = {self.leaf2clade[node.name]}
            else:
                node.clade_set = set()
                for child in node.children:
                    node.clade_set.update(child.clade_set)

    # ------------------------------------------------------------------
    #  Node‑weight factories
    # ------------------------------------------------------------------
    def _node_weights_uniform(self, node) -> Dict[str, float]:
        """Uniform‑per‑node: weight = 1 for every clade covered by the node."""
        return {c: 1.0 for c in node.clade_set}

    def _node_weights_inverse_leaf(self, node) -> Dict[str, float]:
        """Inverse‑leaf: weight = 1 / (# descendant leaves)."""
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
        """LR × uniform (equivalent to pure LR because uniform = 1)."""
        return {c: node.lr_ratio for c in node.clade_set}

    def _node_weights_lr_plus_inv(self, node) -> Dict[str, float]:
        """LR × inverse‑leaf."""
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
    
        for node in self.tree.traverse("postorder"):
            if node.is_leaf():
                # ----- leaf emission (one‑hot) -----
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

            # ----- 1️⃣ lift each child through Φ0 and Φ1 and sum the logs ----------
            logL0 = np.zeros(K)
            logL1 = np.zeros(K)
            for child in node.children:
                prob0 = self.Phi0.T @ child.lik      # shape (K,)
                prob1 = self.Phi1.T @ child.lik
                # avoid log(0); eps is tiny enough not to affect the result
                prob0[prob0 == 0] = eps
                prob1[prob1 == 0] = eps
                logL0 += np.log(prob0)
                logL1 += np.log(prob1)
    
            # ----- 2️⃣ remove a common offset (scaling) -----------------------------
            max0 = logL0.max()
            node.logL0 = logL0 - max0
            node.logL1 = logL1 - max0          # same offset for both components
            node.log_factor = max0
    
            # ----- 3️⃣ LR scalar (for “lr” weighting) -----------------------------
            node.lr_ratio = np.exp(node.logL1.max()) / (np.exp(node.logL0.max()) + eps)
    
            # ----- 4️⃣ clade‑specific mixture probability p_u ----------------------
            clade_weights = defaultdict(float)
            for leaf in node.iter_leaves():
                clade_weights[self.leaf2clade[leaf.name]] += 1.0
            total = sum(clade_weights.values())
            if total == 0:
                p_u = 0.0
            else:
                p_u = sum((w / total) * self.p_frac[c]
                          for c, w in clade_weights.items())
    
            # ----- 5️⃣ mixture upward vector (still a probability vector) ----------
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
    
        # ----- 6️⃣ root contribution ------------------------------------------------
        root = self.tree.get_tree_root()
        root_prior = np.zeros(K)
        root_prior[self.root_type] = 1.0
    
        logL = np.log(np.dot(root_prior, root.lik) + eps) + total_log_factor
    
        # ----- 7️⃣ diagnostics (only once) -----------------------------------------
        if not hasattr(self, "_diag_printed"):
            internal = sum(1 for n in self.tree.traverse() if not n.is_leaf())
            mixed_nodes = sum(1 for n in self.tree.traverse()
                              if not n.is_leaf() and getattr(n, "has_mixture", False))
            self._diag_printed = True
    
        return logL

    
   

    def node_responsibility(self, node, parent_prior_vec=None):
        """
        Public version of the responsibility calculation.
        Returns the posterior probability γ_u that the division represented by ``node``
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
            Posterior responsibility γ_u ∈ [0, 1].
        """
        eps = 1e-12
        
        # ------------------------------------------------------------------
        # 1️⃣  Get the parent prior vector
        # ------------------------------------------------------------------
        if parent_prior_vec is None:
            if node.up is None:
                parent_prior_vec = np.zeros(self.K)
                parent_prior_vec[self.root_type] = 1.0
            else:
                parent_prior_vec = node.up.lik
        
        # ------------------------------------------------------------------
        # 2️⃣  Compute the clade‑weighted prior mixture
        # ------------------------------------------------------------------
        clade_weights = defaultdict(float)
        for leaf in node.iter_leaves():
            clade_weights[self.leaf2clade[leaf.name]] += 1.0
        total = sum(clade_weights.values())
        if total == 0:
            prior_mix = 0.0
        else:
            prior_mix = sum((w / total) * self.p_frac[c]
                            for c, w in clade_weights.items()) 
        
        # ------------------------------------------------------------------
        # 3️⃣  Log‑likelihood pieces that were stored by pruninglikelihood
        # ------------------------------------------------------------------
        log_factor = node.log_factor
        
        # 🔑 CRITICAL FIX: weight by parent_prior_vec BEFORE taking logsumexp
        # We want log(parent_prior · L_vec), not just log(sum(L_vec))
        parent_prior_log = np.log(parent_prior_vec + eps)  # shape (K,)
        
        # log(parent_prior · exp(logL0 + log_factor))
        #  = log(parent_prior) + logL0 + log_factor
        log_term0 = logsumexp(parent_prior_log + node.logL0 + log_factor)
        
        # log(parent_prior · exp(logL1 + log_factor))
        log_term1 = logsumexp(parent_prior_log + node.logL1 + log_factor)
        
        # ------------------------------------------------------------------
        # 4️⃣  Return the posterior γ_u using Bayes rule
        # ------------------------------------------------------------------
        if prior_mix <= 0.0:
            return 0.0
        if prior_mix >= 1.0:
            return 1.0
        
        # γ_u = P(Z_u=1 | data)
        #     = prior_mix · P(data|Z_u=1) / [prior_mix · P(data|Z_u=1) + (1-prior_mix) · P(data|Z_u=0)]
        #     = prior_mix · exp(log_term1) / [prior_mix · exp(log_term1) + (1-prior_mix) · exp(log_term0)]
        
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
        Convert the clade‑level mixture weights (self.p_frac) into leaf‑level
        predictions.
    
        Parameters
        ----------
        threshold : float, default 0.5
            Decision threshold for the hard label.  Leaves with
            ``prob_induced >= threshold`` are labelled “induced”.
        use_path_responsibilities : bool, default False
            If True, incorporate the posterior responsibilities (γ_u) of every
            internal node on the path from the leaf to the root.  The leaf
            probability is then
                P(induced | leaf) = ∏_{u on path} [ p_parent·γ_u + (1-p_parent)·(1-γ_u) ].
            The product is evaluated in log‑space for numerical stability.
            If False (default) the leaf probability is simply the clade mixture
            weight p_c.
    
        Returns
        -------
        pandas.DataFrame with columns:
            - leaf_name      : identifier of the leaf / cell
            - clade          : clade label taken from self.leaf2clade
            - p_clade        : posterior mixture weight for that clade (p_c)
            - prob_induced   : soft probability that the leaf is induced
            - pred_label     : hard label (“induced” or “un‑induced”)
        """
        records = []
    
        for leaf in self.tree.iter_leaves():
            leaf_name = leaf.name
            clade = self.leaf2clade[leaf_name]
            p_c   = float(self.p_frac[clade])      # clade‑level posterior weight
    
            # --------------------------------------------------------------
            # A) Simple clade‑only probability (fast)
            # --------------------------------------------------------------
            prob = p_c
    
            # --------------------------------------------------------------
            # B) Path‑aware probability (if requested)
            # --------------------------------------------------------------
            if use_path_responsibilities:
                # Walk up from leaf to root, accumulating log‑probability.
                node = leaf
                log_prob = np.log(p_c + 1e-12)      # start with the clade weight
                while not node.is_root():
                    parent = node.up
                    # γ_u for the division that produced `node`
                    gamma_u = self.node_responsibility(node)
    
                    # mixture weight for the parent clade (if the parent has a clade)
                    parent_clade = self.leaf2clade.get(parent.name, None)
                    if parent_clade is None:
                        p_parent = 0.5               # neutral prior when unknown
                    else:
                        p_parent = float(self.p_frac[parent_clade])
    
                    # Edge‑level probability that this division was induced:
                    edge_prob = p_parent * gamma_u + (1.0 - p_parent) * (1.0 - gamma_u)
    
                    log_prob += np.log(edge_prob + 1e-12)
                    node = parent
                prob = np.exp(log_prob)          # back to probability space
    
            # --------------------------------------------------------------
            # Hard label
            # --------------------------------------------------------------
            label = "induced" if prob >= threshold else "un‑induced"
    
            records.append({
                "leaf_name": leaf_name,
                "clade": clade,
                "p_clade": p_c,
                "prob_induced": prob,
                "pred_label": label
            })
    
        return pd.DataFrame.from_records(records)
    # ------------------------------------------------------------------
    #  EM algorithm
    # ------------------------------------------------------------------
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
        per‑node contributions are weighted in the M‑step.
        """
        logL_history = []

     

        # ------------------------------------------------------------------
        #  Choose the weight‑lookup function once (so we do not test strings
        #  inside the inner loop).
        # ------------------------------------------------------------------
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
    
        for it in tqdm.trange(max_iter, disable=not verbose):
            # ------------------- E‑step -------------------
            # **IMPORTANT** – call the *correct* pruning routine
            logL = self.pruninglikelihood()          # <-- fixed line
            logL_history.append(logL)
    
            # (optional) sanity check that vectors exist
            # for n in self.tree.traverse():
            #     if not n.is_leaf() and not n.is_root():
            #         assert n.L0_vec is not None and n.L1_vec is not None
    
            if verbose:
                logging.info(f"EM iteration {it+1:02d} – logL = {logL:.6f}")
    
            # ------------------- Gather sufficient statistics -------------------
            clade_stats = {c: [0.0, 0.0] for c in self.clades}   # [numer, denom]
    
            # Root prior vector (used as parent prior for the root's children)
            root_prior = np.zeros(self.K)
            root_prior[self.root_type] = 1.0


            for node in self.tree.traverse("postorder"):
                if node.is_root() or node.is_leaf(): # added is_leaf() check
                    continue                       # root has no Z_u

                clade_weights = defaultdict(float)
                for leaf in node.iter_leaves():
                    clade_weights[self.leaf2clade[leaf.name]] += 1.0
                # if len(clade_weights) > 1:
                #     print(f"  Node {node.name}: split across {len(clade_weights)} clades: {dict(clade_weights)}")
                # elif len(clade_weights) == 1:
                #     print(f"  Node {node.name}: has clade_lengths == 1")
    
                # parent prior comes from the already‑computed node.up.lik
                parent_vec = node.up.lik if node.up else root_prior
    
                # ------------------- posterior responsibility -------------------
                gamma_u = self.node_responsibility(node, parent_vec)   # uses L0_vec/L1_vec
    
                # ------------------- node weight according to strategy -------
                w_dict = weight_func(node)        # e.g. {'cladeA': 5.0, ...}
                # print(f'{w_dict = }')
    
                # accumulate numerator / denominator for each clade
                for c, w in w_dict.items():
                    clade_stats[c][0] += w * gamma_u    # Σ w·γ
                    clade_stats[c][1] += w             # Σ w
    
            # ------------------- M‑step (MAP update) -------------------
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
    
            # ------------------- convergence check -------------------
            logL_change = abs(logL_history[-1] - logL_history[-2]) if it > 0 else np.inf

            if verbose:
                if it > 0 and logL_change < tol:
                    logging.info(f"EM converged (log‑likelihood change {logL_change:.2e} < {tol})")
                    break
            if it > 0 and abs(logL_history[-1] - logL_history[-2]) < tol:
                if verbose:
                    logging.info("EM converged (log‑likelihood change < tol).")
                break


        # final_logL = logL_history[-1]
        # improvement = final_logL - baseline_logL
        # print(f"\nFinal logL:     {final_logL:.8f}")
        # print(f"Improvement:    {improvement:.8f}")
        # if improvement > 0:
        #     print("✓ Model improved over baseline")
        # else:
        #     print("⚠️  Model did NOT improve over baseline")
        self.logL_history_ = logL_history
        return logL_history
    # ------------------------------------------------------------------
    #  Post‑processing utilities
    # ------------------------------------------------------------------
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
        print(df.to_string(index=False))
        if hasattr(self, "logL_history_"):
            pass

    def diagnose_gamma_by_clade(self):
        """
        Inspect per-node γ_u values grouped by clade.
        Should be called AFTER fit() to see which clades prefer which matrix.
        """
        
        # First, recompute the pruning pass to populate L0_vec, L1_vec
        self.pruninglikelihood()
        
        root_prior = np.zeros(self.K)
        # root_prior[self.type_encoder[self.root_type]] = 1.0
        root_prior[self.root_type] = 1.0
        
        # Collect γ values by clade
        gamma_by_clade = defaultdict(list)
        
        for node in self.tree.traverse("postorder"):
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
                clade = max(clade_weights, key=clade_weights.get)
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
                    for node_name, gamma in sorted(gammas, key=lambda x: x[1]):
                        pass
        
        # Overall statistics
        all_gammas = [g for gammas in gamma_by_clade.values() for _, g in gammas]
        
        return gamma_by_clade


from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from sklearn.metrics import roc_auc_score, roc_curve

def FULL_make_all_em_plots(param_combo_list,
                           full_modelres_dict,
                           urid_timepoint_res_dict,
                           return_plotting_dfs = False):
    
    full_plotting_df_dict = {}
    
    for urid in list(full_modelres_dict.keys()):
        
        full_plotting_df_dict[urid] = {}
    
        
        path_to_output_dir = Path('./output')
        path_to_em_fig_dir = path_to_output_dir / 'em_figs' / str(urid)

        os.makedirs(path_to_em_fig_dir, exist_ok=True)

        plot_matrix_sep(urid, path_to_output_dir = path_to_output_dir,
                        save_path = path_to_em_fig_dir / 'matrix_sep')


        for param_combo in param_combo_list:

            path_to_param_combo_dir = path_to_em_fig_dir / param_combo

            os.makedirs(path_to_param_combo_dir, exist_ok=True)

            timepoints = sorted(list(full_modelres_dict[urid].keys()))
            
            tpt_df_dict = {}
            eval_dict_by_tp = {}
            
            for tp in timepoints:
                tpt_df_dict[tp] = full_modelres_dict[urid][tp]['tpt_induction_dict']
                ed = full_modelres_dict[urid][tp]['eval_dict'].get(param_combo)
                if ed is None:
                    print(f'WARNING: no eval_dict for urid={urid} tp={tp} combo={param_combo}', flush=True)
                    ed = {}
                eval_dict_by_tp[tp] = ed


            annotated_dfs = build_annotated_dfs(
                tpt_df_dict, urid=urid, timepoints=timepoints,
                urid_timepoint_res_dict= urid_timepoint_res_dict,
                param_combo=param_combo,
                path_to_output_dir = path_to_output_dir)

            fig1 = plot_p_induced_per_clade(annotated_dfs, timepoints,
                                            save_path=path_to_param_combo_dir / f'em_p_induced_per_clade_{urid}_{param_combo}.pdf')
            plt.close(fig1)

            fig2 = plot_roc_curves(eval_dict_by_tp, timepoints, save_path=path_to_param_combo_dir / f'em_roc_{urid}_{param_combo}.pdf')
            plt.close(fig2)

            fig3 = plot_pr_curves(eval_dict_by_tp, timepoints, save_path=path_to_param_combo_dir / f'em_pr_{urid}_{param_combo}.pdf')
            plt.close(fig3)

            fig4 = plot_p_induced_distribution(annotated_dfs, timepoints,
                                                save_path=path_to_param_combo_dir / f'em_p_induced_distribution_{urid}_{param_combo}.pdf')
            plt.close(fig4)

            cell_level_dfs = build_cell_level_dfs(
                tpt_df_dict, urid=urid, timepoints=timepoints,
                urid_timepoint_res_dict=urid_timepoint_res_dict,
                param_combo = param_combo,
                path_to_output_dir = path_to_output_dir)

            fig5 = plot_cell_concordance(cell_level_dfs, timepoints,
                                        save_path=path_to_param_combo_dir / f'em_cell_concordance_{urid}_{param_combo}.pdf')
            plt.close(fig5)
            
            full_plotting_df_dict[urid][param_combo] = annotated_dfs
            
    if return_plotting_dfs:
        return full_plotting_df_dict

        



def plot_roc_comparison(
    full_modelres_dict,
    urid_fracdel_map,
    timepoints,
    output_dir_pattern='./output/{urid}',
    figsize_per_cell=(2.8, 2.8),
    gt_combos=('20clade_topological_treedist_gt', '2clade_topological_treedist_gt'),
    save_path=None,
):
    """
    Grid of ROC curves: rows = frac_del (sorted), cols = timepoint.
    Each panel shows one curve per param_combo, AUC in label.

    Colour coding:
      GT combos  — solid thick black/grey
      min5 MT    — dashed, blue shades by RP (darker = higher RP)
      min10 MT   — dotted, orange shades by RP
    """
    rp_values  = [0.7, 0.8, 0.9, 1.0]
    rp_blues   = ['#c6dbef', '#6baed6', '#2171b5', '#08306b']
    rp_oranges = ['#fdd0a2', '#fd8d3c', '#d94801', '#7f2704']
    rp_color_min5  = dict(zip(rp_values, rp_blues))
    rp_color_min10 = dict(zip(rp_values, rp_oranges))
    gt_colors  = {
        '20clade_topological_treedist_gt': '#1a1a1a',
        '2clade_topological_treedist_gt':  '#666666',
    }

    def _style(combo):
        if combo in gt_colors:
            return dict(color=gt_colors[combo], lw=2.0, ls='-', alpha=1.0, zorder=5)
        rp = float(combo.split('RP')[0].split('_')[-1]) if 'RP' in combo else 1.0
        ls = '--' if 'min5' in combo else ':'
        c  = rp_color_min5.get(rp, '#000') if 'min5' in combo else rp_color_min10.get(rp, '#000')
        return dict(color=c, lw=1.2, ls=ls, alpha=0.85, zorder=2)

    def _label(combo):
        if 'gt' in combo:
            return combo.replace('clade_topological_treedist_gt', 'clades-GT')
        minp = 'min5' if 'min5' in combo else 'min10'
        rp   = next((p for p in combo.split('_') if 'RP' in p), '')
        af   = next((p for p in combo.split('_') if 'AF' in p), '')
        return f'{minp} {rp} {af}'

    def _gt_map(urid, tp, combo):
        csv = Path(output_dir_pattern.format(urid=urid)) / \
              'clade_assignment_dfs' / str(urid) / str(tp) / \
              f'clade_df_{urid}_t{tp}.csv'
        if not csv.exists():
            return {}
        cdf = pd.read_csv(csv)
        if combo not in cdf.columns:
            return {}
        return (cdf.groupby(combo)['clade_differentiation_induced']
                   .apply(lambda x: 'induced' if x.astype(float).mean() > 0.5 else 'uninduced')
                   .to_dict())

    # collect all combos present across the dict
    all_combos = []
    for urid in full_modelres_dict:
        for tp in full_modelres_dict[urid]:
            for c in full_modelres_dict[urid][tp]['tpt_induction_dict']['clade_name'].unique():
                if c not in all_combos:
                    all_combos.append(c)

    frac_dels = sorted(urid_fracdel_map.values())
    fracdel_to_urid = {v: k for k, v in urid_fracdel_map.items()}
    nrows, ncols = len(frac_dels), len(timepoints)

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(figsize_per_cell[0] * ncols, figsize_per_cell[1] * nrows),
        squeeze=False,
    )

    for ri, fd in enumerate(frac_dels):
        urid = fracdel_to_urid[fd]
        for ci, tp in enumerate(timepoints):
            ax = axes[ri][ci]
            ax.plot([0, 1], [0, 1], color='lightgrey', lw=0.8, ls='--', zorder=0)

            tpt = full_modelres_dict[urid][tp]['tpt_induction_dict']

            for combo in all_combos:
                grp = tpt[tpt['clade_name'] == combo].copy()
                if len(grp) == 0:
                    continue
                gt = _gt_map(urid, tp, combo)
                grp['gt_label'] = grp['clade'].map(gt)
                valid = grp.dropna(subset=['gt_label', 'p_induced'])
                if len(valid) < 2 or valid['gt_label'].nunique() < 2:
                    continue
                y_true  = (valid['gt_label'] == 'induced').astype(int).values
                y_score = valid['p_induced'].astype(float).values
                try:
                    fpr, tpr, _ = roc_curve(y_true, y_score)
                    auc = roc_auc_score(y_true, y_score)
                except Exception:
                    continue
                sty = _style(combo)
                ax.plot(fpr, tpr, label=f'{_label(combo)} ({auc:.2f})', **sty)

            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.set_xticks([0, 0.5, 1])
            ax.set_yticks([0, 0.5, 1])
            ax.tick_params(labelsize=7)
            if ri == 0:
                ax.set_title(f't={tp}', fontsize=9)
            if ci == 0:
                ax.set_ylabel(f'frac_del={fd}\nTPR', fontsize=8)
            if ri == nrows - 1:
                ax.set_xlabel('FPR', fontsize=8)

    legend_handles = [
        mlines.Line2D([], [], label=_label(c), **_style(c))
        for c in all_combos
    ]
    fig.legend(
        handles=legend_handles, loc='lower center',
        ncol=4, fontsize=7, framealpha=0.9,
        bbox_to_anchor=(0.5, -0.02),
    )
    fig.suptitle('ROC curves: frac_del × timepoint', fontsize=11, y=1.01)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
    return fig




# ============================================================
# Proportion-only classifier: does cell-type composition alone
# distinguish induced vs uninduced clades?
# ============================================================
from numpy.linalg import matrix_power
from scipy.special import rel_entr
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import accuracy_score

Phi0 = np.array([[0.9, 0.05, 0.03, 0.02],
                 [0.04, 0.9, 0.04, 0.02],
                 [0.04, 0.04, 0.9, 0.02],
                 [0.02, 0.02, 0.02, 0.94]])

Phi1 = np.array([[0.9, 0.08, 0.0, 0.02],
                 [0.0, 0.9, 0.08, 0.02],
                 [0.08, 0.0, 0.9, 0.02],
                 [0.02, 0.02, 0.02, 0.94]])


# ------------------------------------------------------------------
# 1. Show how separable the two marginal distributions are by step count.
#    If they converge quickly, proportions will never be enough.
# ------------------------------------------------------------------
def matrix_separation(Phi0, Phi1, prior=None, max_steps=60):
    """
    Compare two transition matrices across step counts.

    Metrics per step n:
      marginal_tv  — TV between marginals prior @ Phi0^n vs prior @ Phi1^n.
                     Collapses source type; the original metric.
      expected_tv  — E_src[ TV(Phi0^n[src,:], Phi1^n[src,:]) ] under prior.
                     Keeps co-occurrence with source: how distinguishable are
                     the two matrices starting from a typical cell?
      max_tv       — max over source types of TV(Phi0^n[src,:], Phi1^n[src,:]).
                     Which starting type maximises separation at step n?
      path2_tv     — TV between 2-step path distributions
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
        tv_per_src = 0.5 * np.abs(Pn0 - Pn1).sum(axis=1)   # shape (K,)
        expected_tv = float(prior @ tv_per_src)
        max_tv      = float(tv_per_src.max())

        # 2-step path joint: P(X_0=i, X_1=j, X_n=k) = prior[i]*Phi[i,j]*Phi^(n-1)[j,k]
        if n > 1:
            Pm0 = matrix_power(Phi0, n - 1)
            Pm1 = matrix_power(Phi1, n - 1)
            # broadcasting: (K,1,1)*(K,K,1)*(1,K,K) → (K,K,K)
            path0 = prior[:, None, None] * Phi0[:, :, None] * Pm0[None, :, :]
            path1 = prior[:, None, None] * Phi1[:, :, None] * Pm1[None, :, :]
        else:
            path0 = prior[:, None] * Phi0
            path1 = prior[:, None] * Phi1
        path2_tv = float(0.5 * np.abs(path0 - path1).sum())

        records.append((n, marginal_tv, expected_tv, max_tv, path2_tv))

    return pd.DataFrame(records,
                        columns=["n_steps", "marginal_tv", "expected_tv", "max_tv", "path2_tv"])


def matrix_diff_heatmap(Phi0, Phi1, n_steps_list, cell_type_names=None, prior=None):
    """
    For each n in n_steps_list plot a K×K heatmap of prior[i]*(Phi0^n - Phi1^n)[i,j].
    Red = Phi0 more likely, Blue = Phi1 more likely for that (source, target) pair.
    """
    K = Phi0.shape[0]
    if prior is None:
        prior = np.ones(K) / K
    if cell_type_names is None:
        cell_type_names = [f'ct{i+1}' for i in range(K)]

    diffs = []
    for n in n_steps_list:
        diff = prior[:, None] * (matrix_power(Phi0, n) - matrix_power(Phi1, n))
        diffs.append(diff)
    vmax = max(np.abs(d).max() for d in diffs)

    fig, axes = plt.subplots(1, len(n_steps_list),
                             figsize=(4.5 * len(n_steps_list), 4.2))
    if len(n_steps_list) == 1:
        axes = [axes]
    for ax, n, diff in zip(axes, n_steps_list, diffs):
        im = ax.imshow(diff, cmap='RdBu_r', vmin=-vmax, vmax=vmax, aspect='auto')
        ax.set_xticks(range(K)); ax.set_xticklabels(cell_type_names, rotation=45, ha='right')
        ax.set_yticks(range(K)); ax.set_yticklabels(cell_type_names)
        ax.set_xlabel('Target cell type (step n)')
        ax.set_ylabel('Source cell type (step 0)')
        ax.set_title(f'n = {n}')
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04,
                     label='prior · (Phi0ⁿ − Phi1ⁿ)')
    fig.suptitle('Source→target separation: Phi0ⁿ − Phi1ⁿ  (prior-weighted)',
                 fontweight='bold')
    fig.tight_layout()
    return fig


# sep_df = matrix_separation(Phi0, Phi1)

def plot_matrix_sep(urid, path_to_output_dir=path_to_output_dir, cell_type_names=None,
                    save_path = None):
    phi0, phi1 = map(np.array, get_tms(urid, path_to_output_dir=path_to_output_dir))
    K = phi0.shape[0]
    if cell_type_names is None:
        cell_type_names = [f'ct{i+1}' for i in range(K)]

    sep_df = matrix_separation(phi0, phi1)

    # ── Row 0: TV metrics over step count ─────────────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))

    cols   = ["marginal_tv", "expected_tv", "max_tv", "path2_tv"]
    colors = ["#95a5a6",     "#2980b9",     "#c0392b", "#27ae60"]
    labels = ["marginal TV (old)", "expected TV", "max TV", "2-step path TV"]

    ax = axes[0]
    for col, color, label in zip(cols, colors, labels):
        ax.plot(sep_df["n_steps"], sep_df[col], color=color, linewidth=1.6, label=label)
    peak_n = int(sep_df.loc[sep_df["expected_tv"].idxmax(), "n_steps"])
    ax.axvline(peak_n, color="black", linestyle=":", linewidth=0.9, label=f"peak n={peak_n}")
    ax.set_xlabel("Steps (n)"); ax.set_ylabel("TV distance")
    ax.set_title("Separation metrics vs step count")
    ax.legend(fontsize=7, frameon=False)

    # ── Row 1 & 2: source-level TV at peak n ──────────────────────────────────
    from numpy.linalg import matrix_power as _mp
    Pn0 = _mp(phi0, peak_n)
    Pn1 = _mp(phi1, peak_n)
    tv_per_src = 0.5 * np.abs(Pn0 - Pn1).sum(axis=1)

    ax = axes[1]
    bar_colors = plt.cm.RdYlGn(tv_per_src / tv_per_src.max())
    ax.bar(cell_type_names, tv_per_src, color=bar_colors, edgecolor="white")
    ax.set_xlabel("Source cell type"); ax.set_ylabel(f"TV at n={peak_n}")
    ax.set_title(f"Per-source separation at peak n={peak_n}")

    # ── Row 3: diff heatmap at peak n ─────────────────────────────────────────
    ax = axes[2]
    diff = phi0 - phi1   # 1-step difference; most interpretable
    vmax = np.abs(diff).max()
    im = ax.imshow(diff, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(K)); ax.set_xticklabels(cell_type_names, rotation=45, ha="right")
    ax.set_yticks(range(K)); ax.set_yticklabels(cell_type_names)
    ax.set_xlabel("Target type"); ax.set_ylabel("Source type")
    ax.set_title("Phi0 − Phi1  (1-step, red = Phi0 higher)")
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    fig.suptitle(f"Matrix separation — URID {urid}", fontweight="bold")
    plt.tight_layout()
    
    if save_path:
        fig.savefig(save_path, bbox_inches='tight')
    plt.close(fig)

    return sep_df, peak_n

# print(sep_df[sep_df["n_steps"].isin([1,2,3,5,10,20,50])].to_string(index=False))


# ------------------------------------------------------------------
# 2. Derive ground-truth clade labels from full_df using
#    clade_differentiation_induced (majority vote per 20-clade).
# ------------------------------------------------------------------
def gt_clade_labels(full_df, param_combo):
    """Return {clade_name: 'induced'|'uninduced'} via majority vote."""
    target = full_df[full_df["param_combo"] == param_combo][["cell", "clade"]]
    gt_raw = full_df[full_df["param_combo"] == "clade_differentiation_induced"][["cell", "clade"]].rename(
        columns={"clade": "induced"}
    )
    merged = target.merge(gt_raw, on="cell", how="inner")
    return (
        merged.groupby("clade")["induced"]
              .apply(lambda x: "induced" if np.mean(x.astype(float)) > 0.5 else "uninduced")
              .to_dict()
    )


# ------------------------------------------------------------------
# 3. Proportion-only classifiers.
# ------------------------------------------------------------------
def proportion_classifiers(full_df, param_combo, Phi0, Phi1,
                            type_encoder, root_type="ct1",
                            n_steps_list=(1, 2, 3, 5, 10, 20)):
    eps = 1e-10
    all_types = sorted(type_encoder, key=type_encoder.get)
    root_idx  = type_encoder[root_type]
    K = Phi0.shape[0]
    root_vec  = np.zeros(K); root_vec[root_idx] = 1.0

    filt = full_df[full_df["param_combo"] == param_combo].copy()
    counts = (filt.groupby(["clade", "celltype"]).size()
                  .unstack(fill_value=0)
                  .reindex(columns=all_types, fill_value=0))
    props = counts.div(counts.sum(axis=1), axis=0)

    gt = gt_clade_labels(full_df, param_combo)
    clades = [c for c in props.index if c in gt]
    props  = props.loc[clades]
    y      = np.array([1 if gt[c] == "induced" else 0 for c in clades])

    results = {}

    # Method A: log-likelihood ratio at various step counts
    for n in n_steps_list:
        exp0 = np.clip(root_vec @ matrix_power(Phi0, n), eps, 1)
        exp1 = np.clip(root_vec @ matrix_power(Phi1, n), eps, 1)
        preds = [1 if (props.loc[c].values @ np.log(exp1)) > (props.loc[c].values @ np.log(exp0)) else 0
                 for c in clades]
        results[f"LLR n={n}"] = accuracy_score(y, preds)

    # Method B: KL divergence from stationary distribution
    stat0 = np.clip(root_vec @ matrix_power(Phi0, 200), eps, 1)
    stat1 = np.clip(root_vec @ matrix_power(Phi1, 200), eps, 1)
    preds_kl = []
    for c in clades:
        obs = np.clip(props.loc[c].values, eps, 1); obs /= obs.sum()
        preds_kl.append(1 if rel_entr(obs, stat1).sum() < rel_entr(obs, stat0).sum() else 0)
    results["KL stationary"] = accuracy_score(y, preds_kl)

    # Method C: leave-one-out logistic regression (upper bound with learned weights)
    X = props.values
    loo_preds = []
    for tr, te in LeaveOneOut().split(X):
        if len(np.unique(y[tr])) < 2:
            loo_preds.append(int(y[tr].mean() > 0.5)); continue
        clf = LogisticRegression(max_iter=1000, C=1.0).fit(X[tr], y[tr])
        loo_preds.append(int(clf.predict(X[te])[0]))
    results["Logistic (LOO)"] = accuracy_score(y, loo_preds)

    # Baseline
    results["Majority class"] = max(y.mean(), 1 - y.mean())

    return (pd.DataFrame.from_dict(results, orient="index", columns=["accuracy"])
              .sort_values("accuracy", ascending=False)), props, y, clades




# ============================================================
# EM plotting functions
# ============================================================
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from sklearn.metrics import auc, confusion_matrix

_INDUCED_COLOR   = '#c0392b'
_UNINDUCED_COLOR = '#2980b9'

# ── data preparation ─────────────────────────────────────────────────────────

def build_annotated_dfs(tpt_df_dict, urid, timepoints, urid_timepoint_res_dict,
                        param_combo='20clade_topological_treedist_gt',
                        path_to_output_dir = path_to_output_dir):
    """
    Per-timepoint clade df augmented with majority-vote GT label.
    Returns {tp: DataFrame(clade, p_induced, gt_label)}.
    """
    annotated = {}
    for tp in timepoints:
        full_df = make_clade_celltype_df(urid=urid, timepoint=tp,
                                         urid_timepoint_res_dict=urid_timepoint_res_dict,
                                         path_to_output_dir = path_to_output_dir)
        target = full_df[full_df['param_combo'] == param_combo][['cell', 'clade']]
        gt_raw = (full_df[full_df['param_combo'] == 'clade_differentiation_induced']
                  [['cell', 'clade']].rename(columns={'clade': 'induced'}))
        merged = target.merge(gt_raw, on='cell', how='inner')
        gt_map = (merged.groupby('clade')['induced']
                  .apply(lambda x: 'induced' if np.mean(x.astype(float)) > 0.5 else 'uninduced')
                  .to_dict())
        df = tpt_df_dict[tp][tpt_df_dict[tp]['clade_name'] == param_combo].copy()
        df['gt_label'] = df['clade'].map(gt_map)
        annotated[tp] = df
    return annotated


def build_cell_level_dfs(tpt_df_dict, urid, timepoints, urid_timepoint_res_dict,
                         param_combo='20clade_topological_treedist_gt',
                         path_to_output_dir = path_to_output_dir):
    """
    Expand clade-level EM predictions to individual cells.
    Returns {tp: DataFrame(clade, p_induced, y_true)} where y_true in {0, 1}.
    """
    cell_dfs = {}
    for tp in timepoints:
        full_df = make_clade_celltype_df(urid=urid, timepoint=tp,
                                          urid_timepoint_res_dict=urid_timepoint_res_dict,
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
        }).dropna(subset=['p_induced']).reset_index(drop=True)
        cell_dfs[tp] = cdf
    return cell_dfs


# ── plots ─────────────────────────────────────────────────────────────────────

def plot_p_induced_per_clade(annotated_dfs, timepoints,
                              induced_color=_INDUCED_COLOR, uninduced_color=_UNINDUCED_COLOR,
                              title='EM p(induced) per clade', save_path=None):
    """Sorted scatter of clade p_induced coloured by GT label."""
    fig, axes = plt.subplots(1, len(timepoints), figsize=(10, 3.2), sharey=True)
    if len(timepoints) == 1:
        axes = [axes]
    for ax, tp in zip(axes, timepoints):
        df = annotated_dfs[tp].sort_values('p_induced')
        colors = [induced_color if g == 'induced' else uninduced_color for g in df['gt_label']]
        ax.scatter(range(len(df)), df['p_induced'], c=colors, s=40, zorder=3, edgecolors='none')
        ax.axhline(0.5, color='#7f8c8d', linewidth=0.8, linestyle='--', zorder=2)
        ax.set_title(f'Timepoint {tp}')
        ax.set_xlabel('Clade (sorted by p)')
        ax.set_ylim(0, 1); ax.set_xticks([])
    axes[0].set_ylabel('p(induced)')
    legend_elements = [
        Line2D([0],[0], marker='o', color='w', markerfacecolor=induced_color,  markersize=7, label='Induced (GT)'),
        Line2D([0],[0], marker='o', color='w', markerfacecolor=uninduced_color, markersize=7, label='Uninduced (GT)'),
    ]
    axes[-1].legend(handles=legend_elements, loc='upper left', frameon=False)
    if title:
        fig.suptitle(title, y=1.02)
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, bbox_inches='tight')
    return fig


def plot_roc_curves(eval_dict_by_tp, timepoints, title='ROC — EM classifier', save_path=None):
    """ROC curve per timepoint on a single axes."""
    valid_tps = [tp for tp in timepoints if eval_dict_by_tp.get(tp, {}).get('fpr') is not None]
    tp_colors = plt.cm.viridis(np.linspace(0.1, 0.9, max(len(valid_tps), 1)))
    fig, ax = plt.subplots(figsize=(4, 4))
    for tp, col in zip(valid_tps, tp_colors):
        ed = eval_dict_by_tp[tp]
        ax.plot(ed['fpr'], ed['tpr'], color=col, linewidth=1.5,
                label=f'tp {tp}  (AUC={auc(ed["fpr"], ed["tpr"]):.2f})')
    ax.plot([0,1],[0,1], color='#bdc3c7', linewidth=0.8, linestyle='--')
    ax.set_xlabel('False positive rate'); ax.set_ylabel('True positive rate')
    ax.set_title(title); ax.legend(frameon=False)
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, bbox_inches='tight')
    return fig


def plot_pr_curves(eval_dict_by_tp, timepoints, baseline_tp=None,
                   title='Precision-Recall — EM classifier', save_path=None):
    """Precision-Recall curve per timepoint with no-skill baseline."""
    valid_tps = [tp for tp in timepoints if eval_dict_by_tp.get(tp, {}).get('recall') is not None]
    tp_colors = plt.cm.viridis(np.linspace(0.1, 0.9, max(len(valid_tps), 1)))
    fig, ax = plt.subplots(figsize=(4, 4))
    for tp, col in zip(valid_tps, tp_colors):
        ed = eval_dict_by_tp[tp]
        ax.plot(ed['recall'], ed['precision'], color=col, linewidth=1.5,
                label=f'tp {tp}  (AUC={auc(ed["recall"], ed["precision"]):.2f})')
    ref_tp  = baseline_tp if baseline_tp is not None else (valid_tps[0] if valid_tps else None)
    gt      = eval_dict_by_tp[ref_tp]['ground_truth'] if ref_tp is not None else []
    baseline = sum(1 for g in gt if g == 1) / len(gt) if len(gt) > 0 else 0.5
    ax.axhline(baseline, color='#bdc3c7', linewidth=0.8, linestyle='--',
               label=f'Baseline ({baseline:.2f})')
    ax.set_xlabel('Recall'); ax.set_ylabel('Precision')
    ax.set_title(title); ax.set_xlim(0,1); ax.set_ylim(0,1.05)
    ax.legend(frameon=False)
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, bbox_inches='tight')
    return fig


def plot_p_induced_distribution(annotated_dfs, timepoints,
                                 induced_color=_INDUCED_COLOR, uninduced_color=_UNINDUCED_COLOR,
                                 title='EM p(induced) distribution by ground-truth label',
                                 save_path=None):
    """Violin + jitter of clade p_induced grouped by GT label."""
    fig, axes = plt.subplots(1, len(timepoints), figsize=(10, 3.2), sharey=True)
    if len(timepoints) == 1:
        axes = [axes]
    for ax, tp in zip(axes, timepoints):
        df = annotated_dfs[tp]
        for label, col, xpos in [('induced', induced_color, 1), ('uninduced', uninduced_color, 0)]:
            vals = df[df['gt_label'] == label]['p_induced'].values
            if len(vals) < 2 or np.std(vals) == 0:
                continue
            parts = ax.violinplot(vals, positions=[xpos], widths=0.6, showmedians=True, showextrema=False)
            for pc in parts['bodies']:
                pc.set_facecolor(col); pc.set_alpha(0.6); pc.set_edgecolor('none')
            parts['cmedians'].set_color(col); parts['cmedians'].set_linewidth(2)
            ax.scatter(np.random.normal(xpos, 0.05, len(vals)), vals,
                       color=col, s=18, alpha=0.7, zorder=3, edgecolors='none')
        ax.axhline(0.5, color='#7f8c8d', linewidth=0.8, linestyle='--')
        ax.set_title(f'Timepoint {tp}')
        ax.set_xticks([0,1]); ax.set_xticklabels(['Uninduced','Induced'], rotation=20, ha='right')
        ax.set_ylim(0, 1)
    axes[0].set_ylabel('p(induced)')
    if title:
        fig.suptitle(title, y=1.02)
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, bbox_inches='tight')
    return fig


def plot_cell_concordance(cell_level_dfs, timepoints, threshold=0.5, n_bins=5,
                          induced_color=_INDUCED_COLOR, uninduced_color=_UNINDUCED_COLOR,
                          title='Cell-level concordance: EM p(induced) vs ground truth',
                          save_path=None):
    """
    3-row x len(timepoints) figure:
      Row 0 — Waterfall bars sorted by p_induced, coloured by GT, wrong predictions outlined
      Row 1 — Calibration reliability diagram
      Row 2 — Confusion matrix at threshold with raw counts and fraction-of-total
    """
    fig, axes = plt.subplots(3, len(timepoints), figsize=(13, 8.5),
                             gridspec_kw={'height_ratios':[2.5,1.5,1.5], 'hspace':0.55, 'wspace':0.3})
    if title:
        fig.suptitle(title, fontsize=10, fontweight='bold', y=1.01)

    for col, tp in enumerate(timepoints):
        cdf    = cell_level_dfs[tp].sort_values('p_induced').reset_index(drop=True)
        y_pred = (cdf['p_induced'] >= threshold).astype(int)
        n      = len(cdf)
        bar_colors = [induced_color if y else uninduced_color for y in cdf['y_true']]

        # row 0: waterfall
        ax = axes[0, col]
        ax.bar(np.arange(n), cdf['p_induced'], color=bar_colors, width=1.0, linewidth=0)
        wrong = cdf.index[y_pred != cdf['y_true']].tolist()
        if wrong:
            ax.bar(wrong, cdf.loc[wrong, 'p_induced'],
                   color=[bar_colors[i] for i in wrong],
                   width=1.0, linewidth=0.6, edgecolor='#111111')
        ax.axhline(threshold, color='#7f8c8d', linewidth=0.8, linestyle='--')
        ax.set_xlim(-1, n); ax.set_ylim(0, 1)
        ax.set_xticks([]); ax.set_title(f'Timepoint {tp}', fontsize=9)
        if col == 0:
            ax.set_ylabel('p(induced)')
        acc = (y_pred == cdf['y_true']).mean()
        ax.text(0.97, 0.95, f'acc={acc:.0%}', transform=ax.transAxes,
                ha='right', va='top', fontsize=7.5, color='#333333')

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
        ax.plot([0,1],[0,1], color='#bdc3c7', linewidth=0.8, linestyle='--', zorder=1)
        ax.scatter(bin_centers, frac_pos, s=[10 + 3*c for c in bin_counts], color='#27ae60', zorder=3)
        ax.plot(bin_centers, frac_pos, color='#27ae60', linewidth=1.3, zorder=2)
        ax.set_xlim(0,1); ax.set_ylim(0,1)
        ax.set_xticks([0,0.5,1]); ax.set_yticks([0,0.5,1])
        ax.tick_params(labelsize=7)
        if col == 0:
            ax.set_ylabel('Frac. truly induced', fontsize=8)
        ax.set_xlabel('Mean p(induced)', fontsize=8)

        # row 2: confusion matrix
        ax = axes[2, col]
        cm      = confusion_matrix(cdf['y_true'], y_pred, labels=[1, 0])
        cm_norm = cm.astype(float) / cm.sum()
        ax.imshow(cm_norm, cmap='Blues', vmin=0, vmax=cm_norm.max()*1.2, aspect='auto')
        for i, j in itertools.product(range(2), range(2)):
            ax.text(j, i, f'{cm[i,j]}\n({cm_norm[i,j]:.0%})',
                    ha='center', va='center', fontsize=7.5,
                    color='white' if cm_norm[i,j] > 0.45 else 'black')
        ax.set_xticks([0,1]); ax.set_xticklabels(['Pred +','Pred −'], fontsize=7.5)
        ax.set_yticks([0,1]); ax.set_yticklabels(['GT +','GT −'], fontsize=7.5)
        if col == 0:
            ax.set_ylabel('Ground truth', fontsize=8)

    legend_handles = [
        Patch(facecolor=induced_color,   label='Induced (GT)'),
        Patch(facecolor=uninduced_color, label='Uninduced (GT)'),
        Patch(facecolor='#aaaaaa', edgecolor='#111111', linewidth=0.8, label='Wrong prediction'),
    ]
    axes[0,-1].legend(handles=legend_handles, loc='upper left', frameon=False, fontsize=7)
    if save_path:
        fig.savefig(save_path, bbox_inches='tight')
    return fig


# ── Param dicts (cell 2) ─────────────────────────────────────────────────────
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



# newdict2 = {'name': '8clade_topological_treedist',
#                  'scheme': 'topological',
#                  'use_branch_length': True,
#                  'num_clades': 8,
#                  'distance_metric': None,
#                  'dist_cutoff': None}

# newdict3 = {'name': '8clade_hierarch_treedist',
#                  'scheme': 'hierarchical',
#                  'use_branch_length': True,
#                  'num_clades': 8,
#                  'distance_metric': 'tree_dist',
#                  'dist_cutoff': None,
#                  'include_in_plot': True}

# newdict4 = {'name': '12clade_hierarch_treedist',
#                  'scheme': 'hierarchical',
#                  'use_branch_length': True,
#                  'num_clades': 12,
#                  'distance_metric': 'tree_dist',
#                  'dist_cutoff': None,
#                  'include_in_plot': True}

# ex_param_dict_gt = {'name': '12clade_topological_treedist',
#                  'scheme': 'topological',
#                  'use_branch_length': True,
#                  'num_clades': 12,
#                  'distance_metric': None,
#                  'dist_cutoff': None}


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



# param_dict_list = [ex_param_dict1, ex_param_dict2, ex_param_dict3,
#                   ex_param_dict4, ex_param_dict5, ex_param_dict6]'

# this is what i was using
# param_dict_list = [ex_param_dict1, ex_param_dict2, ex_param_dict4, ex_param_dict5]

# scoremat_param_dict_list = [newdict1, newdict2, newdict3, newdict4, ex_param_dict_gt]

# scoremat_param_dict_list = [newdict1, ex_param_dict_gt, ex_param_dict_gt2]
# scoremat_param_dict_list = [newdict1, ex_param_dict_gt2]
# scoremat_param_dict_list = [newdict1]
# scoremat_param_dict_list = [param_dict_gt20]
scoremat_param_dict_list = [scoremat_min5, scoremat_min10, scoremat_min20]

only_topological_clades = [param_dict_gt5, param_dict_gt10, param_dict_gt20]

# longmito_36_first_pass_list = [scoremat_min10, param_dict_gt10]
# longmito_36_first_pass_list = [scoremat_min5, scoremat_min10, scoremat_min20, param_dict_gt5, param_dict_gt20]

longmito_36_first_pass_list = [scoremat_min5, scoremat_min10, scoremat_min20, param_dict_gt5, param_dict_gt20]



param_dict_list_100clades = [param_dict_gt100]

param_dict_list_50clades = [param_dict_gt50]

param_dict_list_2clades = [param_dict_gt2]


param_dict_list_hsc = [param_dict_gt2, param_dict_gt20, param_dict_gt50]

param_dict_list_rand10 = [random_dict_10clades, param_dict_gt2, param_dict_gt20, param_dict_gt50, scoremat_min10, scoremat_min20, scoremat_min5, scoremat_min1]

param_dict_score_only = [scoremat_min10, scoremat_min20, scoremat_min5, scoremat_min1, scoremat_2clades, scoremat_10clades, scoremat_20clades]

# score_and_20top = 


# ── Parallel worker functions ─────────────────────────────────────────────────
# These must be module-level (not nested) so ProcessPoolExecutor can pickle them.

def _worker_process_urid(urid):
    """Phase-1 worker: load & cluster one urid."""
    sys.setrecursionlimit(100000)
    global path_to_output_dir, path_to_em_induction_prob
    path_to_output_dir = Path('./output')
    path_to_em_induction_prob = path_to_output_dir / 'em_induction_prob'
    timepoint_res_dict, _ = process_urid(
        urid=urid,
        param_dict_list=longmito_36_first_pass_list,
        cluster_vmin=-0.5,
        cluster_vmax=1,
        timepoints=[8, 9, 10, 11, 12],
    )
    return urid, timepoint_res_dict


def _worker_inference(args):
    """Phase-3 worker: run EM inference for one (urid, timepoint) pair."""
    sys.setrecursionlimit(100000)
    urid, timepoint = args
    global path_to_output_dir, path_to_em_induction_prob
    path_to_output_dir = Path('./output')
    path_to_em_induction_prob = path_to_output_dir / 'em_induction_prob'

    print(f'[EM start] urid={urid} tp={timepoint}', flush=True)

    _hsc_enc = {'LT': 0, 'ST': 1, 'MyMPP': 2, 'LyMPP': 3, 'Mye': 4, 'Lym': 5}

    tpt_induction_dict, opt_mod_res, eval_dict = urid_to_induction_inference(
        [urid], [timepoint],
        param_dict_list=longmito_36_first_pass_list,
        weight_strategy='uniform',
        clade_priors=None,
        upper_rarity_bound=1.5,
        include_param_combos=None,
        scoremat_filters=scoremat_filters,
        verbose=False,
        root_type='LT',
        num_restarts=1,
        return_opt_model_dict=True,
        type_encoder=_hsc_enc,
        enrichment_dict=None,
        strength=10,
        use_gmm=True,
        use_cp_class=True,
        path_to_output_dir=path_to_output_dir,
    )
    open(f'./progress_txt_files/progress_{urid}_{timepoint}.txt', 'w').close()
    # Strip tree references from EM models before pickling across the queue.
    # FULL_make_all_em_plots only uses tpt_induction_dict and eval_dict;
    # opt_mod_res is saved for inspection but the tree is not accessed again.
    for tp_dict in opt_mod_res.values():
        for em in tp_dict.values():
            if hasattr(em, 'tree'):
                em.tree = None
    return urid, timepoint, tpt_induction_dict, opt_mod_res, eval_dict


# ─────────────────────────────────────────────────────────────────────────────
# Alpha-sweep helpers and ProcessPoolExecutor worker
# ─────────────────────────────────────────────────────────────────────────────

_ALPHA_GT_COMBO = {
    "gt5":  "5clade_topological_treedist_gt",
    "gt20": "20clade_topological_treedist_gt",
}
_ALPHA_SCOREMAT_COMBO = {
    "scoremat_min5":  "scoremat_min5_mt_250ints_1.0RP_0.0AF_TB_1-1-1-1-1-1",
    "scoremat_min10": "scoremat_min10_mt_250ints_1.0RP_0.0AF_TB_1-1-1-1-1-1",
    "scoremat_min20": "scoremat_min20_mt_250ints_1.0RP_0.0AF_TB_1-1-1-1-1-1",
}


def _alpha_get_tms(urid, output_dir=None):
    """Return (Phi0, Phi1_true, founder_type, type_order) from run_specs JSON."""
    od = Path(output_dir) if output_dir else Path("./output")
    candidates = [
        od / str(urid) / "run_specs" / str(urid),
        od / str(urid) / "run_specs",
        od / "run_specs" / str(urid),
    ]
    jsons = []
    for p in candidates:
        jsons = list(p.glob("*.json"))
        if jsons:
            break
    if not jsons:
        raise FileNotFoundError(f"No run_specs JSON for urid={urid}; tried {candidates}")
    d  = json.loads(jsons[0].read_text())
    ct = d["cell_type_dict"]
    founder    = ct.get("founder_cell_type", "ct1")
    type_order = list(ct["cell_type_params"].keys())
    return (np.array(ct["uninduced_transition_matrix"], dtype=float),
            np.array(ct["induced_transition_matrix"],   dtype=float),
            founder, type_order)


def _alpha_load_clade_df(urid, timepoint, output_dir=None):
    od = Path(output_dir) if output_dir else Path("./output")
    path = od / f"clade_assignment_dfs/{urid}/{timepoint}/clade_df_{urid}_t{timepoint}.csv"
    df = pd.read_csv(path, index_col=0)
    df.index.name = "cell"
    return df.reset_index()


def _alpha_get_gt_tree_path(urid, timepoint, output_dir=None):
    od = Path(output_dir) if output_dir else Path("./output")
    urid_dir = od / f"processed_newicks/{urid}"
    matches = list(urid_dir.glob(f"ground_truth_tree_*_time_{timepoint}.newick"))
    if not matches:
        raise FileNotFoundError(f"GT tree not found for urid={urid} tp={timepoint}")
    return str(matches[0])


def _alpha_get_recon_tree_path(urid, timepoint, output_dir=None):
    od = Path(output_dir) if output_dir else Path("./output")
    candidates = [
        od / "recon_trees" / str(urid),                 # cluster: output/recon_trees/{urid}/
        od / str(urid) / "recon_trees" / str(urid),     # local: output/{urid}/recon_trees/{urid}/
    ]
    glob_pat = f"*RP_1_samp*_time_{timepoint}_CD_T_AF_0_B_T/*.treefile"
    for recon_dir in candidates:
        matches = list(recon_dir.glob(glob_pat))
        if matches:
            return str(matches[0])
    return None


def _alpha_build_leaf_maps(clade_df, param_combo):
    if param_combo not in clade_df.columns:
        return None, None, None
    sub = clade_df[
        ["cell", param_combo, "clade_celltype", "clade_differentiation_induced"]
    ].dropna(subset=[param_combo, "clade_celltype"])
    if len(sub) < 2:
        return None, None, None
    leaf2clade      = dict(zip(sub["cell"], sub[param_combo].astype(str)))
    leaf2type       = dict(zip(sub["cell"], sub["clade_celltype"].astype(str)))
    leaf2induced_gt = {
        cell: ("induced" if bool(val) else "un-induced")
        for cell, val in zip(sub["cell"], sub["clade_differentiation_induced"])
    }
    return leaf2clade, leaf2type, leaf2induced_gt


def _worker_alpha_sweep(args):
    """ProcessPoolExecutor worker: alpha sweep for one (urid, timepoint, param_combo)."""
    from sklearn.metrics import roc_auc_score as _roc_auc, auc as _auc, precision_recall_curve as _prc
    (urid, timepoint, param_combo, tree_path,
     leaf2clade, leaf2type, leaf2induced_gt,
     Phi0, Phi1_true, alpha_values, type_encoder,
     n_restarts, root_type, weight_strategy) = args

    # load and prune tree (done once; all alphas reuse it)
    raw_newick = "".join(open(tree_path).read().split())
    clean_nwk  = re.sub(r"\)([^:);]+):", r"):", raw_newick)
    tree = Tree(clean_nwk, format=1)

    keep = set(leaf2clade.keys()) & {n.name for n in tree.iter_leaves()}
    if len(keep) < 4:
        print(f"[alpha_sweep] {urid} tp={timepoint} {param_combo}: {len(keep)} leaves, skip", flush=True)
        return []
    tree.prune(list(keep), preserve_branch_length=True)

    leaf2clade      = {c: v for c, v in leaf2clade.items()      if c in keep}
    leaf2type       = {c: v for c, v in leaf2type.items()       if c in keep}
    leaf2induced_gt = {c: v for c, v in leaf2induced_gt.items() if c in keep}

    clades      = sorted(set(leaf2clade.values()))
    clade_alpha = {c: 1.0 for c in clades}
    clade_beta  = {c: 1.0 for c in clades}

    cell_names = list(leaf2clade.keys())
    gt_arr = np.array([1 if leaf2induced_gt[c] == "induced" else 0 for c in cell_names])

    print(f"[alpha_sweep] {urid} tp={timepoint} {param_combo}: "
          f"{len(cell_names)} cells, {len(alpha_values)} α values", flush=True)

    try:
        alpha_results = run_em_with_restarts(
            Phi0=Phi0, Phi1=Phi0,           # Phi1 unused when alpha_values provided
            root_type=root_type, tree=tree,
            leaf2clade=leaf2clade, leaf2type=leaf2type,
            type_encoder=type_encoder, weight_strategy=weight_strategy,
            clade_alpha=clade_alpha, clade_beta=clade_beta,
            n_restarts=n_restarts, max_iter=200, tol=1e-3,
            seed_start=42, verbose=False, use_cp_class=True,
            alpha_values=alpha_values, Phi1_true=Phi1_true,
        )
    except Exception as exc:
        print(f"[alpha_sweep] EM error {urid} tp={timepoint} {param_combo}: {exc}", flush=True)
        return []

    rows = []
    for alpha, best in alpha_results.items():
        em      = best["model"]
        df_post = em.clade_posteriors(use_gmm=False)

        clade_to_p  = dict(zip(df_post["clade"].astype(str), df_post["p_induced"]))
        clade_to_hl = dict(zip(df_post["clade"].astype(str), df_post["hard_label"]))

        pred_prob = np.array([clade_to_p.get(str(leaf2clade[c]), 0.5) for c in cell_names])
        pred_hard = np.array([
            1 if clade_to_hl.get(str(leaf2clade[c]), "un-induced") == "induced" else 0
            for c in cell_names
        ])

        if len(np.unique(gt_arr)) < 2:
            roc_auc_val = float("nan")
            pr_auc_val  = float("nan")
        else:
            roc_auc_val = _roc_auc(gt_arr, pred_prob)
            if roc_auc_val < 0.5:
                pred_prob   = 1.0 - pred_prob
                pred_hard   = 1   - pred_hard
                roc_auc_val = 1.0 - roc_auc_val
            prec_c, rec_c, _ = _prc(gt_arr, pred_prob, pos_label=1)
            pr_auc_val = _auc(rec_c, prec_c)

        tp_n = int(np.sum((pred_hard == 1) & (gt_arr == 1)))
        fp_n = int(np.sum((pred_hard == 1) & (gt_arr == 0)))
        tn_n = int(np.sum((pred_hard == 0) & (gt_arr == 0)))
        fn_n = int(np.sum((pred_hard == 0) & (gt_arr == 1)))

        prec_h = tp_n / (tp_n + fp_n) if (tp_n + fp_n) > 0 else 0.0
        rec_h  = tp_n / (tp_n + fn_n) if (tp_n + fn_n) > 0 else 0.0
        f1_h   = 2 * prec_h * rec_h / (prec_h + rec_h) if (prec_h + rec_h) > 0 else 0.0

        best_logL = best.get("final_logL", best.get("logL", float("nan")))

        rows.append(dict(
            urid=urid, timepoint=timepoint, param_combo=param_combo, alpha=alpha,
            roc_auc=round(roc_auc_val, 6),
            pr_auc=round(pr_auc_val, 6),
            f1_hard=round(f1_h, 6),
            precision=round(prec_h, 6),
            recall=round(rec_h, 6),
            tp=tp_n, fp=fp_n, tn=tn_n, fn=fn_n,
            pred_pos_frac=round((tp_n + fp_n) / max(tp_n + fp_n + tn_n + fn_n, 1), 4),
            n_cells=len(cell_names),
            n_clades=len(df_post),
            best_logL=(round(best_logL, 4) if not np.isnan(best_logL) else float("nan")),
        ))

    print(f"[alpha_sweep done] {urid} tp={timepoint} {param_combo} → {len(rows)} rows", flush=True)
    return rows


# ── Main (cell 29) ───────────────────────────────────────────────────────────
# i want to export this cell as a script

# long_mito_36_urids = [8581181233557, 1219288324872]

# long_mito_36_urids = [
#     4807145605041,
#     5174853967582,
#     8823915356286,
#     4856020464336,
#     4126344089951,
#     9796194427324,
#     2852130303988,
#     8605119805741,
#     9547253106784,
#     2412930989432,
#     8581181233557,
#     1219288324872,
#     3400966091602,
#     4621438212904,
#     3255685044575,
#     7073674028991,
#     4971152140979,
#     924078409061,
#     6290262748470,
#     7778854110315,
#     3483020099928,
#     5877308073148,
#     6146677390098,
#     6513792620415,
#     5122970731735,
#     9172933209610,
#     3132521060600,
#     8030859708906,
#     9493229458858,
#     3199712236239,
#     9886194576482,
#     1400131950216,
#     9793329286434,
#     4699787243590,
#     4397598422406,
#     1949496656499
# ]

from concurrent.futures import ProcessPoolExecutor, as_completed

# ── Tunable: set to the number of CPU cores you want to use (None = all cores) ──
N_WORKERS = 45   # set to e.g. 8 to cap usage

# Subset of scoremat RP/AF combinations to run EM on.
# Each string is matched as a substring of the scoremat combo name.
# GT combos (ending in _gt) always run regardless of this filter.
# Set to None to run all scoremats generated by longmito_36_first_pass_list.
# scoremat_filters = ['1.0RP_0.0AF']
# scoremat_filters = ['scoremat_min5_mt_250ints_1.0RP_0.0AF']
# scoremat_filters = ['scoremat_min5_mt_250ints_1.0RP_0.0AF', 'clade_topological_treedist_mt']
scoremat_filters = ['scoremat_min5_mt_250ints_1.0RP_0.0AF', 
                    'scoremat_min5_mt_250ints_1.0RP_0.01AF',
                    'scoremat_min5_mt_250ints_1.0RP_0.05AF',
                    'clade_topological_treedist_mt_250ints_1.0RP_0.0AF']


gt_combos = ['20clade_topological_treedist_gt', '5clade_topological_treedist_gt']
mt_scoremat_clades = ['5clade_topological_treedist_mt_250ints_1.0RP_0.0AF_TB_1-1-1-1-1-1',
             '20clade_topological_treedist_mt_250ints_1.0RP_0.0AF_TB_1-1-1-1-1-1',
             'scoremat_min5_mt_250ints_1.0RP_0.0AF_TB_1-1-1-1-1-1',
             'scoremat_min5_mt_250ints_1.0RP_0.01AF_TB_1-1-1-1-1-1',
             'scoremat_min5_mt_250ints_1.0RP_0.05AF_TB_1-1-1-1-1-1']
# # gt_combos = []
# mt_scoremat_clades = [
#     'scoremat_min5_mt_250ints_1.0RP_0.0AF_TB_1-1-1-1-1-1',
#     # 'scoremat_min10_mt_250ints_1.0RP_0.0AF_TB_1-1-1-1-1-1',
#     # 'scoremat_min20_mt_250ints_1.0RP_0.0AF_TB_1-1-1-1-1-1',
# ]

# mt_scoremat_clades = []
# all_param_combos = gt_combos + mt_scoremat_clades

TIMEPOINTS = [8, 9, 10, 11, 12]

if __name__ == '__main__':

    # ── CLI ──────────────────────────────────────────────────────────────────
    parser = argparse.ArgumentParser(
        prog='longmito_36.py',
        description='Run clade assignment and/or EM induction inference.',
    )
    parser.add_argument(
        '--mode',
        choices=['full', 'scoremat_sweep', 'no_em', 'alpha_sweep'],
        default='full',
        help=(
            'full          — run all phases (process_urid → EM → plots) [default]\n'
            'scoremat_sweep — run scoremat_clade_sweep only, no trees, no EM\n'
            'no_em         — run phases 1-2 only (process_urid + make_clade_df); stop before EM\n'
            'alpha_sweep   — sweep Phi1 prior strength via α; runs EM in parallel with ProcessPoolExecutor'
        ),
    )
    parser.add_argument(
        '--urids',
        nargs='+',
        type=int,
        default=None,
        metavar='URID',
        help='Space-separated list of urids to process (default: long_mito_36_urids)',
    )
    parser.add_argument(
        '--timepoints',
        nargs='+',
        type=int,
        default=None,
        metavar='T',
        help='Timepoints to include (default: TIMEPOINTS constant = 8 9 10 11 12)',
    )
    # scoremat_sweep-specific
    parser.add_argument(
        '--min-clade-sizes',
        nargs='+',
        type=int,
        default=[5, 10, 20],
        metavar='N',
        help='(scoremat_sweep only) min_clade_size values to sweep (default: 5 10 20)',
    )
    parser.add_argument(
        '--output',
        default='scoremat_clade_sweep_df.csv',
        metavar='PATH',
        help='(scoremat_sweep only) path for output CSV (default: scoremat_clade_sweep_df.csv)',
    )
    parser.add_argument(
        '--clade-output',
        default=None,
        metavar='PATH',
        help='Path for combined clade-assignments CSV (written after Phase 2, before EM, '
             'in any mode). If omitted, no CSV is saved unless --mode no_em also sets --output.',
    )
    parser.add_argument(
        '--joblib-dir',
        default='.',
        metavar='DIR',
        help='Directory to write joblib checkpoint and final result files (default: current dir)',
    )
    parser.add_argument(
        '--joblib-prefix',
        default='longmito_results',
        metavar='PREFIX',
        help='Filename stem for joblib outputs (default: longmito_results). '
             'Files are written as <DIR>/<PREFIX>_model_dict_ckpt_N.joblib, '
             '<DIR>/<PREFIX>_model_dict_FINAL.joblib, <DIR>/<PREFIX>_roc_dict_FINAL.joblib',
    )
    # alpha_sweep-specific arguments
    parser.add_argument(
        '--alpha-values',
        nargs='+',
        type=float,
        default=[0.0, 0.1, 0.25, 0.5, 0.75, 1.0],
        metavar='A',
        help='(alpha_sweep only) α values to sweep (default: 0 0.1 0.25 0.5 0.75 1.0). '
             'α=0 means uninformative prior, α=1 means exact prior.',
    )
    parser.add_argument(
        '--urid-file',
        default=None,
        metavar='PATH',
        help='(alpha_sweep only) text file with one urid per line; '
             'lines starting with # are ignored (alternative to --urids).',
    )
    parser.add_argument(
        '--param-combos',
        nargs='+',
        default=['gt5', 'scoremat_min5'],
        choices=['gt5', 'gt20', 'scoremat_min5', 'scoremat_min10', 'scoremat_min20'],
        metavar='COMBO',
        help='(alpha_sweep only) clade-definition schemes to run '
             '(default: gt scoremat_min5).',
    )
    parser.add_argument(
        '--n-restarts',
        type=int,
        default=3,
        metavar='N',
        help='(alpha_sweep only) EM random restarts per α value (default: 3).',
    )
    parser.add_argument(
        '--alpha-out',
        default='alpha_sweep_results.csv',
        metavar='PATH',
        help='(alpha_sweep only) output CSV path (default: alpha_sweep_results.csv).',
    )
    args = parser.parse_args()

    # resolve urids — alpha_sweep also accepts --urid-file
    if args.urid_file:
        with open(args.urid_file) as _fh:
            _file_urids = [int(l.strip()) for l in _fh
                           if l.strip() and not l.strip().startswith('#')]
        _urids = _file_urids
    else:
        _urids = args.urids if args.urids is not None else long_mito_36_urids

    _timepoints = args.timepoints if args.timepoints is not None else TIMEPOINTS

    import os as _os
    _joblib_dir = args.joblib_dir
    _os.makedirs(_joblib_dir, exist_ok=True)
    _joblib_prefix = args.joblib_prefix

    # ── scoremat_sweep mode: no trees, no EM ─────────────────────────────────
    if args.mode == 'scoremat_sweep':
        print(
            f'=== scoremat_sweep  urids={len(_urids)}  '
            f'timepoints={_timepoints}  min_clade_sizes={args.min_clade_sizes} ==='
        )
        sweep_df = scoremat_clade_sweep(
            urid_list       = _urids,
            min_clade_sizes = args.min_clade_sizes,
            timepoints      = _timepoints,
        )
        sweep_df.to_csv(args.output, index=False)
        print(f'Saved {len(sweep_df):,} rows to {args.output}')
        sys.exit(0)

    # ── alpha_sweep mode: parallel EM over α values ──────────────────────────
    if args.mode == 'alpha_sweep':
        import gc as _gc
        from pathlib import Path as _Path

        _alpha_out = _Path(args.alpha_out)
        _alpha_values = args.alpha_values
        _param_combos = args.param_combos
        _n_restarts   = args.n_restarts

        print(
            f'=== alpha_sweep  urids={len(_urids)}  tp={_timepoints}  '
            f'α={_alpha_values}  combos={_param_combos}  '
            f'n_restarts={_n_restarts}  workers={N_WORKERS} ==='
        )

        # build task list by loading TMs and leaf maps in the main process
        # (tree loading happens inside each worker to avoid pickling ete3 trees)
        tasks      = []
        skip_count = 0

        for _urid in _urids:
            try:
                _Phi0, _Phi1_true, _founder, _type_order = _alpha_get_tms(str(_urid))
            except Exception as _e:
                print(f'[alpha_sweep] Could not load TMs for {_urid}: {_e}')
                skip_count += len(_timepoints) * len(_param_combos)
                continue

            _type_encoder = {t: i for i, t in enumerate(_type_order)}

            for _tp in _timepoints:
                try:
                    _clade_df = _alpha_load_clade_df(str(_urid), _tp)
                except FileNotFoundError as _e:
                    print(f'[alpha_sweep] clade_df not found {_urid} tp={_tp}: {_e}')
                    skip_count += len(_param_combos)
                    continue

                for _pc_key in _param_combos:
                    if _pc_key in ('gt5', 'gt20'):
                        _param_combo_full = _ALPHA_GT_COMBO[_pc_key]
                        try:
                            _tree_path = _alpha_get_gt_tree_path(str(_urid), _tp)
                        except FileNotFoundError as _e:
                            print(f'[alpha_sweep] GT tree not found {_urid} tp={_tp}: {_e}')
                            skip_count += 1
                            continue
                    else:
                        _param_combo_full = _ALPHA_SCOREMAT_COMBO.get(_pc_key)
                        _tree_path = _alpha_get_recon_tree_path(str(_urid), _tp)
                        if _tree_path is None:
                            print(f'[alpha_sweep] recon tree not found {_urid} tp={_tp} {_pc_key}')
                            skip_count += 1
                            continue

                    _l2c, _l2t, _l2igt = _alpha_build_leaf_maps(_clade_df, _param_combo_full)
                    if _l2c is None:
                        print(f'[alpha_sweep] column absent {_urid} tp={_tp} {_pc_key}')
                        skip_count += 1
                        continue
                    if len(_l2c) < 4:
                        print(f'[alpha_sweep] too few cells {_urid} tp={_tp} {_pc_key}')
                        skip_count += 1
                        continue

                    tasks.append((
                        str(_urid), _tp, _param_combo_full, _tree_path,
                        _l2c, _l2t, _l2igt,
                        _Phi0, _Phi1_true, _alpha_values, _type_encoder,
                        _n_restarts, _founder, 'leaf_prop',
                    ))

                del _clade_df
                _gc.collect()

        print(f'[alpha_sweep] {len(tasks)} tasks queued, {skip_count} skipped')

        _n_rows = 0
        with ProcessPoolExecutor(max_workers=N_WORKERS) as _pool:
            _futures = {_pool.submit(_worker_alpha_sweep, t): t for t in tasks}
            for _fut in tqdm.tqdm(as_completed(_futures), total=len(tasks), desc='alpha_sweep'):
                try:
                    _rows = _fut.result()
                except Exception as _exc:
                    _task = _futures[_fut]
                    print(f'ERROR for {_task[:3]}: {_exc}')
                    continue
                if _rows:
                    _chunk = pd.DataFrame(_rows)
                    _chunk.to_csv(
                        _alpha_out, mode='a', index=False,
                        header=not _alpha_out.exists(),
                    )
                    _n_rows += len(_rows)

        print(f'[alpha_sweep] done — {_n_rows} rows written to {_alpha_out}')
        sys.exit(0)

    # ── Phase 1: load & cluster all urids (serial) ───────────────────────────
    # Parallelism is not viable here: TimepointRes objects contain deep ete3 Tree
    # objects whose pickle recursion depth exceeds any practical sys limit, causing
    # _ForkingPickler to crash when the worker tries to send results back over the
    # process queue.  Phase 3 (EM inference) is where the real compute lives.
    print('=== Phase 1: process_urid (serial) ===')
    longmito_36_resdict = {}
    for urid in tqdm.tqdm(_urids):
        path_to_output_dir = Path('./output')
        path_to_em_induction_prob = path_to_output_dir / 'em_induction_prob'
        timepoint_res_dict, _ = process_urid(
            urid=urid,
            param_dict_list=longmito_36_first_pass_list,
            cluster_vmin=-0.5,
            cluster_vmax=1,
            timepoints=_timepoints,
        )
        longmito_36_resdict[urid] = timepoint_res_dict

    # ── Phase 2: make_clade_df (fast; keep serial) ───────────────────────────
    print('=== Phase 2: make_clade_df ===')
    for urid, timepoint_res_dict in longmito_36_resdict.items():
        for tp, tpr in timepoint_res_dict.items():
            make_clade_df(timepoint_res=tpr, urid=urid, timepoint=tp)

    # Save clade assignments CSV if requested (works in any mode)
    _clade_csv_path = args.clade_output or (args.output if args.mode == 'no_em' else None)
    if _clade_csv_path:
        print('=== Building combined clade DataFrame ===')
        records = []
        for urid, tpr_dict in longmito_36_resdict.items():
            for tp, tpr in tpr_dict.items():
                for linstring, celldict in tpr.cell_pop.items():
                    for param_combo, clade in celldict.get('clade_assignments', {}).items():
                        records.append((urid, tp, param_combo, linstring, clade))
        combined_df = pd.DataFrame(
            records,
            columns=['urid', 'timepoint', 'param_combo', 'linstring', 'clade']
        )
        combined_df.to_csv(_clade_csv_path, index=False)
        print(f'Saved {len(combined_df):,} rows to {_clade_csv_path}')

    if args.mode == 'no_em':
        sys.exit(0)

    # ── Phase 3: EM inference over all (urid, timepoint) pairs in parallel ───
    print('=== Phase 3: urid_to_induction_inference (parallel) ===')
    tasks = [
        (urid, tp)
        for urid in _urids
        for tp in _timepoints
    ]

    model_longmito_36_dict  = {urid: {} for urid in _urids}
    roc_longmito_36_dict    = {urid: {} for urid in _urids}
    completed = 0

    with ProcessPoolExecutor(max_workers=N_WORKERS) as pool:
        futures = {pool.submit(_worker_inference, task): task for task in tasks}
        for future in tqdm.tqdm(as_completed(futures), total=len(tasks)):
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
            print(f'[EM done]  urid={urid} tp={timepoint} roc={roc_summary}', flush=True)

            completed += 1
            # checkpoint every 50 completed tasks
            if completed % 50 == 0:
                joblib.dump(model_longmito_36_dict, f'{_joblib_dir}/{_joblib_prefix}_model_dict_ckpt_{completed}.joblib')
                # with open(f'./no_gt_combos_model_longmito_36_dict_ckpt_{completed}.pkl', 'wb') as f:
                #     pickle.dump(model_longmito_36_dict, f)

    # ── Save final results ────────────────────────────────────────────────────
    joblib.dump(model_longmito_36_dict, f'{_joblib_dir}/{_joblib_prefix}_model_dict_FINAL.joblib')
    joblib.dump(roc_longmito_36_dict,   f'{_joblib_dir}/{_joblib_prefix}_roc_dict_FINAL.joblib')
    # with open('./no_gt_combos_model_longmito_36_dict_FINAL.pkl', 'wb') as f:
    #     pickle.dump(model_longmito_36_dict, f)
    # with open('./no_gt_combos_roc_longmito_36_dict_FINAL.pkl', 'wb') as f:
    #     pickle.dump(roc_longmito_36_dict, f)

    # ── Phase 4: plots (serial; needs all results) ────────────────────────────
    print('=== Phase 4: FULL_make_all_em_plots ===')
    plotting_dfs = FULL_make_all_em_plots(
        param_combo_list=gt_combos + mt_scoremat_clades,
        full_modelres_dict=model_longmito_36_dict,
        urid_timepoint_res_dict=longmito_36_resdict,
        return_plotting_dfs=True,
    )



