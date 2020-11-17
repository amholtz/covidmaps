from collections import defaultdict

import numpy as np
import pandas as pd
from pastml.tree import read_tree, remove_certain_leaves
import random

if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--input_tree_mutations', default='C:/Users/Andrew/Dropbox/masters_project/NY_alignment/tree_reroot_collapsed.treefile', type=str) # you need to put the path to your own tree for input and output
    parser.add_argument('--input_locs', default='C:/Users/Andrew/Dropbox/masters_project/initial_maps/metadata_0.tab', type=str) # metadata which includes IDs and locations
    parser.add_argument('--input_counts', default='C:/Users/Andrew/PycharmProjects/masters_project/25April2020_perstate.csv', type=str)  # state case counts since April 25, 2020
    parser.add_argument('--output_stats', required=False, type=str, default='C:/Users/Andrew/Dropbox/masters_project/mutations.table')

    params = parser.parse_args()

    # phylogenetic tree read from the newick file
    tree_mutations = read_tree(params.input_tree_mutations)

    def tip_mutations(tip, length):
        mutations = tip.dist * length
        parent = tip.up
        while parent.up: #if there is no parent, none will be treated as false, any non value is treated as true
            mutations += parent.dist * length
            parent = parent.up
        return mutations

    def num_mutations(tree, length):
        mutations = {}
        for tip in tree:
            mutations[tip.name] = tip_mutations(tip,length)
        return mutations
    print("This is mutation dictionary")


    mut_dic = num_mutations(tree_mutations, 29891)

    df_mut = pd.DataFrame.from_dict(mut_dic, orient='index', columns= ["mutations"])

    df_mut.to_csv(params.output_stats, sep='\t', index_label='index')

