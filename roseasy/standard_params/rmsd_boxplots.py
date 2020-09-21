#!/usr/bin/env python2
# encoding: utf-8

"""\
Visualize the results from the loop modeling simulations in PIP and identify 
promising designs.

Usage:
    mean_boxplots.py <ref_dir> <pdb_directories>... [options]

Options:
    -F, --fork
        Do not fork into a background process.

    -f, --force
        Force the cache to be regenerated.

    -q, --quiet
        Build the cache, but don't launch the GUI.


"""

import os, glob, numpy as np, random
from scipy.stats import ttest_ind, ks_2samp, anderson_ksamp
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from roseasy import pipeline
from roseasy import structures
from roseasy import gui


def select_random_subset(pdb_paths, sample_size):
    subset = random.choices(pdb_paths, k=sample_size)

    return subset


def main():
    import docopt
    args = docopt.docopt(__doc__)

    #print(args)

    pdb_dirs = args['<pdb_directories>']
    ref_dir = args['<ref_dir>']
    #num_samples = int(args['--samples'])
    #sample_size = int(args['--sample_size'])

    #num_samples = 50
    #sample_sizes = [1, 5, 10, 20, 30, 50, 100, 150, 200,]

    
    #prop = 'ca_rmsd_no_loop'
    #prop = 'ca_rmsd'
    #ref_prop = 'ca_rmsd'

    models = {}

    ref_records, ref_metadata = structures.load(ref_dir, use_cache=not args['--force'])
    ref_records['label_name'] = 'WT'
    ref_records['model_name'] = 'WT'
    models[ref_dir] = ref_records

    for pdb_dir in pdb_dirs:
        #print(pdb_dir)

        try:
            all_records, metadata = structures.load(pdb_dir, use_cache=not args['--force'])
        except:
            print("No structures in:", pdb_dir, "... continuing")
            continue

        #print(pdb_dir, len(all_records))
        if 'ca_rmsd_no_loop' in all_records.keys():
            prop = 'ca_rmsd_no_loop'
        else:
            all_records = all_records.rename(columns={"ca_rmsd": "ca_rmsd_no_loop"})
            prop = 'ca_rmsd_no_loop' 

        p_val= ttest_ind(all_records[prop], ref_records[prop], equal_var = False)[1]
        AD_p = anderson_ksamp([ref_records[prop],all_records[prop]])[2]
        dir_split = pdb_dir.split("/")
        model_name = dir_split[0]
        label_name = "\n".join([dir_split[0], dir_split[-1], "T p= " + "{:.2e}".format(p_val), "AD p= " + "{:.2e}".format(AD_p)])

        #print(name)
        all_records['name'] = pdb_dir
        all_records['model_name'] = model_name
        all_records['label_name'] = label_name

        print(pdb_dir,":   ", ttest_ind(all_records[prop], ref_records[prop])[1])

        #print(list(all_records.columns))

        models[pdb_dir] = all_records


    models = pd.concat(models)
    print(models)
    #order = [""]
    #fig, axes = plt.subplots(nrows=1, ncols=1)
    #bp = sns.violinplot(x=ref_records['name'], y=ref_records['ca_rmsd_no_loop'],  )
    bp = sns.violinplot(y=models['label_name'], x=models[prop], inner= 'point'  )
    #bp = sns.catplot(y=models['label_name'], x=models['ca_rmsd_no_loop'], data = models, kind = "violin",   )
    #bp.set_xticklabels(bp.get_xticklabels(), wrap = True)
    #bp = models.boxplot(column='ca_rmsd_no_loop', by='name',)


    plt.show()


    

    # fig, axes = plt.subplots(nrows=1, ncols=2)

    # #all_p_vals = all_p_vals.transpose()
    # all_p_vals.plot(kind = 'line', ax=axes[0], marker = 'o', legend = False)
    # all_diff_props.plot(kind='bar', ax = axes[1], legend = False)

    # #plt.legend()
    # plt.show()


    #     ax1.plot(sample_sizes, p_vals,'--o', label = pdb_dir)
    #     bins = [x + width for x in prev]
    #     ax2.bar(bins , diff_counts, label = pdb_dir)
    #     prev = bins
    # ax1.set_ylabel('K-S p value')
    # ax1.set_xlabel('relax sample size')
    # ax2.set_ylabel('Rejection Proportion')
    # ax2.set_xticks(x)
    # ax2.set_xticklabels(sample_sizes)
    # ax2.set_xlabel('relax sample size')
    # ax1.legend()
    # ax2.legend()
    # #plt.xlabel('# relax structures sampled')
    # #plt.ylabel('K-S p value')
    # plt.show()

if __name__=='__main__':
    main()
