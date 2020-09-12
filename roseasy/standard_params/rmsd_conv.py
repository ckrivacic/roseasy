#!/usr/bin/env python2
# encoding: utf-8

"""\
Visualize the results from the loop modeling simulations in PIP and identify 
promising designs.

Usage:
    rmsd_conv.py <ref_dir> <pdb_directories>... [options]

Options:
    -F, --fork
        Do not fork into a background process.

    -f, --force
        Force the cache to be regenerated.

    -q, --quiet
        Build the cache, but don't launch the GUI.

    --samples SAMPLES, -s SAMPLES [default: 10]
        number of pdb samples to randomly select

This command launches a GUI designed to visualize the results for the loop 
modeling simulations in PIP and to help you identify promising designs.  To 
this end, the following features are supported:

1. Extract quality metrics from forward-folded models and plot them against 
   each other in any combination.

2. Easily visualize specific models by right-clicking on plotted points.  
   Add your own visualizations by writing `*.sho' scripts.

3. Plot multiple designs at once, for comparison purposes.

4. Keep notes on each design, and search your notes to find the designs you 
   want to visualize.

Generally, the only arguments you need are the names of one or more directories 
containing the PDB files you want to look at.  For example:

    $ ls -R
    .:
    design_1  design_2 ...

    ./design_1:
    model_1.pdb  model_2.pdb ...

    ./design_2:
    model_1.pdb  model_2.pdb ...

    $ pull_into_place plot_funnels design_*

This last command will launch the GUI.  If you specified more than one design 
on the command line, the GUI will have a panel on the left listing all the 
designs being compared.  You can control what is plotted by selecting one or 
more designs from this list.  The search bar at the top of this panel can be 
used to filter the list for designs that have the search term in their 
descriptions.  The buttons at the bottom can be used to save information about 
whatever designs are selected.  The "Save selected paths" button will save a 
text file listing the path to the lowest scoring model for each selected 
design.  The "Save selected funnels" button will save a PDF with the plot for 
each selected design on a separate page.

The upper right area of the GUI will contain a plot with different metrics on 
the two axes where each point represents a single model.  You can right-click 
on any point to take an action on the model represented by that point.  Usually 
this means visualizing the model in an external program, like pymol or chimera.  
You can also run custom code by writing a script with the extension *.sho that 
takes the path of a model as its only argument.  ``plot_funnels`` will search 
for scripts with this extension in every directory starting with the directory 
containing the model in question and going down all the way to the root of the 
file system.  Any scripts that are found are added to the menu you get by 
right-clicking on a point, using simple rules (the first letter is capitalized 
and underscores are converted to spaces) to convert the file name into a menu 
item name.

The tool bar below the plot can be used to pan around, zoom in or out, save an 
image of the plot, or change the axes.  If the mouse is over the plot, its 
coordinates will be shown just to the right of these controls.  Below the plot 
is a text form which can be used to enter a description of the design.  These 
descriptions can be searched.  I like using the '+', '++', ... convention to 
rank designs so I can easily search for increasingly good designs.

Hotkeys:
    j,f,down: Select the next design, if there is one.
    k,d,up: Select the previous design, if there is one.
    i,a: Focus on the description form.
    z: Use the mouse to zoom on a rectangle.
    x: Use the mouse to pan (left-click) or zoom (right-click).
    c: Return to the original plot view.
    slash: Focus on the search bar.
    tab: Change the y-axis metric.
    space: Change the x-axis metric.
    escape: Unfocus the search and description forms.
"""

import os, glob, numpy as np, random
from scipy.stats import ks_2samp
import matplotlib.pyplot as plt
import pandas as pd
from roseasy import pipeline
from roseasy import structures
from roseasy import gui


def select_random_subset(pdb_paths, sample_size):
    subset = random.choices(pdb_paths, k=sample_size)

    return subset


def main():
    import docopt
    args = docopt.docopt(__doc__)

    print(args)

    pdb_dirs = args['<pdb_directories>']
    ref_dir = args['<ref_dir>']
    num_samples = int(args['--samples'])
    #sample_size = int(args['--sample_size'])

    #num_samples = 50
    sample_sizes = [5, 10, 20, 30, 50, 100, 150, 200]
    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    ref_records, ref_metadata = structures.load(ref_dir, use_cache=not args['--force'])
    #prop = 'ca_rmsd_no_loop'
    prop = 'ca_rmsd'
    ref_prop = 'ca_rmsd'

    width = 0.35
    x = np.arange(len(sample_sizes))
    prev = np.zeros(len(sample_sizes))
    for pdb_dir in pdb_dirs:
        p_vals = []
        all_records, metadata = structures.load(pdb_dir, use_cache=not args['--force'])
        for sample_size in sample_sizes:
            #print(sample_size)
            p_val = []
            diff_count = 0
            diff_counts = []
            for sample in range(num_samples):
                records = all_records.sample(n=sample_size)
                val = ks_2samp(records[prop],ref_records[ref_prop])[1]
                if val < 0.05:
                    diff_count += 1
                p_val.append(val)
            diff_counts.append(diff_count/num_samples)
            p_vals.append(np.average(p_val))
        ax1.plot(sample_sizes, p_vals,'--o', label = pdb_dir)
        bins = [x + width for x in prev]
        ax2.bar(x , diff_counts, label = pdb_dir)
        prev = bins
    ax1.set_ylabel('K-S p value')
    ax1.set_xlabel('relax sample size')
    ax2.set_ylabel('Rejection Proportion')
    ax2.set_xticks(x)
    ax2.set_xticklabels(sample_sizes)
    ax2.set_xlabel('relax sample size')
    ax1.legend()
    ax2.legend()
    #plt.xlabel('# relax structures sampled')
    #plt.ylabel('K-S p value')
    plt.show()

if __name__=='__main__':
    main()
