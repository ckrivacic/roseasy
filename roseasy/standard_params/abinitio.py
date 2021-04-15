"""\
Carry out ab initio folding using the contents of an output folder from 
roseasy generate_fragments (namely the fasta, 3mers.gz, 9mers.gz, and 
psipred_ss2 files therein).

Usage:
    roseasy submit <workspace> <path_to_this_script> [options]

Options:
    -m, --mem-free=MEM  [default: 2]
        The amount of memory (GB) to request from the cluster.  Bigger systems
        may need more memory, but making large memory requests can make jobs
        take much longer to come off the queue (since there may only be a few
        nodes with enough memory to meet the request).

    -d, --dry-run
        Print out the command-line that would be used to generate fragments,
        but don't actually run it.

    -i, --input_model=MODEL
        Carry out ab initio relaxation for a model other than the workspace 
        input modelg

    -n, --nstruct=NSTRUCT
        The number of ab initio folded structures to generate in total.
"""

from pyrosetta.rosetta.core.scoring import CA_rmsd
from pyrosetta.rosetta.core.scoring import all_atom_rmsd
from pyrosetta import pose_from_file
from pyrosetta import init
from pyrosetta.rosetta.core.pose import setPoseExtraScore
import os
import sys
import glob
import subprocess
from klab import scripting, cluster
from roseasy import pipeline, big_jobs
import time

def get_workspace(root_dir, step):
    return pipeline.ValidationWorkspace(root_dir, step)

def main():
    cluster.require_qsub()

    start_time = time.time()

    workspace, job_info = big_jobs.initiate()
    pdbpath = workspace.input_path(job_info)
    if not os.path.exists(workspace.output_prefix(job_info)):
        os.mkdir(workspace.output_prefix(job_info))
    outpath = workspace.output_path(job_info)
    test_run = job_info.get('test_run', False)

    # append location of Rosetta binaries to path
    abinitio = '/wynton/home/kortemme/krivacic/rosetta/source/bin/AbinitioRelax.linuxgccrelease'

    # find necessary files
    print('fdir :', workspace.fragments_dir)
    print('pdb :', pdbpath)
    print(os.path.join(workspace.fragments_dir, 
                       workspace.fragments_tag(pdbpath.lower())+'?', 
                       '*3mers.gz'))
    tmers = glob.glob(
        os.path.join(workspace.fragments_dir, 
                     workspace.fragments_tag(pdbpath.lower())+'?', 
                     '*3mers.gz'))[0]
    nmers = glob.glob(
        os.path.join(workspace.fragments_dir, 
                     workspace.fragments_tag(pdbpath.lower())+'?', 
                     '*9mers.gz'))[0]
    ss2 = glob.glob(
        os.path.join(workspace.fragments_dir, 
                     workspace.fragments_tag(pdbpath.lower())+'?', 
                     '*psipred_ss2'))[0]
    fasta = glob.glob(
        os.path.join(workspace.fragments_dir, 
                     workspace.fragments_tag(pdbpath.lower())+'?', 
                     '*fasta'))[0]

        # Run the ab initio relax script.

    relax_abinitio = [
            abinitio,
            '-abinitio:relax', 
            '-use_filters', 'true', 
            '-abinitio::increase_cycles', '10', 
            '-abinitio::rg_reweight', '0.5', 
            '-abinitio::rsd_wt_helix', '0.5', 
            '-abinitio::rsd_wt_loop', '0.5', 
            '-relax::fast', 
            '-in:file:fasta', fasta, 
            '-in:file:frag3', tmers, 
            '-in:file:frag9', nmers, 
            '-in:file:psipred_ss2', ss2, 
            '-nstruct', '10', 
            # '-out:pdb_gz',
            '-out:file:silent', workspace.output_prefix(job_info) + 'silent_{}.out'.format(workspace.output_suffix(job_info)),
            # '-out:prefix', workspace.output_prefix(job_info),
            # '-out:suffix', workspace.output_suffix(job_info),
            '-out:no_nstruct_label'
            # '--outdir', workspace.fragments_dir,
            # '--memfree', args['--mem-free']
            ]
    
    print('Running ROSETTA command:')
    print(' '.join(relax_abinitio))

    subprocess.call(relax_abinitio)


if __name__ == "__main__":
    main()
