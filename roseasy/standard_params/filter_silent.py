'''
USAGE:
    This script runs filters on the output of silent files.
    Silent files should be part of a ValidationWorkspace. There is some
    prep work before you can run this script:
    First, you must create a new ValidationWorkspace (roseasy
    submit . project_params/filter_silent.py --make-dirs).
    Next, symlink all of the silent files to the inputs directory 
    such that the link is named after the FOLDER name of the silent file and 
    has a PDB extension.
    E.g., 02_validated_designs/outputs/0000/silent.out would have a link
    at 03_validated_designs/inputs/0000.pdb.

    There is a link_silent.py script that can be used for this; simply
    python link_silent.py <previous_focusdir> <current_focusdir>.

    Finally, you can run this script with RosEasy:
    roseasy submit <workspace> <script_path> [options]
    '''
import os, sys, subprocess, gzip
from klab.rosetta import input_files
from roseasy.utils import mover_utils
from pyrosetta import init
from pyrosetta.rosetta.core.pose import setPoseExtraScore
from pyrosetta import pose_from_file
from pyrosetta import poses_from_silent
from pyrosetta.rosetta.core.scoring import CA_rmsd
from pyrosetta.rosetta.core.scoring import all_atom_rmsd
from roseasy import pipeline
from roseasy import big_jobs
from pyrosetta import create_score_function
from roseasy.movers import fastdesign
#from roseasy.standard_params.filters import FilterContainer

def get_workspace(root_dir, step):
    return pipeline.ValidationWorkspace(root_dir, step)

if __name__=='__main__':
    workspace, job_info = big_jobs.initiate()
    test_run = job_info.get('test_run', False)

    # Figure out input pdb and create a pose
    silentpath = workspace.input_path(job_info)

    # Input pdb path for RMSD comparisons
    predecessor = pipeline.workspace_from_dir(workspace.predecessor)
    pdbpath = os.path.join(
            workspace.root_dir,
            predecessor.input_dir,
            os.path.basename(silentpath)
            )
    print('CURRENT INPUT: {}'.format(pdbpath))
    dalphaball_path = os.path.join(workspace.rosetta_dir, 'source',
            'external', 'DAlpahBall', 'DAlphaBall.gcc')
    init('-total_threads 1 -holes:dalphaball {} -in:file:s {}'.format(dalphaball_path, pdbpath))
    poses = poses_from_silent(silentpath)
    i = 0
    tasknum = int(job_info['task_id'])
    posedict = {}
    for pose in poses:
        if i%10 == (tasknum - 1)// len(job_info['inputs']):
            posedict[i] = pose
        i += 1
    sfxn = create_score_function('ref2015')

    # Create new pose from input file for comparison
    input_pose = pose_from_file(pdbpath)

    for i in posedict:
        pose = posedict[i]
        sfxn(input_pose)
        sfxn(pose)

        ca_rmsd = CA_rmsd(pose, input_pose)
        aa_rmsd = all_atom_rmsd(pose, input_pose)

        # Add RMSDs to pose for parsing
        setPoseExtraScore(pose, 'EXTRA_METRIC_CA_RMSD', ca_rmsd)
        setPoseExtraScore(pose, 'EXTRA_METRIC_AllAtom_RMSD', aa_rmsd)
        #input_pose = pose_from_file(workspace.input_pdb_path)
        # Calculate several different types of RMSD

        filters = workspace.get_filters(pose,
                task_id=str(job_info['task_id']) + '_' + str(i), 
                score_fragments=False,
                test_run=test_run,
                fragment_full_chain=1)
        filters.run_filters()

        # Save final pose as a pdb file.
        input_name = os.path.basename(pdbpath).split('.')[0]
        out = workspace.output_prefix(job_info) + \
                input_name + '_{}'.format(i) + '.pdb.gz'
        print('Dumping PDB to {}'.format(out))
        os.makedirs(workspace.output_prefix(job_info), exist_ok=True)
        pose.dump_pdb(out)
