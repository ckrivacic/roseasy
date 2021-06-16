import os, sys, subprocess, gzip
from klab.rosetta import input_files
from roseasy.utils import mover_utils
from pyrosetta import init
from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects
from pyrosetta.rosetta.core.pose import setPoseExtraScore
from pyrosetta import pose_from_file
from pyrosetta.rosetta.core.scoring import CA_rmsd
from pyrosetta.rosetta.core.scoring import all_atom_rmsd
from roseasy import pipeline
from roseasy import big_jobs
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import ReadResfile
from pyrosetta import create_score_function
from roseasy.movers import fastdesign
#from roseasy.standard_params.filters import FilterContainer
import time

def get_workspace(root_dir, step):
    return pipeline.DesignWorkspace(root_dir, step)

if __name__=='__main__':
    start_time = time.time()
    workspace, job_info = big_jobs.initiate()
    test_run = job_info.get('test_run', False)
    init()

    # Figure out input pdb and create a pose
    pdbpath = workspace.input_path(job_info)
    pose = pose_from_file(pdbpath)

    # Create FastDesign object
    fd = fastdesign.FastDesign()
    fd.pose = pose
    dalphaball_path = os.path.join(workspace.rosetta_dir, 'source',
     'external', 'DAlpahBall', 'DAlphaBall.gcc')
    fd.add_init_arg('-holes:dalphaball {} -in:file:s {}'.format(dalphaball_path, pdbpath))
    fd.add_init_arg('-ex1 -ex2 -use_input_sc -ex1aro')
    fd.add_init_arg('-total_threads 1')

    # Create task factory and read the resfile
    taskfactory = TaskFactory()
    readresfile = ReadResfile(workspace.resfile_path)
    taskfactory.push_back(readresfile)
    fd.task_factory = taskfactory

    # Parse resfile & create movemap
    resfile_parser = input_files.Resfile(input_resfile=workspace.resfile_path)
    # chain = 'A'
    designable = []
    repackable = []
    for chain in resfile_parser.design:
        designable.extend([pose.pdb_info().pdb2pose(chain, int(key)) for key in
            resfile_parser.design[chain]])
    for chain in resfile_parser.repack:
        repackable.extend([pose.pdb_info().pdb2pose(chain, int(key)) for key in
            resfile_parser.repack[chain]])
    fd.setup_default_movemap(bb=designable.extend(repackable),
            chi=designable.extend(repackable))

    if test_run:
        fd.rounds = 1

    print(fd.movemap)
    print(fd.task_factory)
    fd.apply()

    # Create new pose from input file for comparison
    input_pose = pose_from_file(pdbpath)
    #input_pose = pose_from_file(workspace.input_pdb_path)
    # Calculate several different types of RMSD
    ca_rmsd = CA_rmsd(fd.pose, input_pose)
    all_atom_rmsd = all_atom_rmsd(fd.pose, input_pose)
    score_fragments = os.path.exists(workspace.loops_path)

    filters = workspace.get_filters(fd.pose,
            task_id=job_info['task_id'], score_fragments=score_fragments,
            test_run=test_run)
    filters.run_filters()

    # Add RMSDs as extra metrics, which will be printed at the end of
    # the PDB file.
    setPoseExtraScore(fd.pose, 'EXTRA_METRIC_CA_RMSD', ca_rmsd)
    setPoseExtraScore(fd.pose, 'EXTRA_METRIC_AllAtom_RMSD', all_atom_rmsd)

    total_time = time.time() - start_time
    setPoseExtraScore(fd.pose, 'EXTRA_METRIC_Run time', total_time)

    # Save final pose as a pdb file.
    input_name = os.path.basename(pdbpath).split(".")[0]
    out = workspace.output_prefix(job_info) + input_name + workspace.output_suffix(job_info) + '.pdb.gz'
    pose.dump_pdb(out)
