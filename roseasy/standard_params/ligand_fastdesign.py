import os, sys, subprocess, gzip
from klab.rosetta import input_files
from roseasy.utils import mover_utils
from roseasy import pipeline
from roseasy import big_jobs
from roseasy.movers import fastdesign

from pyrosetta import init
from pyrosetta import Pose
from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects
from pyrosetta.rosetta.core.pose import setPoseExtraScore
from pyrosetta import pose_from_file
from pyrosetta import pose_from_pdb
from pyrosetta import generate_nonstandard_residue_set
from pyrosetta.rosetta.core.scoring import CA_rmsd
from pyrosetta.rosetta.core.scoring import all_atom_rmsd
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import ReadResfile
from pyrosetta import create_score_function
#from roseasy.standard_params.filters import FilterContainer

def get_workspace(root_dir, step):
    return pipeline.DesignWorkspace(root_dir, step)

if __name__=='__main__':
    workspace, job_info = big_jobs.initiate()
    test_run = job_info.get('test_run', False)
    init()

    # Create residue typeset
    pose = Pose()
    typeset = generate_nonstandard_residue_set(pose,
            workspace.ligand_params_paths)
    # Figure out input pdb and create a pose
    pdbpath = workspace.input_path(job_info)
    pose_from_file(pose, typeset, pdbpath)


    # Create FastDesign object
    fd = fastdesign.FastDesign()
    fd.pose = pose

    # Create task factory and read the resfile
    taskfactory = TaskFactory()
    readresfile = ReadResfile(workspace.resfile_path)
    taskfactory.push_back(readresfile)
    fd.task_factory = taskfactory

    # Parse resfile & create movemap
    resfile_parser = input_files.Resfile(input_resfile=workspace.resfile_path)
    chain = 'A'
    designable = [int(key) for key in resfile_parser.design[chain]]
    if chain in resfile_parser.repack:
        repackable = [int(key) for key in resfile_parser.repack[chain]]
    else:
        repackable = []
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

    filters = workspace.get_filters(fd.pose,
            task_id=job_info['task_id'], score_fragments=False,
            test_run=test_run)
    filters.run_filters()

    # Add RMSDs as extra metrics, which will be printed at the end of
    # the PDB file.
    setPoseExtraScore(fd.pose, 'EXTRA_METRIC_CA_RMSD', ca_rmsd)
    setPoseExtraScore(fd.pose, 'EXTRA_METRIC_AllAtom_RMSD', all_atom_rmsd)

    # Save final pose as a pdb file.
    input_name = os.path.basename(pdbpath).split(".")[0]
    out = workspace.output_prefix(job_info) + input_name + workspace.output_suffix(job_info) + '.pdb.gz'
    pose.dump_pdb(out)
