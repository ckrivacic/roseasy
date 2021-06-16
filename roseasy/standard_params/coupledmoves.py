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

    # Figure out input pdb and create a pose
    pdbpath = workspace.input_path(job_info)
    temp_outdir = os.path.join(
            workspace.focus_dir, 'cm_outputs'
            )
    os.makedirs(temp_outdir, exist_ok=True)

    os.chdir(temp_outdir)

    # Set up init args
    init_args = []
    init_args.append('-ex1 -ex2 -use_input_sc -ex1aro -extrachi_cutoff 0')
    init_args.append('-s {}'.format(pdbpath))
    print('PARAM PATHS:')
    print(workspace.ligand_params_paths)
    for path in workspace.ligand_params_paths:
        init_args.append('-extra_res_fa {}'.format(path))
        print('Appended extra res')
    if len(workspace.ligand_params_paths) > 0:
        init_args.append('-coupled_moves:ligand_mode true')
    init_args.append('-coupled_moves:output_prefix {}'.format(
        os.path.basename(pdbpath).split('.')[0] +
        workspace.output_suffix(job_info)))


    dalphaball_path = os.path.join(workspace.rosetta_dir, 'source',
            'external', 'DAlpahBall', 'DAlphaBall.gcc')
    init_args.append('-holes:dalphaball {}'.format(dalphaball_path))
    init_args.append('-total_threads 1')
    init_args.append('-resfile {}'.format(workspace.resfile_path))

    init(' '.join(init_args))
    pose = pose_from_file(pdbpath)

    # Create task factory and read the resfile
    # taskfactory = TaskFactory()
    # readresfile = ReadResfile(workspace.resfile_path)
    # taskfactory.push_back(readresfile)

    # Set up coupled moves
    from pyrosetta.rosetta.protocols.coupled_moves import CoupledMovesProtocol
    coupledmoves = CoupledMovesProtocol()

    coupledmoves.apply(pose)

    output_prefix = os.path.join(
            temp_outdir,
            'EMPTY_JOB_use_jd2_0000' + 
            os.path.basename(pdbpath).split('.')[0] + 
            workspace.output_suffix(job_info)
            )
    poselast = pose_from_file(output_prefix + '_last.pdb')
    poselow = pose_from_file(output_prefix + '_low.pdb')
    sfxn = create_score_function('ref2015')

    for tup in [(poselast, '_last'), (poselow, '_low')]:
        # Create new pose from input file for comparison
        pose = tup[0]
        sfxn(pose)
        suffix = tup[1]
        input_pose = pose_from_file(pdbpath)
        # Calculate several different types of RMSD
        ca_rmsd = CA_rmsd(pose, input_pose)
        aa_rmsd = all_atom_rmsd(pose, input_pose)

        filters = workspace.get_filters(pose,
                task_id=job_info['task_id'], score_fragments=False,
                test_run=test_run)
        results = filters.run_filters()

        # Add RMSDs as extra metrics, which will be printed at the end of
        # the PDB file.
        setPoseExtraScore(pose, 'EXTRA_METRIC_CA_RMSD', ca_rmsd)
        setPoseExtraScore(pose, 'EXTRA_METRIC_AllAtom_RMSD', aa_rmsd)

        total_time = time.time() - start_time
        setPoseExtraScore(pose, 'EXTRA_METRIC_Run time', total_time)

        # Save final pose as a pdb file.
        input_name = os.path.basename(pdbpath).split(".")[0]
        out = workspace.output_prefix(job_info) + input_name + \
                workspace.output_suffix(job_info) + suffix + '.pdb.gz'
        pose.dump_pdb(out)
