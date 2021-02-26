from pyrosetta import init
from pyrosetta import pose_from_file
from pyrosetta.rosetta.core.scoring import CA_rmsd
from pyrosetta.rosetta.core.scoring import all_atom_rmsd
from pyrosetta.rosetta.core.pose import setPoseExtraScore
from roseasy import pipeline
from roseasy import big_jobs
from roseasy.standard_params.filters import FilterContainer
import os, sys, subprocess


def get_workspace(root_dir, step):
    return pipeline.ValidationWorkspace(root_dir, step)


if __name__=="__main__":
    from roseasy.movers.loopmodeler import LoopModeler
    # Get basic job information
    workspace, job_info = big_jobs.initiate()
    test_run = job_info.get('test_run', False)
    init()
    pdbpath = workspace.input_path(job_info)
    print('LOADING PDB FROM {}'.format(pdbpath))

    # Setup the pose
    pose = pose_from_file(pdbpath)
    pdbinfo = pose.pdb_info()
    input_name = os.path.basename(workspace.input_path(job_info)).split('.')[0]
    # Setup the mover
    lm = LoopModeler()
    # Start and end of loop
    start = workspace.largest_loop.start
    end = workspace.largest_loop.end

    # Don't mutate the first two and last two residues of loop
    nomutate = ''
    nomutate += str(start)
    nomutate += ','
    nomutate += str(int(start) + 1)
    nomutate += ','
    nomutate += str(int(end) - 1)
    nomutate += ','
    nomutate += str(end)

    lm.seqpose_no_mutate=nomutate

    # Setup loop
    # lm.loop_from_range(start, end)
    lm.loops_from_file(workspace.loops_path)

    # Loophash database path (can leave as is if on Wynton)
    lh_db_path = '/wynton/home/kortemme/krivacic/rosetta/database/loophash_db'
    # Loophash needs additional initiation arguments
    init('-lh:db_path {} -lh:loopsizes 6 8 9'.format(lh_db_path))
    lm.add_init_arg('-lh:db_path {}'.format(lh_db_path))
    lm.add_init_arg('-lh:loopsizes 6 8 9')
    # Configure roseasy mover now that we have init args
    lm.config = 'lhk'
    # Dalphaball needed for BUNS calculations
    dalphaball_path = os.path.join(workspace.rosetta_dir, 'source',
            'external', 'DAlphaBall', 'DAlphaBall.gcc')
    lm.add_init_arg('-holes:dalphaball {}'.format(dalphaball_path))
    lm.add_init_arg('-in:file:s {}'.format(pdbpath))

    if test_run:
        lm.mover.centroid_stage().mark_as_test_run()
        lm.mover.fullatom_stage().mark_as_test_run()

    # Apply mover to pose
    lm.pose = pose
    lm.apply()

    # Re-load the input pose to compare loop RMSD
    input_pose = pose_from_file(workspace.input_path(job_info))
    # These are all-atom RMSDs. Can comment out if you don't care about
    # that.
    ca_rmsd = CA_rmsd(lm.pose, input_pose)
    all_atom_rmsd = all_atom_rmsd(lm.pose, input_pose)
    # This is the loop  RMSD. Usually a useful metric for loop modeling
    # simulations.
    loop = workspace.largest_loop
    rmsd = CA_rmsd(input_pose, lm.pose, loop.start,
            loop.end)

    out = workspace.output_prefix(job_info) + input_name + workspace.output_suffix(job_info) + '.pdb.gz'

    filters = FilterContainer(workspace, lm.pose,
            task_id=job_info['task_id'], score_fragments=True,
            test_run=test_run)
    filters.run_filters()

    setPoseExtraScore(lm.pose, 'EXTRA_METRIC_CA_RMSD', ca_rmsd)
    setPoseExtraScore(lm.pose, 'EXTRA_METRIC_AllAtom_RMSD', all_atom_rmsd)
    setPoseExtraScore(lm.pose, 'EXTRA_METRIC_Loop_RMSD', rmsd)

    if not os.path.exists(workspace.output_prefix(job_info)):
        os.mkdir(workspace.output_prefix(job_info))
    pose.dump_pdb(out)
