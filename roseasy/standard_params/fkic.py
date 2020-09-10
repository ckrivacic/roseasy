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
    '''
    Provides information to the submission script about the workspace
    '''
    return pipeline.ValidationWorkspace(root_dir, step)


if __name__=="__main__":
    from roseasy.movers.loopmodeler import LoopModeler
    '''Get information about the task to be run'''
    workspace, job_info = big_jobs.initiate()
    test_run = job_info.get('test_run', False)
    pdbpath = workspace.input_path(job_info)

    init()
    pose = pose_from_file(pdbpath)
    # RosEasy loopmodeler object
    lm = LoopModeler()
    # Fragment flags taken care of by workspace
    lm.fragments_flags = workspace.fragments_flags(pdbpath)
    # Config loop modeler for fkic
    lm.config = 'fkic'
    # Get loops from loops file; can also get loops from a range of
    # residues
    lm.loops_from_file(workspace.loops_path)
    # Argument needed for FragmentScoreFilter
    lm.add_init_arg('-in:file:s {}'.format(pdbpath))
    if test_run:
        # What happens when --test-run is supplied?
        lm.mark_as_test_run()

    # Add pose and apply
    lm.pose = pose
    lm.apply()

    # Get and run default filters
    filters = FilterContainer(workspace, lm.pose,
            task_id=job_info['task_id'], score_fragments=True,
            test_run=test_run)
    filters.run_filters()

    # Calculate RMSD from input
    input_pose = pose_from_file(workspace.input_path(job_info))
    ca_rmsd = CA_rmsd(lm.pose, input_pose)
    all_atom_rmsd = all_atom_rmsd(lm.pose, input_pose)

    # Add RMSDs to pose for parsing
    setPoseExtraScore(lm.pose, 'EXTRA_METRIC_CA_RMSD', ca_rmsd)
    setPoseExtraScore(lm.pose, 'EXTRA_METRIC_AllAtom_RMSD', all_atom_rmsd)

    if not os.path.exists(workspace.output_prefix(job_info)):
        os.mkdir(workspace.output_prefix(job_info))
    pose.dump_pdb(workspace.output_path(job_info))
