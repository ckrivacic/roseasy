from pyrosetta import init
from pyrosetta import pose_from_file
from pyrosetta.rosetta.core.scoring import CA_rmsd
from pyrosetta.rosetta.core.scoring import all_atom_rmsd
from roseasy import pipeline
from roseasy import big_jobs
import os, sys, subprocess


def get_workspace(root_dir, step):
    return pipeline.ValidationWorkspace(root_dir, step)


if __name__=="__main__":
    from roseasy.movers.loopmodeler import LoopModeler
    workspace, job_info = big_jobs.initiate()
    test_run = job_info.get('test_run', False)
    init()
    pdbpath = workspace.input_path(job_info)
    pose = pose_from_file(pdbpath)
    lm = LoopModeler()
    lm.config = 'fkic'
    lm.fragments_flags = workspace.fragments_flags(pdbpath)
    lm.loops_from_file(workspace.loops_path)
    if test_run:
        lm.mover.centroid_stage().mark_as_test_run()
        lm.mover.fullatom_stage().mark_as_test_run()
    lm.pose = pose
    lm.apply()

    input_pose = pose_from_file(workspace.input_path(job_info))
    ca_rmsd = CA_rmsd(lm.pose, input_pose)
    all_atom_rmsd = all_atom_rmsd(lm.pose, input_pose)

    input_name = os.path.basename(workspace.input_path(job_info)).split('.')[0]
    out = workspace.output_prefix(job_info) + input_name + workspace.output_suffix(job_info) + '.pdb'
    if not os.path.exists(workspace.output_prefix(job_info)):
        os.mkdir(workspace.output_prefix(job_info))
    pose.dump_pdb(out)
    with open(out, 'a') as f:
        f.write('\nCA_RMSD {}'.format(ca_rmsd))
        f.write('\nAllAtom_RMSD {}'.format(all_atom_rmsd))
