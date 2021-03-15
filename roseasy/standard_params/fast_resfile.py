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

def get_workspace(root_dir, step):
    return pipeline.DesignWorkspace(root_dir, step)

if __name__=='__main__':
    local=False
    if local:
        workspace = pipeline.workspace_from_dir(sys.argv[1])
        job_info = {}
        job_info['task_id'] = 1
        job_info['inputs'] = ['0000.pdb.gz']
        test_run = True
    else:
        workspace, job_info = big_jobs.initiate()
        test_run = job_info.get('test_run', False)
    init('-total_threads 1')
    pdbpath = workspace.input_path(job_info)
    pose = pose_from_file(pdbpath)
    fd = fastdesign.FastDesign()

    dalphaball_path = os.path.join(workspace.rosetta_dir, 'source',
            'external', 'DAlpahBall', 'DAlphaBall.gcc')
    if os.path.exists(dalphaball_path):
        fd.add_init_arg('-holes:dalphaball {} -in:file:s {}'.format(dalphaball_path, pdbpath))
    fd.add_init_arg('-total_threads 1')

    fd.pose = pose
    taskfactory = TaskFactory()
    readresfile = ReadResfile(workspace.resfile_path)
    taskfactory.push_back(readresfile)
    fd.task_factory = taskfactory
    resfile_parser = input_files.Resfile(input_resfile=workspace.resfile_path)
    designable = [int(key) for key in resfile_parser.design['A']]
    repackable = [int(key) for key in resfile_parser.repack['A']]
    # fd.movemap = mover_utils.setup_movemap_from_resfile(designable,
            # repackable, pdbinfo=fd.pose.pdb_info(), chain='A')
    # fd.movemap = mover_utils.setup_movemap(designable,
            # repackable)
    fd.setup_default_movemap(bb=designable.extend(repackable),
            chi=designable.extend(repackable))
    favornative = '''
    <FavorNativeResidue name="favornative" bonus="1.0"/>

    '''
    favornative = XmlObjects.static_get_mover(favornative)
    favornative.apply(fd.pose)
    if test_run:
        fd.rounds = 1
    # Warning: This is an all-atom movemap. Constrain to input coords if
    # you don't want things to move around a lot.
    # loop = workspace.largest_loop
    # fd.task_factory_from_range(loop.start, loop.end)
    # Get a task factory from the resfile
    print(fd.movemap)
    print(fd.task_factory)
    fd.apply()

    # This will compare it to the input to the step
    input_pose = pose_from_file(pdbpath)
    # But you can uncomment this to compare it to the input to the
    # project
    #input_pose = pose_from_file(workspace.input_pdb_path)
    ca_rmsd = CA_rmsd(fd.pose, input_pose)
    all_atom_rmsd = all_atom_rmsd(fd.pose, input_pose)
    score_fragments = os.path.exists(workspace.loops_path)
    sfxn = create_score_function('ref2015')

    filters = workspace.get_filters(fd.pose,
            task_id=job_info['task_id'], score_fragments=score_fragments,
            test_run=test_run)
    filters.run_filters()

    input_name = os.path.basename(pdbpath).split(".")[0]
    out = workspace.output_prefix(job_info) + input_name + workspace.output_suffix(job_info) + '.pdb.gz'
    # Re-score, just in case the native mover impacts the final score
    # (it shouldn't though).
    final_score = sfxn(pose)
    rmsd = CA_rmsd(input_pose, fd.pose, 997,
            1006)

    setPoseExtraScore(fd.pose, 'EXTRA_METRIC_CA_RMSD', ca_rmsd)
    setPoseExtraScore(fd.pose, 'EXTRA_METRIC_AllAtom_RMSD', all_atom_rmsd)
    setPoseExtraScore(fd.pose, 'EXTRA_METRIC_Final_Score_(REU)',
            final_score)
    setPoseExtraScore(fd.pose, 'EXTRA_METRIC_Loop_RMSD', rmsd)

    pose.dump_pdb(out)
