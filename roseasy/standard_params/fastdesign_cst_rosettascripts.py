from pyrosetta import init, create_score_function
from pyrosetta import rosetta
from pyrosetta.rosetta.core.pose import setPoseExtraScore
from pyrosetta import pose_from_file
from pyrosetta.rosetta.core.scoring import CA_rmsd
from pyrosetta.rosetta.protocols import rosetta_scripts
from pyrosetta.rosetta.core.scoring import all_atom_rmsd
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task import operation
from pyrosetta.rosetta.core.pack.task import residue_selector
from pyrosetta.rosetta.core import select
from pyrosetta.rosetta.protocols.minimization_packing import RotamerTrialsMover
from roseasy import pipeline
from roseasy import big_jobs
from roseasy.utils.mover_utils import setup_movemap_from_resselectors
import os, sys, subprocess, gzip
import json
from roseasy.movers import fastdesign
from roseasy.utils import numeric
#from roseasy.standard_params.filters import FilterContainer


def get_workspace(root_dir, step):
    return pipeline.DesignWorkspace(root_dir, step)

def boollist_to_vector1_bool(a):
    # Should also go in utils eventually
    from pyrosetta.rosetta.utility import vector1_bool
    vector = vector1_bool()
    for item in a:
        vector.append(item)
    return vector

def size_list_to_res_selector(list_of_residues, pose):
    # Should go in utils eventually
    bool_list = []
    pdbinfo = pose.pdb_info()
    for res in range(0, pose.size()):
        bool_list.append(res in list_of_residues)
    return boollist_to_vector1_bool(bool_list)

def strlist_to_vector1_str(strlist):
    # Should also go in utils eventually
    from pyrosetta.rosetta.utility import vector1_std_string
    vector = vector1_std_string()
    for string in strlist:
        vector.append(string)
    return vector

def get_insertion(pdbpath):
    # hardcoding for now
    insertion_path = os.path.join(
            os.environ['HOME'],
            'cas', 'dels', 'rec2_helix',
            'lucs', 'data', 'compatible'
            )
    model_number = os.path.basename(pdbpath).split('.')[0].split('_')[-1]
    insertion_json = os.path.join(
            insertion_path,
            'insertion_points_{}.json'.format(model_number)
            )
    with open(insertion_json, 'r') as f:
        insertion = json.load(f)
    # Return only first insertion for this script because there is only
    # one
    return insertion[0]

def clash_based_taskfactory(pdbpath, pose):
    insertion = get_insertion(pdbpath)
    # How much to add to insertion['stop'] to get the residue number to
    # design
    lucs_relative_resfile = [4, 5, 6, 8, 12, 16, 20, 23, 24, 26, 27, 30]
    # subtract 1 because LUCS loop definitions stop on the residue after the
    # insertion, whereas normal loop definition stops at the end of the
    # loop.
    lucs_relative_resfile = [i + insertion['stop'] - 1 for i in lucs_relative_resfile]
    # This should give us all designed positions
    loop_list = [j for j in range(insertion['start'], insertion['stop'])]
    lucs_relative_resfile += loop_list
    design_mask = size_list_to_res_selector(lucs_relative_resfile, pose)
    design_selector = select.residue_selector.ResidueIndexSelector(
            numeric.intlist_to_vector1_size(lucs_relative_resfile)
            )

    # Takes a resfile and returns a task factory with a clash-based
    # repack shell.
    tf = TaskFactory()
    cl = operation.InitializeFromCommandline()
    notaa = operation.ProhibitSpecifiedBaseResidueTypes(
            strlist_to_vector1_str(['C', 'H']),
            design_selector)
    # read = operation.ReadResfile(resfile)
    tf.push_back(cl)
    tf.push_back(notaa)
    # tf.push_back(read)

    # Select repackable residues from designed residues
    repack_only = operation.RestrictToRepackingRLT()
    repack = operation.OperateOnResidueSubset(repack_only,
            design_selector, False)
    repack.flip_subset(True)
    tf.push_back(repack)

    all_selector = residue_selector.ClashBasedShellSelector(design_mask)
    all_selector.set_num_shells(2)
    all_selector.set_include_focus(True)
    all_selector.invert(True)

    no_packing = operation.PreventRepackingRLT()
    static = operation.OperateOnResidueSubset(no_packing, all_selector,
            False)

    packertask = tf.create_task_and_apply_taskoperations(pose)

    ld = rosetta_scripts.XmlObjects.static_get_task_operation(
        '''<LayerDesign name="layer_all" layer="core_boundary_surface_Nterm_Cterm" use_sidechain_neighbors="True">
            <Nterm>
                    <all append="DEGHKNQRST" />
                    <all exclude="CAFILMPVWY" />
            </Nterm>
            <Cterm>
                    <all append="DEGHKNQRST" />
                    <all exclude="CAFILMPVWY" />
            </Cterm>
    </LayerDesign>''')
    tf.push_back(static)
    tf.push_back(ld)
    packertask = tf.create_task_and_apply_taskoperations(pose)
    print('PRINTING PACKERTASK')
    print(packertask)

    repack_mask = packertask.repacking_residues()
    design_mask = packertask.designing_residues()
    movemap = setup_movemap_from_resselectors(design_mask, repack_mask)

    return tf, movemap



if __name__=='__main__':
    test=False
    # test = True
    if not test:
        workspace, job_info = big_jobs.initiate()
        test_run = job_info.get('test_run', False)
        pdbpath = workspace.input_path(job_info)
    else:
        workspace = pipeline.workspace_from_dir(sys.argv[1])
        pdbpath = '02_designs/inputs/model_0.pdb.gz'
        test_run = True
        task_id = 1
    init('-total_threads 1 -packing:ex1 -packing:ex2 -packing:ex1aro '\
            '-use_input_sc')
    pose = pose_from_file(pdbpath)
    task_factory, movemap = clash_based_taskfactory(pdbpath,
            pose)
    ref = create_score_function('ref2015')
    rot = RotamerTrialsMover(ref, task_factory)
    print('APPLYING ROTAMERTRIALS')
    rot.apply(pose)

    # cst_gen_str = '''
    # <CoordinateConstraintGenerator name="backbone"
    # ca_only="true"
    # '''
    # backbone_cst = rosetta_scripts.XmlObjects.static_get_mover(cst_gen_str)
    # csts = rosetta.protocols.constraint_generator.CoordinateConstraintGenerator()
    # csts.set_ca_only(True)

    # fd = fastdesign.FastDesign()
    xml = rosetta_scripts.XmlObjects.create_from_file(os.path.join(
        workspace.find_path('fastdesign_cst.xml')
        ))
    protocol = xml.get_mover('ParsedProtocol')

    init_args = []
    dalphaball_path = os.path.join(workspace.rosetta_dir, 'source',
            'external', 'DAlpahBall', 'DAlphaBall.gcc')
    init_args.append('-holes:dalphaball {} -in:file:s {}'.format(dalphaball_path, pdbpath))
    init_args.append('-total_threads 1')
    init_args.append('-packing:ex1 -packing:ex2')
    init_args.append('-packing:ex1aro')
    init_args.append('-use_input_sc')
    init_args.append('-relax:constrain_relax_to_start_coords')
    if test_run:
        protocol.get_mover(2).rounds = 1
        # fd.rounds = 1
    # fd.pose = pose
    # Warning: This is an all-atom movemap. Constrain to input coords if
    # you don't want things to move around a lot.
    # loop = workspace.largest_loop
    protocol.get_mover(2).set_task_factory(task_factory)
    protocol.get_mover(2).set_movemap(movemap)
    # print('PRINTING MOVEMAP AND TASK FACTORY')
    # print(fd.movemap)
    # print(fd.task_factory)
    # fd.mover.add_constraint_generator(csts)
    # fd.mover.constrain_relax_to_start_coords(True)
    # fd.mover.ramp_down_constraints(False)
    # Before we apply FastDesign, also setup and run RotamerTrials
    # fd.apply()
    protocol.apply(pose)

    # This will compare it to the input to the step
    input_pose = pose_from_file(pdbpath)
    # But you can uncomment this to compare it to the input to the
    # project
    #input_pose = pose_from_file(workspace.input_pdb_path)
    ca_rmsd = CA_rmsd(pose, input_pose)
    all_atom_rmsd = all_atom_rmsd(pose, input_pose)
    # score_fragments = os.path.exists(workspace.loops_path)
    score_fragments=False

    filters = workspace.get_filters(pose,
            task_id=job_info['task_id'], score_fragments=score_fragments,
            test_run=test_run)
    filters.run_filters()

    input_name = os.path.basename(pdbpath).split(".")[0]
    out = workspace.output_prefix(job_info) + input_name + workspace.output_suffix(job_info) + '.pdb.gz'

    setPoseExtraScore(pose, 'EXTRA_METRIC_CA_RMSD', ca_rmsd)
    setPoseExtraScore(pose, 'EXTRA_METRIC_AllAtom_RMSD', all_atom_rmsd)

    pose.dump_pdb(out)
