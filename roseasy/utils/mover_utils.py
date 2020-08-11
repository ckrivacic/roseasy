from pyrosetta import *
from roseasy.utils.utils import list_to_str
import random
"""
Functions to help set up movers.
"""


def setup_movemap(residues_bb_movable, residues_sc_movable):
    mm = rosetta.core.kinematics.MoveMap()
    mm.set_chi(False)
    mm.set_bb(False)

    for i in residues_bb_movable:
        mm.set_bb(i, True)
        mm.set_chi(i, True)

    for i in residues_sc_movable:
        mm.set_chi(i, True)

    return mm


def setup_movemap_from_resselectors(designable_selector, repackable_selector):
    """
    Temporary function. Ultimately we want a more sophisticated movemap,
    probably a neighbor movemap or something using clash-based shell selector.
    """
    mm = rosetta.core.kinematics.MoveMap()
    mm.set_chi(False)
    mm.set_bb(False)

    for i in range(1, len(repackable_selector) + 1):
        if designable_selector[i] or repackable_selector[i]:
            mm.set_bb(i, True)
            mm.set_chi(i, True)

    # for i in residues_sc_movable:
    #    mm.set_chi(i, True)

    return mm


def setup_restrained_sfxn(restraint_types, weights):
    sfxn = create_score_function("ref2015_cst")
    score_manager = rosetta.core.scoring.ScoreTypeManager()
    for i in range(0, len(restraint_types)):
        score_term = score_manager.score_type_from_name(restraint_types[i])
        sfxn.set_weight(score_term, weights[i])

    return sfxn


def setup_task_factory(pose, designable_residue_selector,
        repackable_residue_selector,
        motif_dict=None,
        extra_rotamers=True, limit_aro_chi2=True, layered_design=True,
        designable_aa_types=None, prepare_focus=False):
    """
    Adapted from XingJie Pan's code at
    git@github.com:Kortemme-Lab/local_protein_sequence_design.git:
    local_protein_sequence_design/basic.py

    motif_dict should have the resnum as the key for the single-letter restype,
    ex. {38:'E'}

    If you want to set up a task factory to only design the 'focus' residue,
    pass it a motif dict and prepare_focus=True
    """

    def list_to_str(l):
        return ','.join(list(str(i) for i in l))

    task_factory = rosetta.core.pack.task.TaskFactory()

    if len(designable_residue_selector) > 0 and not prepare_focus:
        racaa = rosetta.core.pack.task.operation.RestrictAbsentCanonicalAASRLT()

        if designable_aa_types is None or\
                len(list(compress(xrange(len(designable_residue_selector)),\
                    designable_residue_selector))) != len(designable_aa_types):
            racaa.aas_to_keep('GAPVILMFYWSTKRDENQ')  # No Cys or His
        else:
            racaa.aas_to_keep(designable_aa_types[i])

        designable_operation = rosetta.core.pack.task.operation.OperateOnResidueSubset(
                racaa, designable_residue_selector)
        task_factory.push_back(designable_operation)

    if motif_dict:
        for resnum in motif_dict:
            selector = rosetta.core.select.residue_selector.ResidueIndexSelector(str(resnum))
            print(selector)
            racaa = rosetta.core.pack.task.operation.RestrictAbsentCanonicalAASRLT()
            racaa.aas_to_keep(motif_dict[resnum])

            motif_operation = rosetta.core.pack.task.operation.OperateOnResidueSubset(
                    racaa, selector
                    )
            task_factory.push_back(motif_operation)

    if len(repackable_residue_selector) > 0: # This is always going to be true
        # given a res selector; should alter this logic at some point
        repackable_operation = rosetta.core.pack.task.operation.OperateOnResidueSubset(
                rosetta.core.pack.task.operation.RestrictToRepackingRLT(),
                repackable_residue_selector)
        task_factory.push_back(repackable_operation)
    if prepare_focus and len(designable_residue_selector) > 0:
        repackable_operation = rosetta.core.pack.task.operation.OperateOnResidueSubset(
                rosetta.core.pack.task.operation.RestrictToRepackingRLT(),
                designable_residue_selector
                )
        task_factory.push_back(repackable_operation)

    natro_residues = [i for i in range(1, pose.size() + 1) if (not
            (designable_residue_selector[i] or repackable_residue_selector[i])
            and i not in motif_dict)]
    if len(natro_residues) > 0:
        natro_selector =\
            rosetta.core.select.residue_selector.ResidueIndexSelector(list_to_str(natro_residues))
        natro_operation = rosetta.core.pack.task.operation.OperateOnResidueSubset(
                rosetta.core.pack.task.operation.PreventRepackingRLT(),
                natro_selector)
        task_factory.push_back(natro_operation)

    if extra_rotamers:
        ers = rosetta.core.pack.task.operation.ExtraRotamersGeneric()
        ers.ex1(True)
        ers.ex2(True)
        ers.extrachi_cutoff(18)
        task_factory.push_back(ers)

    if limit_aro_chi2:
        lac = rosetta.protocols.task_operations.LimitAromaChi2Operation()
        task_factory.push_back(lac)

    if layered_design:
        ld = rosetta.protocols.rosetta_scripts.XmlObjects.static_get_task_operation(
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
        task_factory.push_back(ld)

    return task_factory

    #task_design.restrict_to_residues(residue_selector_output)


def generate_loop_simple(pose, focus_residue, resbuffer=3):
    '''Function to get a loops object for LoopModeler. For now, simply
    do focus residue +/- 3 residues, but we can make this more
    sophisticated later (the constructor for Loops can also take a
    residue selector).'''
    loop = rosetta.protocols.loops.Loop(max(focus_residue - resbuffer, 1),
            min(focus_residue + resbuffer, pose.size()))
    loop.set_cut(focus_residue-1)
    return loop


def generate_loops_simple(pose, focus_residue, resbuffer=3):
    '''Function to get a loops object for LoopModeler. For now, simply
    do focus residue +/- 3 residues, but we can make this more
    sophisticated later (the constructor for Loops can also take a
    residue selector).'''
    loop = rosetta.protocols.loops.Loop(max(focus_residue - resbuffer, 1),
            min(focus_residue + resbuffer, pose.size()))
    loop.set_cut(focus_residue)
    loops = rosetta.protocols.loops.Loops()
    loops.add_loop(loop)
    print(loops)
    return loops


def generate_loops_from_range(start, end, cut=None):
    '''Create a loops object from a start and finish position.'''
    loop = rosetta.protocols.loops.Loop(start, end)
    if cut:
        loop.set_cut(cut)
    else:
        loop.set_cut(start + ((end - start)//2))
    loops = rosetta.protocol.loops.Loops()
    loops.add_loop(loop)
    return loops


def generate_loop_from_range(start, end, cut=None, skip=0.0,
        extended=False):
    '''Create a loops object from a start and finish position.'''
    try:
        # Can only set skip rate in constructor.
        loop = rosetta.protocols.loops.Loop(start, end, cut, skip, extended)
    except:
        loop = rosetta.protocols.loops.Loop(start, end)
        if cut:
            loop.set_cut(cut)
        else:
            loop.set_cut(start + ((end - start)//2))

        loop.set_extended(extended)
    return loop


def generate_loops_from_res_selector(pose, designable_selector, focus_residue,
        resbuffer=3, randomize_cutpoints=False):
    """
    Generate loops from designable residues selector. For each stretch of 3 or
    more designable residues, generate a loop object for it.

    This is adopted from the constructor for the Loops class in Rosetta, except
    that here I make sure that all loops are >= 3 residues long.
    """
    loops = rosetta.protocols.loops.Loops()

    prev = False
    for i in range(1, len(designable_selector) + 1):
        if designable_selector[i] and not prev:
            start = i
        elif not designable_selector[i] and prev:
            assert(start != 0)
            if i - start >= 3 and not (start < focus_residue < i):
                if randomize_cutpoints:
                    loops.add_loop(rosetta.protocols.loops.Loop(start, i-1,
                        random.randint(start+1, i-1)))
                else:
                    loops.add_loop(rosetta.protocols.loops.Loop(start, i-1,
                        int(start + (i-start)/2)))
            start = 0
        prev = designable_selector[i]

    """
    Add terminal loops.
    """
    if start:
        end = len(designable_selector)
        if end - start >= 4 and not (start < focus_residue < end):
            if randomize_cutpoints:
                loops.add_loop(rosetta.protocols.loops.Loop(start, end,
                    random.randint(start + 1, end)))
            else:
                loops.add_loop(rosetta.protocols.loops.Loop(start, end,
                    int(start + (end -
                    start + 1)/2)))

    """
    Add loop main loop surrounding focus residue. Done separately to ensure
    focus residue gets added, and we can have different logic for focus
    residue (can have bigger loop or loop size depend on how far the bb should
    move).
    """
    loopstart = max(focus_residue - resbuffer, 2)
    loopend = min(focus_residue + resbuffer, len(designable_selector))
    print('LOOPEND')
    if randomize_cutpoints:
        loops.add_loop(rosetta.protocols.loops.Loop(loopstart,
            loopend, random.randint(loopstart+1, loopend)))
    else:
        loops.add_loop(rosetta.protocols.loops.Loop(loopstart,
            loopend,focus_residue))

    print(loops)
    return loops



def choose_designable_residues(pose, focus_residues, include_focus=True,
        dshell=8.0, pshell=12.0):
    """
    Chooses a shell (for now, might make more sophisticated later) of residues
    to design around the motif residue.
    """

    focus_residue_selector =\
            rosetta.core.select.residue_selector.ResidueIndexSelector(list_to_str(focus_residues))
    designable_selector =\
            rosetta.core.select.residue_selector.NeighborhoodResidueSelector(
                    focus_residue_selector, dshell,
                    include_focus_in_subset=include_focus
            )
    designable_not_selector =\
            rosetta.core.select.residue_selector.NotResidueSelector(
                    designable_selector
                    )
    packable_selector =\
            rosetta.core.select.residue_selector.NeighborhoodResidueSelector(
                    focus_residue_selector, pshell,
                    include_focus_in_subset=include_focus
            )
    repack_only_selector =\
            rosetta.core.select.residue_selector.AndResidueSelector(
                    designable_not_selector, packable_selector
                    )

    design_residues = designable_selector.apply(pose)
    repack_residues = repack_only_selector.apply(pose)

    return design_residues, repack_residues
