from pyrosetta import *


class Mover(object):
    '''
    A super-clas to handle data commonly associated with movers. Unlike
    Rosetta movers, these are capable of holding a Pose; this is to make
    it possible to alter movemaps, loops, etc. from information about
    the pose.
    '''
    def __init__(self):
        self.edited_movemap = False
        self.init_args = []

    def setup_default_movemap(self, bb=None, chi=None):
        """
        Set up the movemap. Optionally specify backbone and sidechain
        degrees of freedom, using either a list of residues or a residue
        selector.
        Should probably get rid of "default" in the function name
        because this can be used as non-default movemap.
        """
        mm = rosetta.core.kinematics.MoveMap()
        if bb is not None:
            mm.set_bb(False)
            if max(bb) > 1:
                # This means it's a list of residues
                for i in bb:
                    mm.set_bb(i, True)
            else:
                # Otherwise it's a res selector
                for i in range(1, len(bb) + 1):
                    if bb[i]:
                        mm.set_bb(i, True)
        else:
            mm.set_bb(True)

        if chi is not None:
            mm.set_chi(False)
            if max(chi) > 1:
                # This means it's a list of residues
                for i in chi:
                    mm.set_chi(i, True)
            else:
                # Otherwise it's a res selector
                for i in range(1, len(chi) + 1):
                    if chi[i]:
                        mm.set_chi(i, True)
        else:
            mm.set_chi(True)

        self.edited_movemap = True

        return mm
    
    '''
    I initially wanted workspace to be a part of the mover class to
    automatically get information about things like resfiles, loops
    files, fragments files, etc. but ultimately decided this would make
    design scripts harder to follow, and saving a few lines of code here
    and there wasn't worth it.
    @property
    def workspace(self):
        return self._workspace

    @workspace.setter
    def workspace(self, workspace):
        self._workspace = workspace
    '''

    @property
    def init_args(self):
        return self._init_args

    @init_args.setter
    def init_args(self, args):
        self._init_args = []
        for arg in args:
            self._init_args.append(arg)

    def add_init_args(self, args):
        for arg in args:
            self._init_args.append(arg)

    def add_init_arg(self, arg):
        self._init_args.append(arg)

    def init(self):
        args = ''
        for arg in self._init_args:
            args += arg
            args += ' '
        init(args)

    @property
    def pose(self):
        if hasattr(self, '_pose'):
            return self._pose
        else:
            return None

    @pose.setter
    def pose(self, pose):
        self._pose = pose

    @property
    def movemap(self):
        if hasattr(self, '_mm'):
            return self._mm
        else:
            # Don't want to accidentally call setup_default_movemap(),
            # so raise an error.
            raise
            #self._mm = self.setup_default_movemap()
            #return self._mm

    @movemap.setter
    def movemap(self, movemap):
        self._mm = movemap
        # Some movers will have their own defaults for the movemap, so
        # we want to keep track of if we change it so that if not, those
        # movers can use their own defaults.
        self.edited_movemap = True

    def reset_movemap(self):
        setup_movemap()
        self.edited_movemap = False

    @property
    def sfxn(self):
        if hasattr(self, '_sfxn'):
            return self._sfxn
        else:
            return create_score_function('ref2015')

    @sfxn.setter
    def sfxn(self, scorefunction):
        '''
        Select a scorefunction. Can either create a custom scorefunction or
        pass the string of a default scorefunction.
        '''
        if type(scorefunction) == type('string'):
            self._sfxn = create_score_function(scorefunction)
        else:
            self._sfxn = scorefunction
        print(self._sfxn)

    '''
    Let's see if LoopModeler will automatically find a centroid sfxn
    first.
    @property
    def sfxn_cen(self):
        if hasattr(self, '_sfxn_cen'):
            return self._sfxn_cen
        else:
            return create_score_function('???')
    '''

    @property
    def apply(self):
        '''Meant to be overriden in child classes.'''
        raise NotImplementedError


def insert_alas(pose, position, length, insert_after=True, reset_fold_tree=True, fold_tree_root=1):
    '''Insert a poly-ALA peptide before or after a given position.,
    Set the fold tree to have a cutpoint before or after inserted residues.
    Author: XingJie Pan
    '''
    assert(1 <= position <= pose.size())

    # Set the fold tree with a single cutpoint

    def sub_fold_tree_add_edges_no_jump(ft, root, start, stop):
        '''Add edges to a sub-fold-tree that does not have
        and jumps.'''
        if start < root:
            ft.add_edge(root, start, -1)
        if stop > root:
            ft.add_edge(root, stop, -1)

    if reset_fold_tree:
        cutpoint = position if insert_after else position - 1
        ft = rosetta.core.kinematics.FoldTree()
        
        if fold_tree_root <= cutpoint and cutpoint < pose.size():
            sub_root = pose.size()
            ft.add_edge(fold_tree_root, sub_root, 1)
            sub_fold_tree_add_edges_no_jump(ft, sub_root, cutpoint + 1, pose.size())
            sub_fold_tree_add_edges_no_jump(ft, fold_tree_root, 1, cutpoint)
        
        elif fold_tree_root > cutpoint and cutpoint > 0:
            sub_root = 1
            ft.add_edge(fold_tree_root, sub_root, 1)
            sub_fold_tree_add_edges_no_jump(ft, sub_root, 1, cutpoint)
            sub_fold_tree_add_edges_no_jump(ft, fold_tree_root, cutpoint + 1, pose.size())

        else:
            sub_fold_tree_add_edges_no_jump(ft, fold_tree_root,  1, pose.size())
        
        pose.fold_tree(ft)

    # Append the residues

    residue_type_set = pose.residue_type_set_for_pose()
    new_rsd = rosetta.core.conformation.ResidueFactory.create_residue( residue_type_set.name_map("ALA") )
   
    for i in range(length):
        if insert_after:
            pose.conformation().safely_append_polymer_residue_after_seqpos(new_rsd, position + i, True)
            pose.set_omega(position + i, 180)
        else:
            pose.conformation().safely_prepend_polymer_residue_before_seqpos(new_rsd, position, True)
            pose.set_omega(position, 180)

    if insert_after:
        rosetta.core.conformation.idealize_position(position + length, pose.conformation())
        
        if position + length + 1 <= pose.size():
            rosetta.core.conformation.idealize_position(position + length + 1, pose.conformation())
    else:
        if position - 1 > 0:
            rosetta.core.conformation.idealize_position(position - 1, pose.conformation())
        rosetta.core.conformation.idealize_position(position, pose.conformation())


def mutate_residues(pose, res_list, aa_list, protein_only=True):
    '''Mutate a list of residues. The list of AAs could
    either be 1 letter code or 3 letter code.
    Author: XingJie Pan
    '''
    aa_name_map = {'A':'ALA', 'P':'PRO', 'V':'VAL', 'L':'LEU', 'I':'ILE', 'M':'MET',
                   'F':'PHE', 'Y':'TYR', 'W':'TRP', 'S':'SER', 'T':'THR', 'C':'CYS',
                   'K':'LYS', 'R':'ARG', 'H':'HIS', 'D':'ASP', 'E':'GLU', 'N':'ASN',
                   'Q':'GLN', 'G':'GLY'}


    mutater = rosetta.protocols.simple_moves.MutateResidue()
    for i in range(len(res_list)):
        if protein_only and (not pose.residue(res_list[i]).is_protein()):
            continue

        name = aa_list[i] if len(aa_list[i]) == 3 else aa_name_map[aa_list[i]]
        mutater.set_res_name(name)
        mutater.set_target(res_list[i])
        mutater.apply(pose)


def add_aas(pose, position, sequence, pdbnum=False, chain='A'):
    if pdbnum:
        position = pose.pdb_info().pdb2pose(chain, position)

    insert_alas(pose, position, len(sequence))
    mutate_residues(pose, list(range(position + 1, position + 1 +
        len(sequence))), list(sequence), True)

def close_helix_by_minimization(pose, movable_region_start, movable_region_end, helix_start, helix_end):
    '''Close a gap inside a helix by minimization.
    Return true if the gap could be closed.
    '''
    # Make a clone of poly ALA pose for minimization

    #simple_pose_moves.mutate_pose_to_single_AA(pose, 'ALA')
    rosetta.core.pose.correctly_add_cutpoint_variants(pose)

    # Set hydrogen bond constraints for the linkers and helix

    linker_residues = list(range(movable_region_start, helix_start + 1)) + list(range(helix_end, movable_region_end + 1))
    linker_hbonds = find_bb_hbonds_involving_residues(pose, linker_residues)

    pose.constraint_set().clear()
    helix_hbs = [(i + 4, i) for i in range(helix_start, helix_end - 3)]
    constraint.add_constraints_to_pose(pose, constraint.get_bb_hbond_constraint(linker_hbonds + helix_hbs))

    # Set score function

    sfxn = rosetta.core.scoring.get_score_function()
    sfxn.set_weight(rosetta.core.scoring.base_pair_constraint, 1) #H-bond constraint

    # Set movemap

    mm = rosetta.core.kinematics.MoveMap()
    
    for i in range(movable_region_start, movable_region_end + 1):
        mm.set_bb(i, True)
   
    # Set the minimization mover

    min_opts = rosetta.core.optimization.MinimizerOptions( "lbfgs_armijo_nonmonotone", 0.01, True )
    min_mover = rosetta.protocols.minimization_packing.MinMover()
    min_mover.movemap(mm)
    min_mover.min_options(min_opts)

    # Close the chain

    for chainbreak_weight in [0.5, 1, 5, 10]:
        sfxn.set_weight(rosetta.core.scoring.chainbreak, chainbreak_weight)
        min_mover.score_function(sfxn)
        min_mover.apply(pose)

    chainbreak_energy = pose.energies().total_energies()[rosetta.core.scoring.chainbreak] 
    if chainbreak_energy > 0.2:
        return False

    # Minimize without constraints
    
    sfxn.set_weight(rosetta.core.scoring.base_pair_constraint, 0)
    min_mover.score_function(sfxn)
    min_mover.apply(pose)
    
    return True

def find_bb_hbonds_involving_residues(pose, residues):
    '''Find backbone hbonds involving a given set of residues.
    An Hbond is defined as (donor_res, acceptor_res).
    Ignore the terminal residues.
    '''
    hbset = rosetta.core.scoring.hbonds.HBondSet(pose, bb_only=True)
    hbonds = []

    for i in range(1, hbset.nhbonds() + 1):
        acc = hbset.hbond(i).acc_res()
        don = hbset.hbond(i).don_res()

        # Ignore terminal residues
        if acc in [1, pose.size()] or don in [1, pose.size()]:
            continue

        if acc in residues or don in residues:
            hbonds.append((don, acc))

    return hbonds
