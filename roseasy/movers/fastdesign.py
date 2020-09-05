from roseasy.movers.mover import Mover
from roseasy.utils.mover_utils import setup_task_factory
from roseasy.utils.mover_utils import choose_designable_residues
from roseasy.utils.mover_utils import setup_movemap_from_resselectors
from pyrosetta import init
from pyrosetta.rosetta.protocols.denovo_design.movers import FastDesign as rFastDesign


class FastDesign(Mover):
    '''Relax your pose'''
    def __init__(self):
        self.edited_movemap =  False
        super().__init__()

    @property
    def mover(self):
        if hasattr(self, '_mvr'):
            return self._mvr
        else:
            return rFastDesign(self.sfxn, self.rounds)
    
    @mover.setter
    def mover(self, mover):
        self._mvr = mover

    @property
    def rounds(self):
        if hasattr(self, '_relax_rounds'):
            return self._relax_rounds
        else:
            return 5

    @rounds.setter
    def rounds(self, relax_rounds):
        self._relax_rounds = relax_rounds

    def task_factory_from_range(self, start, end):
        focus = [i for i in range(start, end+1)]
        design_residues, repack_residues =\
                choose_designable_residues(self.pose, focus)
        task_factory = setup_task_factory(self.pose, design_residues,
                repack_residues)
        movemap = setup_movemap_from_resselectors(design_residues, repack_residues)
        self.movemap = movemap
        self.task_factory = task_factory

    def update_mover(self):
        self.mover = rFastDesign(self.sfxn, self.rounds)
        self.mover.set_up_default_task_factory()
        if hasattr(self, '_task_factor'):
            self.mover.set_task_factory(self.task_factory)
        # Only update movemap if edited (includes calling
        # self.setup_default_movemap())
        if self.edited_movemap:
            self.mover.set_movemap(self.movemap)
        init(' '.join(self.init_args))

    def apply(self):
        self.update_mover()
        self.mover.apply(self.pose)
