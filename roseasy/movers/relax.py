from roseasy.movers.mover import Mover
from pyrosetta import init
from pyrosetta.rosetta.protocols.relax import FastRelax


class Relax(Mover):
    '''Relax your pose'''
    def __init__(self):
        self.edited_movemap =  False
        super().__init__()

    @property
    def mover(self):
        if hasattr(self, '_mvr'):
            return self._mvr
        else:
            return FastRelax(self.sfxn, self.rounds)
    
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

    def update_mover(self):
        self.mover = FastRelax(self.sfxn, self.rounds)
        # Only update movemap if edited (includes calling
        # self.setup_default_movemap())
        if self.edited_movemap:
            self.mover.set_movemap(self.movemap)
        init(' '.join(self.init_args))

    def apply(self):
        self.update_mover()
        self.mover.apply(self.pose)
