from roseasy.movers.mover import Mover
from pyrosetta.rosetta.protocols.relax import FastRelax


class Relax(Mover):
    '''Relax your pose'''
    def __init__(self):
        self.edited_movemap =  False

    @property
    def mvr(self):
        if hasattr(self, '_mvr'):
            return self._mvr
        else:
            return FastRelax(self.sfxn, self.rounds)
    
    @mvr.setter
    def mvr(self, mover):
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
        self.mvr = FastRelax(self.sfxn, self.rounds)
        self.mvr.set_movemap(self.movemap)

    def apply(self):
        self.update_mover()
        self.mvr.apply(self.pose)
