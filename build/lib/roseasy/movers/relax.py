from movers import Mover
from pyrosetta.rosetta.protocols.relax import FastRelax


class Relax(Mover):
    '''Relax your pose'''
    def __init(self):
        self.edited_movemap =  False

    @property
    def mvr(self):
        return FastRelax(self.sfxn, self.rounds)

    @property
    def rounds(self):
        if self.relax_rounds:
            return self.relax_rounds
        else:
            return 5

    @rounds.setter
    def rounds(self, relax_rounds):
        self.relax_rounds = relax_rounds

    def update_mover(self):
        self.mvr = FastRelax(self.sfxn, self.rounds)
        self.mvr.set_movemap(self.movemap)

    def apply(self):
        self.update_mover()
        self.mvr.apply(self.pose)
