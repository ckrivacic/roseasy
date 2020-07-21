from roseasy.movers.mover import Mover
from pyrosetta.rosetta.protocols.loop_modeler import LoopModeler


class LoopModeler(Mover):
    def __init__(self, config):
        self.mvr = LoopModeler()
        self.config = config
        return

    @property
    def config(self):
        return self._config

    @config.setter
    def config(self, config):
        self._config = config
        if config=='ngk':
            self.mvr.setup_kic_config()
        elif config=='lhk':
            self.mvr.setup_loophash_kic_config()
    
    @property
    def mvr(self):
        if hasattr(self, '_mvr'):
            return self._mvr
        else:
            return LoopModeler()
    
    @mvr.setter
    def mvr(self, mover):
        self._mvr = mover

    @property
    def loops(self):
        return self._loops

    @loops.setter
    def loops(self, loops):
        self._loops = loops

    def update_mover(self):
        self.mvr = FastRelax(self.sfxn, self.rounds)
        self.mvr.set_movemap(self.movemap)

    def apply(self):
        self.update_mover()
        self.mvr.apply(self.pose)
