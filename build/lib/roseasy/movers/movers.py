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
        return

    def setup_default_movemap(self, bb=None, chi=None):
        """
        Set up the movemap. Optionally specify backbone and sidechain
        degrees of freedom, using either a list of residues or a residue
        selector.
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

        self.movemap(mm)
        return mm

    @property
    def pose(self):
        return self.pose

    @pose.setter
    def pose(self, pose):
        self.pose = pose

    @property
    def movemap(self):
        if self.mm:
            return self.mm
        else:
            return setup_default_movemap()

    @movemap.setter
    def movemap(self, movemap):
        self.mm = movemap
        # Some movers will have their own defaults for the movemap, so
        # we want to keep track of if we change it so that if not, those
        # movers can use their own defaults.
        self.edited_movemap=True

    def reset_movemap(self):
        setup_movemap()
        self.edited_movemap = False

    @property
    def sfxn(self):
        if self.scorefunction:
            return self.scorefunction
        else:
            return create_score_function('ref2015')

    @sfxn.setter
    def sfxn(self, scorefunction):
        '''
        Select a scorefunction. Can either create a custom scorefunction or
        pass the string of a default scorefunction.
        '''
        if type(scorefunction) == type('string'):
            self.scorefunction = create_score_function(scorefunction)
        else:
            self.scorefunction = scorefunction

    @property
    def apply(self):
        raise NotImplementedError
