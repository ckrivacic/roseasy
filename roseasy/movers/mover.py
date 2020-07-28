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
            self._mm = self.setup_default_movemap()
            return self._mm

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
