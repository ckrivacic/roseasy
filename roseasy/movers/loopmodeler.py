import sys
from roseasy.movers.mover import Mover
from roseasy.parsers import parse_loops
from roseasy.utils.mover_utils import generate_loops_from_range
from roseasy.utils.mover_utils import generate_loop_from_range
from pyrosetta.rosetta.protocols.loop_modeler import LoopModeler as rLoopModeler
from pyrosetta import init


class FragmentError(Exception):
    pass

class LoopModeler(Mover):
    def __init__(self):
        self.mover = rLoopModeler()
        self.edited_fa_cycles = False
        self.edited_cen_cycles = False
        self.perturb_seq=True
        self.seqpose_no_mutate=''
        super().__init__()

    def configure(self):
        if self.config=='ngk':
            self.mover.setup_kic_config()
        elif self.config=='lhk':
            self.mover.setup_loophash_kic_config(self.perturb_seq,
                    self.seqpose_no_mutate)
        elif self.config=='fkic':
            '''
            Commented out for now b/c I don't think I want to go with
            workspace as attributes for movers.
            if not hasattr(self.workspace, 'fragments_flags'):
                print('Error: Incorrect workspace type for fragment KIC. '\
                        ' Make sure you are initiating the mover with a '\
                        'WithFragmentLibs workspace.')
                sys.exit()
            Instead...
            '''
            self.mover.setup_kic_with_fragments_config()

    @property
    def config(self):
        return self._config

    def mark_as_test_run(self):
        self.mover.centroid_stage().mark_as_test_run()
        self.mover.fullatom_stage().mark_as_test_run()

    @config.setter
    def config(self, config):
        self._config = config
        self.configure()

    @property
    def fa_temp_cycles(self):
        return self._fa_temp_cycles

    @fa_temp_cycles.setter
    def fa_temp_cycles(self, cycles):
        self.edited_fa_cycles = True
        self._fa_temp_cycles = cycles

    @property
    def cen_temp_cycles(self):
        return self._cen_temp_cycles

    @cen_temp_cycles.setter
    def cen_temp_cycles(self, cycles):
        self.edited_cen_cycles = True
        self._cen_temp_cycles = cycles

    @property
    def frag_sizes(self):
        '''Set the fragment sizes.
        Note:  Recommended that you use workspace.fragments_flags to
        set flags for this mover, but these are typical defaults.'''
        if hasattr(self, '_frag_sizes'):
            return self._frag_sizes
        else:
            return [3, 9]

    @property
    def frag_files(self):
        try:
            return self._frag_files
        except AttributeError:
            print('No definition found for fragment files. Pass '\
                    'fragment files before trying to run KIC with '\
                    'fragments using LoopModeler.frag_files = <path>')
            raise

    @frag_files.setter
    def frag_files(self, paths):
        '''Set path to fragment files for Fragment KIC. Takes a list of
        paths!'''
        '''
        Not doing workspaces here anymore.
        if not hasattr(self.workspace, 'fragments_flags'):
            print('Error: Incorrect workspace type for fragment KIC. '\
                    ' Make sure you are initiating the mover with a '\
                    'WithFragmentLibs workspace.')
            sys.exit()
        '''
        '''
        TODO: Automatically get sizes from paths
        '''
        self._frag_files = paths

    @frag_sizes.setter
    def frag_sizes(self, sizes):
        '''Takes a list of ints'''
        self._frag_sizes = sizes

    @property
    def fragments_flags(self):
        if not hasattr(self, '_fragments_flags'):
            flags = []
            # Make flags based on other defined attributes
            if hasattr(self, '_frag_files') and hasattr(self,
                    '_frag_sizes'):
                flags.append('-loops:frag_sizes')
                flags.extend(list(map(str, self.frag_sizes)))
                flags.append('-loops:frag_files')
                flags.extend(self.frag_files)
                self._fragments_flags = flags
            else:
                raise FragmentError('Error: Either fragment files or '\
                        'fragment sizes were undefined, and no '\
                        'fragment flags were provided.')
        else:
            return self._fragments_flags
    
    @fragments_flags.setter
    def fragments_flags(self, flags):
        '''Set up custom fragments flags. Takes a list of command-line
        arguments for Rosetta to be used with Init. Recommended that you
        use a workspace to set this flag.'''
        self._fragments_flags = flags
        init(' '.join(self._fragments_flags))
    
    @property
    def mover(self):
        if hasattr(self, '_mvr'):
            return self._mvr
        else:
            return LoopModeler()
    
    @mover.setter
    def mover(self, mover):
        self._mvr = mover

    @property
    def loops(self):
        if hasattr(self, '_loops'):
            return self._loops
        else:
            print('WARNING: Loops have not been set!')
            return rosetta.protocols.loops.Loops()

    @loops.setter
    def loops(self, loops):
        self._loops = loops
        self.mover.set_loops(self._loops)

    def add_loop_obj(self, loop):
        self.loops.add_loop(loop)

    def loop_from_range(self, start, end):
        '''Set loop from a residue range'''
        loops = generate_loops_from_range(start, end)
        self.loops = loops

    def loops_from_file(self, loops_path):
        '''Set loop from a loops file'''
        self.loops = parse_loops(loops_path)

    def add_loop_range(self, start, end):
        loop = generate_loop_from_range(start, end)
        self.add_loop_obj(loop)

    def update_mover(self):
        if self.edited_movemap:
            self.mover.set_movemap(self.movemap)
        if self.edited_fa_cycles:
            self.mover.fullatom_stage().set_temp_cycles(self.fa_temp_cycles)
        if self.edited_cen_cycles:
            self.mover.centroid_stage().set_temp_cycles(self.cen_temp_cycles)
        self.mover.set_loops(self.loops)
        self.mover.set_fa_scorefxn(self.sfxn)

        if self.config == 'fkic':
            # Need to call fragment_flags in case it needs to be
            # constructed, otherwise it won't have the attribute and
            # will throw an error.
            self.fragments_flags
            init_args = ' '.join(self.init_args) + ' ' + ' '.join(self.fragments_flags)
            print('Running init with args {}'.format(init_args))
            init(init_args)
        else:
            init(' '.join(self.init_args))
        self.configure()

    def apply(self):
        self.update_mover()
        self.mover.apply(self.pose)
