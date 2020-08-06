"""\
Generate fragments for the initial model building simulations.  Note that it's
a little bit weird to use fragments even though the models are allowed to
design in these simulations.  Conformations that are common for the current
sequence but rare for the original one might not get sampled.  However, we
believe that the improved sampling that fragments offer outweighs this
potential drawback.

Usage:
    roseasy generate_fragments <workspace> [options]

Options:
    -L, --ignore-loop-file
        Generate fragments for the entire input structure, not just for the
        region that will be remodeled as specified in the loop file.  This is
        currently necessary only if multiple loops are being remodeled.

    -m, --mem-free=MEM  [default: 2]
        The amount of memory (GB) to request from the cluster.  Bigger systems
        may need more memory, but making large memory requests can make jobs
        take much longer to come off the queue (since there may only be a few
        nodes with enough memory to meet the request).

    -d, --dry-run
        Print out the command-line that would be used to generate fragments,
        but don't actually run it.

    -i, --input_model=MODEL
        Generate fragments for a model other than the workspace input
        model.
"""

import subprocess
from klab import scripting, cluster
import docopt
from roseasy import pipeline

def main():
    args = docopt.docopt(__doc__)
    print(args)
    cluster.require_qsub()

    workspace = pipeline.workspace_from_dir(args['<workspace>'])
    # If not a fragment workspace, make a new one
    if not hasattr(workspace, 'fasta_path'):
        step = workspace.get_next_step()
        workspace = pipeline.FKICModels(workspace.root_dir, step)
    workspace.check_paths()
    workspace.make_dirs()
    workspace.clear_fragments()

    if not '--input_model' in args:
        model = workspace.input_pdb_path
    else:
        model = args['--input_model']

    # Run the fragment generation script.

    generate_fragments = [
            'klab_generate_fragments',
            workspace.input_pdb_path,
            '--outdir', workspace.fragments_dir,
            '--memfree', args['--mem-free'],
            ]
    if not args['--ignore-loop-file']:
        generate_fragments += [
                '--loops_file', workspace.loops_path,
                ]

    if args['--dry-run']:
        print(' '.join(generate_fragments))
    else:
        subprocess.call(generate_fragments)
