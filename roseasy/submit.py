"""
Usage:
    submit.py <workspace> [options]
    submit.py <workspace> <script> [options]

Options:
    --nstruct NUM, -n NUM  [default: 10]
        The number of jobs to run.

    --max-runtime TIME     [default: 12:00:00]
        The runtime limit for each job.

    --max-memory MEM       [default: 2G]
        The memory limit for each job.

    --local, -l
        Run locally

    --clear
        Clear existing results before submitting new job.

    --test-run
        Do a test run

"""

from klab import scripting, cluster
import docopt
import sys, os, importlib, shutil
from roseasy.workspace import pipeline
from roseasy import big_jobs

@scripting.catch_and_print_errors()
def main():
    args = docopt.docopt(__doc__)
    if not args['--local']:
        cluster.require_qsub()

    workspace = pipeline.workspace_from_dir(args['<workspace>'])
    step = workspace.get_next_step()
    if '<script>' not in args:
        script = os.path.join(workspace.focus_dir, 'run.py')
    else:
        script = args['<script>']

    print(script)

    if not os.path.exists(script):
        raise pipeline.PathNotFound(script)


    # Workspace type can be defined in the run script.
    script_path = os.path.dirname(script)
    sys.path.insert(0, script_path)
    script_name = os.path.basename(script)[:-3]
    imp = importlib.import_module(script_name)

    workspace = imp.Workspace(workspace.root_dir, step)
    workspace.check_paths()
    workspace.check_rosetta()
    workspace.make_dirs()
    shutil.copyfile(script, workspace.script_path)

    if args['--clear'] or args['--test-run']:
        workspace.clear_outputs()

    # Submit the job

    big_jobs.submit(
            script, workspace,
            nstruct=args['--nstruct'],
            max_runtime=args['--max-runtime'],
            max_memory=args['--max-memory'],
            test_run=args['--test-run']
            )

if __name__=='__main__':
    main()
