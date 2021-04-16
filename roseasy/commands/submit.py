"""
Usage:
    roseasy submit <workspace> [options]
    roseasy submit <workspace> <script> [options]

Options:
    --step, -s
        Which step are you on? Otherwise, this will increment the step
        by 1 from the last step folder (sometimes this is undesirable if
        there is an empty step; to avoid this, make sure your workspace
        is the folder for the step you want run).

    --nstruct NUM, -n NUM  [default: 10]
        The number of jobs to run.

    --max-runtime TIME     [default: 12:00:00]
        The runtime limit for each job.

    --max-memory MEM       [default: 2G]
        The memory limit for each job.

    --local, -l
        Run locally

    --slurm, -u  
        Submit job on a system with a Slurm task scheduler (otherwise
        will default to SGE)

    --clear
        Clear existing results before submitting new job.

    --test-run
        Do a test run

    --make-dirs  
        Just make the directories for the new workspace without
        submitting any jobs.

"""

from klab import scripting, cluster
import docopt
import sys, os, importlib, shutil
from roseasy import pipeline
from roseasy import big_jobs

import subprocess
from io import StringIO

import asyncio
import sys
from asyncio.subprocess import PIPE

@asyncio.coroutine
def read_and_display(*cmd):
    """Read cmd's stdout, stderr while displaying them as they arrive."""
    # start process
    process = yield from asyncio.create_subprocess_exec(*cmd,
            stdout=PIPE, stderr=PIPE)

    # read child's stdout/stderr concurrently
    stdout, stderr = [], [] # stderr, stdout buffers
    tasks = {
        asyncio.Task(process.stdout.readline()): (
            stdout, process.stdout, sys.stdout.buffer),
        asyncio.Task(process.stderr.readline()): (
            stderr, process.stderr, sys.stderr.buffer)}
    while tasks:
        done, pending = yield from asyncio.wait(tasks,
                return_when=asyncio.FIRST_COMPLETED)
        assert done
        for future in done:
            buf, stream, display = tasks.pop(future)
            line = future.result()
            if line: # not EOF
                buf.append(line)    # save for later
                display.write(line) # display in terminal
                # schedule to read the next line
                tasks[asyncio.Task(stream.readline())] = buf, stream, display

    # wait for the process to exit
    rc = yield from process.wait()
    return rc, b''.join(stdout), b''.join(stderr)

def execute(cmd):
    print('Running command:')
    print(' '.join(cmd))
    import time

    with subprocess.Popen(cmd, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, bufsize=1,
           universal_newlines=True) as p, StringIO() as buf:
        for line in p.stdout:
            print(line, end='')
            buf.write(line)
        output = buf.getvalue()
    rc = p.returncode

@scripting.catch_and_print_errors()
def main():
    args = docopt.docopt(__doc__)
    if not args['--local'] and not args['--slurm'] and not args['--make-dirs']:
        cluster.require_qsub()

    workspace = pipeline.workspace_from_dir(args['<workspace>'])
    if args['--step']:
        step = args['--step']
    elif hasattr(workspace, 'step'):
        step = workspace.step
    else:
        step = workspace.get_next_step()

    if not args['<script>']:
        script = os.path.join(workspace.focus_dir, 'run.py')
    else:
        script = args['<script>']

    if not os.path.exists(script):
        raise pipeline.PathNotFound(script)


    # Workspace type is defined in the run script, so we first need to
    # import that.
    script_path = os.path.dirname(script)
    sys.path.insert(1, script_path)
    script_name = os.path.basename(script)[:-3]
    imp = importlib.import_module(script_name)

    workspace = imp.get_workspace(workspace.root_dir, step)
    workspace.check_paths()
    # workspace.check_rosetta()
    workspace.make_dirs()
    if args['--make-dirs']:
        sys.exit()
    # Copying the script to the focus directory helps track exactly what
    # we did at each step.
    shutil.copyfile(script, workspace.script_path)

    if args['--clear'] or args['--test-run']:
        workspace.clear_outputs()

    inputs = [
            x for x in workspace.unclaimed_inputs
            ]

    if len(inputs)==0:
        num_inputs = 1
    else:
        num_inputs = len(inputs)

    if args['--test-run']:
        nstruct = num_inputs * 10
    else:
        nstruct = num_inputs * int(args['--nstruct'])

    if workspace.subdirs:
        for inp in inputs:
            subdir = workspace.output_subdir(inp)
            # scripting.clear_directory(subdir)

    # Submit the job

    if args['--local']:
        print('Running locally.')
        for n in range(1,nstruct + 1):
            cmd = [workspace.python_path]
            cmd.append(workspace.script_path)
            cmd.append(workspace.focus_dir)
            cmd.append(str(n))
            execute(cmd)
            # read_and_display(cmd)

    elif args['--slurm']:
        big_jobs.submit_slurm(
                workspace, 
                nstruct=nstruct,
                max_runtime=args['--max-runtime'],
                max_memory=args['--max-memory'],
                test_run=args['--test-run'],
                job_name=script_name,
                inputs=inputs
                )

    else:
        big_jobs.submit(
                workspace.script_path, workspace,
                nstruct=nstruct,
                max_runtime=args['--max-runtime'],
                max_memory=args['--max-memory'],
                test_run=args['--test-run'],
                job_name=script_name,
                inputs=inputs
                )

if __name__=='__main__':
    main()
