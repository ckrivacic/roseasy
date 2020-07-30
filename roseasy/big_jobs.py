#!/usr/bin/env python2

import sys, os, re, json, subprocess
from roseasy.workspace import pipeline

def submit(script, workspace, **params):
    """Submit a job with the given parameters."""
    from klab import cluster, process

    # Make sure the rosetta symlink has been created.

    if not os.path.exists(workspace.rosetta_dir):
        raise pipeline.RosettaNotFound(workspace)

    # Parse some job parameters for the keyword arguments.

    params = dict((k, v) for k, v in list(params.items()) if v is not None)
    test_run = params.get('test_run', False)
    nstruct = params.get('nstruct')
    max_runtime = params.get('max_runtime', '6:00:00')
    max_memory = params.get('max_memory', '1G')

    if test_run:
        nstruct = 50
        max_runtime = '0:30:00'

    if nstruct is None:
        raise TypeError("sumbit() requires the keyword argument 'nstruct' for production runs.")

    # Submit the job and put it immediately into the hold state.

    qsub_command = 'qsub', '-h', '-cwd'
    qsub_command += '-o', workspace.log_dir
    qsub_command += '-e', workspace.log_dir
    qsub_command += '-t', '1-{0}'.format(nstruct),
    qsub_command += '-l', 'h_rt={0}'.format(max_runtime),
    qsub_command += '-l', 'mem_free={0}'.format(max_memory),
    qsub_command += '-b', 'y',
    qsub_command += workspace.python_path,
    qsub_command += script,
    qsub_command += workspace.focus_dir,

    status = process.check_output(qsub_command).decode('utf-8')
    status_pattern = re.compile(r'Your job-array (\d+).[0-9:-]+ \(".*"\) has been submitted')
    status_match = status_pattern.match(status)

    if not status_match:
        print(status)
        sys.exit()

    # Figure out the job id, then make a params file specifically for it.

    job_id = status_match.group(1)

    with open(workspace.job_info_path(job_id), 'w') as file:
        json.dump(params, file)

    # Release the hold on the job.

    qrls_command = 'qrls', job_id
    process.check_output(qrls_command)
    print(status, end=' ')

def initiate():
    """Return some relevant information about the currently running job."""
    print_debug_header()

    workspace = pipeline.workspace_from_dir(sys.argv[1])
    workspace.cd_to_root()

    job_info = read_job_info(workspace.job_info_path(os.environ['JOB_ID']))
    job_info['job_id'] = int(os.environ['JOB_ID'])
    job_info['task_id'] = int(os.environ['SGE_TASK_ID']) - 1

    return workspace, job_info

def read_job_info(params_path):
    with open(params_path) as file:
        return json.load(file)

def print_debug_info():
    from datetime import datetime
    from socket import gethostname

    print("Date:", datetime.now())
    print("Host:", gethostname())
    print("Command: JOB_ID={0[JOB_ID]} SGE_TASK_ID={0[SGE_TASK_ID]} {1}".format(
            os.environ, ' '.join(sys.argv)))
    print()
    sys.stdout.flush()

def run_command(command):
    print("Working directory:", os.getcwd())
    print("Command:", ' '.join(command))
    sys.stdout.flush()

    process = subprocess.Popen(command)

    print("Process ID:", process.pid)
    print()
    sys.stdout.flush()

    process.wait()

def print_debug_header():
    from datetime import datetime
    from socket import gethostname

    print("Date:", datetime.now())
    print("Host:", gethostname())
    print("Python:", sys.executable or 'unknown!')
    print("Command: JOB_ID={0[JOB_ID]} SGE_TASK_ID={0[SGE_TASK_ID]} {1}".format(
            os.environ, ' '.join(sys.argv)))
    print()
    sys.stdout.flush()
