#!/usr/bin/env python2

import sys, os, re, json, subprocess
from roseasy import pipeline


def submit_slurm(workspace, **params):
    from klab import process
    """Submit a job to a workstation with a Slurm job scheduler."""
    
    # Parse some job parameters for the keyword arguments.

    params = dict((k, v) for k, v in list(params.items()) if v is not None)
    test_run = params.get('test_run', False)
    nstruct = params.get('nstruct')
    max_runtime = params.get('max_runtime', '12:00:00')
    max_memory = params.get('max_memory', '3G')

    if test_run:
        max_runtime = '0:30:00'

    if nstruct is None:
        raise TypeError("sumbit() requires the keyword argument 'nstruct' for production runs.")

    # Create the job file

    workspace_jobno = workspace.slurm_custom_jobno + 1
    job_file = os.path.join(workspace.focus_dir,
            'slurm_{}'.format(workspace_jobno))

    nstruct = params.get('nstruct')
    f = open(job_file, 'w')
    for n in range(1, nstruct + 1):
        cmd = '{0} {1} {2} {3}'.format(
                workspace.python_path,
                workspace.script_path,
                workspace.focus_dir,
                n
                )
        f.write(cmd + '\n')
    f.close()


    submission_script = '''#!/bin/bash
#
# Simple SLURM script for submitting multiple serial
# commands (e.g. parametric studies) using a script wrapper
# to launch the commands.
#
# To use, change this job script to accommodate
# running your serial application(s) in your WORKDIR
# directory (usually the directory of submission).
# Edit the commands file to specify the executions
# for the launcher (paramrun) to perform.
#-------------------------------------------------------
#-------------------------------------------------------
#
#         <------ Setup Parameters ------>
#
#SBATCH -J {name}
#SBATCH -N 1
#SBATCH -n 1             #use site recommended # of cores
#SBATCH -p skx-normal
#SBATCH -o {logs}.o%j
#SBATCH -e {logs}.e%j
#SBATCH -t {runtime}
##SBATCH -A <acct_name>   #uncomment and insert acct name if necessary
#------------------------------------------------------

#                         # USING SLURM; plugins defines SLURM env. vars.
export LAUNCHER_RMI=SLURM
export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins


#                         # JOB_FILE is a list of executions to run

export LAUNCHER_JOB_FILE={commands}
export LAUNCHER_SCHED=dynamic
export LAUNCHER_WORKDIR={focus_dir}

$LAUNCHER_DIR/paramrun    # will run the executions in the LAUNCHER_JOB_FILE file
                          # "JOB" is a misnomer--these are not slurm jobs
                          # Each line in the commands file is an
                          # execution.'''.format(
          name=params.get('job_name', 'roseasy_job'),
          logs=os.path.join(workspace.focus_dir,
              'logs', params.get('job_name',
                  'roseasy_job')) + '.',
          runtime=max_runtime,
          commands=job_file,
          focus_dir=workspace.focus_dir
          )
    
    with open(workspace.slurm_submit_file, 'w') as f:
        f.write(submission_script)

    status = process.check_output(('bash',workspace.slurm_submit_file)).decode('utf-8')
    print('Job submission status:')
    print(status)

    #submission = 

    with open(workspace.job_info_path(workspace_jobno), 'w') as file:
        json.dump(params, file)


def submit(script, workspace, **params):
    """Submit a job with the given parameters."""
    from klab import cluster, process

    # Make sure the rosetta symlink has been created.

    # if not os.path.exists(workspace.rosetta_dir):
        # raise pipeline.RosettaNotFound(workspace)

    # Parse some job parameters for the keyword arguments.

    params = dict((k, v) for k, v in list(params.items()) if v is not None)
    test_run = params.get('test_run', False)
    nstruct = params.get('nstruct')
    max_runtime = params.get('max_runtime', '6:00:00')
    max_memory = params.get('max_memory', '1G')

    if test_run:
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
    qsub_command += '-N', params.get('job_name'),
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

    workspace = pipeline.workspace_from_dir(sys.argv[1])
    workspace.cd_to_root()

    try:
        print('Trying qsub')
        """Return some relevant information about the currently running job."""
        print_debug_header()
        job_info = read_job_info(workspace.job_info_path(os.environ['JOB_ID']))
        job_info['job_id'] = int(os.environ['JOB_ID'])
        job_info['task_id'] = int(os.environ['SGE_TASK_ID']) - 1
    except:
        try:
            print('Trying slurm')
            # If not qsub, slurm?
            job_info = read_job_info(workspace.slurm_cmd_file)
            print('Read job info')
            job_info['task_id'] = int(sys.argv[2])
            print('Assigned task id')
        except:
            print('Trying local')
            # Apparently this is a local job.
            # TODO: Need a better way to get job info for local jobs.
            job_info = {
                    'inputs': [x for x in workspace.unclaimed_inputs],
                    'nstruct': 1,
                    'test_run': False,
                    'task_id': int(sys.argv[2])
                    }

    return workspace, job_info

def read_job_info(params_path):
    print('reading job info')
    with open(params_path) as f:
        print('opened {}'.format(params_path))
        return json.load(f)
    print('loaded')

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
