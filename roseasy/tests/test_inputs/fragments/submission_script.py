#!/usr/bin/python
#$ -S /usr/bin/python
#$ -N fragment_generation
#$ -o /wynton/home/kortemme/krivacic/cas/roseasy/roseasy/tests/test_inputs/fragments
#$ -e /wynton/home/kortemme/krivacic/cas/roseasy/roseasy/tests/test_inputs/fragments
#$ -cwd
#$ -r y
#$ -t 1-1
#$ -l mem_free=40G
#$ -l scratch=10G
#$ -l h_rt=6:00:00
#$ -q long.q

import sys
from time import strftime
import socket
import os
import platform
import subprocess
import tempfile
import shutil
import glob
import re
import shlex
import traceback


# Utility functions

class ProcessOutput(object):

    def __init__(self, stdout, stderr, errorcode):
        self.stdout = stdout
        self.stderr = stderr
        self.errorcode = errorcode

    def getError(self):
        if self.errorcode != 0:
            return("Errorcode: %d\n%s" % (self.errorcode, self.stderr))
        return None

def Popen(outdir, args):
    subp = subprocess.Popen(shlex.split(" ".join([str(arg) for arg in args])), stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=outdir, env={'SPARKSXDIR' : '/netapp/home/klabqb3backrub/tools/sparks-x'})
    output = subp.communicate()
    return ProcessOutput(output[0], output[1], subp.returncode) # 0 is stdout, 1 is stderr

def shell_execute(command_line):
    subp = subprocess.Popen(command_line, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    output = subp.communicate()
    return ProcessOutput(output[0], output[1], subp.returncode) # 0 is stdout, 1 is stderr

def create_scratch_path():
    path = tempfile.mkdtemp(dir = '/scratch')
    if not os.path.isdir(path):
        raise os.error
    return path

def print_tag(tag_name, content):
    print('<%s>%s</%s>' % (tag_name, content, tag_name))

def print_subprocess_output(subp):
    '''Prints the stdout and stderr output.'''
    if subp:
        if subp.errorcode != 0:
            print('<error errorcode="%s">' % str(subp.errorcode))
            print(subp.stderr)
            print("</error>")
            print_tag('stdout', '\n%s\n' % subp.stdout)
        else:
            print_tag('success', '\n%s\n' % subp.stdout)
            print_tag('warnings', '\n%s\n' % subp.stderr)


# Job/task parameters

task_id = os.environ.get('SGE_TASK_ID')
job_id = os.environ.get('JOB_ID')
array_idx = int(task_id) - 1              # this can be used to index Python arrays (0-indexed) rather than task_id (based on 1-indexing)
task_root_dir = None
subp = None
errorcode = 0 # failed jobs should set errorcode


# Markup opener
print('<task type="fragment_generation" job_id="%s" task_id="%s">' % (job_id, task_id))


# Standard task properties - start time, host, architecture
print_tag("start_time", strftime("%Y-%m-%d %H:%M:%S"))
print_tag("host", socket.gethostname())
print_tag("architecture", platform.machine() + ', ' + platform.processor() + ', ' + platform.platform())


# Set up a scratch directory on the node
scratch_path = create_scratch_path()
print_tag("cwd", scratch_path)


# Job data arrays. This section defines arrays that are used for tasks.
chains = ['A']
pdb_ids = ['4cmp']
fasta_files = ['/wynton/home/kortemme/krivacic/cas/roseasy/roseasy/tests/test_inputs/fragments/4cmpA/4cmpA.fasta']


# Job setup. The job's root directory must be specified inside this block.
chain = chains[array_idx]
pdb_id = pdb_ids[array_idx]
fasta_file = fasta_files[array_idx]
task_root_dir = os.path.split(fasta_file)[0]
job_root_dir = os.path.split(task_root_dir)[0]
print_tag('job_root_dir', job_root_dir)

sys.path.insert(0, job_root_dir)
from post_processing import post_process

# Copy resources
shutil.copy(fasta_file, scratch_path)
if not os.path.exists(task_root_dir):
    raise Exception("You must set the task's root directory so that the script can clean up the job.")


# Job execution block. Note: failed jobs should set errorcode.
print("<output>")
cmd_args = [c for c in ['/wynton/home/kortemme/krivacic/software/klab/klab/bio/fragments/make_fragments_QB3_cluster.pl', '-verbose', '-id', pdb_id + chain, '', '-frag_sizes 3,9', '-n_frags 200', '-n_candidates 1000', fasta_file] if c]
print_tag('cmd', ' '.join(cmd_args))

subp = Popen(scratch_path, cmd_args)
sys.stdout.write(subp.stdout)

if True:
    print("<gzip>")
    for f in glob.glob(os.path.join(scratch_path, "*mers")) + [os.path.join(scratch_path, 'ss_blast')]:
        if os.path.exists(f):
            subpzip = Popen(scratch_path, ['gzip', f])
            print(f)
    print("</gzip>")

os.remove(fasta_file)
print("</output>")


# Post-processing. Copy files from scratch back to /netapp.
shutil.rmtree(task_root_dir, ignore_errors=True)
shutil.copytree(scratch_path, task_root_dir)
shutil.rmtree(scratch_path)
# Run post-processing script
task_dirname = os.path.split(task_root_dir)[1]
post_process(task_dirname)

# Print task run details. The full path to qstat seems necessary on the QB3 cluster if you are not using a bash shell.
task_usages = shell_execute('/usr/local/sge/bin/linux-x64/qstat -j %s' % job_id)
if task_usages.errorcode == 0:
  try:
    print("<qstat>")
    print(task_usages.stdout)
    print("</qstat>")
    mtchs = re.match('.*?usage\s*(\d+):(.*?)\n.*', task_usages.stdout, re.DOTALL)
    print(mtchs)
    if mtchs and str(mtchs.group(1)) == str(task_id):
      task_properties = [s.strip() for s in mtchs.group(2).strip().split(",")]
      for tp in task_properties:
         if tp:
           prp=tp.split('=')[0]
           v=tp.split('=')[1]
           print('<task_%s>%s</task_%s>' % (prp, v, prp))
  except Exception, e:
    print('<qstat_parse_error>')
    print(str(e))
    print(traceback.format_exc())
    print('</qstat_parse_error>')
else:
    print_tag('qstat_error', task_usages.stderr)


# Print the end walltime and close the outer tag
print_tag("end_time", strftime("%Y-%m-%d %H:%M:%S"))
print("</task>")


# Exit the job with the errorcode set in the execution block
if errorcode != 0:
    sys.exit(subp.errorcode)