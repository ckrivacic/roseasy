from roseasy import pipeline
from roseasy import big_jobs
import sys, os, glob
import subprocess

def get_workspace(root_dir, step):
    return pipeline.ValidationWorkspace(root_dir, step)

def main():

    workspace, job_info = big_jobs.initiate()
    pdbpath = workspace.input_path(job_info)
    if not os.path.exists(workspace.output_prefix(job_info)):
        os.mkdir(workspace.output_prefix(job_info))
    outpath = workspace.output_path(job_info)
    folder = workspace.output_prefix(job_info)
    test_run = job_info.get('test_run', False)
    combiner = '/wynton/home/kortemme/krivacic/rosetta/source/bin/combine_silent.linuxgccrelease'
    # for folder in workspace.output_subdirs:
    cmd = [
            combiner,
            '-database', '/wynton/home/kortemme/krivacic/rosetta/database',
            '-out:file:silent', os.path.join(folder, 'silent.out')
            ]

    for f in glob.glob(folder + '/*.out'):
            cmd.extend(['-in:file:silent', f,])

    print('Running ROSETTA command:')
    print(' '.join(cmd))

    subprocess.call(cmd)

if __name__=='__main__':
    main()
