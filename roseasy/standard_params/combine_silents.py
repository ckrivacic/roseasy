from roseasy import pipeline
import sys, os, glob
import subprocess

def main():
    combiner = '/wynton/home/kortemme/krivacic/rosetta/source/bin/combine_silent.linuxgccrelease'
    workspace = pipeline.workspace_from_dir(sys.argv[1])
    for folder in workspace.output_subdirs:
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
