from roseasy import pipeline
import sys, os
import subprocess

def main():
    combiner = '/wynton/home/kortemme/krivacic/rosetta/source/bin/combine_silent.linuxgccrelease'
    workspace = pipeline.workspace_from_dir(sys.argv[1])
    for folder in workspace.output_subdirs:
        cmd = [
                combiner,
                '-database', '/wynton/home/kortemme/krivacic/rosetta/database',
                '-in:file:silent', os.path.join(folder, '*.out'),
                '-out:file:silent', os.path.join(folder, 'silent.out')
                ]

        print('Running ROSETTA command:')
        print(' '.join(cmd))

        subprocess.call(cmd)

if __name__=='__main__':
    main()
