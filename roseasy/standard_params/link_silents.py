from roseasy.pipeline import workspace_from_dir
import sys, os, glob

workspace = workspace_from_dir(sys.argv[1])
outpath = os.path.join(
        workspace.focus_dir,
        'outputs'
        )

for silentfile in glob.glob(outpath + '/*/silent.out'):
    os.symlink(silentfile,
            os.path.join(
                sys.argv[2],
                'inputs',
                os.path.basename(os.path.dirname(silentfile)) + '.pdb'
            ))
