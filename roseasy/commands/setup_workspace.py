"""\
Query the user for all the input data needed for a design.  This includes a 
starting PDB file, the backbone regions that will be remodeled, the residues 
that will be allowed to design, and more.  A brief description of each field is 
given below.  This information is used to build a workspace for this design 
that will be used by the rest of the scripts in this pipeline.  

Usage:
    roseasy setup_workspace <workspace> [--remote] [--overwrite]

Options:
    --remote, -r
        Setup a link to a design directory on a remote machine, to help with 
        transferring data between a workstation and a cluster.  Note: the 
        remote and local design directories must have the same name.

    --overwrite, -o
        If a design with the given name already exists, remove it and replace 
        it with the new design created by this script.
"""

import os, re, shutil, subprocess, glob

def ensure_path_exists(path):
    path = os.path.abspath(os.path.expanduser(path))
    if not os.path.exists(path):
        raise ValueError("'{0}' does not exist.".format(path))
    return path

class PythonPath:
    prompt = "Path to python binary where PyRosetta is installed: "
    description = """\
Python path: Path to your python binary. Make sure PyRosetta package is
installed. You can skip this step in the future by adding a ROSEASY_PYTHON
variable to your environment.
    """

    if 'ROSEASY_PYTHON' in os.environ:
        print('ROSEASY_PYTHON found in environment: {}'.format(
            os.environ.get('ROSEASY_PYTHON')))
        setting_env_var = os.environ.get('ROSEASY_PYTHON')

    @staticmethod
    def install(workspace, python_path):
        python_path = ensure_path_exists(python_path)

        os.symlink(python_path, workspace.python_path)

class RosettaDir:
    prompt = "Path to rosetta: "
    description = """\
    Rosetta checkout: Path to the main directory of a Rosetta source code checkout.  
    This is the directory called 'main' in a normal rosetta checkout.  Rosetta is 
    used both locally and on the cluster, but the path you specify here probably 
    won't apply to both machines.  You can manually correct the path by changing 
    the symlink called 'rosetta' in the workspace directory."""

    if 'ROSETTA_DIR' in os.environ:
        setting_env_var = os.environ.get('ROSETTA_DIR')

    @staticmethod
    def install(workspace, rosetta_dir):
        rosetta_dir = ensure_path_exists(rosetta_dir)
        rosetta_subdirs = [
                os.path.join(rosetta_dir, 'database'),
                os.path.join(rosetta_dir, 'tests'),
                os.path.join(rosetta_dir, 'source'),
                os.path.join(rosetta_dir, 'source', 'bin'),
        ]
        rosetta_subdirs_exist = list(map(os.path.exists, rosetta_subdirs))

        if not all(rosetta_subdirs_exist):
            message = [
                    "'{0}' does not appear to be the main rosetta directory.".format(rosetta_dir),
                    "The following subdirectories are missing:"
            ]
            for path in rosetta_subdirs:
                if not os.path.exists(path):
                    message.append('    ' + path)
            raise ValueError('\n'.join(message))

        os.symlink(rosetta_dir, workspace.rosetta_dir)


class InputPdb:
    prompt = "Path to the input PDB file: "
    description = """\
Input PDB file: A structure containing the functional groups to be positioned.  
This file should already be parse-able by rosetta, which often means it must be 
stripped of waters and extraneous ligands."""

    @staticmethod
    def install(workspace, pdb_path):
        pdb_path = ensure_path_exists(pdb_path)
        destination = os.path.join(workspace.root_dir, 'input.pdb.gz')
        if pdb_path.endswith('.pdb.gz'):
            shutil.copyfile(pdb_path, destination)
        elif pdb_path.endswith('.pdb'):
            subprocess.call('gzip -c {0} > {1}'.format(
                    pdb_path, destination), shell=True)
        else:
            raise ValueError("'{0}' is not a PDB file.".format(pdb_path))


class LoopsFile:
    prompt = "Path to the loops file: "
    description = """\
Loops file: A file specifying which backbone regions will be allowed to move.  
These backbone regions do not have to be contiguous, but each region must span 
at least 4 residues. Leave blank if you don't have a default loops file
you want to install or if you want to add one later. Loops files can be
added to subdirectories to apply them to only specific steps."""

    @staticmethod
    def install(workspace, loops_path):
        if loops_path:
            loops_path = ensure_path_exists(loops_path)
            shutil.copyfile(loops_path, workspace.loops_path)
        else:
            pass


class ParamsFile:
    prompt = "Path to any ligand params files: "
    description = """\
Loops file: A file specifying scoring parameters for ligands not
recognized by Rosetta by default. Leave blank if you don't have any params files
you want to install or if you want to add one later. """

    @staticmethod
    def install(workspace, params_path):
        if params_path:
            params_path = ensure_path_exists(params_path)
            shutil.copyfile(params_path,
                    os.path.join(workspace.root_dir, 'project_params',
                        os.path.basename(params_path)))
        else:
            pass


class Resfile:
    prompt = "Path to resfile: "
    description = """\
Resfile: A file specifying which positions to design and which positions to 
repack.  I recommend designing as few residues as possible outside the loops."""

    @staticmethod
    def install(workspace, resfile_path):
        if resfile_path:
            resfile_path = ensure_path_exists(resfile_path)
            shutil.copyfile(resfile_path, workspace.resfile_path)
        else:
            pass


class RestraintsFile:
    prompt = "Path to restraints file: "
    description = """\
Restraints file: A file describing the geometry you're trying to design.  In 
rosetta parlance, this is more often (inaccurately) called a constraint file.  
Note that restraints are not used during the validation step."""

    @staticmethod
    def install(workspace, restraints_path):
        restraints_path = ensure_path_exists(restraints_path)
        shutil.copyfile(restraints_path, workspace.restraints_path)


class ScoreFunction:
    prompt = "Path to weights file [optional]: "
    description = """\
Score function: A file that specifies weights for all the terms in the score 
function, or the name of a standard rosetta score function.  The default is 
talaris2014.  That should be ok unless you have some particular interaction 
(e.g. ligand, DNA, etc.) that you want to score in a particular way."""

    @staticmethod
    def install(workspace, scorefxn_path):

        # If the user didn't specify a score function, use talaris2014 by 
        # default.

        if not scorefxn_path:
            scorefxn_path = 'talaris2014'

        # Figure out if the user is specifying the name of a standard score 
        # function.  If so, get the path to the real score file.

        if not os.path.exists(scorefxn_path):
            builtin_scorefxn_path = workspace.rosetta_subpath(
                    'database', 'scoring', 'weights', scorefxn_path + '.wts')
            if os.path.exists(builtin_scorefxn_path):
                scorefxn_path = builtin_scorefxn_path

        # Copy the score function into the workspace.

        if scorefxn_path:
            scorefxn_path = ensure_path_exists(scorefxn_path)
            shutil.copyfile(scorefxn_path, workspace.scorefxn_path)


class DefaultScripts:
    prompt = None
    Description = """\
Installing default scripts."""
    
    @staticmethod
    def install(workspace):
        script_dir = os.path.join(os.path.dirname(__file__), '..',
                'standard_params')
        python = glob.glob(script_dir + '/*.py')
        yaml = glob.glob(script_dir + '/*.yml')
        wts = glob.glob(script_dir + '/*.wts')
        sho = glob.glob(script_dir + '/*.sho')
        for script in python + yaml + wts:
            script_path = os.path.join(script_dir, script)
            workspace_path = os.path.join(workspace.standard_params_dir,
                    os.path.basename(script))
            shutil.copyfile(script_path, workspace_path)
        for script in sho:
            script_path = os.path.join(script_dir, script)
            workspace_path = os.path.join(workspace.root_dir,
                    os.path.basename(script))
            shutil.copyfile(script_path, workspace_path)


class FlagsFile:
    prompt = "Path to flags file [optional]: "
    description = """\
Flags file: A file containing command line flags that should be passed to every 
invocation of rosetta for this design.  For example, if your design involves a 
ligand, put flags related to the ligand parameter files in this file."""

    @staticmethod
    def install(workspace, flags_path):
        if flags_path:
            flags_path = ensure_path_exists(flags_path)
            shutil.copyfile(flags_path, workspace.flags_path)
        else:
            scripting.touch(workspace.flags_path)


class RsyncUrl:
    prompt = "Path to project on remote host: "
    description = """\
Rsync URL: An ssh-style path to the directory that contains (i.e. is one level 
above) the remote workspace.  This workspace must have the same name as the 
remote one.  For example, to link to "~/path/to/my_design" on chef, name this 
workspace "my_design" and set its rsync URL to "chef:path/to"."""

    @staticmethod
    def install(workspace, rsync_url):
        with open(workspace.rsync_url_path, 'w') as file:
            file.write(rsync_url.strip() + '\n')



from klab import scripting
import docopt
from roseasy import pipeline

@scripting.catch_and_print_errors()
def main():
    arguments = docopt.docopt(__doc__)
    workspace = pipeline.Workspace(arguments['<workspace>'])

    # Make a new workspace directory.

    if workspace.incompatible_with_fragments_script:
        scripting.print_error_and_die("""\
Illegal character(s) found in workspace path:

  {}

The full path to a workspace must contain only characters that are alphanumeric
or '.' or '_'.  The reason for this ridiculous rule is the fragment generation
script, which will silently fail if the full path to its input file contains 
any characters but those.""", workspace.abs_root_dir)

    if workspace.exists():
        if arguments['--overwrite']:
            shutil.rmtree(workspace.root_dir)
        else:
            scripting.print_error_and_die("""\
Design '{0}' already exists.  Use '-o' to overwrite.""", workspace.root_dir)

    workspace.make_dirs()

    # Decide which settings to ask for.

    if arguments['--remote']:
        installers = (
                RosettaDir,
                RsyncUrl,
                PythonPath,
        )
    else:
        installers = (
                RosettaDir,
                InputPdb,
                PythonPath,
                DefaultScripts,
                LoopsFile,
                Resfile,
                ParamsFile,
        )

    # Get the necessary settings from the user and use them to fill in the 
    # workspace.

    print("Please provide the following pieces of information:")
    print()

    scripting.use_path_completion()

    for installer in installers:

        # If the installer doesn't have a prompt, just install it without 
        # asking any questions.

        if installer.prompt is None:
            installer.install(workspace)
            continue

        # Check if an environment variable is defined for the installer
        if hasattr(installer, 'setting_env_var'):
            installer.install(workspace, installer.setting_env_var)
            continue

        # Otherwise, print a description of the setting being installed and 
        # prompt the user for a value.


        print(installer.description)
        print()

        while True:
            try:
                setting = input(installer.prompt)
                installer.install(workspace, setting)
            except (ValueError, IOError) as problem:
                print(problem)
                continue
            except (KeyboardInterrupt, EOFError):
                shutil.rmtree(workspace.root_dir)
                scripting.print_error_and_die("\nReceived exit command, no workspace created.")
            else:
                break

        print()

    # If we made a link to a remote workspace, immediately try to synchronize 
    # with it.  Rsync will say whether or not it succeeded.  Otherwise just 
    # print a success message.

    if arguments['--remote']:
        pipeline.fetch_data(workspace.root_dir)
    else:
        print("Setup successful for design '{0}'.".format(workspace.root_dir))

if __name__=='__main__':
    main()
