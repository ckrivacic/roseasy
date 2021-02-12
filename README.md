The purpose of this repository will be to help track designs and make Rosetta 
job submissions easier. 


### DEPENDENCIES

matplotlib version 3.1.3  
pygobject (try `conda install -c conda-forge pygobject`)  
pango (easiest to install via `conda install -c anaconda pango`)  
~~pygtk (easiest to install via `conda install -c vgauthier pygtk`)~~  
pygtk (conda install -c pkgw-roge gtk3)
pandas  
numpy  
docopt  
yaml  
klab (python3_compatibility branch)  


### INSTALLATION


`git clone https://github.com/ckrivacic/roseasy.git`  
`python setup.py develop`  

### USAGE  

`roseasy <command> <arguents> [options]`

  Note that you do not have the type the full command name; RosEasy is perfectly 
  happy with the first word or even few letters of a command unless there is 
  ambiguity.


### COMMANDS  

`setup_workspace <name> [-r]`  
  Set up a new or remote workspace. This will walk you through the workspace
  setup process. If setting up a remote workspace (i.e., on your own 
  computer for analysis), pass the -r option. This will set up a link between 
  the new workspace and the workspace on the cluster (make sure they are named 
  the same).


`submit <workspace> <script path> [options]`  
  Submits a job to be run on the cluster. This command takes a 
  python script as an argument. The only hard requirement is that the 
  python script defines which type of workspace Roseasy should create 
  for the task. The rest of the script may run a PyRosetta job, or simply 
  call RosettaScripts. It is recommended that you take an existing script 
  (found in standard_params) as a template, as there are a few lines at 
  the beginning of each script that handle figuring out which input file 
  goes with which task number and whether the script should be run in test 
  mode or not. You can also pass the --make-dirs option, which 
  causes the workspace to set up the sub-workspace without submitting any jobs.
  This is useful for example when you want to handle step inputs manually.
  Note that when the root directory is given as the <workspace> argument, 
  RosEasy will automatically create a NEW sub-workspace, even if the previous 
  step is empty, so it is recommended that whenever possible you create the new 
  sub-workspace FIRST, either by passing --make-dirs or by running the 
  pick_designs command, and then pass that sub-workspace as the <workspace> 
  argument.


`generate_fragments <workspace> <step_number> [options]`  
  Generates fragments for an FKIC run.


`add_residues <pdb file or folder> <residue_string> <pdb_position>`  
  Inserts residues into a PDB file and does a brief minimization 
  to close the chain break. You may specify any sequence to insert 
  using one-letter AA name formatting (i.e. 'SALTY'). The residues 
  will be inserted immediately after the <pdb_position> argument. 
  Defaults to chain A, pass the --chain option to specify a different 
  chain.


`pick_designs_to_validate <step> [<picks_file>] [options]`  
  Pick the designs that are at or near the Pareto front of the given metrics to
  validate in the next step. You may also specify thresholds that designs much pass 
  in order to be chosen for the next step. This will create the workspace for the 
  next step and symlink the chosen input files, and is reccommended before any job.
  For more information type 'roseasy pick_designs --help'
  

`fetch_data <workspace>`  
  On a "remote" workspace, fetch data from the main workspace.


`push_data <workspace>`  
  Push data from a "remote" workspace to the main workspace.


`plot <directory> [options]`  
  Generates a GUI for viewing designs.


### EXAMPLE USAGE (basic)

  First, set up a workspace on the cluster. It will ask you for a path to your 
  Rosetta directory (the 'main' folder, rosetta/main in most installations), an 
  input PDB file, and your Python executable where PyRosetta is installed. You 
  will also be asked to provide a loop file, but you can ignore this if it is not 
  applicable to your usage, or provide one later.


*Cluster:*  
`roseasy setup test`


  On your local workstation, set up the remote workspace.

*Local:*  
`roseasy setup test -r`


  Now, on the cluster, run your first script:

*Cluster:*  
`roseasy submit test test/standard_params/relax.py`


  When it's done running (check job status with the 'qstat' command), 
  pull the data to your computer for analysis and to pick designs 
  for the next step, then push your changes to the cluster.

*Local:*  
`roseasy fetch test`  
`roseasy pick 2 test/standard_params/picks.yml`  
`roseasy push test`  


  Now you're ready to run the next step, in this example FKIC.

*Cluster:*

`roseasy generate_fragments test/02`  
(wait for fragment generation to finish)  
`roseasy submit 02 standard_params/fkic.py`  


  Pull your decoys to your local workstation and view them in the GUI.

*Local:*  
`roseasy fetch test`  
`roseasy plot test/02/outputs/*`  
