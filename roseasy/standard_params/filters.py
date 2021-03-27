from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects
from pyrosetta.rosetta.core.pose import setPoseExtraScore
import os

class FilterContainer(object):
    def __init__(self, workspace, pose, task_id='0000', score_fragments=False,
            test_run=False):
        self.workspace = workspace
        self.task_id = task_id
        self.pose = pose
        self.filters = self.get_default_filters(score_fragments=score_fragments, test_run=test_run)

    def get_default_filters(self, score_fragments=False,
            test_run=False):
        writer = '''
        <WriteFiltersToPose name="writer" prefix="EXTRA_SCORE_"/>
        '''

        buns_all = '''
        <BuriedUnsatHbonds name="Buried Unsat [[-]]"
        report_all_heavy_atom_unsats="true" scorefxn="ref2015" cutoff="4" residue_surface_cutoff="20.0" ignore_surface_res="true" print_out_info_to_pdb="true" confidence="0" />

        '''

        buns_sc = '''
        <BuriedUnsatHbonds name="Buried Unsat Sidechains [[-]]"
        report_sc_heavy_atom_unsats="true" scorefxn="ref2015" cutoff="4" residue_surface_cutoff="20.0" ignore_surface_res="true" print_out_info_to_pdb="true" confidence="0" />

        '''

        packstat = '''
          <PackStat
            name="PackStat Score [[+]]"
            threshold="0"
          />
        '''

        prepro = '''
          <PreProline
            name="Pre-Proline Potential [[-]]"
            use_statistical_potential="true"
          />
        '''

        exposed_hydrophobics = '''
          <ExposedHydrophobics
            name="Exposed Hydrophobic Residue SASA [[-]]"
            sasa_cutoff="20"
            threshold="-1"
          />
          '''

        filters = [buns_all, packstat, prepro, exposed_hydrophobics]
        # filters = [packstat, prepro, exposed_hydrophobics]
        if os.path.exists(os.path.join(
            self.workspace.rosetta_dir,
            'source',
            'external',
            'DAlpahBall',
            'DAlphaBall.gcc'
                )
            ):
            filters.append(buns_all)
            filters.append(buns_sc)
        filter_objs = []
        for filt in filters:
            filter_objs.append(XmlObjects.static_get_filter(filt))
        if score_fragments:

            fsf = '''
              <FragmentScoreFilter
                name="Avg. 9-Residue Fragment Crmsd[[-]]"
                scoretype="FragmentCrmsd"
                sort_by="FragmentCrmsd"
                threshold="9999" 
                direction="-"
                start_res="{largest_loop_start}"
                end_res="{largest_loop_end}"
                compute="average"
                outputs_folder="{seqprof_dir}"
                outputs_name="{task_id}" 
                csblast="/wynton/home/kortemme/krivacic/software/fragments/csblast-2.2.3_linux64"  
                blast_pgp="/wynton/home/kortemme/krivacic/software/fragments/blast/blast-2.2.26/bin/blastpgp" 
                psipred="/wynton/home/kortemme/krivacic/software/fragments/psipred/runpsipred_single" 
                sparks-x="/wynton/home/kortemme/krivacic/software/fragments/sparks-x" 
                sparks-x_query="/wynton/home/kortemme/krivacic/software/fragments/sparks-x/bin/buildinp_query.sh" 
                frags_scoring_config="{fragment_weights_path}"
                placeholder_seqs="/wynton/home/database/blast/blastdb/pdbaa"
                print_to_pdb="true"
                n_frags="200"
                n_candidates="1000" 
                fragment_size="9"
                vall_path="{vall_path}"
              />
            '''.format(largest_loop_start=self.workspace.largest_loop.start,
                    largest_loop_end=self.workspace.largest_loop.end,
                    seqprof_dir=self.workspace.seqprof_dir, task_id=self.task_id,
                    fragment_weights_path=self.workspace.fragment_weights_path,
                    vall_path=self.workspace.rosetta_vall_path(test_run))
                #placeholder_seqs="/wynton/home/kortemme/krivacic/software/fragments/derived_data/pdb_seqres.txt"
            filter_objs.append(XmlObjects.static_get_filter(fsf))

        return filter_objs

    #placeholder_seqs="/netapp/home/xingjiepan/Databases/BLAST/placeholder/placeholder_seqs" 

    def run_filter(self, f):
        if f.name() == 'FragmentScoreFilter':
            # FragmentScoreFilter can be finicky and sometimes requires you
            # to run sparks-x outside of C++. I have no idea why, just roll
            # with it.
            try:
                # First try to run FSF because it needs to generate the
                # other sequence profile files
                score = f.report_sm(self.pose)
            except:
                import subprocess
                # Now get rid of the empty .fasta.pssm file it creates
                # and try again in Python
                rempath = os.path.join(self.workspace.seqprof_dir,
                        '{}.fasta.phipsi'.format(self.task_id))
                print('REMOVING {}'.format(rempath))
                os.remove(rempath)
                print('CWD: {}'.format(os.getcwd()))
                print('TASK: {}'.format(self.task_id))
                print('FILES IN CWD:')
                print(os.listdir(os.getcwd()))
                cmd = [
                        '/wynton/home/kortemme/krivacic/software/fragments/sparks-x/bin/buildinp_query.sh',
                        '{}.fasta'.format(self.task_id),
                        ]

                process = subprocess.run(cmd,
                        env=dict(SPARKSXDIR='/wynton/home/kortemme/krivacic/software/fragments/sparks-x',
                            **os.environ))
                score = f.report_sm(self.pose)

        elif f.name() == 'BuriedUnsatHbondFilter':
            try:
                score = f.report_sm(self.pose)
            except:
                print('Buried Unsats FAILED; make sure DAlphaBall is '\
                        'installed and working. You may need to install '\
                        'gfortran \n(conda install -c anaconda '\
                        'gfortran_<OS>-64 gmp)')
                print('Once installed, add the library path to your '\
                        'environment.\n'\
                        'Ex. export '\
                        'LD_LIBRARY_PATH=anaconda/lib:LD_LIBRARY_PATH')
        else:
            score = f.report_sm(self.pose)

        fname = 'EXTRA_METRIC_' + f.get_user_defined_name()
        setPoseExtraScore(self.pose, fname, score)

    def run_filters(self):
        for f in self.filters:
            self.run_filter(f)
