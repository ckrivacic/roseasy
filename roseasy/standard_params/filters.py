from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects

def get_filters(workspace, job_id='0000', score_fragments=False,
        test_run=False):
    writer = '''
    <WriteFiltersToPose name="writer" prefix="EXTRA_SCORE_"/>
    '''

    buns = '''
    <BuriedUnsatHbonds name="Buried Unsat [[-]]"
    report_all_heavy_atom_unsats="true" scorefxn="ref2015" cutoff="4" residue_surface_cutoff="20.0" ignore_surface_res="true" print_out_info_to_pdb="true" dalphaball_sasa="1" probe_radius="1.1" confidence="0" />

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

    fsf = '''
      <FragmentScoreFilter
        name="Max 9-Residue Fragment Crmsd[[-]]"
        scoretype="FragmentCrmsd"
        sort_by="FragmentCrmsd"
        threshold="9999" 
        direction="-"
        start_res="{largest_loop_start}"
        end_res="{largest_loop_end}"
        compute="maximum"
        outputs_folder="{seqprof_dir}"
        outputs_name="{job_id}" 
        csblast="/wynton/home/kortemme/krivacic/software/fragments/csblast-2.2.3_linux64"  
        blast_pgp="/wynton/home/kortemme/krivacic/software/fragments/blast/blast-2.2.26/bin/blastpgp" 
        placeholder_seqs="/wynton/home/kortemme/krivacic/software/fragments/derived_data/pdb_seqres.txt"
        psipred="/wynton/home/kortemme/krivacic/software/fragments/psipred/runpsipred_single" 
        sparks-x="/wynton/home/kortemme/krivacic/software/fragments/sparks-x" 
        sparks-x_query="/wynton/home/kortemme/krivacic/software/fragments/sparks-x/bin/buildinp_query.sh" 
        frags_scoring_config="{fragment_weights_path}"
        n_frags="200"
        n_candidates="1000" 
        fragment_size="9"
        vall_path="{vall_path}"
      />
    '''.format(largest_loop_start=workspace.largest_loop.start,
            largest_loop_end=workspace.largest_loop.end,
            seqprof_dir=workspace.seqprof_dir, job_id=job_id,
            fragment_weights_path=workspace.fragment_weights_path,
            vall_path=workspace.rosetta_vall_path(test_run))

    filters = [buns, packstat, prepro, exposed_hydrophobics]
    filter_objs = []
    for filt in filters:
        filter_objs.append(XmlObjects.static_get_filter(filt))
    if score_fragments:
        filter_objs.append(XmlObjects.static_get_filter(fsf))

    return filter_objs

#placeholder_seqs="/netapp/home/xingjiepan/Databases/BLAST/placeholder/placeholder_seqs" 
