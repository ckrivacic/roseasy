#!/usr/bin/env python2

"""\

Simple wrapper around a klab script to make it more user-friendly (for me anyway).

Generate a well-behaved DNA sequence from the given protein sequence.  If a
template DNA sequence is specified, the returned DNA sequence will be as
similar to it as possible.  All of the given restriction sites will be removed
from the generated sequence.  And finally, the given leading and trailing
sequences will be appropriately concatenated.

If you provide the --from-pdb-folder flag, you can also point it to a
folder of pdb files and it will get the design sequences from there. 

Usage:
    roseasy reverse_translate [options] <input_fasta_or_folder> <output_file>

Options:
    --template-dna DNA
        A (presumably natural) DNA sequence that should be used as a template
        for the reverse translation process.  In other words, when applicable
        codons from this template are used in favor of "optimized" codons.

    --leading-dna DNA
        A fixed DNA sequence to prepend to the generated sequence.

    --trailing-dna DNA
        A fixed DNA sequence to append to the generated sequence.

    --restriction-sites SITES
        Comma-separated list of restriction sites that may not appear in the
        generated sequence.  You may specify either case-sensitive restriction
        site names (must be listed in cloning.py) or DNA sequences.

    --using-gen9
        Indicate that this gene is going to be ordered from Gen9.  This enables
        a number of checks and filters that are specific to Gen9.

    --make-part3
        Adds restriction sites for part 3s at either end of the DNA sequence.

    --make-part3a
        Adds restriction sites/overhangs for part 3a at either end of the DNA
        sequence.

    --from-pdb-folder
        Get design protein sequences from a folder of pdb files.

    --print-protein-seq
        print the protein sequence before writing dna seq

    --combine-chains FASTA
        If you have mutations spread across two identical chains, use
        this flag to combine them into one amino acid sequence. Provide
        a FASTA file after this flag for use as a template amino acid
        sequence.

    --reversion-mutations EXCEL_FILE

    --first-res=<offset>
        If your protein sequence does not start at residue 1, enter the
        number of the first residue here.

    --chain=CHAR
        Only import a certain chain. Only works for PDB folder input.


"""

import docopt
from roseasy import pipeline
from klab import cloning
from copy import deepcopy
import re, os
import pandas as pd
from Bio import pairwise2

def import_wt_protein_sequence(wt_sequence_file):
    wt_sequence_str = ''
    with open(wt_sequence_file,'r') as file:
        carot_counter = 0
        for line in file:
            if line.startswith('>'):
                carot_counter += 1
                if carot_counter > 1:
                    sys.exit("Aborted: FASTA file contained too many\
 sequences. Please provide a FASTA file with only WT protein sequence.")
            else:
                wt_sequence_str += line.rstrip()
    return wt_sequence_str

def combine_chains(chain1, chain2, wtchain):
    final_seq = deepcopy(wtchain)
    for index,aa in enumerate(chain1):
        if aa != wtchain[index]:
            final_seq[index] = aa
    for index,aa in enumerate(chain2):
        if aa != wtchain[index]:
            final_seq[index] = aa

    # final_seq[37] = 'E'
    # final_seq[82] = 'T'

    return final_seq
        
def import_protein_fasta(inputs):
    protein_sequences = {}
    with open(inputs) as f:
        current_sequence = ''
        for line in f:
            if line[0] == '>':
                current_sequence = line[0:].rstrip()
            else:
                protein_sequences[current_sequence] += line.rstrip()

    return protein_sequences

def import_protein_structure(inputs, wt_protein_fasta_file, chain):
    
    from Bio import PDB
    import gzip
    file_list = []
    for subdir, dirs, files in os.walk(inputs):
        for file in files:
            if file.endswith('.pdb') or file.endswith('.pdb.gz'):
                file_list.append(os.path.join(subdir, file))
    
    protein_sequences = {}

    for file in file_list:
        subdir = file.split('/')[-2]
        name = subdir + os.path.basename(file).split('.')[0]

        parser = PDB.PDBParser()
        if file.endswith('.gz'):
            pdb = gzip.open(file,'rt')
        else:
            pdb = open(file,'r')
        io = PDB.PDBIO
        struct = parser.get_structure(name,pdb)
        ppb = PDB.PPBuilder()
        
        chains = []
        for pp in ppb.build_peptides(struct):
            chains.append(list(pp.get_sequence()))

        if wt_protein_fasta_file:
            print('Template being used')
            wt_seq_list = \
            list(import_wt_protein_sequence(wt_protein_fasta_file))
            final_design_seq = combine_chains(chains[0],chains[1],wt_seq_list)
            protein_sequences[name] = "".join(final_design_seq)
        else:
            if len(chains) > 1:
                print("Warning: Multiple chains found. Splitting "\
"sequence into ", len(chains), " DNA sequences for ordering.")
            for index,chain in enumerate(chains):
                protein_sequences[name + "_chain_" + str(index)] = "".join(chain)

    return protein_sequences


def pairwise_align(wt_aa, mut_aa):
    '''Align wt and mutant AA strings and return modified (aligned)
    strings'''
    from Bio.Seq import Seq
    wt_aa = Seq(wt_aa)
    mut_aa = Seq(mut_aa)
    alignments = pairwise2.align.globalxd(wt_aa, mut_aa, -100, -100, -5,
            -1)
    alignment = alignments[0]
    wt = alignment[0]
    mut = alignment[1]
    return wt, mut


def make_reversion_mutations(reversion_file,protein_sequences):
    import difflib
  
    if arguments['--first-res']:
        first_res = int(arguments['--first-res'])
    else:
        first_res = 1

    reversions = pd.read_excel(reversion_file)
    reversions = reversions.fillna(value='')
    aa_pattern = re.compile(r'([A-Z]|[a-z]){1,1}[0-9]{1,4}([A-Z]|[a-z]){1,1}')

    reverted_sequences = {}

    for reversion in reversions:
        original_design = \
        difflib.get_close_matches(reversion, [design for
            design in protein_sequences])[0]
        print("Reversion mutations: matched ", reversion, " with ",
                original_design)
        reversion_seq = list(deepcopy(protein_sequences[original_design]))
        for mutation in reversions[reversion]:
            if aa_pattern.match(mutation):
                mut_split = re.split(r'([0-9]{1,4})', mutation)
                wt_aa = mut_split[0]
                aa_num = int(mut_split[1])
                mut_aa = mut_split[2]

                while aa_num > len(reversion_seq):
                    aa_num -=len(reversion_seq)

                assert reversion_seq[aa_num - first_res] == mut_aa,\
                "Your reversion mutation did not match the design \
sequence. Check that (a) you don't have any duplicates in your mutation \
list, (b) that your WT sequence is correct, and that (c) the numbering \
of your protein is correctly set in the options. \n Reversion aa: \
" + reversion_seq[aa_num - first_res] + "\n Design aa: " + mut_aa + \
"\n aa number: " + str(aa_num)

                reversion_seq[aa_num - first_res] = wt_aa

        name = original_design + "_reversion"
        num_reversion_designs = 1
        while name in reverted_sequences:
            num_reversion_designs += 1
            name = name + str(num_reversion_designs)
        reverted_sequences[name] = "".join(reversion_seq)
    
    for design in reverted_sequences:
        protein_sequences[design] = reverted_sequences[design]

    return protein_sequences

def main():

    arguments = docopt.docopt(__doc__)

    inputs = arguments['<input_fasta_or_folder>']

    workspace = \
            pipeline.workspace_from_dir(arguments['<input_fasta_or_folder>'])
    input_pdb = workspace.input_pdb_path
    import prody
    atoms = prody.parsePDB(input_pdb).select('chain A and name CA')
    wt_aa_sequence = atoms.getSequence()

    template_dna = arguments['--template-dna']
    wt_sequence_file = arguments['--combine-chains']
    wt_sequence_str = ''

    reversion_mutations = arguments['--reversion-mutations']

    if wt_sequence_file:
        wt_sequence_str = import_wt_protein_sequence(wt_sequence_file)

    from_files = arguments['--from-pdb-folder']

    if from_files:
        protein_sequences = \
                import_protein_structure(inputs,arguments['--combine-chains'],
                arguments['--chain'])

    else:
        protein_sequences = import_protein_fasta(inputs)

    if reversion_mutations:
        protein_sequences = make_reversion_mutations(reversion_mutations, protein_sequences)

    #protein_seq = arguments['<protein_seq>']
    leading_dna = arguments['--leading-dna'] or ""
    trailing_dna = arguments['--trailing-dna'] or ""
    part3 = arguments['--make-part3']
    part3a = arguments['--make-part3a']

    if part3:
        leading_dna = 'gcatCGTCTCaGGTCTCaT'
        trailing_dna = 'ggATCCTGAGACCTGAGACGGCAT'

    if part3a:
        leading_dna = 'gcatCGTCTCaGGTCTCaT'
        trailing_dna = 'ggTTCTTGAGACCTGAGACGGCAT'

    try: restriction_sites = arguments['--restriction-sites'].split(',')
    except: restriction_sites = ()
    using_gen9 = arguments['--using-gen9']

    output_sequences = {}
    if template_dna.endswith('.fasta'):
        with open(template_dna, 'r') as f:
            for line in f:
                if not line.startswith('>'):
                    template_dna = line

    for key in protein_sequences:
        protein_seq = protein_sequences[key]
        protein_seq = protein_seq.replace('*','.')
        wt_aligned, protein_aligned = pairwise_align(wt_aa_sequence,
                protein_seq)
        assert('-' not in wt_aligned)
        output_sequences[key] = cloning.reverse_translate(
                protein_aligned,
                template_dna,
                leading_dna,
                trailing_dna,
                restriction_sites)

    output_file = arguments['<output_file>']
    with open(output_file, 'w') as f:
        for key in sorted(output_sequences.items()):
            if not str(key).startswith('>'):
                f.write('>')
            f.write(key[0].replace(":","_") + '\n')
            if arguments['--print-protein-seq']:
                print('PROTEIN SEQ FOR {}'.format(key[0]))
                print(protein_sequences[key[0]])
            f.write(output_sequences[key[0]] + '\n')

