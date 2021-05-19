#!/usr/bin/env python2

"""\
Create a web logos for sequences generated by the design pipeline.

Usage:
    pull_into_place web_logo <workspace> <round> [options]
    pull_into_place web_logo <directory> [options]

Options:
    --output PATH, -o PATH
        The path where the logo should be saved.  The desired file format is 
        deduced from the extension, with '.pdf', '.svg', '.eps', '.png', 
        '.jpeg', and '.txt' all supported.  By default, the logo is not saved 
        and is instead opened using the program specified in the $PDF 
        environment variable (or evince if no such variable is defined).

If you specify a workspace and a round number, the weblogo will be calculated 
for all the validated designs in that round.  If you specify a directory, the 
weblogo will just be calculated for the models in that directory.

It would be nice to pass all unparsed options through to weblogo.  I'll have to 
think a bit about how to do that.
"""

import os, tempfile, subprocess
# import weblogolib as weblogo, corebio
import weblogo

from klab import scripting
import docopt
from .. import pipeline, structures

@scripting.catch_and_print_errors()
def main():
    args = docopt.docopt(__doc__)
    root = args['<workspace>']
    round = args['<round>']
    directory = args['<directory>']

    if directory:
        models, filters = structures.load(directory)
        resfile = pipeline.load_resfile(directory)
        resis = sorted(int(i) for i in resfile.designable)
        print(resis)
        title = directory
        sequences = [
                ''.join(models['sequence'][i][j-1] for j in resis)
                for i in models.index
        ]

    else:
        workspace = pipeline.ValidatedDesigns(root, round)
        workspace.check_paths()
        title = workspace.focus_dir
        designs = [structures.Design(x) for x in workspace.output_subdirs]
        sequences = [x.resfile_sequence for x in designs]

    sequences = weblogo.seq.SeqList(
            [weblogo.seq.Seq(x) for x in sequences],
            alphabet=weblogo.seq.unambiguous_protein_alphabet,
    )

    logo_data = weblogo.LogoData.from_seqs(sequences)
    logo_options = weblogo.LogoOptions()
    logo_options.title = title
    logo_format = weblogo.LogoFormat(logo_data, logo_options)

    if args['--output']:
        preview = False
        logo_file = open(args['--output'], 'wb')
        with open(args['--output'][:-len('pdf')] + 'txt', 'w') as f:
            for line in str(logo_data): 
                f.write(line)
    else:
        preview = True
        logo_file = tempfile.NamedTemporaryFile(
                'wb', prefix='weblogo_', suffix='.pdf')

    ext = os.path.splitext(logo_file.name)[-1]
    formatters = {
            '.pdf': weblogo.pdf_formatter,
            '.svg': weblogo.svg_formatter,
            '.eps': weblogo.eps_formatter,
            '.png': weblogo.png_formatter,
            '.jpeg': weblogo.jpeg_formatter,
            '.txt': weblogo.txt_formatter,
    }
    if ext not in formatters:
        scripting.print_error_and_die("'{0}' is not a supported file format".format(ext))

    document = formatters[ext](logo_data, logo_format)
    logo_file.write(document)
    logo_file.flush()

    if preview:
        pdf = os.environ.get('PDF', 'evince'), logo_file.name
        subprocess.call(pdf)

    logo_file.close()



