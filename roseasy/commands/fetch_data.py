
"""\
Copy design files from a remote source.  A common application is to copy 
simulation results from the cluster to a workstation for analysis.  The given 
directory must be contained within a workspace created by 01_setup_workspace.

Usage:
    roseasy fetch_data <directory> [options]

Options:
    --remote URL, -r URL
        Specify the URL to fetch data from.  You can put this value in a file 
        called "rsync_url" in the local workspace if you don't want to specify 
        it on the command-line every time.

    --include-logs, -i
        Fetch log files (i.e. stdout and stderr) in addition to everything 
        else.  Note that these files are often quite large, so this may take 
        significantly longer.

    --dry-run, -d
        Output the rsync command that would be used to fetch data.
        
"""

from klab import scripting
import docopt
from roseasy.workspace import pipeline

@scripting.catch_and_print_errors()
def main():
    args = docopt.docopt(__doc__)

    pipeline.fetch_data(
            args['<directory>'],
            remote_url=args['--remote'],
            include_logs=args['--include-logs'],
            dry_run=args['--dry-run'],
    )


if __name__=='__main__':
    main()
