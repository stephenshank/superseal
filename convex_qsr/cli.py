import os
import argparse

from .io import full_pipeline_io
from .viz import show_viz


def full_pipeline():
    parser = argparse.ArgumentParser(
        description='''
            Run the entire quasispecies reconstruction pipeline.
        '''
    )

    parser.add_argument(
        '-b', '--bam',
        metavar="BAM",
        type=str,
        help="Input BAM file",
        required=False,
        default=None
    )
    parser.add_argument(
        '-f', '--fasta',
        metavar="FASTA",
        type=str,
        help="Output FASTA file of reconstructed quasispecies",
        required=False,
        default=None
    )
    parser.add_argument(
        '-o', '--output',
        metavar="OUTPUT",
        type=str,
        help="Directory to output all files",
        required=False,
        default=None
    )

    args = parser.parse_args()
    full_pipeline_io(args.bam, args.output, args.fasta)


visualization_description = '''
    View an interactive dashboard of various steps of the pipeline.
'''


def visualization():
    parser = argparse.ArgumentParser(
        description=visualization_description
    )

    parser.add_argument(
        '-o', '--output',
        metavar="OUTPUT",
        type=str,
        help="Output directory where pipeline has been run",
        required=True,
        default=None
    )
    parser.add_argument(
        '-H', '--host',
        metavar="HOST",
        type=str,
        help="Host to run webserver on",
        required=False,
        default='0.0.0.0'
    )
    parser.add_argument(
        '-p', '--port',
        metavar="PORT",
        type=int,
        help="Port to run webserver on",
        required=False,
        default=16180
    )
    browser_parser = parser.add_mutually_exclusive_group(required=False)
    browser_parser.add_argument(
        '--browser',
        dest='browser',
        action='store_true',
        help="Open a web browser"
    )
    browser_parser.add_argument(
        '--no-browser',
        dest='browser',
        action='store_false',
        help="Do not open a web browser"
    )
    parser.set_defaults(browser=True)

    args = parser.parse_args()
    show_viz(args.output, args.port, args.host, args.browser)
