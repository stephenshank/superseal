import sys
import argparse

from .reads import covarying_sites_io
from .reads import superread_json_io 
from .reads import superread_fasta_io 
from .reads import resolvable_regions_io


description="""
superseal - SUPERread Seed Expansion, Absorption, and Localization

Written by Stephen D. Shank, Ph. D.
Acme Computational Molecular Evolution Group - http://lab.hyphy.org/
https://github.com/stephenshank/superseal

Further help:
    superseal covariation --help
    superseal superreads --help
    superseal resolve --help
"""


def covarying_sites_cli(args):
    covarying_sites_io(
        args.bam, args.sites, args.fasta, args.csv, args.threshold
    )


def superreads_cli(args):
    superread_json_io(args.bam, args.cov, args.sr)
    if args.fasta:
        superread_fasta_io(args.cov, args.sr, args.fasta)


def resolve_cli(args):
    resolvable_regions_io(args.sr, args.rr)


def command_line_interface():
    if len(sys.argv) == 1:
        print(description)
        print("Usage: superseal -h or superseal --help")
        sys.exit(0)

    main_parser = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.RawTextHelpFormatter
    )
    subparsers = main_parser.add_subparsers()


    cvs_description = "Obtain covarying sites from mapped reads."
    cvs_subparser = subparsers.add_parser(
        "covariation", description=cvs_description
    )
    cvs_subparser.set_defaults(func=covarying_sites_cli)
    cvs_subparser.add_argument(
        "-b", "--bam", help="input BAM file", dest="bam", required=True
    )
    cvs_subparser.add_argument(
        "-s", "--sites", help="covarying site JSON", dest="sites", required=True
    )
    cvs_subparser.add_argument(
        "-f", "--fasta", help="consensus FASTA file", dest="fasta"
    )
    cvs_subparser.add_argument(
        "-c", "--csv", help="covariable site CSV", dest="csv"
    )
    cvs_subparser.add_argument(
        "-t",
        "--threshold",
        help="threshold for discerning error from variation",
        dest="threshold",
        default=.01
    )


    sr_description = "Compress mapped reads into superreads."
    sr_subparser = subparsers.add_parser(
        "superreads", description=sr_description
    )
    sr_subparser.set_defaults(func=superreads_cli)
    sr_subparser.add_argument(
        "-b", "--bam", help="input BAM file", dest="bam", required=True
    )
    sr_subparser.add_argument(
        "-c",
        "--covariation",
        help="input covariation",
        dest="cov",
        required=True
    )
    sr_subparser.add_argument(
        "-s", "--superreads", help="superread json", dest="sr", required=True
    )
    sr_subparser.add_argument(
        "-f",
        "--fasta",
        help="superread FASTA file",
        dest='fasta',
        default=None
    )


    rr_description = "Determine resolvable regions from superreads."
    rr_subparser = subparsers.add_parser(
        "resolve", description=rr_description
    )
    rr_subparser.set_defaults(func=resolve_cli)
    rr_subparser.add_argument(
        "-s", "--superreads", help="input superreads", dest="sr", required=True
    )
    rr_subparser.add_argument(
        "-r",
        "--resolvable",
        help="resolvable regions",
        dest="rr",
        required=True
    )

    args = main_parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    command_line_interface()
