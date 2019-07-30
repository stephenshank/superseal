import os
import argparse

from .io import error_correction_io
from .io import read_graph_io
from .io import regression_io
from .viz import show_viz


ec_description = '''
    Correct errors in aligned NGS reads in a quasi-species aware manner.
'''


def error_correction():
    parser = argparse.ArgumentParser(description=ec_description)
    parser.add_argument(
        '-b', '--bam',
        metavar="BAM",
        type=str,
        help="BAM file of mapped, sorted reads",
        required=True
    )
    parser.add_argument(
        '-o', '--output',
        metavar="OUTPUT",
        type=str,
        help="Directory to output all files",
        required=False,
        default=None
    )
    parser.add_argument(
        '-c', '--corrected',
        metavar="CORRECTED",
        type=str,
        help="BAM file of error-corrected reads",
        required=False,
        default=None
    )
    parser.add_argument(
        '-j', '--json',
        metavar="JSON",
        type=str,
        help="JSON file of mapped, sorted reads",
        required=False,
        default=None
    )
    parser.add_argument(
        '-f', '--fasta',
        metavar="FASTA",
        type=str,
        help="FASTA file of consensus sequence (used downstream)",
        required=False,
        default=None
    )
    parser.add_argument(
        '-e', '--end-correction',
        metavar="ENDCORRECTION",
        type=int,
        help="Correct bases at ends of reference to consensus",
        required=False,
        default=None
    )

    args = parser.parse_args()

    if args.output:
        corrected = os.path.join(args.output, "corrected.bam")
        json = os.path.join(args.output, "covarying.json")
        fasta = os.path.join(args.output, "consensus.fasta")
    else:
        corrected = args.corrected
        json = args.json
        fasta = args.fasta

    error_correction_io(
        args.bam, corrected, json, fasta, end_correction=args.end_correction
    )


rg_description = '''
    Build a read graph from a BAM file of error corrected reads.
'''


def read_graph():
    parser = argparse.ArgumentParser(
        description=rg_description
    )

    parser.add_argument(
        '-i', '--input',
        metavar="INPUT",
        type=str,
        help="Input directory, containing results of error correction step",
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
    parser.add_argument(
        '-b', '--bam',
        metavar="BAM",
        type=str,
        help="BAM file of corrected, mapped reads",
        required=False,
        default=None
    )
    parser.add_argument(
        '-s', '--sites',
        metavar="SITES",
        type=str,
        help="JSON file of covarying sites",
        required=False,
        default=None
    )
    parser.add_argument(
        '-c', '--consensus',
        metavar="CONSENSUS",
        type=str,
        help="FASTA file of consensus sequence (full superreads/candidates)",
        required=False,
        default=None
    )
    parser.add_argument(
        '-f', '--full',
        metavar="FULL",
        type=str,
        help="Output fasta file of full superreads",
        required=False,
        default=None
    )
    parser.add_argument(
        '-r', '--restricted',
        metavar="RESTRICTED",
        type=str,
        help="Output file of superreads restricted to covarying sites",
        required=False,
        default=None
    )
    parser.add_argument(
        '-d', '--describing-superreads',
        metavar="DESCRIBING",
        type=str,
        help="Output JSON file of describing superreads",
        required=False,
        default=None
    )
    parser.add_argument(
        '-g', '--graph',
        metavar="GRAPH",
        type=str,
        help="Output JSON file of superread graph",
        required=False,
        default=None
    )
    parser.add_argument(
        '-C', '--candidates',
        metavar="CANDIDATES",
        type=str,
        help="Output FASTA file of candidate quasispecies",
        required=False,
        default=None
    )
    parser.add_argument(
        '-m', '--minimum-weight',
        metavar="MINIMUMWEIGHT",
        type=int,
        help="Minimum weight to discard from superreads",
        required=False,
        default=3
    )

    args = parser.parse_args()

    if args.input:
        bam = os.path.join(args.input, "corrected.bam")
        sites = os.path.join(args.input, "covarying.json")
        consensus = os.path.join(args.input, "consensus.fasta")
    else:
        bam = args.bam
        sites = args.sites
        consensus = args.consensus

    if args.output:
        full = os.path.join(args.output, "full.fasta")
        restricted = os.path.join(args.output, "restricted.fasta")
        describing = os.path.join(args.output, "describing.json")
        graph = os.path.join(args.output, "graph.json")
        candidates = os.path.join(args.output, "candidates.fasta")
    else:
        full = args.full
        restricted = args.restricted
        describing = args.describing_superreads
        graph = args.graph
        candidates = args.candidates

    read_graph_io(
        bam, sites, consensus,
        full, restricted, describing, graph, candidates,
        args.minimum_weight
    )


qr_description = '''
    Reconstruct quasispecies from read graph information.
'''


def quasispecies_reconstruction():
    parser = argparse.ArgumentParser(
        description=qr_description
    )

    parser.add_argument(
        '-i', '--input',
        metavar="INPUT",
        type=str,
        help="Input directory, containing results of error correction step",
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
    parser.add_argument(
        '-g', '--graph',
        metavar="GRAPH",
        type=str,
        help="JSON file describing the superread graph",
        required=False,
        default=None
    )
    parser.add_argument(
        '-d', '--describing-superreads',
        metavar="DESCRIBING",
        type=str,
        help="Input JSON file of describing superreads",
        required=False,
        default=None
    )
    parser.add_argument(
        '-C', '--candidates',
        metavar="CANDIDATES",
        type=str,
        help="Input FASTA file of candidate quasispecies",
        required=False,
        default=None
    )
    parser.add_argument(
        '-q', '--quasispecies',
        metavar="QUASISPECIES",
        type=str,
        help="Output FASTA file of reconstructed quasispecies",
        required=False,
        default=None
    )

    args = parser.parse_args()

    if args.input:
        describing = os.path.join(args.output, "describing.json")
        graph = os.path.join(args.output, "graph.json")
        candidates = os.path.join(args.output, "candidates.fasta")
    else:
        describing = args.describing_superreads
        graph = args.graph
        candidates = args.candidates

    if args.output:
        quasispecies = os.path.join(args.output, "quasispecies.fasta")
    else:
        quasispecies = args.quasispecies

    regression_io(
        graph, describing, candidates, quasispecies
    )


full_pipeline_description = '''
    Run the entire quasispecies reconstruction pipeline.
'''


def full_pipeline():
    parser = argparse.ArgumentParser(
        description=full_pipeline_description
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

    print('TODO: Implement full pipeline for:', args.bam, args.fasta)


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
