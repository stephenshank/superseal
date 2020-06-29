"Command line interface."
import os
import sys
import argparse
from http.server import BaseHTTPRequestHandler
from http.server import HTTPServer
import webbrowser


from .reads import covarying_sites_io
from .reads import superread_json_io
from .reads import superread_fasta_io
from .reads import resolvable_regions_io

from .assembly import assemble_io
from .assembly import local_reconstruction_io


DESCRIPTION = """
superseal - SUPERread Seed Expansion, Absorption, and Localization

Written by Stephen D. Shank, Ph. D.
Acme Computational Molecular Evolution Group - http://lab.hyphy.org/
https://github.com/stephenshank/superseal

Further help:
    superseal covariation --help
    superseal superreads --help
    superseal resolve --help
    superseal assemble --help
    superseal localreconstruct --help
    superseal view --help
"""


class Server(BaseHTTPRequestHandler):
    "Basic HTTP server for CLI."
    # pylint: disable=invalid-name
    def do_GET(self):
        "Handle GET request."
        print(os.getcwd())

        base_dir = os.path.dirname(__file__)
        static_resources = [
            'index.html',
            'main.js',
            'style.css',
            'favicon.ico',
            'logo.png'
        ]
        if self.path == '/':
            path = os.path.join(base_dir, 'viz', 'index.html')
            header = "text/html"
        elif self.path[1:] in static_resources:
            path = os.path.join(base_dir, 'viz', self.path[1:])
            if self.path[1:] == 'style.css':
                header = "text/css"
            elif self.path[-3:] == 'png':
                header = "image/png"
            else:
                header = "text/html"
        else:
            path = os.path.join(os.getcwd(), 'superreads.json')
            header = "text/json"

        self.send_response(200)
        self.send_header("Content-type", header)
        self.end_headers()
        with open(path, 'rb') as f:
            self.wfile.write(f.read())


def covarying_sites_cli(args):
    "Command line interface for covarying sites."
    covarying_sites_io(
        args.bam, args.sites, args.fasta, args.csv, args.threshold
    )


def superreads_cli(args):
    "Command line interface for superreads."
    superread_json_io(args.bam, args.cov, args.sr)
    if args.fasta:
        superread_fasta_io(args.cov, args.sr, args.fasta)


def resolve_cli(args):
    "Command line interface for resolving regions."
    resolvable_regions_io(args.sr, args.rr)


def assemble_cli(args):
    "Command line interface for performing assembly."
    assemble_io(args.sr, args.rr, args.assembly, max_qs=args.mq)


def local_reconstruction_cli(args):
    "Command line interface for local reconstruction."
    local_reconstruction_io(
        args.sr, args.assem, args.con, args.var, args.reg, args.fasta
    )


def view_cli():
    "Command line interface for interactive result viewer."
    server = HTTPServer(("localhost", 8749), Server)

    try:
        webbrowser.open("http://localhost:8749", autoraise=True)
        server.serve_forever()
    except KeyboardInterrupt:
        pass

    server.server_close()


def command_line_interface():
    "Full command line interface function."
    if len(sys.argv) == 1:
        print(DESCRIPTION)
        print("Usage: superseal -h or superseal --help")
        sys.exit(0)

    main_parser = argparse.ArgumentParser(
        description=DESCRIPTION,
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
        default=.01,
        type=float
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

    as_description = "Assemble regions."
    as_subparser = subparsers.add_parser(
        "assemble", description=as_description
    )
    as_subparser.set_defaults(func=assemble_cli)
    as_subparser.add_argument(
        "-s", "--superreads", help="input superreads", dest="sr", required=True
    )
    as_subparser.add_argument(
        "-r",
        "--resolvable",
        help="resolvable regions",
        dest="rr",
        required=True
    )
    as_subparser.add_argument(
        "-a",
        "--assembly",
        help="assembled superreads",
        dest="assembly",
        required=True
    )
    as_subparser.add_argument(
        "-m",
        "--maximum_quasispecies",
        help="maximum number of quasispecies",
        dest="mq",
        default=2,
        type=int
    )

    lr_description = "Perform local reconstruction."
    lr_subparser = subparsers.add_parser(
        "localreconstruct", description=lr_description
    )
    lr_subparser.set_defaults(func=local_reconstruction_cli)
    lr_subparser.add_argument(
        "-s", "--superreads", help="input superreads", dest="sr", required=True
    )
    lr_subparser.add_argument(
        "-a",
        "--assembly",
        help="assembled superreads",
        dest="assem",
        required=True
    )
    lr_subparser.add_argument(
        "-c",
        "--consensus",
        help="consensus sequence",
        dest="con",
        required=True
    )
    lr_subparser.add_argument(
        "-v",
        "--variation",
        help="covarying sites",
        dest="var",
        required=True
    )
    lr_subparser.add_argument(
        "-r",
        "--region",
        help="region to resolve",
        dest="reg",
        required=True,
        type=int
    )
    lr_subparser.add_argument(
        "-f",
        "--f",
        help="quasispecies fasta",
        dest="fasta",
        required=True
    )

    view_description = "Visualize data."
    view_subparser = subparsers.add_parser(
        "view", description=view_description
    )
    view_subparser.set_defaults(func=view_cli)

    args = main_parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    command_line_interface()
