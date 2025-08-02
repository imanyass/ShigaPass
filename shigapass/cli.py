"""Command line interface setup"""

import argparse

from .__version__ import VERSION

def get_args():
    parser = argparse.ArgumentParser(description="Tool for predicting Shigella serotypes in silico")
    parser.add_argument("-l", "--list_file", dest="list_file", required=True, help="List of input file(s) (FASTA) with their path(s)")
    parser.add_argument("-o", "--outdir", dest="outdir", required=True, help="Output directory")
    parser.add_argument("-d", "--db", dest="db", required=True, help="Path to databases directory")
    parser.add_argument("-t", "--threads", dest="threads", type=int, default=2, help="Number of threads (default: 2)")
    parser.add_argument("-u", "--mkdb", dest="mkdb", action="store_true", help="Call the makeblastdb utility for databases initialisation")
    parser.add_argument("-k", "--keep", dest="keep_files", action="store_true", help="Do not remove subdirectories")
    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {VERSION}", help="Display version")
    return parser
