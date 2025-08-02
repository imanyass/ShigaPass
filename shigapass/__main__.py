import sys
import logging

from .__version__ import VERSION
from .cli import get_args
from .main import OptionsParser


def print_help():
    """Print help output"""
    print(f"""
                    ...::: ShigaPass v{VERSION} :::...
Author(s): Iman Yassine & Ryan Kennedy

Description: This tool is used to predict Shigella serotypes.

Usage: shigapass [options]
Options:
    -l, --list_file         List of input file(s) (FASTA) with their path(s) (mandatory)
    -o, --outdir       Output directory (mandatory)
    -p, --db           Path to databases directory (mandatory)
    -t, --threads      Number of threads (optional, default: 2)
    -u, --mkdb         Call the makeblastdb utility for databases initialisation (optional, but required when running the script for the first time)
    -k, --keep         Do not remove subdirectories (optional)
    -v, --version      Display the version and exit
    -h, --help         Display this help and exit

Example: shigapass -l list_of_fasta.txt -o ShigaPass_Results -p ShigaPass/ShigaPass_DataBases -t 4 -u -k

NOTE: The -u option should be used when running the script for the first time and after databases updates.
""")

def main():
    """Main function for initiating software"""
    if len(sys.argv) == 1:
        print_help()
        sys.exit(0)
    elif sys.argv[1] in {"-V", "--version"}:
        print(f"ShigaPass version {VERSION}")
        sys.exit(0)
    elif sys.argv[1] in {"-h", "--h", "-help", "--help"}:
        print_help()
        sys.exit(0)

    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s] %(levelname)s in %(module)s: %(message)s"
    )
    logger = logging.getLogger(__name__)

    try:
        logger.info("Starting analysis")
        args = get_args().parse_args()
        sp_parser = OptionsParser(VERSION)
        sp_parser.parse_options(args)
        logger.info("Analysis completed successfully")
    except SystemExit:
        logger.warning("Controlled exit resulting from early termination.")
        sys.exit(1)
    except KeyboardInterrupt:
        logger.warning("Controlled exit resulting from interrupt signal.")
        sys.exit(1)
    except Exception as e:
        error_message = "Uncontrolled exit resulting from an unexpected error.\n\n"
        error_message += "-" * 80 + "\n"
        error_message += f"EXCEPTION: {type(e).__name__}\n"
        error_message += f"MESSAGE: {e}\n"
        error_message += "-" * 80 + "\n\n"
        logger.error(error_message)
        sys.exit(1)

if __name__ == "__main__":
    main()
