import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description="PSI-esg: PSI by Event Splicing Graphs weighting"
    )
    parser.add_argument(
        "FA",
        type=str,
        help="Path to FASTA reference",
    )
    parser.add_argument(
        "GTF",
        type=str,
        help="Path to GTF annotation",
    )
    parser.add_argument(
        "FQs",
        type=str,
        help="Path to FASTQ RNA-Seq sample(s)",
    )
    parser.add_argument(
        "-e",
        "--events",
        dest="events",
        default="SE,SS,RI,MX,FL",
        type=str,
        help="Comma separated list of events (default: SE,SS,RI,MX,FL)",
    )
    parser.add_argument(
        "--galig",
        dest="galig",
        default=".",
        type=str,
        help="Path to galig folder (default: .)",
    )
    parser.add_argument(
        "--chroms",
        dest="chroms",
        default=None,
        type=str,
        help="Path to chromosome folder (default: None)",
    )
    parser.add_argument(
        "-l",
        dest="l",
        default=10,
        type=int,
        help="MEMs length (default: 10)",
    )
    parser.add_argument(
        "-t",
        "--threads",
        dest="threads",
        default=1,
        type=int,
        help="Number of threads (default: 1)",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        default=False,
        help="Verbose mode",
    )
    parser.add_argument(
        "--version",
        dest="version",
        action="store_true",
        default=False,
        help="Print version",
    )
    parser.add_argument(
        "--wd",
        dest="wd",
        default="psiesg-wd",
        type=str,
        help="Output directory (default: eflux-wd)",
    )

    return parser.parse_args()
