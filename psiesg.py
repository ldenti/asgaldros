import sys
import os
import glob
import subprocess
import multiprocessing as mp
import logging

from Bio import SeqIO
import pysam

from rich.console import Console
from rich.logging import RichHandler
from rich.progress import (
    Progress,
    SpinnerColumn,
    TextColumn,
    BarColumn,
    TaskProgressColumn,
)

from psiesg.gtf import open_gtf, extract_exons, run_gffread
from psiesg.events import run_suppa, build_esg
from psiesg import cli, shark, asgal, psi

console = Console(stderr=True)

logging.basicConfig(
    format="%(message)s",
    level=logging.INFO,
    datefmt="[%x %X]",
    handlers=[RichHandler(console=console)],
)


def main_build(args):
    pass


def main_reads(args):
    pass


def main_bam(args):
    pass


def main(args):
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    args.wd = os.path.join(os.getcwd(), args.wd)

    logging.info("\nCMD: " + " ".join(sys.argv) + "\n")
    try:
        os.makedirs(args.wd)
    except FileExistsError:
        logging.critical("Output folder already exits.")
        logging.critical("Halting..\n")
        sys.exit(1)

    args.events = args.events.split(",")
    FQ1, FQ2 = args.FQs, None
    if "," in args.FQs:
        FQ1, FQ2 = args.FQs.split(",")

    logging.info("Running SUPPA2 to extract events..")
    suppa2_wd = os.path.join(args.wd, "suppa2")
    os.makedirs(suppa2_wd, exist_ok=True)
    retcode, suppa2log_path = run_suppa(args.GTF, suppa2_wd, args.events)
    if retcode != 0:
        logging.critical(f"SUPPA2 did not run succesfully (return code {retcode}).")
        logging.critical(f"See {suppa2log_path} for more details.")
        logging.critical("Halting..\n")
        sys.exit(1)
    suppa2_ioe_path = os.path.join(args.wd, "events.ioe")
    print_header = True
    with open(suppa2_ioe_path, "w") as outioe:
        for inioe_path in glob.glob(os.path.join(suppa2_wd, "*.ioe")):
            with open(inioe_path) as inioe:
                if not print_header:
                    next(inioe)
                print_header = False
                for line in inioe:
                    outioe.write(line)

    logging.info("Opening input annotation..")
    gtf = open_gtf(args.GTF)
    logging.info("Extracting exons from input GTF..")
    exons = extract_exons(gtf)
    esg_gtf_path = os.path.join(args.wd, "esgraphs.gtf")
    build_esg(exons, suppa2_wd, esg_gtf_path, args.fr)

    logging.info("Extracting local isoforms..")
    esg_fa_path = os.path.join(args.wd, "esgraphs.fa")
    gffread_log_path = os.path.join(args.wd, "gffread.log")
    retcode = run_gffread(args.FA, esg_gtf_path, esg_fa_path, gffread_log_path)
    if retcode != 0:
        logging.critical(f"gffread did not run succesfully (return code {retcode}).")
        logging.critical(f"See {gffread_log_path} for more details.")
        logging.critical("Halting..\n")
        sys.exit(1)

    logging.info("Filtering input RNA-Seq sample..")
    shark_wd = os.path.join(args.wd, "shark")
    os.makedirs(shark_wd, exist_ok=True)
    shark_ssv = os.path.join(shark_wd, "associations.ssv")
    sharked_1 = os.path.join(shark_wd, "sample_1.fq")
    sharked_2 = None
    if FQ2 is not None:
        sharked_2 = os.path.join(shark_wd, "sample_2.fq")
    shark_log = os.path.join(args.wd, "shark.log")
    retcode = shark.run(
        esg_fa_path, FQ1, FQ2, args.threads, shark_ssv, sharked_1, sharked_2, shark_log
    )
    if retcode != 0:
        logging.critical(f"shark did not run succesfully (return code {retcode}).")
        logging.critical(f"See {shark_log} for more details.")
        logging.critical("Halting..\n")
        sys.exit(1)

    logging.info("Splitting sharked sample..")
    esg_wd = os.path.join(args.wd, "esg")
    os.makedirs(esg_wd, exist_ok=True)
    shark.split(shark_ssv, sharked_1, sharked_2, esg_wd)

    chrom_dir = args.chroms
    if chrom_dir is None:
        logging.info("Splitting reference genome..")
        chrom_dir = os.path.join(args.wd, "chroms")
        os.makedirs(chrom_dir, exist_ok=True)
        for record in SeqIO.parse(args.FA, "fasta"):
            out_fa = open(os.path.join(chrom_dir, "{}.fa".format(record.id)), "w")
            SeqIO.write(record, out_fa, "fasta")
            out_fa.close()

    logging.info("Splitting event splicing graphs annotation..")
    # get events from splitted fqs (some events may not have reads, so we filter them here)
    events = [
        os.path.basename(fq)[:-3] for fq in glob.glob(os.path.join(esg_wd, "*.fq"))
    ]
    esg_gtf = open_gtf(esg_gtf_path)
    with console.status(f"[bold green]Splitting {len(events)} ESGs...") as _:
        for gene in esg_gtf.features_of_type("gene"):
            gidx = gene.id
            if '"' in gene.id:
                gidx = gidx.split('"')[1]
            if gidx not in events:
                continue
            gene.id = gidx
            with open(os.path.join(esg_wd, f"{gidx}.gtf"), "w") as ogtf:
                # FIXME don't know why but gffutils put additional "" around gene_id - and only for gene record
                gene_str = str(gene).replace('""', '"').replace('";"', '"')
                print(gene_str, file=ogtf)
                for transcript in esg_gtf.children(
                    gene, featuretype="transcript", order_by="start"
                ):
                    print(transcript, file=ogtf)
                    for exon in esg_gtf.children(
                        transcript, featuretype="exon", order_by="start"
                    ):
                        print(exon, file=ogtf)

    logging.info("Mapping reads to event splicing graphs..")
    # TODO: move what follows in the asgal.py module vvv
    spliceawarealigner = (
        "SpliceAwareAligner"
        if args.galig == "."
        else f"{args.galig}/bin/SpliceAwareAligner"
    )
    formatsam = (
        "formatSAM.py"
        if args.galig == "."
        else f"python3 {args.galig}/scripts/formatSAM.py"
    )
    asgal_args = [
        (chrom_dir, esg_wd, event, spliceawarealigner, formatsam, args.l)
        for event in events
    ]
    with console.status(f"[bold green]Mapping against {len(events)} ESGs...") as _:
        with mp.Pool(processes=args.threads) as pool:
            pool.starmap(asgal.run, asgal_args, chunksize=1000)

    logging.info("Merging BAMs..")
    bam_list = os.path.join(args.wd, "bams.list")
    bigbam_1 = os.path.join(args.wd, "asgal.merged.bam")
    with open(bam_list, "w") as txt:
        for bam in glob.glob(os.path.join(esg_wd, "*.bam")):
            txt.write(bam + "\n")
    subprocess.run(["samtools", "merge", bigbam_1, "-b", bam_list])
    subprocess.run(["samtools", "index", bigbam_1])

    logging.info("Cleaning BAM..")
    bigbam_2 = os.path.join(args.wd, "asgal.merged.nodups.bam")
    bigbam_final_unsrt = os.path.join(args.wd, "asgal.unsorted.bam")
    bigbam_final = os.path.join(args.wd, "asgal.bam")
    # ^^^
    asgal.remove_duplicates(bigbam_1, bigbam_2)
    pysam.index(bigbam_2)
    asgal.flag_secondary(bigbam_2, bigbam_final_unsrt)
    pysam.sort("-o", bigbam_final, bigbam_final_unsrt)
    pysam.index(bigbam_final)

    logging.info("Computing PSI..")
    psi_path = os.path.join(args.wd, "events.psi")
    psi.psi(bigbam_final, suppa2_ioe_path, psi_path)

    logging.info(f"Done! Please check {psi_path} :*\n")


if __name__ == "__main__":
    main(cli.parse_args())
