import subprocess

from Bio import SeqIO


def run(fa, fq1, fq2, threads, out, outfq1, outfq2, log):
    cmd = [
        "shark",
        "-r",
        fa,
        "-1",
        fq1,
        "-2",
        fq2,
        "--out1",
        outfq1,
        "--out2",
        outfq2,
        "-t",
        str(threads),
    ]
    if fq2 is None:
        cmd = [
            "shark",
            "-r",
            fa,
            "-1",
            fq1,
            "--out1",
            outfq1,
            "-t",
            str(threads),
        ]
    cmd_str = " ".join(cmd) + f" > {out} 2> {log}"
    # CHECKME: if I use shell=False, subprocess doesn't wait for content to be written to stdout file
    p = subprocess.run(cmd_str, shell=True)
    # p = subprocess.run(cmd_str, stdout=open(out, "w"), stderr=open(log, "w"))
    return p.returncode


def split(ssv, fq1, fq2, odir):
    associations = {}
    for line in open(ssv):
        ridx, gidx = line.strip("\n").split(" ")
        gidx = "_".join(gidx.split("_")[:-1]) + "_G"
        if gidx not in associations:
            associations[gidx] = set()
        associations[gidx].add(ridx)

    # FIXME this may use too much memory if sharked reads set is huge
    # e.g., 1 272 662 reads -> 2.3GB
    reads = {}
    for record in SeqIO.parse(fq1, "fastq"):
        reads[record.id] = record

    for gidx in associations:
        ofile = open(odir + "/" + gidx + ".fq", "w")
        for idx in associations[gidx]:
            SeqIO.write(reads[idx], ofile, "fastq")
        ofile.close()

    if fq2 is not None:
        reads = {}
        for record in SeqIO.parse(fq2, "fastq"):
            reads[record.id] = record

        for gidx in associations:
            ofile = open(odir + "/" + gidx + ".fq", "a")
            for idx in associations[gidx]:
                SeqIO.write(reads[idx], ofile, "fastq")
            ofile.close()
