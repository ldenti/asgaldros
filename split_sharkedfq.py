import sys, os
from Bio import SeqIO


def main():
    ssv_path = sys.argv[1]
    fq_path = sys.argv[2]
    odir = sys.argv[3]

    os.mkdir(odir)

    associations = {}
    for line in open(ssv_path):
        ridx, gidx = line.strip("\n").split(" ")
        gidx = "_".join(gidx.split("_")[:-1]) + "_G"
        if gidx not in associations:
            associations[gidx] = set()
        associations[gidx].add(ridx)

    # FIXME this may use too much memory if sharked reads set is huge
    # e.g., 1 272 662 reads -> 2.3GB
    reads = {}
    for record in SeqIO.parse(fq_path, "fastq"):
        reads[record.id] = record

    for gidx in associations:
        ofile = open(odir + "/" + gidx + ".fq", "w")
        for idx in associations[gidx]:
            SeqIO.write(reads[idx], ofile, "fastq")
        ofile.close()


if __name__ == "__main__":
    main()
