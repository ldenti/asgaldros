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

    for gidx in associations:
        print(f"Parsing {gidx}..")
        ofile = open(odir + "/" + gidx + ".fq", "w")
        for record in SeqIO.parse(fq_path, "fastq"):
            if record.id in associations[gidx]:
                SeqIO.write(record, ofile, "fastq")
        ofile.close()


if __name__ == "__main__":
    main()
