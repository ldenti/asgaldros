import sys
import os
import glob
import pandas as pd

def main():
    indir = sys.argv[1]
    data = []
    for csv in glob.glob(os.path.join(indir, "*", "*", "ASGAL", "events.wpsi.csv")):
        etype = csv.split("/")[-3].split("_")[1]

        sample = csv.split("/")[-4]
        if etype == "SE":
            chrom, info, w1, w2, w3, psi = open(csv).readline().strip("\n").split(",")
        else: # A3, A5, RI
            chrom, info, w1, w2, psi = open(csv).readline().strip("\n").split(",")
            w3 = -1
        if etype == "RI":
            gene, etype, _chrom, E1, I, E2, strand = info.split("_")
            I1 = I.split("-")[0]
            I2 = I.split("-")[1]
        else:
            gene, etype, _chrom, I1, I2, strand = info.split("_")
        data.append([sample, chrom, gene, etype, f"{I1}-{I2}", strand, w1, w2, w3, psi])

    df = pd.DataFrame(data, columns = ["Sample", "Chromosome", "Gene", "Event", "Introns", "Strand", "W1", "W2", "W3", "PSI"])
    df.to_csv(sys.stdout, index=False)

if __name__ == "__main__":
    main()
