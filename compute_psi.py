import sys
import gffutils
import pysam


def open_gtf(gtf_path):
    try:
        gtf = gffutils.FeatureDB("{}.db".format(gtf_path), keep_order=True)
    except ValueError:
        gtf = gffutils.create_db(
            gtf_path,
            dbfn="{}.db".format(gtf_path),
            force=True,
            keep_order=True,
            disable_infer_genes=True,
            disable_infer_transcripts=True,
            merge_strategy="merge",
            sort_attribute_values=True,
        )
    return gtf


def i2s(I):
    return f"{I[0]}-{I[1]}"


def se_psi(introns, bam):
    introns = sorted(introns, key=lambda x: len(x), reverse=True)
    I1, I2, I3 = introns[0][0], introns[0][1], introns[1][0]
    w1, w2, w3 = 0, 0, 0
    bamintrons = bam.find_introns((read for read in bam.fetch()))
    for (s, e), w in bamintrons.items():
        if (s, e) == I1:
            w1 = w
        elif (s, e) == I2:
            w2 = w
        elif (s, e) == I3:
            w3 = w
    PSI = 0
    if w1 + w2 + w3 == 0:
        PSI = -1
    else:
        PSI = ((w1 + w2) / 2) / ((w1 + w2) / 2 + w3)
    PSI = round(PSI, 3)
    positions = f"{i2s(I1)}/{i2s(I2)}|{i2s(I3)}"
    return positions, w1, w2, w3, PSI


def ss_psi(introns, bam):
    introns = sorted(introns, key=lambda x: x[0][1] - x[0][0])
    I1, I2 = introns[0][0], introns[1][0]
    w1, w2 = 0, 0
    introns = bam.find_introns((read for read in bam.fetch()))
    for (s, e), w in introns.items():
        if (s, e) == I1:
            w1 = w
        elif (s, e) == I2:
            w2 = w
    PSI = 0
    if w1 + w2 == 0:
        PSI = -1
    else:
        PSI = w1 / (w1 + w2)
    PSI = round(PSI, 3)
    positions = f"{i2s(I1)}|{i2s(I2)}"
    return positions, w1, w2, PSI


def ri_psi(introns, bam, chrom):
    I1 = introns[0][0]
    w1 = 0
    # Spliced reads
    introns = bam.find_introns((read for read in bam.fetch()))
    for (s, e), w in introns.items():
        if (s, e) == I1:
            w1 = w
    w2 = 0
    # Reads fully aligned to intron
    for al in bam.fetch(chrom, I1[0], I1[1]):
        if I1[0] <= al.reference_start and al.reference_end <= I1[1]:
            w2 += 1
    PSI = 0
    if w1 + w2 == 0:
        PSI = -1
    else:
        PSI = w1 / (w1 + w2)
    PSI = round(PSI, 3)
    positions = f"{i2s(I1)}"
    return positions, w1, w2, PSI


def mx_psi(introns, bam):
    I1, I2, I3, I4 = introns[0][0], introns[0][1], introns[1][0], introns[1][1]
    w1, w2, w3, w4 = 0, 0, 0, 0
    bam_introns = bam.find_introns((read for read in bam.fetch()))
    for (s, e), w in bam_introns.items():
        if (s, e) == I1:
            w1 = w
        elif (s, e) == I2:
            w2 = w
        elif (s, e) == I3:
            w3 = w
        elif (s, e) == I4:
            w4 = w
    PSI = 0
    if w1 + w2 + w3 + w4 == 0:
        PSI = -1
    else:
        PSI = (w1 + w2) / (w1 + w2 + w3 + w4)
    PSI = round(PSI, 3)
    positions = f"{i2s(I1)}/{i2s(I2)}|{i2s(I3)}/{i2s(I4)}"
    return positions, w1, w2, w3, w4, PSI


def fl_psi(introns, bam):
    introns = sorted(introns, key=lambda x: x[0][1] - x[0][0], reverse=True)
    I1, I2 = introns[0][0], introns[1][0]
    w1, w2 = 0, 0
    bam_introns = bam.find_introns((read for read in bam.fetch()))
    for (s, e), w in bam_introns.items():
        if (s, e) == I1:
            w1 = w
        elif (s, e) == I2:
            w2 = w
    PSI = 0
    if w1 + w2 == 0:
        PSI = -1
    else:
        PSI = w1 / (w1 + w2)
    PSI = round(PSI, 3)
    positions = f"{i2s(I1)}|{i2s(I2)}"
    return positions, w1, w2, PSI


def main():
    gtf_path = sys.argv[1]
    bam_path = sys.argv[2]

    gtf = open_gtf(gtf_path)
    bam = pysam.AlignmentFile(bam_path, "rb")
    etype = bam_path.split("/")[-3].split("_")[1]

    chrom, idx = None, None
    introns = []
    for transcript in gtf.features_of_type("transcript"):
        chrom = transcript.chrom
        idx = transcript.attributes["gene_id"][0][:-2]
        _exons = [
            (exon.start, exon.end)
            for exon in gtf.children(transcript, featuretype="exon", order_by="start")
        ]
        _introns = [(_exons[i - 1][1], _exons[i][0] - 1) for i in range(1, len(_exons))]
        introns.append(_introns)
    gene = idx.split("_")[0]
    strand = idx.split("_")[-1]

    w1, w2, w3, w4 = -1, -1, -1, -1
    PSI = -1
    if etype == "SE":
        positions, w1, w2, w3, PSI = se_psi(introns, bam)
    elif etype == "A3" or etype == "A5":
        positions, w1, w2, PSI = ss_psi(introns, bam)
    elif etype == "RI":
        positions, w1, w2, PSI = ri_psi(introns, bam, chrom)
    elif etype == "MX":
        positions, w1, w2, w3, w4, PSI = mx_psi(introns, bam)
    elif etype == "AF" or etype == "AL":
        positions, w1, w2, PSI = fl_psi(introns, bam)

    w1 = w1 if w1 != -1 else "."
    w2 = w2 if w2 != -1 else "."
    w3 = w3 if w3 != -1 else "."
    w4 = w4 if w4 != -1 else "."
    PSI = PSI if PSI != -1 else "NaN"
    print("Type,Chrom,Gene,Strand,Event,W1,W2,W3,W4,PSI")
    print(etype, chrom, gene, strand, positions, w1, w2, w3, w4, PSI, sep=",")


if __name__ == "__main__":
    main()
