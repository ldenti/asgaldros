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


def main():
    gtf_path = sys.argv[1]
    bam_path = sys.argv[2]
    MINW = 2

    chrom, idx = None, None
    introns = set()
    gtf = open_gtf(gtf_path)
    for transcript in gtf.features_of_type("transcript"):
        chrom = transcript.chrom
        idx = transcript.attributes["gene_id"][0][:-2]
        exons = [
            (exon.start, exon.end)
            for exon in gtf.children(transcript, featuretype="exon", order_by="start")
        ]
        _introns = set((exons[i - 1][1], exons[i][0] - 1) for i in range(1, len(exons)))
        introns |= _introns

    I1 = min(introns, key=lambda x: x[1])
    I2 = max(introns, key=lambda x: x[0])
    I3 = max(introns, key=lambda x: x[1] - x[0])
    w1, w2, w3 = 0, 0, 0

    bam = pysam.AlignmentFile(bam_path, "rb")
    introns = bam.find_introns((read for read in bam.fetch()))
    for (s, e), w in introns.items():
        if (s, e) == I1:
            w1 = w
        elif (s, e) == I2:
            w2 = w
        elif (s, e) == I3:
            w3 = w
    if w1 > MINW and w2 > MINW and w3 > MINW:
        print(chrom, idx, w1, w2, w3, (w1 + w2) / (w1 + w2 + 2 * w3))


if __name__ == "__main__":
    main()
