import argparse
import pysam


def main():
    parser = argparse.ArgumentParser(description="Compute PSI from .bam file")
    parser.add_argument("IOE", type=str, help="Path to SUPPA2 .ioe file")
    parser.add_argument("BAM", type=str, help="Path to .bam alignments")
    args = parser.parse_args()

    bam = pysam.AlignmentFile(args.BAM, "rb")

    print(
        "Gene",
        "Chrom",
        "Event",
        "Intron1",
        "Intron2",
        "Intron3",
        "w1",
        "w2",
        "w3",
        "PSI",
        "OtherIntrons",
        sep=",",
    )
    for line in open(args.IOE):
        if line.startswith("seqname"):
            continue
        idx = line.strip("\n").split("\t")[2]
        gene, rest = idx.split(";")
        _etype, chrom, intron1, intron2, strand = rest.split(":")

        intron1 = (int(intron1.split("-")[0]), int(intron1.split("-")[1]))
        intron2 = (int(intron2.split("-")[0]), int(intron2.split("-")[1]))
        intron3 = (0, 0)
        idx = idx.replace(";", "_")
        idx = idx.replace(":", "_")
        mode = (strand == "+" and _etype == "A5") or (strand == "-" and _etype == "A3")
        if _etype == "SE":
            assert (
                intron1[0] < intron1[1]
                and intron1[1] < intron2[0]
                and intron2[0] < intron2[1]
            )
            intron3 = (intron1[0], intron2[1])
        else:  # A3, A5
            if mode:
                assert intron1[1] == intron2[1]
            else:
                assert intron1[0] == intron2[0]
            intron1_tmp = (
                intron1
                if intron1[1] - intron1[0] > intron2[1] - intron2[0]
                else intron2
            )
            intron2_tmp = (
                intron2
                if intron1[1] - intron1[0] > intron2[1] - intron2[0]
                else intron1
            )
            intron1 = intron1_tmp
            intron2 = intron2_tmp

        bam_introns = bam.find_introns(
            (
                read
                for read in bam.fetch(
                    chrom,
                    min(intron1[0], intron2[0]) - 50,
                    max(intron1[1], intron2[1]) + 50,
                )
            )
        )
        w1, w2, w3 = 0, 0, 0
        new_introns = {}
        for (s, e), w in bam_introns.items():
            e += 1  # CHECKME why do we need this?
            if (s, e) == intron1:
                w1 = w
            elif (s, e) == intron2:
                w2 = w
            elif (s, e) == intron3:
                w3 = w
            else:
                new_introns[(s, e)] = w

        PSI = 0
        if _etype == "SE":
            if w1 == 0 and w2 == 0 and w3 == 0:
                PSI = -1
            else:
                PSI = ((w1 + w2) / 2) / ((w1 + w2) / 2 + w3)
        else:  # A3, A5
            if w1 == 0 and w2 == 0:
                PSI = -1
            else:
                PSI = w1 / (w1 + w2)
        PSI = round(PSI, 3)

        print(
            gene,
            chrom,
            _etype,
            f"{intron1[0]}-{intron1[1]}",
            f"{intron2[0]}-{intron2[1]}",
            f"{intron3[0]}-{intron3[1]}" if _etype == "SE" else ".",
            w1,
            w2,
            w3 if _etype == "SE" else -1,
            PSI,
            "/".join(f"{s}-{e}:{w}" for (s, e), w in new_introns.items()),
            sep=",",
        )


if __name__ == "__main__":
    main()
