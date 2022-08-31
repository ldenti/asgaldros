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
        "Intron4",
        "w1",
        "w2",
        "w3",
        "w4",
        "PSI",
        "OtherIntrons",
        sep=",",
    )
    for line in open(args.IOE):
        if line.startswith("seqname"):
            continue
        idx = line.strip("\n").split("\t")[2]
        gene, rest = idx.split(";")
        idx = idx.replace(";", "_")
        idx = idx.replace(":", "_")
        _etype, chrom, *positions, strand = rest.split(":")

        intron1 = (0, 0)
        intron2 = (0, 0)
        intron3 = (0, 0)
        intron4 = (0, 0)

        if _etype == "SE":
            ab, cd = positions
            intron1 = ab.split("-")
            intron2 = cd.split("-")
            intron3 = (intron1[0], intron2[1])
        elif _etype == "A3" or _etype == "A5":
            ab, cd = positions
            intron1 = ab.split("-")
            intron2 = cd.split("-")
        elif _etype == "RI":
            a, bc, d = positions
            intron1 = (b, c)
        elif _etype == "AF" or _etype == "AL":
            mode = (strand == "+" and _etype == "AF") or (
                strand == "-" and _etype == "AL"
            )
            if mode:
                a, bc, d, ef = positions
                b, c = bc.split("-")
                e, f = ef.split("-")
                intron1 = (b, c)
                intron2 = (e, f)
            else:
                ab, c, de, f = positions
                a, b = ab.split("-")
                d, e = de.split("-")
                intron1 = (a, b)
                intron2 = (d, e)
        elif _etype == "MX":
            ab, cd, ef, gh = positions
            intron1 = ab.split("-")
            intron2 = cd.split("-")
            intron3 = ef.split("-")
            intron4 = gh.split("-")

        intron1 = [int(p) for p in intron1]
        intron2 = [int(p) for p in intron2]
        intron3 = [int(p) for p in intron3]
        intron4 = [int(p) for p in intron4]

        bam_introns = bam.find_introns(
            (
                read
                for read in bam.fetch(
                    chrom,
                    min(intron1[0] if intron1[0] != 0 else float("inf"), intron2[0] if intron2[0] != 0 else float("inf"), intron3[0] if intron3[0] != 0 else float("inf"), intron4[0] if intron4[0] != 0 else float("inf")) - 50,
                    max(intron1[1], intron2[1], intron3[1], intron4[1]) + 50,
                )
            )
        )

        w1, w2, w3, w4 = 0, 0, 0, 0
        new_introns = {}
        for (s, e), w in bam_introns.items():
            e += 1  # CHECKME why do we need this?
            if [s, e] == intron1:
                w1 = w
            elif [s, e] == intron2:
                w2 = w
            elif [s, e] == intron3:
                w3 = w
            elif [s, e] == intron4:
                w4 = w
            else:
                new_introns[(s, e)] = w
        if _etype == "RI":
            for al in bam.fetch(chrom, intron1[0], intron1[1]):
                if intron1[0] <= al.reference_start and al.reference_end <= intron1[1]:
                    w2 += 1

        PSI = 0
        if _etype == "SE":
            if w1 + w2 + w3 == 0:
                PSI = -1
            else:
                PSI = ((w1 + w2) / 2) / ((w1 + w2) / 2 + w3)
        elif _etype == "A3" or _etype == "A5":
            if w1 + w2 == 0:
                PSI = -1
            else:
                PSI = w1 / (w1 + w2)
        elif _etype == "RI":
            if w1 + w2 == 0:
                PSI = -1
            else:
                PSI = w1 / (w1 + w2)
        elif _etype == "AF" or _etype == "AL":
            if w1 + w2 == 0:
                PSI = -1
            else:
                PSI = w1 / (w1 + w2)
        elif _etype == "MX":
            if w1 + w2 + w3 + w4 == 0:
                PSI = -1
            else:
                PSI = (w1 + w2) / (w1 + w2 + w3 + w4)
        PSI = round(PSI, 3)

        print(
            gene,
            chrom,
            _etype,
            f"{intron1[0]}-{intron1[1]}",
            f"{intron2[0]}-{intron2[1]}" if intron2 != [0, 0] else ".",
            f"{intron3[0]}-{intron3[1]}" if intron3 != [0, 0] else ".",
            f"{intron4[0]}-{intron4[1]}" if intron4 != [0, 0] else ".",
            w1,
            w2,
            w3,
            w4,
            PSI if PSI != -1 else "NaN",
            "/".join(f"{s}-{e}:{w}" for (s, e), w in new_introns.items()),
            sep=",",
        )


if __name__ == "__main__":
    main()
