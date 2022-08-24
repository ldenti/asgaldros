import sys
import os
import glob
import gffutils


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


def analyze_ES(exons, fpath):
    for line in open(fpath):
        if line.startswith("seqname"):
            continue
        idx = line.strip("\n").split("\t")[2]
        gene, rest = idx.split(";")
        _etype, chrom, intron1, intron2, strand = rest.split(":")
        intron1 = (int(intron1.split("-")[0]), int(intron1.split("-")[1]))
        intron2 = (int(intron2.split("-")[0]), int(intron2.split("-")[1]))
        idx = idx.replace(";", "_")
        idx = idx.replace(":", "_")
        assert (
            intron1[0] < intron1[1]
            and intron1[1] < intron2[0]
            and intron2[0] < intron2[1]
        )
        pre, mid, post = (0, 0), (0, 0), (0, 0)
        for (start, end) in exons[chrom][gene]:
            if end == intron1[0]:
                if pre != (0, 0):
                    if start < pre[0]:
                        pre = (start, end)
                else:
                    pre = (start, end)
            elif start == intron1[1] and end == intron2[0]:
                mid = (start, end)
            elif start == intron2[1]:
                if post != (0, 0):
                    if end > post[1]:
                        post = (start, end)
                else:
                    post = (start, end)
        print(
            chrom,
            ".",
            "gene",
            pre[0],
            post[1],
            ".",
            "+",
            ".",
            f'gene_id "{idx}_G";',
            sep="\t",
        )
        print(
            chrom,
            ".",
            "transcript",
            pre[0],
            post[1],
            ".",
            "+",
            ".",
            f'gene_id "{idx}_G"; transcript_id "{idx}_T1"',
            sep="\t",
        )
        for i, (s, e) in enumerate([pre, mid, post], 1):
            print(
                chrom,
                ".",
                "exon",
                s,
                e,
                ".",
                "+",
                ".",
                f'gene_id "{idx}_G"; transcript_id "{idx}_T1"; exon_id "{idx}_E{i}"',
                sep="\t",
            )
        print(
            chrom,
            ".",
            "transcript",
            pre[0],
            post[1],
            ".",
            "+",
            ".",
            f'gene_id "{idx}_G"; transcript_id "{idx}_T2"',
            sep="\t",
        )
        for i, (s, e) in enumerate([pre, (0, 0), post], 1):
            if s == 0 and e == 0:
                continue
            print(
                chrom,
                ".",
                "exon",
                s,
                e,
                ".",
                "+",
                ".",
                f'gene_id "{idx}_G"; transcript_id "{idx}_T2"; exon_id "{idx}_E{i}"',
                sep="\t",
            )


def analyze_SS(exons, fpath):
    for line in open(fpath):
        if line.startswith("seqname"):
            continue
        idx = line.strip("\n").split("\t")[2]
        gene, rest = idx.split(";")
        _etype, chrom, intron1, intron2, strand = rest.split(":")
        intron1 = (int(intron1.split("-")[0]), int(intron1.split("-")[1]))
        intron2 = (int(intron2.split("-")[0]), int(intron2.split("-")[1]))
        idx = idx.replace(";", "_")
        idx = idx.replace(":", "_")
        mode = (strand == "+" and _etype == "A5") or (strand == "-" and _etype == "A3")

        if mode:
            assert intron1[1] == intron2[1]
        else:
            assert intron1[0] == intron2[0]

        alt1, alt2, const = (0, 0), (0, 0), (0, 0)
        for (start, end) in exons[chrom][gene]:
            if mode:
                if start == intron1[1]:
                    if const != (0, 0):
                        if end > const[1]:
                            const = (start, end)
                    else:
                        const = (start, end)
                elif end == intron1[0]:
                    if alt1 != (0, 0):
                        if start < alt1[0]:
                            alt1 = (start, end)
                    else:
                        alt1 = (start, end)
                elif end == intron2[0]:
                    if alt2 != (0, 0):
                        if start < alt2[0]:
                            alt2 = (start, end)
                    else:
                        alt2 = (start, end)
            else:
                if end == intron1[0]:
                    if const != (0, 0):
                        if start < const[0]:
                            const = (start, end)
                    else:
                        const = (start, end)
                elif start == intron1[1]:
                    if alt1 != (0, 0):
                        if end > alt1[1]:
                            alt1 = (start, end)
                    else:
                        alt1 = (start, end)
                elif start == intron2[1]:
                    if alt2 != (0, 0):
                        if end > alt2[1]:
                            alt2 = (start, end)
                    else:
                        alt2 = (start, end)

        Ts = {}
        begin, end = 0, 0
        if mode:
            begin = min(alt1[0], alt2[0])
            end = const[1]
            Ts = {1: [(1, alt1), (2, const)], 2: [(3, alt2), (2, const)]}
        else:
            begin = const[0]
            end = max(alt1[1], alt2[1])
            Ts = {1: [(1, const), (2, alt1)], 2: [(1, const), (3, alt2)]}

        print(
            chrom,
            ".",
            "gene",
            begin,
            end,
            ".",
            "+",
            ".",
            f'gene_id "{idx}_G";',
            sep="\t",
        )
        for T, Es in Ts.items():
            print(
                chrom,
                ".",
                "transcript",
                Es[0][1][0],
                Es[-1][1][-1],
                ".",
                "+",
                ".",
                f'gene_id "{idx}_G"; transcript_id "{idx}_T{T}"',
                sep="\t",
            )
            for (i, (s, e)) in Es:
                print(
                    chrom,
                    ".",
                    "exon",
                    s,
                    e,
                    ".",
                    "+",
                    ".",
                    f'gene_id "{idx}_G"; transcript_id "{idx}_T{T}"; exon_id "{idx}_E{i}"',
                    sep="\t",
                )


def main():
    gtf_path = sys.argv[1]
    indir = sys.argv[2]

    gtf = open_gtf(gtf_path)
    print("Extracting exons..", file=sys.stderr)
    exons = {}
    for exon in gtf.features_of_type("exon"):
        chrom = exon.chrom
        gidx = exon.attributes["gene_id"][0]
        if chrom not in exons:
            exons[chrom] = {}
        if gidx not in exons[chrom]:
            exons[chrom][gidx] = set()
        exons[chrom][gidx].add((exon.start, exon.end))

    print("Analyzing ES..", file=sys.stderr)
    analyze_ES(exons, os.path.join(indir, "_SE_strict.ioe"))
    print("Analyzing A3..", file=sys.stderr)
    analyze_SS(exons, os.path.join(indir, "_A3_strict.ioe"))
    print("Analyzing A5..", file=sys.stderr)
    analyze_SS(exons, os.path.join(indir, "_A5_strict.ioe"))


if __name__ == "__main__":
    main()
