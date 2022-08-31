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


def print_gene(chrom, idx, begin, end, Ts):
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


def analyze_ES(exons, fpath):
    if not os.path.exists(fpath):
        return
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
        pre, mid, post = (
            (float("inf"), intron1[0]),
            (intron1[1], intron2[0]),
            (intron2[1], -1),
        )
        for (start, end) in exons[chrom][gene]:
            if end == pre[1]:
                if start < pre[0]:
                    pre = (start, end)
            elif start == post[0]:
                if end > post[1]:
                    post = (start, end)
        begin = pre[0]
        end = post[1]
        Ts = {1: [(1, pre), (2, mid), (3, post)], 2: [(1, pre), (3, post)]}
        print_gene(chrom, idx, begin, end, Ts)


def analyze_SS(exons, fpath):
    if not os.path.exists(fpath):
        return
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
            alt1, alt2, const = (
                (float("inf"), intron1[0]),
                (float("inf"), intron2[0]),
                (intron1[1], -1),
            )
        else:
            alt1, alt2, const = (
                (intron1[1], -1),
                (intron2[1], -1),
                (float("inf"), intron1[0]),
            )
        for (start, end) in exons[chrom][gene]:
            if mode:
                if start == const[0]:
                    if end > const[1]:
                        const = (start, end)
                elif end == alt1[1]:
                    if start < alt1[0]:
                        alt1 = (start, end)
                elif end == alt2[1]:
                    if start < alt2[0]:
                        alt2 = (start, end)
            else:
                if end == const[1]:
                    if start < const[0]:
                        const = (start, end)
                elif start == alt1[0]:
                    if end > alt1[1]:
                        alt1 = (start, end)
                elif start == alt2[0]:
                    if end > alt2[1]:
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

        print_gene(chrom, idx, begin, end, Ts)


def analyze_IR(fpath):
    if not os.path.exists(fpath):
        return
    for line in open(fpath):
        if line.startswith("seqname"):
            continue
        idx = line.strip("\n").split("\t")[2]
        gene, rest = idx.split(";")
        _etype, chrom, exon1_s, intron, exon2_e, strand = rest.split(":")
        exon1 = (int(exon1_s), int(intron.split("-")[0]))
        exon2 = (int(intron.split("-")[1]), int(exon2_e))
        idx = idx.replace(";", "_")
        idx = idx.replace(":", "_")
        begin = exon1[0]
        end = exon2[1]
        Ts = {1: [(1, exon1), (2, exon2)], 2: [(3, (exon1[1] + 1, exon2[0] - 1))]}
        print_gene(chrom, idx, begin, end, Ts)


def analyze_FL(exons, fpath):
    if not os.path.exists(fpath):
        return
    for line in open(fpath):
        if line.startswith("seqname"):
            continue
        idx = line.strip("\n").split("\t")[2]
        gene, rest = idx.split(";")
        _etype, chrom, *positions, strand = rest.split(":")
        mode = (strand == "+" and _etype == "AF") or (strand == "-" and _etype == "AL")
        idx = idx.replace(";", "_")
        idx = idx.replace(":", "_")

        if mode:
            a, bc, d, ef = positions
            b, c = bc.split("-")
            e, f = ef.split("-")
            a, b, c, d, e, f = int(a), int(b), int(c), int(d), int(e), int(f)
            exon1 = (a, b)
            exon2 = (d, e)
            const = (c, -1)
        else:
            ab, c, de, f = positions
            a, b = ab.split("-")
            d, e = de.split("-")
            a, b, c, d, e, f = int(a), int(b), int(c), int(d), int(e), int(f)
            exon1 = (b, c)
            exon2 = (e, f)
            const = (float("inf"), a)

        for (start, end) in exons[chrom][gene]:
            if mode:
                # I'm looking for the end
                if start == const[0]:
                    if end > const[1]:
                        const = (start, end)
            else:
                # I'm looking for the start
                if end == const[1]:
                    if start < const[0]:
                        const = (start, end)

        Ts = {}
        begin, end = 0, 0
        if mode:
            begin = exon1[0]
            end = const[1]
            Ts = {1: [(1, exon1), (2, const)], 2: [(3, exon2), (2, const)]}
        else:
            begin = const[0]
            end = exon2[1]
            Ts = {1: [(1, const), (2, exon1)], 2: [(1, const), (3, exon2)]}

        print_gene(chrom, idx, begin, end, Ts)


def analyze_MX(exons, fpath):
    if not os.path.exists(fpath):
        return
    for line in open(fpath):
        if line.startswith("seqname"):
            continue
        idx = line.strip("\n").split("\t")[2]
        gene, rest = idx.split(";")
        _etype, chrom, intron1, intron2, intron3, intron4, strand = rest.split(":")
        intron1 = intron1.split("-")
        intron2 = intron2.split("-")
        intron3 = intron3.split("-")
        intron4 = intron4.split("-")
        const1 = (float("inf"), int(intron1[0]))
        exon1 = (int(intron1[1]), int(intron2[0]))
        exon2 = (int(intron3[1]), int(intron4[0]))
        const2 = (int(intron4[1]), -1)
        idx = idx.replace(";", "_")
        idx = idx.replace(":", "_")

        for (start, end) in exons[chrom][gene]:
            if end == const1[1]:
                if start < const1[0]:
                    const1 = (start, end)
            elif start == const2[0]:
                if end > const2[1]:
                    const2 = (start, end)
        assert const1[0] != float("inf")
        assert const2[1] != -1
        begin = const1[0]
        end = const2[1]
        Ts = {
            1: [(1, const1), (2, exon1), (3, const2)],
            2: [(1, const1), (4, exon2), (3, const2)],
        }
        print_gene(chrom, idx, begin, end, Ts)


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
    print("Analyzing IR..", file=sys.stderr)
    analyze_IR(os.path.join(indir, "_RI_strict.ioe"))
    print("Analyzing AF..", file=sys.stderr)
    analyze_FL(exons, os.path.join(indir, "_AF_strict.ioe"))
    print("Analyzing AL..", file=sys.stderr)
    analyze_FL(exons, os.path.join(indir, "_AL_strict.ioe"))
    print("Analyzing MX..", file=sys.stderr)
    analyze_MX(exons, os.path.join(indir, "_MX_strict.ioe"))


if __name__ == "__main__":
    main()
