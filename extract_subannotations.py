import sys
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


def main():
    gtf_path = sys.argv[1]
    suppa_events = sys.argv[2]

    gtf = open_gtf(gtf_path)

    for line in open(suppa_events):
        if line.startswith("seqname"):
            continue
        idx = line.strip("\n").split("\t")[2]
        gene, rest = idx.split(";")
        etype, chrom, intron1, intron2, strand = rest.split(":")
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
        for exon in gtf.features_of_type("exon"):
            if exon.chrom == chrom and exon.attributes["gene_id"][0] == gene:
                if exon.end == intron1[0]:
                    if pre != (0, 0):
                        if exon.start < pre[0]:
                            pre = (exon.start, exon.end)
                    else:
                        pre = (exon.start, exon.end)
                elif exon.start == intron1[1] and exon.end == intron2[0]:
                    mid = (exon.start, exon.end)
                elif exon.start == intron2[1]:
                    if post != (0, 0):
                        if exon.end > post[1]:
                            post = (exon.start, exon.end)
                    else:
                        post = (exon.start, exon.end)
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


if __name__ == "__main__":
    main()
