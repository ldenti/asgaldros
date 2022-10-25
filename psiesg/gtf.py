import subprocess
import gffutils


def open_gtf(gtf_path):
    try:
        gtf = gffutils.FeatureDB(f"{gtf_path}.db", keep_order=True)
    except ValueError:
        gtf = gffutils.create_db(
            gtf_path,
            dbfn=f"{gtf_path}.db",
            force=True,
            keep_order=True,
            disable_infer_genes=True,
            disable_infer_transcripts=True,
            merge_strategy="merge",
            sort_attribute_values=True,
        )
    return gtf


def extract_exons(gtf):
    exons = {}
    for exon in gtf.features_of_type("exon"):
        chrom = exon.chrom
        gidx = exon.attributes["gene_id"][0]
        if chrom not in exons:
            exons[chrom] = {}
        if gidx not in exons[chrom]:
            exons[chrom][gidx] = set()
        exons[chrom][gidx].add((exon.start, exon.end))
    return exons


def run_gffread(FA, GTF, fa, log):
    # TODO: we can do this manually
    cmd = ["gffread", "-g", FA, GTF, "-w", fa]
    p = subprocess.run(cmd, stderr=open(log, "w"))
    return p.returncode
