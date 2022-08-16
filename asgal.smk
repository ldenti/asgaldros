configfile: "config.yaml"

import re
from os.path import join as pjoin

FA = config["fa"]
ODIR = config["odir"]
GTF = pjoin(ODIR, "_SE_strict.ioe.local.gtf")
FQDIR = pjoin(ODIR, "sharked")
ASGAL_DIR = config["galig"]

genes = {}
for line in open(GTF):
    if line.startswith("#"):
        continue
    line = line.strip("\n").split("\t")
    chrom = line[0]
    if line[2] == "gene":
        gene = line[-1].split("\"")[1]
        # gene = re.match("gene_id \"([A-Za-z0-9\.]+\|\-)\";", line[-1]).group(1)
        genes[gene] = chrom
print(genes)
print(f"{len(genes)} genes")

rule run:
    input:
        expand(pjoin(ODIR, "{gene}", "ASGAL", "events.wpsi.csv"),
                gene = genes.keys())

rule split_reference:
    input:
        fa = FA
    output:
        fa = pjoin(ODIR, "chroms", "{chrom}.fa")
    threads: 1
    shell:
        """
        samtools faidx {input.fa} {wildcards.chrom} > {output.fa}
        """

rule split_annotation:
    input:
        gtf = GTF
    output:
        gtf = pjoin(ODIR, "{gene}", "annotation.gtf")
    threads: 1
    shell:
        """
        grep {wildcards.gene} {input.gtf} > {output.gtf}
        """

rule asgal:
    input:
        fa = lambda wildcards: pjoin(ODIR, "chroms", f"{genes[wildcards.gene]}.fa"),
        gtf = pjoin(ODIR, "{gene}", "annotation.gtf"),
        fq = lambda wildcards: pjoin(FQDIR, "{gene}.fq")
    output:
        sam = pjoin(ODIR, "{gene}", "ASGAL", "aligns.sam")
    params:
        mem = pjoin(ODIR, "{gene}", "ASGAL", "aligns.mem"),
    threads: 1
    shell:
        """
        {ASGAL_DIR}/bin/SpliceAwareAligner -g {input.fa} -a {input.gtf} -s {input.fq} -o {params.mem} -l 7
        python3 {ASGAL_DIR}/scripts/formatSAM.py -m {params.mem} -g {input.fa} -a {input.gtf} -o {output.sam}
        """

rule sam2bam:
    input:
        "{f}.sam"
    output:
        "{f}.bam"
    threads: 1
    shell:
        """
        samtools view -bS {input} | samtools sort > {output}
        samtools index {output}
        """

rule compute_psi:
    input:
        gtf = pjoin(ODIR, "{gene}", "annotation.gtf"),
        bam = pjoin(ODIR, "{gene}", "ASGAL", "aligns.bam")
    output:
        csv = pjoin(ODIR, "{gene}", "ASGAL", "events.wpsi.csv")
    shell:
        """
        python3 compute_psi.py {input.gtf} {input.bam} > {output.csv}
        """
