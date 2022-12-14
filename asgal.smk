configfile: "config.yaml"

import re
from os.path import join as pjoin
import glob

FA = config["fa"]
ODIR = config["odir"]
GTF = pjoin(ODIR, "strict.ioe.local.gtf")
FQDIR = pjoin(ODIR, "sharked")
ASGAL_DIR = config["galig"]

genes = {}
# here we take genes from the sharked set of samples since some gene in the annotation
# may not have reads (accordingly to shark output)
for fq in glob.glob(pjoin(FQDIR, "*.fq")):
    gene = os.path.basename(fq)[:-3]
    genes[gene] = ""
for line in open(GTF):
    if line.startswith("#"):
        continue
    line = line.strip("\n").split("\t")
    chrom = line[0]
    if line[2] == "gene":
        gene = line[-1].split("\"")[1]
        if gene not in genes:
            continue
        genes[gene] = chrom

rule run:
    input:
        pjoin(ODIR, "events.csv")

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

rule list_bams:
    input:
        expand(pjoin(ODIR, "{gene}", "ASGAL", "aligns.bam"),
               gene = genes.keys())
    output:
        pjoin(ODIR, "bams.list")
    threads: 1
    shell:
        """
        ls {input} > {output}
        """

rule merge_bams:
    input:
        pjoin(ODIR, "bams.list")
    output:
        pjoin(ODIR, "asgal.merge.bam")
    threads: 1
    shell:
        """
        samtools merge {output} -b {input}
        samtools index {output}
        """

rule clean_asgal_bam:
    input:
        pjoin(ODIR, "asgal.merge.bam")
    output:
        pjoin(ODIR, "asgal.bam")
    params:
        pjoin(ODIR, "asgal.tmp.bam")
    shell:
        """
        python3 clean_asgal_bam.py rdup {input} | samtools view -bS | samtools sort > {params}
        samtools index {params}
        python3 clean_asgal_bam.py fsec {params} | samtools view -bS | samtools sort > {output}
        samtools index {output}
        rm {params} {params}.bai
        """

rule compute_psi:
    input:
        ioe = pjoin(ODIR, "strict.ioe"),
        bam = pjoin(ODIR, "asgal.bam")
    output:
        pjoin(ODIR, "events.csv")
    shell:
        """
        python3 compute_psi.py {input.ioe} {input.bam} > {output}
        """
