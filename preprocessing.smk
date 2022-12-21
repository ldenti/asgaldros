configfile: "config.yaml"

import re
from os.path import join as pjoin

FA = config["fa"]
GTF = config["gtf"]
FQ = config["fq"]
ODIR = config["odir"]
THREADS = config["threads"]

genes = {}
for line in open(GTF):
    if line.startswith("#"):
        continue
    line = line.strip("\n").split("\t")
    chrom = line[0]
    if line[2] == "gene":
        gene = re.match("gene_id \"([A-Za-z0-9\.]+)\";", line[-1]).group(1)
        genes[gene] = chrom

rule run:
    input:
        pjoin(ODIR, "strict.ioe"),
        pjoin(ODIR, "sharked")

rule trim:
    input:
        FQ
    output:
        FQ + ".trimmed.fq"
    conda:
        "envs/preprocessing.yaml"
    threads: THREADS
    shell:
        """
        cutadapt -j {threads} -m 40 -a "A{{20}}" -g "T{{20}}" {input} > {output}
        """

rule suppa:
    input:
        GTF
    output:
        ioe_se = pjoin(ODIR, "_SE_strict.ioe"),
        ioe_a3 = pjoin(ODIR, "_A3_strict.ioe"),
        ioe_a5 = pjoin(ODIR, "_A5_strict.ioe"),
        ioe_ri = pjoin(ODIR, "_RI_strict.ioe"),
        ioe_af = pjoin(ODIR, "_AF_strict.ioe"),
        ioe_al = pjoin(ODIR, "_AL_strict.ioe"),
        ioe_mx = pjoin(ODIR, "_MX_strict.ioe")
    params:
        oprefix = ODIR + "/"
    conda:
        "envs/preprocessing.yaml"
    threads: 1
    shell:
        """
        mkdir -p {ODIR}
        suppa.py generateEvents -i {input} -o {params.oprefix} -f ioe -e SE SS RI FL MX
        """

rule merge_ioes:
    input:
        ioe_se = pjoin(ODIR, "_SE_strict.ioe"),
        ioe_a3 = pjoin(ODIR, "_A3_strict.ioe"),
        ioe_a5 = pjoin(ODIR, "_A5_strict.ioe"),
        ioe_ri = pjoin(ODIR, "_RI_strict.ioe"),
        ioe_af = pjoin(ODIR, "_AF_strict.ioe"),
        ioe_al = pjoin(ODIR, "_AL_strict.ioe"),
        ioe_mx = pjoin(ODIR, "_MX_strict.ioe")
    output:
        pjoin(ODIR, "strict.ioe")
    shell:
        """
        head -1 {input.ioe_se} > {output}
        cat {input} | grep -v "seqname" >> {output}
        """

rule extract_local:
    input:
        gtf = GTF,
        ioe_se = pjoin(ODIR, "_SE_strict.ioe"),
        ioe_a3 = pjoin(ODIR, "_A3_strict.ioe"),
        ioe_a5 = pjoin(ODIR, "_A5_strict.ioe"),
        ioe_ir = pjoin(ODIR, "_RI_strict.ioe"),
        ioe_af = pjoin(ODIR, "_AF_strict.ioe"),
        ioe_al = pjoin(ODIR, "_AL_strict.ioe"),
        ioe_mx = pjoin(ODIR, "_MX_strict.ioe")
    output:
        gtf = pjoin(ODIR, "strict.ioe.local.gtf")
    params:
        odir = ODIR
    conda:
        "envs/preprocessing.yaml"
    threads: 1
    shell:
        """
        python3 extract_subannotations.py {input.gtf} {params.odir} > {output.gtf}
        """

rule get_local_transcripts:
    input:
        fa = FA,
        gtf = pjoin(ODIR, "strict.ioe.local.gtf")
    output:
        fa = pjoin(ODIR, "strict.ioe.local.fa")
    conda:
        "envs/preprocessing.yaml"
    threads: 1
    shell:
        """
        gffread -g {input.fa} {input.gtf} -w {output.fa}
        """

rule shark:
    input:
        fa = pjoin(ODIR, "strict.ioe.local.fa"),
        fq = FQ + ".trimmed.fq"
    output:
        ssv = pjoin(ODIR, "shark.ssv"),
        fq = pjoin(ODIR, "sharked.fq")
    conda:
        "envs/preprocessing.yaml"
    threads: THREADS
    shell:
        """
        shark -r {input.fa} -1 {input.fq} -o {output.fq} -t {threads} > {output.ssv}
        """

rule split_sharked:
    input:
        ssv = pjoin(ODIR, "shark.ssv"),
        fq = pjoin(ODIR, "sharked.fq")
    output:
        directory(pjoin(ODIR, "sharked"))
    conda:
        "envs/preprocessing.yaml"
    threads: 1
    shell:
        """
        python3 split_sharkedfq.py {input.ssv} {input.fq} {output}
        """
