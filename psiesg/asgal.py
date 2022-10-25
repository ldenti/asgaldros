import os
import subprocess
import pysam


def compare(aln1, aln2):
    return aln1.query_name == aln2.query_name and aln1.cigarstring == aln2.cigarstring


def clean(alns):
    clustered_alns = {}
    for aln in alns:
        if aln.query_name not in clustered_alns:
            clustered_alns[aln.query_name] = []
        clustered_alns[aln.query_name].append(aln)
    cleaned_alns = []
    for _, alns in clustered_alns.items():
        curr_alns = []
        for aln1 in alns:
            keep = True
            for aln2 in curr_alns:
                keep = not compare(aln1, aln2)
                if not keep:
                    break
            if keep:
                curr_alns.append(aln1)
        cleaned_alns += curr_alns
    return cleaned_alns


def remove_duplicates(inbam_path, outbam_path):
    inbam = pysam.AlignmentFile(inbam_path, "rb")
    outbam = pysam.AlignmentFile(outbam_path, "wb", header=inbam.header)
    last_rname = ""
    last_pos = -1
    alns = []
    for aln in inbam.fetch():
        rname = aln.reference_name
        pos = aln.pos  # 0-based
        if last_rname != rname or last_pos != pos:
            alns = clean(alns)
            for a in alns:
                outbam.write(a)
            alns = []
        alns.append(aln)
        last_pos = pos
        last_rname = rname
    alns = clean(alns)
    for a in alns:
        outbam.write(a)
    outbam.close()


def flag_secondary(inbam_path, outbam_path):
    inbam = pysam.AlignmentFile(inbam_path, "rb")
    outbam = pysam.AlignmentFile(outbam_path, "wb", header=inbam.header)
    clustered_alns = {}
    for aln in inbam.fetch():
        if aln.query_name not in clustered_alns:
            clustered_alns[aln.query_name] = []
        clustered_alns[aln.query_name].append(aln)

    for _, alns in clustered_alns.items():
        best_aln = None
        best_score = float("inf")
        for aln in alns:
            nm = aln.get_tag("NM")
            if nm < best_score:
                best_score = nm
                best_aln = aln
        for aln in alns:
            if aln != best_aln:
                aln.flag = 256 if aln.is_forward else 272
            outbam.write(aln)
    outbam.close()


def run(chrom_dir, esg_wd, event, spliceawarealigner, formatsam, l):
    prefix = os.path.join(esg_wd, f"{event}")
    chrom = open(prefix + ".gtf").readline().split("\t")[0]
    fa = os.path.join(chrom_dir, f"{chrom}.fa")
    asgal1_cmd = [
        spliceawarealigner,
        "-l",
        str(l),
        "-g",
        fa,
        "-a",
        prefix + ".gtf",
        "-s",
        prefix + ".fq",
        "-o",
        prefix + ".mem",
    ]
    asgal1_cmd_str = " ".join(asgal1_cmd) + " &> " + prefix + ".mem.log"
    asgal2_cmd = formatsam.split(" ") + [
        "-m",
        prefix + ".mem",
        "-g",
        fa,
        "-a",
        prefix + ".gtf",
        "-o",
        prefix + ".sam",
    ]
    asgal2_cmd_str = " ".join(asgal2_cmd) + " &> " + prefix + ".sam.log"
    subprocess.run(
        asgal1_cmd,
        stdout=open(prefix + ".mem.log", "w"),
        stderr=subprocess.STDOUT,
    )
    subprocess.run(
        asgal2_cmd,
        stdout=open(prefix + ".mem.log", "w"),
        stderr=subprocess.STDOUT,
    )
    view_p = subprocess.run(
        ["samtools", "view", "-bS", f"{prefix}.sam"],
        stdout=subprocess.PIPE,
        check=True,
    )
    subprocess.run(["samtools", "sort", "-o", f"{prefix}.bam"], input=view_p.stdout)
    subprocess.run(["samtools", "index", f"{prefix}.bam"])

    with open(prefix + ".cmd.log", "w") as cmdlog:
        cmdlog.write(asgal1_cmd_str + "\n")
        cmdlog.write(asgal2_cmd_str + "\n")
