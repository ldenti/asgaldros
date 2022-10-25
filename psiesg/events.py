import os
import subprocess
from rich.progress import (
    Progress,
    SpinnerColumn,
    TextColumn,
    BarColumn,
    TaskProgressColumn,
)


def run_suppa(GTF, WD, events):
    log_path = os.path.join(WD, "suppa2.log")
    cmd = [
        "suppa.py",
        "generateEvents",
        "-i",
        GTF,
        "-f",
        "ioe",
        "-o",
        WD + "/",
        "-e",
    ] + events
    p = subprocess.run(cmd, stderr=open(log_path, "w"))
    return p.returncode, log_path


def build_esg(exons, suppa2_wd, esg_gtf_path):
    esg_gtf = open(esg_gtf_path, "w")
    events = ["SE", "A3", "A5", "RI", "MX", "AF", "AL"]
    with Progress(
        SpinnerColumn(),
        TextColumn(f"[bold blue]Building ESG..", justify="right"),
        BarColumn(),
        TaskProgressColumn(),
    ) as progress:
        task_id = progress.add_task("Working...", total=len(events))
        for event in events:
            analyze(exons, suppa2_wd, esg_gtf, event)
            progress.advance(task_id)
    esg_gtf.close()


def analyze(exons, suppa2_wd, esg_gtf, event):
    if event == "SE":
        analyze_SE(exons, os.path.join(suppa2_wd, f"_SE_strict.ioe"), esg_gtf)
    elif event == "RI":
        analyze_RI(os.path.join(suppa2_wd, "_RI_strict.ioe"), esg_gtf)
    elif event == "MX":
        analyze_MX(exons, os.path.join(suppa2_wd, "_MX_strict.ioe"), esg_gtf)
    elif event in ["A3", "A5"]:
        analyze_SS(exons, os.path.join(suppa2_wd, f"_{event}_strict.ioe"), esg_gtf)
    elif event in ["AF", "AL"]:
        analyze_FL(exons, os.path.join(suppa2_wd, f"_{event}_strict.ioe"), esg_gtf)


def write_gene(chrom, idx, begin, end, Ts, ogtf):
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
        file=ogtf,
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
            file=ogtf,
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
                file=ogtf,
            )


def analyze_SE(exons, fpath, ogtf):
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
        write_gene(chrom, idx, begin, end, Ts, ogtf)


def analyze_SS(exons, fpath, ogtf):
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

        write_gene(chrom, idx, begin, end, Ts, ogtf)


def analyze_RI(fpath, ogtf):
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
        write_gene(chrom, idx, begin, end, Ts, ogtf)


def analyze_FL(exons, fpath, ogtf):
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

        write_gene(chrom, idx, begin, end, Ts, ogtf)


def analyze_MX(exons, fpath, ogtf):
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
        write_gene(chrom, idx, begin, end, Ts, ogtf)
