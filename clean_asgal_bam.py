import sys
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


def remove_duplicates():
    bam_path = sys.argv[1]
    bam = pysam.AlignmentFile(bam_path, "rb")
    last_rname = ""
    last_pos = -1
    alns = []
    print(bam.header, end="")
    for aln in bam.fetch():
        qname = aln.query_name
        rname = aln.reference_name
        pos = aln.pos  # 0-based
        if last_rname != rname or last_pos != pos:
            alns = clean(alns)
            for a in alns:
                print(a.tostring(bam))
            alns = []
        alns.append(aln)
        last_pos = pos
        last_rname = rname
    alns = clean(alns)
    for a in alns:
        print(a.tostring(bam))


def flag_secondary():
    bam_path = sys.argv[1]
    bam = pysam.AlignmentFile(bam_path, "rb")
    clustered_alns = {}
    for aln in bam.fetch():
        if aln.query_name not in clustered_alns:
            clustered_alns[aln.query_name] = []
        clustered_alns[aln.query_name].append(aln)

    print(bam.header, end="")
    for _, alns in clustered_alns.items():
        flagged_alns = []
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
            print(aln.tostring(bam))


if __name__ == "__main__":
    mode = sys.argv.pop(1)
    if mode == "rdup":
        remove_duplicates()
    elif mode == "fsec":
        flag_secondary()
