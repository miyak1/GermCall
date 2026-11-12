#!/usr/bin/env python3
"""GermCall run summary: reference, raw reads, trimming, mapping, variants, filters.

Called by `germcall.sh --steps summary`; reads what the earlier stages left in the
working directory and writes, into <workdir>/summary/:

    run_summary.txt        the whole report in plain text (read this first)
    reference_summary.tsv  genome size, contigs, N50, GC, N
    sample_summary.tsv     one row per sample: raw reads, trimming, mapping, depth, genotypes
    variant_summary.tsv    SNP / indel counts: raw joint call, after hard filters
    filter_criteria.tsv    every hard-filter threshold and how many variants it removed
    variants_per_contig.tsv  raw variant records per chromosome / contig

Every section is optional: if a stage has not been run, that section says so.
Only the Python standard library (plus PyYAML, already needed by cfg.py) is used.
"""
import argparse, collections, datetime, glob, gzip, os, re, statistics, subprocess, sys, zipfile

try:
    import yaml
except ImportError:
    sys.exit("ERROR: PyYAML is required (conda install pyyaml)")

# directory names: keep identical to germcall.sh
DIRS = dict(fqc="00_fastqc_raw", trim="01_trimmed", bam="02_bam", gvcf="03_gvcf",
            joint="04_joint", filt="05_filtered")


# ------------------------------------------------------------------ helpers
def run(cmd, check=True):
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if check and p.returncode != 0:
        raise RuntimeError("%s failed: %s" % (" ".join(cmd[:3]), p.stderr.strip()[:300]))
    return p.stdout


def first_line(cmd, pat=None):
    try:
        p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, errors="replace", timeout=120)
        for ln in p.stdout.splitlines():
            if ln.strip() and (pat is None or re.search(pat, ln)):
                return ln.strip()
    except Exception:
        pass
    return "unknown"


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def fmt(n, d=0):
    if n is None:
        return "NA"
    return ("{:,.%df}" % d).format(n)


def pct(a, b, d=2):
    return None if not b else 100.0 * a / b


def stats_line(vals, d=2):
    vals = [v for v in vals if v is not None]
    if not vals:
        return "NA"
    return "mean %s | median %s | min %s | max %s" % (
        fmt(statistics.mean(vals), d), fmt(statistics.median(vals), d), fmt(min(vals), d), fmt(max(vals), d))


def write_tsv(path, header, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join("" if v is None else (("%g" % v) if isinstance(v, float) else str(v)) for v in r) + "\n")


def cfg_get(conf, key, default):
    node = conf
    for part in key.split("."):
        if not isinstance(node, dict) or part not in node or node[part] in (None, ""):
            return default
        node = node[part]
    return node


# ------------------------------------------------------------------ reference
def reference_info(ref, fai_path):
    contigs = []
    with open(fai_path) as fh:
        for ln in fh:
            f = ln.split("\t")
            contigs.append((f[0], int(f[1])))
    total = sum(l for _, l in contigs)
    # N50
    acc, n50 = 0, 0
    for _, l in sorted(contigs, key=lambda x: -x[1]):
        acc += l
        if acc >= total / 2:
            n50 = l
            break
    gc = at = nn = 0
    opener = gzip.open if ref.endswith(".gz") else open
    with opener(ref, "rt") as fh:
        for ln in fh:
            if ln.startswith(">"):
                continue
            u = ln.strip().upper()
            gc += u.count("G") + u.count("C")
            at += u.count("A") + u.count("T")
            nn += len(u) - u.count("G") - u.count("C") - u.count("A") - u.count("T")
    acgt = gc + at
    return dict(file=os.path.basename(ref), contigs=len(contigs), total=total, n50=n50,
                longest=max(l for _, l in contigs), gc=pct(gc, acgt), n_pct=pct(nn, gc + at + nn),
                list=contigs)


# ------------------------------------------------------------------ raw reads (FastQC)
def fastqc_info(fqc_dir, s):
    """Total sequences / length / %GC of read 1 from the FastQC zip, plus read 2 length."""
    out = {}
    zips = sorted(glob.glob(os.path.join(fqc_dir, s, "*_fastqc.zip")))
    for i, z in enumerate(zips[:2]):
        try:
            with zipfile.ZipFile(z) as zf:
                name = [n for n in zf.namelist() if n.endswith("fastqc_data.txt")][0]
                txt = zf.read(name).decode("utf-8", "replace")
        except Exception:
            continue
        d = {}
        for ln in txt.splitlines():
            if ln.startswith(("Total Sequences", "Sequence length", "%GC", "Sequences flagged")):
                k, _, v = ln.partition("\t")
                d[k] = v
        out["R%d" % (i + 1)] = d
    return out


# ------------------------------------------------------------------ Trim Galore reports
def parse_trim_report(path):
    d = {}
    with open(path, errors="replace") as fh:
        txt = fh.read()

    def grab(pat, conv=lambda x: int(x.replace(",", ""))):
        m = re.search(pat, txt)
        return conv(m.group(1)) if m else None
    d["reads"] = grab(r"Total reads processed:\s+([\d,]+)")
    d["with_adapter"] = grab(r"Reads with adapters:\s+([\d,]+)")
    d["bp_in"] = grab(r"Total basepairs processed:\s+([\d,]+) bp")
    d["bp_qtrim"] = grab(r"Quality-trimmed:\s+([\d,]+) bp")
    d["bp_out"] = grab(r"Total written \(filtered\):\s+([\d,]+) bp")
    d["pairs_removed"] = grab(r"sequence pairs removed because at least one read was shorter than the length cutoff \(\d+ bp\):\s+([\d,]+)")
    d["reads_removed"] = grab(r"Sequences removed because they became shorter than the length cutoff of \d+ bp:\s+([\d,]+)")
    d["pairs_total"] = grab(r"Total number of sequence pairs processed:\s+([\d,]+)")
    return d


def trim_info(trim_dir, s):
    reps = sorted(glob.glob(os.path.join(trim_dir, s, "*trimming_report.txt")))
    if not reps:
        return None
    parsed = [parse_trim_report(r) for r in reps]
    reads = parsed[0]["reads"]
    bp_in = sum(p["bp_in"] or 0 for p in parsed)
    bp_out = sum(p["bp_out"] or 0 for p in parsed)
    removed = None
    for p in parsed:
        if p["pairs_removed"] is not None:
            removed = p["pairs_removed"]
        elif p["reads_removed"] is not None and removed is None:
            removed = p["reads_removed"]
    return dict(reads=reads, bp_in=bp_in, bp_out=bp_out, removed=removed,
                adapter_pct=pct(sum(p["with_adapter"] or 0 for p in parsed), sum(p["reads"] or 0 for p in parsed)))


# ------------------------------------------------------------------ mapping (samtools)
def mapping_info(bam_dir, s):
    f = os.path.join(bam_dir, "stats", s + ".flagstat.txt")
    c = os.path.join(bam_dir, "stats", s + ".coverage.txt")
    if not (os.path.isfile(f) and os.path.isfile(c)):
        return None
    d = dict(total=None, mapped=None, proper=None, dup=None, mapped_pct=None, proper_pct=None)
    for ln in open(f):
        n = re.match(r"(\d+) \+ \d+ (.*)", ln)
        if not n:
            continue
        cnt, what = int(n.group(1)), n.group(2)
        if what.startswith("in total"):
            d["total"] = cnt
        elif what.startswith("duplicates"):
            d["dup"] = cnt
        elif what.startswith("mapped ("):
            d["mapped"] = cnt
        elif what.startswith("properly paired"):
            d["proper"] = cnt
    d["mapped_pct"] = pct(d["mapped"], d["total"])
    d["proper_pct"] = pct(d["proper"], d["total"])
    d["dup_pct"] = pct(d["dup"], d["total"])
    w = L = cov = 0
    for i, ln in enumerate(open(c)):
        if i == 0 or not ln.strip():
            continue
        p = ln.rstrip("\n").split("\t")
        length = int(p[2]) - int(p[1]) + 1
        w += float(p[6]) * length
        cov += float(p[4])  # covered bases
        L += length
    d["depth"] = w / L if L else None
    d["breadth_pct"] = pct(cov, L)
    return d


# ------------------------------------------------------------------ bcftools stats
def ensure_stats(bcftools, vcf, stats):
    if os.path.isfile(vcf) and (not os.path.isfile(stats) or os.path.getmtime(stats) < os.path.getmtime(vcf)):
        with open(stats, "w") as fh:
            fh.write(run([bcftools, "stats", "-s", "-", vcf]))
    return os.path.isfile(stats)


def parse_stats(stats):
    sn, tstv, psc = {}, None, {}
    psc_cols = []
    with open(stats) as fh:
        for ln in fh:
            if ln.startswith("# PSC"):
                psc_cols = [re.sub(r"^\[\d+\]", "", x) for x in ln.rstrip("\n").split("\t")[1:]]
                continue
            if ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if f[0] == "SN":
                sn[f[2].rstrip(":").replace("number of ", "")] = int(f[3])
            elif f[0] == "TSTV":
                ts_, tv_ = int(f[2]), int(f[3]); tstv = dict(ts=ts_, tv=tv_, ratio=(ts_ / tv_ if tv_ else None))
            elif f[0] == "PSC" and psc_cols:
                row = dict(zip(psc_cols, f[1:]))
                psc[row["sample"]] = row
    return sn, tstv, psc


def filter_counts(bcftools, vcf):
    """How many sites carry each FILTER name, and how many fail at least one."""
    c, failed, total = collections.Counter(), 0, 0
    p = subprocess.Popen([bcftools, "query", "-f", "%FILTER\n", vcf], stdout=subprocess.PIPE, text=True)
    for ln in p.stdout:
        total += 1
        v = ln.strip()
        if v not in ("PASS", "."):
            failed += 1
            for name in v.split(";"):
                c[name] += 1
    p.wait()
    return total, failed, c


def per_contig(bcftools, vcf):
    out = []
    for ln in run([bcftools, "index", "-s", vcf]).splitlines():
        f = ln.split("\t")
        out.append((f[0], int(f[2])))
    return out


def chromosome_names(conf, contigs):
    """Chromosomes as germcall.sh defines them: calling.chromosomes, else the N longest contigs."""
    listed = cfg_get(conf, "calling.chromosomes", [])
    if isinstance(listed, list) and listed:
        return [str(c) for c in listed]
    n = cfg_get(conf, "calling.n_chromosomes", "")
    try:
        n = int(n)
    except (TypeError, ValueError):
        return []
    return [c for c, _ in sorted(contigs, key=lambda x: -x[1])[:n]]


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser()
    for a in ("config", "workdir", "prefix", "ref", "samples", "bcftools", "samtools", "bwa", "gatk", "trimgalore", "out"):
        ap.add_argument("--" + a, required=True)
    a = ap.parse_args()
    conf = yaml.safe_load(open(a.config)) or {}
    W = a.workdir
    d = {k: os.path.join(W, v) for k, v in DIRS.items()}
    os.makedirs(a.out, exist_ok=True)
    samples = [x.strip() for x in open(a.samples) if x.strip()]
    paired = str(cfg_get(conf, "calling.paired", True)).lower() == "true"
    rep, warn = [], []
    P = rep.append

    def head(t):
        P("")
        P("=" * 78)
        P(t)
        P("=" * 78)

    P("GermCall run summary")
    P("Project : %s" % cfg_get(conf, "project.name", a.prefix))
    P("Created : %s" % datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
    P("Samples : %d (%s) | ploidy used for calling: %s" % (len(samples), "paired-end" if paired else "single-end", cfg_get(conf, "calling.ploidy", 2)))
    P("Tools   : %s | %s | %s" % (
        first_line([a.bwa], r"Version").replace("Version:", "bwa"),
        first_line([a.samtools, "--version"]), first_line([a.bcftools, "--version"])))
    P("          %s | %s" % (first_line([a.gatk, "--version"], r"Genome Analysis Toolkit"),
                              re.sub(r"(?i)version", "Trim Galore", first_line([a.trimgalore, "--version"], r"(?i)version"))))

    # ---------------- 1. reference
    head("1. REFERENCE GENOME")
    fai = a.ref + ".fai"
    refrows, chroms = [], []
    if os.path.isfile(a.ref) and os.path.isfile(fai):
        r = reference_info(a.ref, fai)
        chroms = chromosome_names(conf, r["list"])
        n_chr = len(chroms)
        chr_len = sum(l for c, l in r["list"] if c in chroms)
        P("File                : %s" % r["file"])
        P("Total length        : %s bp (%.1f Mb)" % (fmt(r["total"]), r["total"] / 1e6))
        P("Contigs / scaffolds : %s   (longest %s bp, N50 %s bp)" % (fmt(r["contigs"]), fmt(r["longest"]), fmt(r["n50"])))
        P("GC content          : %s %%   |   unknown bases (N): %s %%" % (fmt(r["gc"], 2), fmt(r["n_pct"], 2)))
        if n_chr:
            P("Chromosomes used    : %d contigs = %.1f Mb (%.1f %% of the genome); the other %d are unplaced"
              % (n_chr, chr_len / 1e6, 100.0 * chr_len / r["total"], r["contigs"] - n_chr))
        refrows = [("file", r["file"]), ("total_length_bp", r["total"]), ("contigs", r["contigs"]),
                   ("longest_contig_bp", r["longest"]), ("N50_bp", r["n50"]), ("GC_pct", round(r["gc"], 2)),
                   ("N_pct", round(r["n_pct"], 2)), ("chromosomes_n", n_chr), ("chromosomes_bp", chr_len)]
        write_tsv(os.path.join(a.out, "reference_summary.tsv"), ["item", "value"], refrows)
    else:
        P("(reference or its .fai not found: run stage prepare)")

    # ---------------- per-sample collection
    rows, S = [], {}
    for s in samples:
        S[s] = dict(fq=fastqc_info(d["fqc"], s), tr=trim_info(d["trim"], s), mp=mapping_info(d["bam"], s))

    # ---------------- 2. raw reads
    head("2. RAW READS (FastQC + Trim Galore reports)")
    raw_pairs = {}
    for s in samples:
        fq = S[s]["fq"].get("R1", {})
        n = num(fq.get("Total Sequences"))
        if n is None and S[s]["tr"]:
            n = S[s]["tr"]["reads"]
        raw_pairs[s] = int(n) if n else None
    have = [v for v in raw_pairs.values() if v]
    unit = "read pairs" if paired else "reads"
    if have:
        P("Raw %s per sample : %s" % (unit, stats_line(have, 0)))
        P("Total raw %s      : %s" % (unit, fmt(sum(have))))
        bp = [S[s]["tr"]["bp_in"] for s in samples if S[s]["tr"] and S[s]["tr"]["bp_in"]]
        if bp:
            P("Raw bases per sample : %s Gb (mean) | total %s Gb"
              % (fmt(statistics.mean(bp) / 1e9, 2), fmt(sum(bp) / 1e9, 1)))
        lens = sorted({S[s]["fq"]["R1"].get("Sequence length") for s in samples if "R1" in S[s]["fq"]} - {None})
        if lens:
            P("Read length          : %s" % ", ".join(lens[:6]) + (" ..." if len(lens) > 6 else ""))
        else:  # no FastQC zip: estimate the mean read length from the Trim Galore report
            ml = [t["bp_in"] / (t["reads"] * (2 if paired else 1)) for t in (S[s]["tr"] for s in samples)
                  if t and t["reads"] and t["bp_in"]]
            if ml:
                P("Read length (mean)   : %s bp" % stats_line(ml, 0))
        gcs = [num(S[s]["fq"]["R1"].get("%GC")) for s in samples if "R1" in S[s]["fq"]]
        gcs = [g for g in gcs if g is not None]
        if gcs:
            P("Read GC %%            : %s" % stats_line(gcs, 1))
        P("Samples with < 5 M %s: %d" % (unit, sum(1 for v in have if v < 5e6)))
        if len(have) < len(samples):
            warn.append("raw read counts missing for %d sample(s) (run stage qc / trim_align)" % (len(samples) - len(have)))
    else:
        P("(no FastQC / Trim Galore results found: run stages qc and trim_align)")

    # ---------------- 3. trimming
    head("3. TRIMMING (Trim Galore: quality %s, minimum length %s)" % (
        cfg_get(conf, "calling.trim_quality", 20), cfg_get(conf, "calling.trim_min_length", 50)))
    trs = [S[s]["tr"] for s in samples if S[s]["tr"]]
    if trs:
        keep = [pct(t["bp_out"], t["bp_in"]) for t in trs if t["bp_in"]]
        ad = [t["adapter_pct"] for t in trs if t["adapter_pct"] is not None]
        rm = [pct(t["removed"], t["reads"]) for t in trs if t["removed"] is not None and t["reads"]]
        P("Bases kept after trimming      : %s %%" % stats_line(keep, 1))
        P("Reads with adapter found       : %s %%" % stats_line(ad, 1))
        P("%s dropped (too short)   : %s %%" % ("Pairs" if paired else "Reads", stats_line(rm, 2)))
        bad = [s_ for s_ in samples if S[s_]["tr"] and S[s_]["tr"]["removed"] is not None and S[s_]["tr"]["reads"]
               and 100.0 * S[s_]["tr"]["removed"] / S[s_]["tr"]["reads"] > 30]
        if bad:
            warn.append("%d sample(s) lost > 30 %% of their %s as too short after trimming (very short fragments or a poor library): %s"
                        % (len(bad), "pairs" if paired else "reads", ", ".join(bad[:8]) + (" ..." if len(bad) > 8 else "")))
    else:
        P("(no Trim Galore reports found)")

    # ---------------- 4. mapping
    head("4. MAPPING (BWA-MEM, duplicates marked by GATK MarkDuplicates)")
    mps = {s: S[s]["mp"] for s in samples if S[s]["mp"]}
    if mps:
        dep = [m["depth"] for m in mps.values()]
        P("Mapped reads         : %s %%" % stats_line([m["mapped_pct"] for m in mps.values()], 2))
        if paired:
            P("Properly paired      : %s %%" % stats_line([m["proper_pct"] for m in mps.values()], 2))
        P("Duplicates           : %s %%" % stats_line([m["dup_pct"] for m in mps.values()], 2))
        P("Mean depth           : %s x" % stats_line(dep, 2))
        P("Genome covered >= 1x : %s %%" % stats_line([m["breadth_pct"] for m in mps.values()], 1))
        for t in (5, 10, 20):
            P("Samples with depth < %2dx : %d of %d" % (t, sum(1 for x in dep if x < t), len(dep)))
        low = sorted(mps, key=lambda k: mps[k]["depth"])[:5]
        P("Lowest depth         : " + ", ".join("%s (%.1fx)" % (k, mps[k]["depth"]) for k in low))
        n5 = sum(1 for x in dep if x < 5)
        if n5:
            warn.append("%d sample(s) have mean depth < 5x: many uncertain / missing genotypes, "
                        "check their missing-data rate in section 7 before GWAS" % n5)
        lowpp = [k for k, m in mps.items() if paired and m["proper_pct"] is not None and m["proper_pct"] < 50]
        if lowpp:
            warn.append("%d sample(s) have < 50 %% properly paired reads (short or degraded fragments, or a library problem): %s"
                        % (len(lowpp), ", ".join(lowpp[:8]) + (" ..." if len(lowpp) > 8 else "")))
        lowmap = [k for k, m in mps.items() if m["mapped_pct"] is not None and m["mapped_pct"] < 90]
        if lowmap:
            warn.append("%d sample(s) map < 90 %% (contamination or wrong reference?): %s"
                        % (len(lowmap), ", ".join(lowmap[:8]) + (" ..." if len(lowmap) > 8 else "")))
        if len(mps) < len(samples):
            warn.append("mapping statistics missing for %d sample(s)" % (len(samples) - len(mps)))
    else:
        P("(no mapping statistics found: run stage trim_align)")

    # ---------------- 5. raw variants
    raw_vcf = os.path.join(d["joint"], a.prefix + "_joint_raw.vcf.gz")
    raw_stats = os.path.join(d["joint"], a.prefix + "_joint_raw.stats.txt")
    vrows, psc_raw = [], {}
    head("5. VARIANT CALLING (GATK HaplotypeCaller -> joint genotyping, %s)" % cfg_get(conf, "calling.combine_method", "genomicsdb"))
    if ensure_stats(a.bcftools, raw_vcf, raw_stats):
        sn, tstv, psc_raw = parse_stats(raw_stats)
        P("Samples in the VCF            : %s" % fmt(sn.get("samples")))
        P("Variant sites (records)       : %s" % fmt(sn.get("records")))
        P("   SNPs                       : %s" % fmt(sn.get("SNPs")))
        P("   Indels                     : %s" % fmt(sn.get("indels")))
        P("   MNPs / other               : %s / %s" % (fmt(sn.get("MNPs")), fmt(sn.get("others"))))
        P("   Multi-allelic sites        : %s (SNP: %s)" % (fmt(sn.get("multiallelic sites")), fmt(sn.get("multiallelic SNP sites"))))
        if tstv and tstv["ratio"] is not None:
            P("Ts/Tv ratio (SNPs)            : %.2f   (random errors give ~0.5; real SNPs are clearly higher)" % tstv["ratio"])
            if tstv["ratio"] < 1.5:
                warn.append("raw Ts/Tv is %.2f (low): the raw call set contains many false SNPs, filtering is needed" % tstv["ratio"])
        if sn.get("samples") and sn["samples"] != len(samples):
            warn.append("VCF has %d samples but the sample table has %d" % (sn["samples"], len(samples)))
        if refrows and sn.get("SNPs"):
            P("SNP density (raw)             : 1 SNP per %.0f bp of genome" % (refrows[1][1] / sn["SNPs"]))
        vrows += [("raw_joint_call", "records", sn.get("records")), ("raw_joint_call", "SNPs", sn.get("SNPs")),
                  ("raw_joint_call", "indels", sn.get("indels")), ("raw_joint_call", "MNPs", sn.get("MNPs")),
                  ("raw_joint_call", "other", sn.get("others")),
                  ("raw_joint_call", "multiallelic_sites", sn.get("multiallelic sites")),
                  ("raw_joint_call", "Ts/Tv", tstv["ratio"] if tstv else None)]
        try:
            pc = per_contig(a.bcftools, raw_vcf)
            big = [(c, n) for c, n in pc if c in set(chroms)] if refrows else []
            lens = dict(r["list"]) if refrows else {}
            write_tsv(os.path.join(a.out, "variants_per_contig.tsv"), ["contig", "length_bp", "raw_records", "records_per_Mb"],
                      [(c, lens.get(c), n, round(n / (lens[c] / 1e6), 1) if lens.get(c) else None) for c, n in pc])
            if big:
                P("")
                P("Raw variant records per chromosome:")
                for c, n in big:
                    P("   %-14s %12s records   %8.1f per Mb" % (c, fmt(n), n / (lens[c] / 1e6)))
                rest = sum(n for c, n in pc) - sum(n for _, n in big)
                P("   %-14s %12s records   (all unplaced contigs)" % ("other", fmt(rest)))
        except Exception as e:
            P("(per-contig counts unavailable: %s)" % e)
    else:
        P("(raw joint VCF not found: run stage joint)")

    # ---------------- 6. filtering
    head("6. HARD FILTERING (GATK generic recommendations; run with --steps filter)")
    fcrit = []
    snp_pass = os.path.join(d["filt"], a.prefix + "_snps_pass.vcf.gz")
    ind_pass = os.path.join(d["filt"], a.prefix + "_indels_pass.vcf.gz")
    sets = [("SNPs", "snps",
             [("QD", "<", "filter_QD", 2.0), ("QUAL", "<", "filter_QUAL", 30.0), ("SOR", ">", "filter_SOR", 3.0),
              ("FS", ">", "filter_FS", 60.0), ("MQ", "<", "filter_MQ", 40.0),
              ("MQRankSum", "<", "filter_MQRankSum", -12.5), ("ReadPosRankSum", "<", "filter_ReadPosRankSum", -8.0)]),
            ("Indels", "indels",
             [("QD", "<", "filter_indel_QD", 2.0), ("QUAL", "<", "filter_indel_QUAL", 30.0),
              ("FS", ">", "filter_indel_FS", 200.0), ("ReadPosRankSum", "<", "filter_indel_ReadPosRankSum", -20.0)])]
    meaning = {"QD": "variant quality / depth (low = weak evidence per read)", "QUAL": "Phred-scaled call quality",
               "SOR": "strand bias (symmetric odds ratio)", "FS": "strand bias (Fisher's test)",
               "MQ": "mapping quality of reads at the site", "MQRankSum": "ALT reads map worse than REF reads",
               "ReadPosRankSum": "ALT alleles sit near read ends"}
    PASSV, BEFORE = {}, {}
    for label, key, crit in sets:
        flt = os.path.join(d["filt"], "%s_%s_filtered.vcf.gz" % (a.prefix, key))
        passv = snp_pass if key == "snps" else ind_pass
        P("")
        P("%s - a site is REMOVED if any of these is true:" % label)
        for name, op, ck, dflt in crit:
            val = cfg_get(conf, "calling." + ck, dflt)
            P("   %-15s %s %-7s  %s" % (name, op, val, meaning[name]))
            fcrit.append([label, name, "%s %s" % (op, val), meaning[name], None])
        if os.path.isfile(flt) and os.path.isfile(passv):
            tot, fail, cnt = filter_counts(a.bcftools, flt)
            P("   -> %s before filters: %s | removed: %s (%.1f %%) | PASSING: %s"
              % (label, fmt(tot), fmt(fail), 100.0 * fail / tot if tot else 0, fmt(tot - fail)))
            P("   -> sites failing each filter (one site can fail several): " +
              (", ".join("%s %s" % (k, fmt(cnt.get(k, 0))) for k, _, _, _ in crit)))
            for row in fcrit:
                if row[0] == label:
                    row[4] = cnt.get(row[1], 0)
            vrows += [("hard_filter_" + key, "before_filters", tot), ("hard_filter_" + key, "removed", fail),
                      ("hard_filter_" + key, "passing", tot - fail)]
            PASSV[key] = passv
            BEFORE[key] = tot
            ps = os.path.join(d["filt"], "%s_%s_pass.stats.txt" % (a.prefix, key))
            ensure_stats(a.bcftools, passv, ps)
            if key == "snps" and os.path.isfile(ps):
                _, ts2, _ = parse_stats(ps)
                if ts2 and ts2["ratio"] is not None:
                    P("   -> Ts/Tv of passing SNPs: %.2f" % ts2["ratio"])
                    vrows.append(("hard_filter_snps", "Ts/Tv_passing", ts2["ratio"]))
                    if tot - fail and refrows:
                        P("   -> SNP density after filtering: 1 SNP per %.0f bp" % (refrows[1][1] / (tot - fail)))
        else:
            P("   (not run yet: --steps filter)")
    if len(BEFORE) == 2 and psc_raw is not None:
        try:
            nrec = parse_stats(raw_stats)[0].get("records")
            mixed = nrec - BEFORE["snps"] - BEFORE["indels"]
            if mixed > 0:
                P("")
                P("Sites in neither file: %s (%.1f %% of the raw records). These carry both a SNP and an indel allele"
                  % (fmt(mixed), 100.0 * mixed / nrec))
                P("(GATK calls them MIXED), so SelectVariants puts them in neither the SNP nor the indel file.")
        except Exception:
            pass
    P("")
    P("Not used here: genotype-level filters (per-genotype depth / quality, missing data, minor")
    P("allele frequency, HWE). GermGWAS applies those in its QC step (qc.min_dp, qc.min_gq, geno, maf, hwe).")
    write_tsv(os.path.join(a.out, "filter_criteria.tsv"), ["variant_type", "filter", "removes_if", "meaning", "sites_failing"], fcrit)
    write_tsv(os.path.join(a.out, "variant_summary.tsv"), ["stage", "metric", "value"], vrows)

    # ---------------- 7. per-sample genotype statistics
    head("7. PER-SAMPLE GENOTYPES (from the raw joint VCF%s)" % (", missing % also after SNP filtering" if "snps" in PASSV else ""))
    psc_pass = {}
    if "snps" in PASSV:
        ps = os.path.join(d["filt"], "%s_snps_pass.stats.txt" % a.prefix)
        if os.path.isfile(ps):
            psc_pass = parse_stats(ps)[2]
    gen = {}
    for sname, r_ in psc_raw.items():
        tot = sum(int(r_.get(k, 0)) for k in ("nRefHom", "nNonRefHom", "nHets", "nMissing"))
        called = tot - int(r_.get("nMissing", 0))
        gen[sname] = dict(
            missing=pct(int(r_.get("nMissing", 0)), tot), het=pct(int(r_.get("nHets", 0)), called),
            sing=int(r_.get("nSingletons", 0)))
        if sname in psc_pass:
            q = psc_pass[sname]
            t2 = sum(int(q.get(k, 0)) for k in ("nRefHom", "nNonRefHom", "nHets", "nMissing"))
            gen[sname]["missing_pass"] = pct(int(q.get("nMissing", 0)), t2)
    if gen:
        P("Missing genotypes (raw)   : %s %%" % stats_line([g["missing"] for g in gen.values()], 1))
        if psc_pass:
            P("Missing genotypes (SNP pass) : %s %%" % stats_line([g.get("missing_pass") for g in gen.values()], 1))
        P("Heterozygous calls        : %s %%" % stats_line([g["het"] for g in gen.values()], 2))
        worst = sorted(gen, key=lambda k: -(gen[k]["missing"] or 0))[:5]
        P("Highest missing           : " + ", ".join("%s (%.1f %%)" % (k, gen[k]["missing"]) for k in worst))
        bad = [k for k, g in gen.items() if g["missing"] is not None and g["missing"] > 30]
        if bad:
            warn.append("%d sample(s) have > 30 %% missing genotypes (%s%s): consider dropping them or re-sequencing"
                        % (len(bad), ", ".join(bad[:6]), " ..." if len(bad) > 6 else ""))
        hets = [g["het"] for g in gen.values() if g["het"] is not None]
        if hets:
            hm = statistics.median(hets)
            odd = [k for k, g in gen.items() if g["het"] is not None and g["het"] > max(3 * hm, hm + 5)]
            if odd:
                warn.append("unusually high heterozygosity in %s (possible contamination or mixed samples)" % ", ".join(odd[:8]))
    else:
        P("(no genotype statistics: raw joint VCF missing)")

    # ---------------- sample table
    for s in samples:
        m, t, g = S[s]["mp"], S[s]["tr"], gen.get(s, {})
        rows.append([s, raw_pairs.get(s), (round(t["bp_in"] / 1e9, 3) if t and t["bp_in"] else None),
                     round(pct(t["bp_out"], t["bp_in"]), 1) if t and t["bp_in"] else None,
                     round(t["adapter_pct"], 1) if t and t["adapter_pct"] is not None else None,
                     m["total"] if m else None, round(m["mapped_pct"], 2) if m and m["mapped_pct"] is not None else None,
                     round(m["proper_pct"], 2) if m and m["proper_pct"] is not None else None,
                     round(m["dup_pct"], 2) if m and m["dup_pct"] is not None else None,
                     round(m["depth"], 2) if m else None, round(m["breadth_pct"], 1) if m else None,
                     round(g["missing"], 2) if g.get("missing") is not None else None,
                     round(g["missing_pass"], 2) if g.get("missing_pass") is not None else None,
                     round(g["het"], 2) if g.get("het") is not None else None, g.get("sing")])
    write_tsv(os.path.join(a.out, "sample_summary.tsv"),
              ["Sample", "Raw_read_pairs" if paired else "Raw_reads", "Raw_Gb", "Bases_kept_after_trim_pct",
               "Reads_with_adapter_pct", "Reads_after_trim", "Mapped_pct", "Properly_paired_pct",
               "Duplicate_pct", "Mean_depth", "Genome_covered_1x_pct", "Missing_GT_raw_pct", "Missing_GT_SNPpass_pct",
               "Het_pct", "Singletons"], rows)

    head("8. THINGS TO CHECK")
    if warn:
        for w in warn:
            P(" * " + w)
    else:
        P("No warnings.")
    P("")
    P("Tables written next to this report: sample_summary.tsv (one row per sample),")
    P("variant_summary.tsv, filter_criteria.tsv, variants_per_contig.tsv, reference_summary.tsv.")

    txt = "\n".join(rep) + "\n"
    open(os.path.join(a.out, "run_summary.txt"), "w").write(txt)
    print(txt)


if __name__ == "__main__":
    main()
