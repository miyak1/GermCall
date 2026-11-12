#!/usr/bin/env python3
"""Generate a tiny SIMULATED data set for testing GermCall (no real data).

    python3 examples/make_test_data.py OUTDIR

Writes in OUTDIR:  ref.fa  (2 chromosomes + 2 scaffolds, 105 kb)
                   fastq/C01..C04 _R1/_R2.fastq.gz  (about 12x, 100 bp pairs, some adapter read-through)
                   samples.tsv
"""
import gzip, os, random, sys

random.seed(21)
out = sys.argv[1] if len(sys.argv) == 2 else sys.exit(__doc__)
os.makedirs(os.path.join(out, "fastq"), exist_ok=True)
rnd = lambda n: "".join(random.choice("ACGT") for _ in range(n))
contigs = {"chr1": 60000, "chr2": 40000, "scaf1": 3000, "scaf2": 2000}
ref = {k: rnd(v) for k, v in contigs.items()}
with open(os.path.join(out, "ref.fa"), "w") as f:
    for k, s in ref.items():
        f.write(f">{k}\n")
        for i in range(0, len(s), 60):
            f.write(s[i:i + 60] + "\n")
comp = str.maketrans("ACGT", "TGCA")
rc = lambda s: s.translate(comp)[::-1]
ADAPT = "AGATCGGAAGAGC"
samples = ["C01", "C02", "C03", "C04"]
snps = {k: sorted(random.sample(range(200, v - 200), 20)) for k, v in contigs.items()}
flip = {"A": "G", "G": "A", "C": "T", "T": "C"}
for sm in samples:
    hap = {}
    for k in ref:
        s = list(ref[k])
        for p in snps[k]:
            if random.random() < 0.5:
                s[p] = flip[s[p]]
        hap[k] = "".join(s)
    n = 0
    with gzip.open(os.path.join(out, "fastq", f"{sm}_R1.fastq.gz"), "wt") as r1, \
         gzip.open(os.path.join(out, "fastq", f"{sm}_R2.fastq.gz"), "wt") as r2:
        for k, s in hap.items():
            for _ in range(int(len(s) * 12 / 200)):
                ins = random.randint(150, 300); st = random.randint(0, len(s) - ins)
                frag = s[st:st + ins]; a = frag[:100]; b = rc(frag)[:100]
                if random.random() < 0.15:
                    a = a[:80] + ADAPT; b = b[:80] + ADAPT
                n += 1
                r1.write(f"@{sm}_{n}/1\n{a}\n+\n{'I' * len(a)}\n")
                r2.write(f"@{sm}_{n}/2\n{b}\n+\n{'I' * len(b)}\n")
with open(os.path.join(out, "samples.tsv"), "w") as f:
    f.write("# sample\tR1\tR2\n")
    for s in samples:
        f.write(f"{s}\t{os.path.abspath(out)}/fastq/{s}_R1.fastq.gz\t{os.path.abspath(out)}/fastq/{s}_R2.fastq.gz\n")
print("Test data written to", out)
