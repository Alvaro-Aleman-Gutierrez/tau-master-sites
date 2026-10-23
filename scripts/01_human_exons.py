"""Paso 1: péptidos por exón (humano) y presencia de cada exón en todas las isoformas de mamíferos.

Cada isoforma se alinea globalmente (extremos libres) contra una super-isoforma humana con todos
los exones codificantes en orden genómico; la cobertura de un exón es la fracción de sus residuos
alineados con residuos (no gaps) de la isoforma.
"""
import csv, json
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner, substitution_matrices

# exones codificantes humanos (GRCh38, NC_000017.11), nomenclatura canónica
EXONS = [("E1", 45962338, 45962470), ("E2", 45971859, 45971945), ("E3", 45974385, 45974471),
         ("E4", 45978375, 45978440), ("E4aL", 45982866, 45983177), ("E4a", 45983178, 45983930),
         ("E5", 45987040, 45987095), ("E6", 45989878, 45990075), ("E7", 45991460, 45991586),
         ("E8", 45993924, 45993977), ("E9", 45996399, 45996664), ("E10", 46010310, 46010402),
         ("E11", 46014243, 46014324), ("E12", 46018618, 46018730), ("E13", 46023956, 46024171)]
reports = json.load(open("data/raw/mammal_product_reports.json"))
prod = {str(r["product"]["tax_id"]): r["product"] for r in reports}
prot = {r.id: str(r.seq) for r in SeqIO.parse("data/raw/mammal_protein.faa", "fasta")}
meta = {r.id: r.description for r in SeqIO.parse("data/raw/mammal_protein.faa", "fasta")}
chrom = SeqIO.read("data/raw/human_MAPT_region.fasta", "fasta")   # NC_000017.11:45894000-46030000
OFF = 45894000

def exon_peptides():
    """Traduce cada exón en el marco de un transcrito humano que lo contiene."""
    peps = {}
    for t in prod["9606"]["transcripts"]:
        if not t["accession_version"].startswith("NM_"):
            continue
        ex = sorted((int(e["begin"]), int(e["end"])) for e in t["genomic_locations"][0]["exons"])
        c0 = int(t["cds"]["range"][0]["begin"]); c1 = int(t["cds"]["range"][0]["end"])
        pos, coding = 0, []        # posiciones genómicas codificantes en orden
        for b, e in ex:
            for g in range(b, e + 1):
                pos += 1
                if c0 <= pos <= c1:
                    coding.append(g)
        for name, b, e in EXONS:
            if name in peps or not (b in coding and e in coding):
                continue
            # codones cuya primera base cae en el exón
            first = [i for i in range(0, len(coding) - 2, 3) if b <= coding[i] <= e]
            aa = "".join(str(Seq("".join(chrom.seq[coding[i + k] - OFF] for k in range(3))).translate())
                         for i in first)
            peps[name] = aa.rstrip("*")
    return peps

peps = exon_peptides()
with open("results/human_exons.tsv", "w") as f:
    f.write("exon\tstart\tend\tnt\taa\tpeptide\n")
    for name, b, e in EXONS:
        f.write(f"{name}\t{b}\t{e}\t{e - b + 1}\t{len(peps.get(name, ''))}\t{peps.get(name, '')}\n")
super_iso = "".join(peps[n] for n, _, _ in EXONS)
bounds, p = {}, 0
for n, _, _ in EXONS:
    bounds[n] = (p, p + len(peps[n])); p += len(peps[n])

al = PairwiseAligner(mode="global", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                     open_gap_score=-10, extend_gap_score=-0.5)
al.target_end_gap_score = 0; al.query_end_gap_score = 0

names = [n for n, _, _ in EXONS]
with open("results/exon_presence_all_isoforms.tsv", "w") as f:
    f.write("protein\tspecies\tgene_id\tlength\tidentity_to_human\t" + "\t".join(names) + "\n")
    for pid, seq in prot.items():
        a = al.align(super_iso, seq)[0]
        cov = dict.fromkeys(names, 0); ident = 0; alnd = 0
        for (ts, te), (qs, qe) in zip(*a.aligned):
            for k in range(te - ts):
                i = ts + k
                n = next(n for n in names if bounds[n][0] <= i < bounds[n][1])
                cov[n] += 1; alnd += 1
                ident += super_iso[i] == seq[qs + k]
        d = meta[pid]
        sp = d.split("[organism=")[1].split("]")[0]; gid = d.split("[GeneID=")[1].split("]")[0]
        f.write(f"{pid}\t{sp}\t{gid}\t{len(seq)}\t{ident / max(alnd, 1):.3f}\t" +
                "\t".join(f"{cov[n] / (bounds[n][1] - bounds[n][0]):.2f}" for n in names) + "\n")
print(super_iso.__len__(), "aa super-isoforma;", len(prot), "isoformas evaluadas")
