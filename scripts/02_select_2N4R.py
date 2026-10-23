"""Paso 2: una secuencia 2N4R (proteína + CDS) por especie.

1) Isoforma anotada con todos los exones core (cobertura >= 0.5) y sin 4aL/4a/6/8 (<= 0.2).
2) Si no existe: la isoforma más completa en exones core, a la que se le eliminan in silico los exones
   que cubren 4aL/4a/6/8, usando la estructura de exones del propio transcrito.
Preferencias: NP_ antes que XP_, luego mayor identidad con el humano.
"""
import csv, json
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner, substitution_matrices

CORE = ["E1", "E2", "E3", "E4", "E5", "E7", "E9", "E10", "E11", "E12", "E13"]
EXTRA = ["E4aL", "E4a", "E6", "E8"]
PRESENT, ABSENT = 0.5, 0.2

pres = pd.read_csv("results/exon_presence_all_isoforms.tsv", sep="\t")
taxa = pd.read_csv("data/raw/mammal_orthologs.tsv", sep="\t", dtype=str).set_index("gene_id")
prot = {r.id: str(r.seq) for r in SeqIO.parse("data/raw/mammal_protein.faa", "fasta")}
cds = {r.id.split(":")[0]: str(r.seq).upper() for r in SeqIO.parse("data/raw/mammal_cds.fna", "fasta")}
tx_of, tx_info = {}, {}
for rep in json.load(open("data/raw/mammal_product_reports.json")):
    for t in rep["product"]["transcripts"]:
        if "protein" in t:
            tx_of[t["protein"]["accession_version"]] = t["accession_version"]
            tx_info[t["accession_version"]] = t
hx = pd.read_csv("results/human_exons.tsv", sep="\t").set_index("exon")
super_iso = "".join(hx.loc[e, "peptide"] for e in hx.index)
bounds, p = {}, 0
for e in hx.index:
    bounds[e] = (p, p + len(hx.loc[e, "peptide"])); p += len(hx.loc[e, "peptide"])

al = PairwiseAligner(mode="global", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                     open_gap_score=-10, extend_gap_score=-0.5)
al.target_end_gap_score = 0; al.query_end_gap_score = 0

def human_exon_of_residues(seq):
    """Para cada residuo de seq, el exón humano con el que se alinea (o None)."""
    lab = [None] * len(seq)
    a = al.align(super_iso, seq)[0]
    for (ts, te), (qs, qe) in zip(*a.aligned):
        for k in range(te - ts):
            lab[qs + k] = next(e for e in bounds if bounds[e][0] <= ts + k < bounds[e][1])
    return lab

def excise(pid):
    """Elimina de la CDS los exones del transcrito que no corresponden a exones core humanos."""
    t = tx_info[tx_of[pid]]
    ex = sorted((int(e["begin"]), int(e["end"])) for e in t["genomic_locations"][0]["exons"])
    if t["genomic_locations"][0]["genomic_range"]["orientation"] == "minus":
        ex = ex[::-1]
    c0 = int(t["cds"]["range"][0]["begin"]); c1 = int(t["cds"]["range"][0]["end"])
    nt = cds[tx_of[pid]]
    lab = human_exon_of_residues(prot[pid])
    blocks, pos = [], 0                                    # (i0, i1, eliminar?)
    for b, e in ex:
        s, f = max(pos + 1, c0), min(pos + abs(e - b) + 1, c1)
        pos += abs(e - b) + 1
        if s > f:
            continue
        i0, i1 = s - c0, f - c0 + 1                       # rango en la CDS
        res = [lab[i // 3] for i in range(i0, i1, 3) if i // 3 < len(lab)]
        core = sum(r in CORE for r in res)
        first_or_last = i0 == 0 or i1 >= len(nt) - 3
        blocks.append((i0, i1, bool(res) and core / len(res) < 0.5 and not first_or_last))
    keep, removed = [], 0                                 # 4aL/4a/6/8 o exones no homólogos
    for i0, i1, drop in blocks:
        if drop:
            removed += i1 - i0
            continue
        if removed % 3:
            return None                                   # el bloque eliminado rompería el marco
        removed = 0
        keep.append(nt[i0:i1])
    new = "".join(keep)
    aa = str(Seq(new).translate()).rstrip("*")
    return (new, aa) if "*" not in aa else None

pres["np"] = pres.protein.str.startswith("NP_")
pres["n_core"] = (pres[CORE] >= PRESENT).sum(axis=1)
pres["natural"] = (pres[CORE] >= PRESENT).all(axis=1) & (pres[EXTRA] <= ABSENT).all(axis=1)

rows, out_p, out_n = [], [], []
for gid, g in pres.groupby("gene_id"):
    tx = taxa.loc[str(gid)]
    nat = g[g.natural].sort_values(["np", "identity_to_human"], ascending=False)
    seq = None
    if len(nat):
        pid = nat.iloc[0].protein; how = "anotada"
        nt = cds[tx_of[pid]]; aa = prot[pid]
        seq = (nt[:len(aa) * 3], aa)
    else:
        for _, c in g[g.n_core == g.n_core.max()].sort_values(["np", "identity_to_human"], ascending=False).iterrows():
            seq = excise(c.protein)
            if seq:
                pid = c.protein; how = "in_silico"; break
    if not seq:
        rows.append([tx.species, tx.superorder, tx.order, gid, "", "", "excluida", "", ""]); continue
    nt, aa = seq
    lab = human_exon_of_residues(aa)
    # inicio del primer tramo de >= 10 residuos alineados con el humano (evita alineamientos espurios)
    first = next(i for i in range(len(lab)) if all(l is not None for l in lab[i:i + 10]))
    mets = [i for i, c in enumerate(aa[:first + 6]) if c == "M"]
    if first > 10 and mets:                               # extensión N-terminal sin homología: se recorta
        cut = min(mets, key=lambda i: abs(i - first))     # Met más cercana al M1 homólogo humano
        aa, nt = aa[cut:], nt[cut * 3:]
        how += "+recorte_Nter"
        lab = human_exon_of_residues(aa)
    core_cov = sum(l in CORE for l in lab) / sum(bounds[e][1] - bounds[e][0] for e in CORE)
    missing = [e for e in CORE if sum(l == e for l in lab) / (bounds[e][1] - bounds[e][0]) < PRESENT]
    rows.append([tx.species, tx.superorder, tx.order, gid, pid, tx_of[pid], how, len(aa), f"{core_cov:.2f}", ",".join(missing)])
    name = tx.species.replace(" ", "_")
    out_p.append(f">{name} {pid}\n{aa}\n"); out_n.append(f">{name} {tx_of[pid]}\n{nt}\n")

cols = ["species", "superorder", "order", "gene_id", "protein", "transcript", "source", "length", "core_coverage", "missing_core_exons"]
with open("results/selected_2N4R.tsv", "w") as f:
    f.write("\t".join(cols) + "\n")
    for r in rows:
        f.write("\t".join(map(str, r + [""] * (len(cols) - len(r)))) + "\n")
open("results/mapt_2N4R_protein.faa", "w").write("".join(out_p))
open("results/mapt_2N4R_cds.fna", "w").write("".join(out_n))
