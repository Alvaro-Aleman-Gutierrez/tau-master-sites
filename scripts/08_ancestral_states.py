"""Paso 8: estados ancestrales de los sitios maestros sobre el árbol de 259 mamíferos.

IQ-TREE reconstruye la secuencia ancestral (marginal, empírico-bayesiano) de cada
nodo interno con el alineamiento, el árbol con restricción taxonómica y el modelo
elegido por ModelFinder en el paso 5, sin reestimar la topología. Luego se leen las
probabilidades posteriores de T en los nodos que definen la historia de cada sitio.

Uso: python3 scripts/08_ancestral_states.py  (después del paso 7)
"""
import csv
import subprocess
from pathlib import Path

from Bio import Phylo

ROOT = Path(__file__).resolve().parents[1]
ALN = ROOT / "results/mapt_2N4R_final.aln.faa"
TREE = ROOT / "results/phylo_protein/mapt_2N4R_constrained.treefile"
LINEAGES = ROOT / "data/raw/ncbi_lineages.tsv"
OUT = ROOT / "results/asr"
def best_model():
    """El modelo que eligió ModelFinder en el paso 5 (HIVb+F+R4 en la corrida publicada)."""
    for line in open(ROOT / "results/phylo_protein/mapt_2N4R.log"):
        if line.startswith("Best-fit model:"):
            return line.split()[2]
    raise RuntimeError("no se encontró el modelo de ModelFinder")
SITES = [50, 52, 69, 111, 113, 181, 184, 205]

# Clados definidos por un taxón del linaje NCBI.
CLADES = ["Theria", "Metatheria", "Eutheria", "Glires", "Rodentia", "Myomorpha",
          "Muroidea", "Cricetidae+Muridae", "Sciuromorpha", "Hystricomorpha", "Lagomorpha",
          "Primates", "Haplorrhini", "Strepsirrhini", "Chiroptera", "Eulipotyphla", "Afrotheria",
          "Laurasiatheria"]


def run_iqtree():
    OUT.mkdir(parents=True, exist_ok=True)
    if (OUT / "asr.state").exists():
        return
    subprocess.run(["iqtree2", "-s", str(ALN), "-te", str(TREE), "-m", best_model(),
                    "--ancestral", "-T", "4", "--prefix", str(OUT / "asr"), "-redo", "--quiet"],
                   check=True)


def human_columns():
    seq, name = "", None
    for line in open(ALN):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
        elif name == "Homo_sapiens":
            seq += line
    col, i = {}, 0
    for c, a in enumerate(seq, 1):
        if a != "-":
            i += 1
            col[i] = c
    return col


def read_states():
    """nodo -> columna -> (estado, {aa: probabilidad})"""
    states = {}
    with open(OUT / "asr.state") as fh:
        header = None
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if header is None:
                header = f
                aas = [h.replace("p_", "") for h in header[3:]]
                continue
            probs = dict(zip(aas, map(float, f[3:])))
            states.setdefault(f[0], {})[int(f[1])] = (f[2], probs)
    return states


def main():
    run_iqtree()
    col = human_columns()
    states = read_states()
    lin = {r["id"]: r["lineage"].split(";") for r in csv.DictReader(open(LINEAGES), delimiter="\t")}
    tree = Phylo.read(OUT / "asr.treefile", "newick")
    tips = [t.name for t in tree.get_terminals()]
    mono = [t for t in tips if "Monotremata" in lin.get(t, [])]
    tree.root_with_outgroup(*mono)

    rows = []
    for clade in CLADES:
        parts = clade.split("+")  # "Cricetidae+Muridae" es la unión de dos familias
        members = [t for t in tips if any(p in lin.get(t, []) for p in parts)]
        if len(members) < 2:
            continue
        node = tree.common_ancestor(members)
        n_in = len(node.get_terminals())
        row = {"clado": clade, "n_especies": len(members), "monofiletico": n_in == len(members),
               "nodo": node.confidence if isinstance(node.confidence, str) else node.name}
        for s in SITES:
            state, probs = states[row["nodo"]][col[s]]
            row[f"pos{s}_estado"] = state
            row[f"pos{s}_P(T)"] = round(probs.get("T", 0.0), 3)
            row[f"pos{s}_P(S)"] = round(probs.get("S", 0.0), 3)
        rows.append(row)

    with open(OUT / "ancestral_states.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print(r["clado"], r["n_especies"], r["monofiletico"], r["nodo"],
              *(f"{s}:{r[f'pos{s}_estado']}({r[f'pos{s}_P(T)']})" for s in SITES))


if __name__ == "__main__":
    main()
