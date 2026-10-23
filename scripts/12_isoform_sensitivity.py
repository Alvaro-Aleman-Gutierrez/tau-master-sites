"""Paso 12: sensibilidad de la conservación en tetrápodos a la regla de selección de isoformas.

Regla A (la del manuscrito): en cada posición, la isoforma con mayor identidad local (±5) a la
humana, si llega a 0,5. Puede combinar isoformas distintas en posiciones distintas.
Regla B: una sola isoforma por especie, la forma 4R, elegida como la isoforma con mayor cobertura
del exón 10 humano (al menos 0,9) y, entre ellas, la más larga. No mira los sitios estudiados. El
desempate por largo evita elegir isoformas con un primer exón alternativo, que tienen algo más de
identidad global con el humano pero alinean secuencia no homóloga sobre los exones 2 y 3.
Regla C: una sola isoforma por especie, la más larga.

En B y C se informa el residuo alineado y si el entorno del sitio es homólogo (identidad local
≥ 0,5); sin ese filtro, una isoforma sin exón 2 aportaría un residuo no homólogo.

Uso: python3 scripts/12_isoform_sensitivity.py  → results/isoform_sensitivity.tsv
"""
import csv
import importlib.util
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("tm", ROOT / "scripts/07_site_table.py")
tm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tm)

from Bio import Align  # noqa: E402
from Bio.Align import substitution_matrices  # noqa: E402

SITES = [50, 69, 181, 217, 403]
AMPHIBIA = {"Anura", "Caudata", "Gymnophiona"}


def main():
    _, _, _, human = tm.mammal_data()
    groups = {r["species"]: r["group"] for r in csv.DictReader(open(tm.TETRA_SP), delimiter="\t")}
    seqs = defaultdict(dict)
    for header, seq in tm.read_fasta(tm.TETRA_PROT).items():
        sp = header.split("[organism=")[1].split("]")[0]
        if sp in groups:
            seqs[sp][header.split()[0]] = seq.rstrip("*")
    aligner = Align.PairwiseAligner(mode="global", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                                    open_gap_score=-10, extend_gap_score=-0.5, end_gap_score=0)

    # residuos del exón 10 en la numeración de 2N4R
    e10 = [p for p, ex in tm.exon_map(human).items() if ex == "E10"]

    def choose_4r(isoforms):
        best = None
        for acc, seq in isoforms.items():
            if len(seq) < 200:
                continue
            a = aligner.align(human, seq)[0]
            res = {}
            for (h0, h1), (s0, _) in zip(*a.aligned):
                for k in range(h1 - h0):
                    res[h0 + k + 1] = seq[s0 + k]
            cov = sum(p in res for p in e10) / len(e10)
            if best is None or (cov, len(seq)) > best[0]:
                best = ((cov, len(seq)), acc)
        return best[1] if best and best[0][0] >= 0.9 else None

    def read(seq):
        a = aligner.align(human, seq)[0]
        res = {}
        for (h0, h1), (s0, _) in zip(*a.aligned):
            for k in range(h1 - h0):
                res[h0 + k + 1] = seq[s0 + k]
        out = {}
        for pos in SITES:
            win = [p for p in range(pos - tm.WINDOW, pos + tm.WINDOW + 1) if p != pos]
            wid = sum(res.get(p) == human[p - 1] for p in win) / len(win)
            out[pos] = res[pos] if pos in res and wid >= tm.TETRA_MIN_WINDOW_ID else None
        return out

    rule_a, _ = tm.tetrapod_residues(human)
    rows = []
    for sp in sorted(seqs):
        g = "Amphibia" if groups[sp] in AMPHIBIA else "Sauropsida"
        acc4r = choose_4r(seqs[sp])
        b = read(seqs[sp][acc4r]) if acc4r else {p: None for p in SITES}
        longest = max(seqs[sp].values(), key=len)
        c = read(longest)
        for pos in SITES:
            a = next(iter(rule_a.get(sp, {}).get(pos, {None})))
            rows.append({"species": sp, "group": g, "site": pos, "rule_A": a or "", "rule_B": b[pos] or "",
                         "rule_C": c[pos] or ""})

    (ROOT / "tables").mkdir(exist_ok=True)
    for path in (ROOT / "results/isoform_sensitivity.tsv", ROOT / "tables/Table_S2.tsv"):
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
            w.writeheader()
            w.writerows(rows)

    for pos in SITES:
        print(f"\n== posición {pos} (humano {human[pos - 1]})")
        for g in ("Sauropsida", "Amphibia"):
            sub = [r for r in rows if r["site"] == pos and r["group"] == g]
            line = []
            for rule in ("rule_A", "rule_B", "rule_C"):
                cnt = Counter(r[rule] for r in sub if r[rule])
                line.append(f"{rule[-1]}: n={sum(cnt.values())} {dict(cnt.most_common(4))}")
            both = [r for r in sub if r["rule_A"] and r["rule_B"]]
            agree = sum(r["rule_A"] == r["rule_B"] for r in both)
            print(f"  {g:10} " + " | ".join(line) + f" | A=B en {agree}/{len(both)}")


if __name__ == "__main__":
    main()
