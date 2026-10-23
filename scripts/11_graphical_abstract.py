"""Graphical abstract: residuo de los tres sitios maestros en seis linajes de tetrápodos.

Elsevier pide al menos 1328 px de ancho, el lado largo horizontal y que se lea a 13 cm de ancho.
Usa el mismo código de color que la Figura 2 (T, S, P, otro residuo).

Uso: python3 scripts/11_graphical_abstract.py  (después del paso 7)
"""
import csv
import importlib.util
import sys
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

sys.path.insert(0, str(Path(__file__).parent))
from figstyle import REF, MARK, FOCUS, LIGHT, INK, MUTED, MM, save  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("tm", ROOT / "scripts/07_site_table.py")
tm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tm)

SITES = [50, 69, 181]
LINEAGES = ["Amphibians", "Sauropsids", "Monotremes", "Marsupials", "Other placentals",
            "Mice, rats and hamsters"]
COLORS = {"T": REF, "S": MARK, "P": FOCUS}
AMPHIBIA = {"Anura", "Caudata", "Gymnophiona"}


def counts():
    aln, order, col, human = tm.mammal_data()
    tetra, groups = tm.tetrapod_residues(human)
    lin = {r["id"]: r["lineage"] for r in csv.DictReader(open(tm.LINEAGES), delimiter="\t")}

    def mgroup(sp):
        L = lin[sp]
        if "Monotremata" in L:
            return "Monotremes"
        if "Metatheria" in L:
            return "Marsupials"
        if "Cricetidae" in L or "Muridae" in L:
            return "Mice, rats and hamsters"
        return "Other placentals"

    c = defaultdict(Counter)
    for site in SITES:
        for sp, s in aln.items():
            if s[col[site]] != "-":
                c[(mgroup(sp), site)][s[col[site]]] += 1
        for sp, d in tetra.items():
            if site in d:
                g = "Amphibians" if groups[sp] in AMPHIBIA else "Sauropsids"
                c[(g, site)][next(iter(d[site]))] += 1
    return c


def main():
    c = counts()
    fig = plt.figure(figsize=(130 * MM, 62 * MM))
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    ax.set_xlim(0, 130)
    ax.set_ylim(62, 0)
    ax.axis("off")

    ax.text(4, 5, "Tau master phosphorylation sites across tetrapods", fontsize=9, fontweight="bold",
            va="center")

    # cladograma simplificado: anfibios, saurópsidos, monotremas, marsupiales, placentarios
    y = {g: 16 + i * 7.3 for i, g in enumerate(LINEAGES)}
    x_tip, lw = 34, 0.9
    yo, yr = y["Other placentals"], y["Mice, rats and hamsters"]
    nodes = {
        "plac": (30, (yo + yr) / 2),
    }
    nodes["ther"] = (25, (y["Marsupials"] + nodes["plac"][1]) / 2)
    nodes["mam"] = (20, (y["Monotremes"] + nodes["ther"][1]) / 2)
    nodes["amn"] = (14, (y["Sauropsids"] + nodes["mam"][1]) / 2)
    nodes["tet"] = (8, (y["Amphibians"] + nodes["amn"][1]) / 2)

    def join(parent, children):
        px, py = nodes[parent]
        ys = [nodes[c][1] if c in nodes else y[c] for c in children]
        ax.plot([px, px], [min(ys), max(ys)], color=MUTED, lw=lw)
        for ch, cy in zip(children, ys):
            cx = nodes[ch][0] if ch in nodes else x_tip
            ax.plot([px, cx], [cy, cy], color=MUTED, lw=lw)

    join("plac", ["Other placentals", "Mice, rats and hamsters"])
    join("ther", ["Marsupials", "plac"])
    join("mam", ["Monotremes", "ther"])
    join("amn", ["Sauropsids", "mam"])
    join("tet", ["Amphibians", "amn"])
    ax.plot([5, 8], [nodes["tet"][1]] * 2, color=MUTED, lw=lw)

    for g in LINEAGES:
        ax.text(x_tip + 1.5, y[g], g, va="center", fontsize=7.5,
                fontweight="bold" if g.startswith("Mice") else "normal")

    # matriz de residuos
    x0, w, h = 78, 13, 6.2
    for j, site in enumerate(SITES):
        ax.text(x0 + j * (w + 1.5) + w / 2, 10.5, f"T{site}", ha="center", va="center", fontsize=8,
                fontweight="bold")
        for g in LINEAGES:
            cnt = c[(g, site)]
            n = sum(cnt.values())
            cx, cy = x0 + j * (w + 1.5), y[g] - h / 2
            aa, k = cnt.most_common(1)[0]
            # con pocas especies se muestran conteos, no porcentajes
            label = f"{aa}  {k}/{n}" if n < 5 else f"{aa}  {round(100 * k / n)}%"
            ax.add_patch(Rectangle((cx, cy), w, h, fc=COLORS.get(aa, LIGHT), ec="white", lw=1.0))
            ax.text(cx + w / 2, y[g], label, ha="center", va="center",
                    fontsize=7, color="white" if aa in COLORS else INK, fontweight="bold")

    # alineada al borde derecho de la matriz para que no se salga de la figura
    ax.text(x0 + len(SITES) * w + (len(SITES) - 1) * 1.5, 58.5,
            "Most common residue (share or count of species)", fontsize=7, color=MUTED,
            ha="right", va="center")
    save(fig, "Graphical_abstract", ROOT / "figures", dpi=300)


if __name__ == "__main__":
    main()
