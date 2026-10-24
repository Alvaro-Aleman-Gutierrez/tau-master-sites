"""Paso 14: Figura 3, esquema de la historia evolutiva de los sitios estudiados.

Cladograma de nueve linajes de tetrápodos con el residuo más común en las posiciones 50, 69, 181,
113 y 184, y las transiciones inferidas dibujadas sobre las ramas. Dentro de mamíferos las
transiciones vienen de la reconstrucción ancestral (paso 8); la de los saurópsidos, de la
parsimonia sobre la topología de tetrápodos (paso 13). Las relaciones dentro de Placentalia se
dejan como politomía, porque la figura no depende de ellas.

Uso: python3 scripts/14_figure3_schematic.py  (después de los pasos 7, 8 y 13)
"""
import csv
import importlib.util
import sys
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

sys.path.insert(0, str(Path(__file__).parent))
from figstyle import REF, MARK, FOCUS, LIGHT, INK, MUTED, MM, W15, save  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("tm", ROOT / "scripts/07_site_table.py")
tm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tm)

SITES = [50, 69, 181, 113, 184]
LINEAGES = ["Amphibians", "Sauropsids", "Monotremes", "Marsupials", "Other placentals",
            "Strepsirrhine primates", "Haplorrhine primates", "Other rodents", "Cricetidae + Muridae"]
COLORS = {"T": REF, "S": MARK, "P": FOCUS}
AMPHIBIA = {"Anura", "Caudata", "Gymnophiona"}


def lineage_of(L):
    if "Monotremata" in L:
        return "Monotremes"
    if "Metatheria" in L:
        return "Marsupials"
    if "Cricetidae" in L or "Muridae" in L:
        return "Cricetidae + Muridae"
    if "Rodentia" in L:
        return "Other rodents"
    if "Strepsirrhini" in L:
        return "Strepsirrhine primates"
    if "Haplorrhini" in L:
        return "Haplorrhine primates"
    return "Other placentals"


def counts():
    aln, _, col, human = tm.mammal_data()
    tetra, groups = tm.tetrapod_residues(human)
    lin = {r["id"]: r["lineage"] for r in csv.DictReader(open(tm.LINEAGES), delimiter="\t")}
    c = defaultdict(Counter)
    for site in SITES:
        for sp, s in aln.items():
            if s[col[site]] != "-":
                c[(lineage_of(lin[sp]), site)][s[col[site]]] += 1
        for sp, d in tetra.items():
            if site in d:
                g = "Amphibians" if groups[sp] in AMPHIBIA else "Sauropsids"
                c[(g, site)][next(iter(d[site]))] += 1
    return c


def main():
    c = counts()
    fig = plt.figure(figsize=(W15, 92 * MM))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 140)
    ax.set_ylim(92, 0)
    ax.axis("off")

    y = {g: 12 + i * 8.6 for i, g in enumerate(LINEAGES)}
    tip = 54
    lw = 0.9
    # nodos internos (x, y)
    node = {}
    node["rod"] = (44, (y["Other rodents"] + y["Cricetidae + Muridae"]) / 2)
    node["pri"] = (44, (y["Strepsirrhine primates"] + y["Haplorrhine primates"]) / 2)
    node["pla"] = (35, (y["Other placentals"] + node["rod"][1]) / 2)
    node["the"] = (27, (y["Marsupials"] + node["pla"][1]) / 2)
    node["mam"] = (21, (y["Monotremes"] + node["the"][1]) / 2)
    node["amn"] = (14, (y["Sauropsids"] + node["mam"][1]) / 2)
    node["tet"] = (7, (y["Amphibians"] + node["amn"][1]) / 2)

    def yy(k):
        return node[k][1] if k in node else y[k]

    def xx(k):
        return node[k][0] if k in node else tip

    def join(parent, kids):
        px, py = node[parent]
        ys = [yy(k) for k in kids]
        ax.plot([px, px], [min(ys), max(ys)], color=MUTED, lw=lw)
        for k in kids:
            ax.plot([px, xx(k)], [yy(k), yy(k)], color=MUTED, lw=lw)

    join("rod", ["Other rodents", "Cricetidae + Muridae"])
    join("pri", ["Strepsirrhine primates", "Haplorrhine primates"])
    join("pla", ["Other placentals", "pri", "rod"])
    join("the", ["Marsupials", "pla"])
    join("mam", ["Monotremes", "the"])
    join("amn", ["Sauropsids", "mam"])
    join("tet", ["Amphibians", "amn"])
    ax.plot([3, 7], [node["tet"][1]] * 2, color=MUTED, lw=lw)

    for g in LINEAGES:
        ax.text(tip + 1.5, y[g], g, va="center", fontsize=7,
                fontweight="bold" if g == "Cricetidae + Muridae" else "normal")

    # Transiciones sobre las ramas: (x de la marca, y de la rama, rótulo). Dentro de mamíferos vienen
    # de la reconstrucción ancestral; la de los saurópsidos, de la parsimonia del paso 13.
    # El cuarto valor indica si el rótulo va arriba (True) o abajo (False) de la rama.
    events = [((node["amn"][0] + tip) / 2 - 8, y["Sauropsids"], "T50→I", True),
              ((node["the"][0] + node["pla"][0]) / 2, node["pla"][1], "S184→G", False),
              ((node["pri"][0] + tip) / 2, y["Haplorrhine primates"], "N113→S\nG184→S", True),
              ((node["rod"][0] + tip) / 2, y["Cricetidae + Muridae"], "T50→P", True)]
    for xm, ym, label, above in events:
        ax.plot([xm, xm], [ym - 1.6, ym + 1.6], color=FOCUS, lw=2.2, solid_capstyle="butt", zorder=3)
        ax.text(xm, ym - 2.0 if above else ym + 2.0, label, ha="center", va="bottom" if above else "top",
                fontsize=7, color=INK, zorder=4, linespacing=1.1,
                bbox=dict(boxstyle="square,pad=0.1", fc="white", ec="none"))

    # matriz de residuos
    x0, w, h = 86.5, 9.6, 6.6
    for j, site in enumerate(SITES):
        aa_h = {50: "T", 69: "T", 181: "T", 113: "S", 184: "S"}[site]
        ax.text(x0 + j * (w + 1.0) + w / 2, 5.5, f"{aa_h}{site}", ha="center", va="center", fontsize=7,
                fontweight="bold")
        for g in LINEAGES:
            cnt = c[(g, site)]
            n = sum(cnt.values())
            cx, cy = x0 + j * (w + 1.0), y[g] - h / 2
            if not n:
                ax.add_patch(Rectangle((cx, cy), w, h, fc="white", ec=LIGHT, lw=0.6))
                ax.text(cx + w / 2, y[g], "–", ha="center", va="center", fontsize=7, color=MUTED)
                continue
            aa, k = cnt.most_common(1)[0]
            label = f"{aa} {k}/{n}" if n < 5 else f"{aa} {round(100 * k / n)}%"
            ax.add_patch(Rectangle((cx, cy), w, h, fc=COLORS.get(aa, LIGHT), ec="white", lw=1.0))
            ax.text(cx + w / 2, y[g], label, ha="center", va="center", fontsize=7,
                    color="white" if aa in COLORS else INK)
    ax.text(x0 + len(SITES) * (w + 1.0) - 1.0, 88.5, "Most common residue (share or count of species)",
            ha="right", va="center", fontsize=7, color=MUTED)
    save(fig, "Figure_3", ROOT / "figures")


if __name__ == "__main__":
    main()
