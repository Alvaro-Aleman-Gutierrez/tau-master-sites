"""Paso 9: figuras 1 y 2 del manuscrito.

Figura 1. Conservación de los 85 S/T/Y de la tau 2N4R humana en 259 mamíferos (a) y en
tetrápodos no mamíferos (b), sobre el mapa de exones y regiones.
Figura 2. Estado de la posición 50 en las puntas del árbol de 259 mamíferos, con los nodos
reconstruidos (a), y el alineamiento local en especies representativas (b).

Uso: python3 scripts/09_figures.py  (después de los pasos 7 y 8)
"""
import csv
import importlib.util
import sys
from pathlib import Path

from Bio import Phylo
from matplotlib.collections import LineCollection
from matplotlib.patches import Patch, Rectangle
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
from figstyle import (FOCUS, REF, MARK, INK, MUTED, LIGHT, GRID, SHADE, W2, MM,  # noqa: E402
                      panel, despine, save)

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures"
spec = importlib.util.spec_from_file_location("tm", ROOT / "scripts/07_site_table.py")
tm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tm)

MASTER = {50, 69, 181}
# Rótulos inclinados a 45° que salen del punto con una línea guía, hacia la derecha ("r")
# o la izquierda ("l"), con fondo blanco para no quedar sobre los tallos vecinos.
LABEL_DIR = {
    ("a", 50): "l", ("a", 52): "r", ("a", 69): "r", ("a", 113): "r", ("a", 181): "l",
    ("a", 184): "r", ("a", 217): "r",
    ("b", 50): "l", ("b", 52): "r", ("b", 69): "r", ("b", 113): "r", ("b", 181): "l",
    ("b", 184): "r", ("b", 217): "r",
}
LEADER = 5  # largo de la línea guía, en puntos
LABELLED = {111, 205, 199, 212, 217, 231, 235, 396, 404, 422, 113, 184, 52}


def read_table():
    return list(csv.DictReader(open(ROOT / "results/site_table.tsv"), delimiter="\t"))


def site_color(pos):
    if pos in MASTER:
        return FOCUS
    if pos in LABELLED:
        return REF
    return LIGHT


def exon_blocks():
    rows = {r["exon"]: len(r["peptide"]) for r in csv.DictReader(open(tm.EXONS), delimiter="\t")}
    start, out = 1, []
    for ex in tm.EXONS_2N4R:
        out.append((ex, start, start + rows[ex] - 1))
        start += rows[ex]
    return out


def figure1(rows):
    fig = plt.figure(figsize=(W2, 106 * MM))
    # Los paneles a y b se separan más que b y la franja de exones, que describe a ambos.
    outer = fig.add_gridspec(2, 1, height_ratios=[2.5, 0.42], hspace=0.12,
                             left=0.07, right=0.99, top=0.93, bottom=0.16)
    inner = outer[0].subgridspec(2, 1, hspace=0.5)
    axes = [fig.add_subplot(inner[0]), fig.add_subplot(inner[1])]
    track = fig.add_subplot(outer[1], sharex=axes[0])

    series = [
        ("Mammals (n = 259)", lambda r: int(r["mam_identico"]) / int(r["mam_n"])),
        ("Non-mammalian tetrapods", lambda r: (int(r["tetra_identico"]) / int(r["tetra_n_con_posicion"])
                                               if int(r["tetra_n_con_posicion"]) else None)),
    ]
    for ax, (label, fn), letter in zip(axes, series, "ab"):
        for order in (LIGHT, REF, FOCUS):  # el fondo primero, los sitios anotados encima
            for r in rows:
                pos, y = int(r["pos"]), fn(r)
                if y is None or site_color(pos) != order:
                    continue
                ax.vlines(pos, 0, y, color=order, lw=0.8 if order != LIGHT else 0.6)
                ax.plot(pos, y, "o", ms=3.2 if order != LIGHT else 2.4, color=order,
                        mec="white", mew=0.5, zorder=3)
        for r in rows:
            pos, y = int(r["pos"]), fn(r)
            if y is None or (pos not in MASTER and pos not in {52, 113, 184, 217}):
                continue
            right = LABEL_DIR.get((letter, pos), "r") == "r"
            ax.annotate(r["sitio"], (pos, y), xytext=(LEADER if right else -LEADER, LEADER),
                        textcoords="offset points", rotation=45 if right else -45,
                        rotation_mode="anchor", ha="left" if right else "right", va="bottom",
                        fontsize=7, color=INK, fontweight="bold" if pos in MASTER else "normal",
                        bbox=dict(boxstyle="square,pad=0.08", fc="white", ec="none"),
                        arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.5, shrinkA=0, shrinkB=0),
                        annotation_clip=False, zorder=5)
        ax.set_xlim(0, 442)
        ax.set_ylim(0, 1.08)
        ax.set_yticks([0, 0.5, 1])
        ax.set_ylabel("Fraction identical\nto human")
        ax.grid(axis="y", color=GRID, lw=0.5)
        ax.set_axisbelow(True)
        ax.set_title(label, loc="right", fontsize=7, color=MUTED, pad=2)
        despine(ax)
        # el panel a lleva su propia escala; la de b es la del eje inferior, bajo la franja de exones
        ax.tick_params(labelbottom=letter == "a")
        panel(ax, letter, dx=-0.055)

    for i, (ex, a, b) in enumerate(exon_blocks()):
        track.add_patch(Rectangle((a - 0.5, 0.58), b - a + 1, 0.4,
                                  fc="#e9e9e9" if i % 2 else "#d4d4d4", ec="white", lw=0.6))
        track.text((a + b) / 2, 0.78, ex.replace("E", ""), ha="center", va="center", fontsize=7)
    for a, b, name in tm.REGIONS:
        track.plot([a + 1, b - 1], [0.46, 0.46], color=MUTED, lw=0.8, solid_capstyle="butt")
        track.text((a + b) / 2, 0.40, name, ha="center", va="top", fontsize=7, color=MUTED)
    track.set_ylim(0, 1)
    track.set_yticks([])
    track.set_xlim(0, 442)
    track.set_xlabel("Residue (human 2N4R tau)")
    despine(track, keep=("bottom",))
    track.text(-0.01, 0.78, "Exon", transform=track.get_yaxis_transform(), ha="right",
               va="center", fontsize=7, color=MUTED)

    handles = [Patch(color=FOCUS, label="Master sites"),
               Patch(color=REF, label="Other annotated sites"),
               Patch(color=LIGHT, label="Other S/T/Y")]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
               bbox_to_anchor=(0.53, 0.0), handlelength=1, handleheight=0.8)
    save(fig, "Figure_1", FIG)


STATE_COLORS = {"T": REF, "S": MARK, "P": FOCUS}


def tip_state(res):
    if res == "-":
        return None
    return res if res in STATE_COLORS else "other"


def figure2():
    aln, order, col, human = tm.mammal_data()
    lin = {r["id"]: r["lineage"].split(";")
           for r in csv.DictReader(open(tm.LINEAGES), delimiter="\t")}
    tree = Phylo.read(ROOT / "results/asr/asr.treefile", "newick")
    tips = [t.name for t in tree.get_terminals()]
    tree.root_with_outgroup(*[t for t in tips if "Monotremata" in lin[t]])
    tree.ladderize(reverse=True)
    tips = [t.name for t in tree.get_terminals()]
    ypos = {name: i for i, name in enumerate(tips)}
    depths = tree.depths()
    if max(depths.values()) == 0:
        depths = tree.depths(unit_branch_lengths=True)
    xmax = max(depths.values())

    fig = plt.figure(figsize=(W2, 150 * MM))
    ax = fig.add_axes([0.01, 0.03, 0.56, 0.93])
    segs = []

    def y_of(cl):
        if cl.is_terminal():
            return ypos[cl.name]
        ys = [y_of(c) for c in cl.clades]
        cl._y = sum((min(ys), max(ys))) / 2
        for c in cl.clades:
            cy = ypos[c.name] if c.is_terminal() else c._y
            segs.append([(depths[cl], cy), (depths[c], cy)])
        segs.append([(depths[cl], min(ys)), (depths[cl], max(ys))])
        return cl._y

    y_of(tree.root)

    # clado Cricetidae + Muridae resaltado
    cm = [t for t in tips if "Cricetidae" in lin[t] or "Muridae" in lin[t]]
    cm_node = tree.common_ancestor(cm)
    y0, y1 = min(ypos[t] for t in cm) - 0.5, max(ypos[t] for t in cm) + 0.5
    ax.add_patch(Rectangle((depths[cm_node], y0), xmax * 1.16 - depths[cm_node], y1 - y0,
                           fc=SHADE, ec="none", zorder=0))
    ax.add_collection(LineCollection(segs, colors=MUTED, linewidths=0.35, zorder=1))

    strip = xmax * 1.03
    for t in tips:
        st = tip_state(aln[t][col[50]])
        if st is None:
            ax.add_patch(Rectangle((strip, ypos[t] - 0.45), xmax * 0.05, 0.9,
                                   fc=INK, ec=INK, lw=0.3, zorder=2))
        else:
            ax.add_patch(Rectangle((strip, ypos[t] - 0.5), xmax * 0.05, 1.0,
                                   fc=STATE_COLORS.get(st, LIGHT), ec="none", zorder=2))

    # órdenes con 4 o más especies, rotulados a la derecha
    lab_x = xmax * 1.11
    def group(t):
        if "Metatheria" in lin[t]:
            return "Marsupials"
        if "Monotremata" in lin[t]:
            return "Monotremes"
        return order.get(t, "?")

    runs, cur = [], None
    for t in tips:
        o = group(t)
        if cur and cur[0] == o:
            cur[2] = ypos[t]
        else:
            cur = [o, ypos[t], ypos[t]]
            runs.append(cur)
    for o, a, b in runs:
        if b - a + 1 < 4 and o != "Monotremes":
            continue
        ax.plot([lab_x - xmax * 0.01, lab_x, lab_x, lab_x - xmax * 0.01],
                [a + 0.2, a + 0.2, b - 0.2, b - 0.2], color=INK, lw=0.6)
        ax.text(lab_x + xmax * 0.015, (a + b) / 2, f"{o} ({b - a + 1})", va="center",
                fontsize=7)

    # nodos reconstruidos
    states = {}
    asr = list(csv.DictReader(open(ROOT / "results/asr/ancestral_states.tsv"), delimiter="\t"))
    for r in asr:
        states[r["clado"]] = r
    marks = [("Theria", f"T {float(states['Theria']['pos50_P(T)']):.2f}"),
             ("Rodentia", f"T {float(states['Rodentia']['pos50_P(T)']):.2f}"),
             ("Myomorpha", f"T {float(states['Myomorpha']['pos50_P(T)']):.2f}")]
    # Los rótulos de nodos dentro de Eutheria van a la izquierda de su eje principal, en espacio
    # libre, con una línea guía hasta el nodo, para no quedar sobre las ramas.
    eu_x = depths[tree.common_ancestor([t for t in tips if "Eutheria" in lin[t]])]
    for clade, text in marks:
        members = [t for t in tips if clade in lin[t]]
        nd = tree.common_ancestor(members)
        ax.plot(depths[nd], nd._y, "o", ms=3.5, color=REF, mec="white", mew=0.5, zorder=4)
        if clade == "Theria":
            ax.text(depths[nd] - xmax * 0.008, nd._y + 0.8, f"{clade}\n{text}", ha="right",
                    va="top", fontsize=7, color=INK)
            continue
        lx = eu_x - xmax * 0.03
        ax.plot([lx + xmax * 0.004, depths[nd]], [nd._y, nd._y], color=MUTED, lw=0.4, ls=":", zorder=3)
        ax.text(lx, nd._y, f"{clade} {text}", ha="right", va="center", fontsize=7, color=INK)
    ax.plot(depths[cm_node], cm_node._y, "o", ms=3.5, color=FOCUS, mec="white", mew=0.5, zorder=4)
    tip_x = max(depths[tree.find_any(name=t)] for t in cm)
    ax.text(tip_x + xmax * 0.015, y0 + 0.8, "Cricetidae\n+ Muridae\nT\u2192P\n(P 1.00)",
            ha="left", va="top", fontsize=7, color=INK)

    ax.set_xlim(0, xmax * 1.42)
    ax.set_ylim(len(tips) - 0.5, -0.5)
    ax.axis("off")
    ax.plot([0, 0.1], [len(tips) + 1.5] * 2, color=INK, lw=0.6, clip_on=False)
    ax.text(0, len(tips) + 2.5, "0.1 substitutions per site", ha="left", va="top", fontsize=7)
    panel(ax, "a", dx=0.0, dy=1.0)

    handles = [Patch(color=REF, label="T"), Patch(color=MARK, label="S"),
               Patch(color=FOCUS, label="P"), Patch(color=LIGHT, label="Other residue"),
               Patch(color=INK, label="Gap")]
    fig.legend(handles=handles, title="Residue at position 50", loc="upper left",
               bbox_to_anchor=(0.62, 0.32), frameon=False, handlelength=1, handleheight=0.8,
               title_fontsize=7, fontsize=7)

    figure2_alignment(fig, aln, col, human)
    save(fig, "Figure_2", FIG)


ALIGN_SPECIES = [
    ("Homo sapiens", "Homo_sapiens", "Primates"),
    ("Sciurus carolinensis", None, "Rodentia, Sciuridae"),
    ("Cavia porcellus", "Cavia_porcellus", "Rodentia, Caviidae"),
    ("Jaculus jaculus", "Jaculus_jaculus", "Rodentia, Dipodidae"),
    ("Nannospalax galili", "Nannospalax_galili", "Rodentia, Spalacidae"),
    ("Mus musculus", "Mus_musculus", "Rodentia, Muridae"),
    ("Rattus norvegicus", "Rattus_norvegicus", "Rodentia, Muridae"),
    ("Mesocricetus auratus", "Mesocricetus_auratus", "Rodentia, Cricetidae"),
    ("Monodelphis domestica", "Monodelphis_domestica", "Didelphimorphia"),
    ("Sarcophilus harrisii", "Sarcophilus_harrisii", "Dasyuromorphia"),
    ("Ornithorhynchus anatinus", "Ornithorhynchus_anatinus", "Monotremata"),
    ("Gallus gallus", "tetra", "Aves"),
    ("Anolis carolinensis", "tetra", "Lepidosauria"),
    ("Chrysemys picta bellii", "tetra", "Testudines"),
]
SPAN = (43, 57)
ROW = 1.4  # separación vertical entre especies en el panel b


def figure2_alignment(fig, aln, col, human):
    tetra, _ = tm.tetrapod_residues(human)
    ax = fig.add_axes([0.62, 0.36, 0.37, 0.60])
    ax.axis("off")
    panel(ax, "b", dx=0.0, dy=1.0)
    positions = list(range(SPAN[0], SPAN[1] + 1))
    n = len(ALIGN_SPECIES)
    x0, dx = 0.52, 0.032
    for j, p in enumerate(positions):
        if p in (50, 52):
            ax.add_patch(Rectangle((x0 + j * dx - dx / 2, -0.2), dx, n * ROW + 0.2,
                                   transform=ax.transData, fc=SHADE, ec="none"))
        if p % 5 == 0 or p in (50, 52):
            ax.text(x0 + j * dx, -0.8, str(p), ha="center", va="bottom", fontsize=7,
                    color=INK if p in (50, 52) else MUTED, fontweight="bold" if p == 50 else "normal")
    shown = []
    for name, key, group in ALIGN_SPECIES:
        if key is None:
            key = next((k for k in aln if k.startswith(name.split()[0])), None)
            if key is None:
                continue
            name = key.replace("_", " ")
        if key == "tetra":
            res = tetra.get(name, {})
            seq = "".join(next(iter(res[p])) if p in res else "-" for p in positions)
        else:
            seq = "".join(aln[key][col[p]] for p in positions)
        shown.append((name, group, seq))
    for i, (name, group, seq) in enumerate(shown):
        ax.text(0.0, i * ROW + 0.4, name, fontsize=7, style="italic", va="center")
        ax.text(0.0, i * ROW + 1.0, group, fontsize=7, color=MUTED, va="center")
        for j, aa in enumerate(seq):
            same = aa == human[positions[j] - 1]
            color = INK
            if positions[j] == 50:
                color = STATE_COLORS.get(aa, MUTED) if aa != "-" else MUTED
            ax.text(x0 + j * dx, i * ROW + 0.4, aa if (not same or i == 0 or positions[j] == 50) else ".",
                    ha="center", va="center", fontsize=7, family="Courier New", color=color,
                    fontweight="bold" if positions[j] == 50 else "normal")
    ax.set_xlim(0, x0 + len(positions) * dx)
    ax.set_ylim(len(shown) * ROW + 0.2, -1.0)


def main():
    rows = read_table()
    figure1(rows)
    figure2()


if __name__ == "__main__":
    main()
