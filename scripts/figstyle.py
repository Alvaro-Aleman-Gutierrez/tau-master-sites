"""Estilo común de las figuras, según las instrucciones de arte de Elsevier.

Elsevier pide archivos separados por figura, EPS o PDF vectorial con las fuentes incrustadas
(o TIFF a 300 dpi o más para medios tonos, 1000 dpi para línea), ancho de 90 mm (una columna),
140 mm (1,5) o 190 mm (doble), y color accesible.
"""
import glob
import os

import matplotlib as mpl
mpl.use("Agg")
from matplotlib import font_manager as fm
from PIL import Image

MM = 1 / 25.4
W1, W15, W2 = 90 * MM, 140 * MM, 190 * MM

# Paleta validada para daltonismo (separación CVD mínima ΔE 10,9).
FOCUS = "#B4700E"   # sitios maestros
REF = "#2A6F9E"     # otros sitios anotados
MARK = "#8C2D5A"
INK = "#1a1a1a"
MUTED = "#595959"
LIGHT = "#b8b8b8"   # sitios de fondo
GRID = "#dcdcdc"
SHADE = "#f3eadc"   # fondo del clado resaltado (FOCUS muy aclarado, sin alfa)
PAPER = "#ffffff"


def _arial():
    pats = ["/mnt/c/Windows/Fonts/arial*.ttf", "/mnt/c/Windows/Fonts/cour*.ttf",
            "/usr/share/fonts/**/LiberationSans*.ttf",
            "/usr/share/fonts/**/NimbusSans*.otf"]
    for f in (f for p in pats for f in glob.glob(p, recursive=True)):
        try:
            fm.fontManager.addfont(f)
        except Exception:
            pass
    names = {f.name for f in fm.fontManager.ttflist}
    return next((w for w in ("Arial", "Liberation Sans", "Nimbus Sans", "Helvetica") if w in names),
                "DejaVu Sans")


FONT = _arial()
mpl.rcParams.update({
    "font.family": FONT, "font.size": 7,
    "axes.labelsize": 7, "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "grid.linewidth": 0.5, "lines.linewidth": 0.8,
    "figure.facecolor": PAPER, "axes.facecolor": PAPER, "savefig.facecolor": PAPER,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
})


def panel(ax, letter, dx=-0.06, dy=1.02):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=8, fontweight="bold",
            va="bottom", ha="left", color=INK)


def despine(ax, keep=("left", "bottom")):
    for s, sp in ax.spines.items():
        sp.set_visible(s in keep)


def save(fig, name, outdir, dpi=1000):
    """EPS y PDF para el envío, TIFF de respaldo y PNG para revisar.

    Elsevier pide 1000 dpi para dibujos de línea en mapa de bits (7480 px a doble columna).
    """
    os.makedirs(outdir, exist_ok=True)
    base = os.path.join(outdir, name)
    # Sin fecha de creación, para que volver a generar una figura idéntica no cambie el archivo.
    fig.savefig(base + ".eps", format="eps", metadata={"CreationDate": None})
    fig.savefig(base + ".pdf", metadata={"CreationDate": None})
    fig.savefig(base + ".png", dpi=dpi)
    im = Image.open(base + ".png").convert("RGB")
    im.save(base + ".tif", format="TIFF", compression="tiff_lzw", dpi=(dpi, dpi))
    print(f"{name}  {im.size[0]}x{im.size[1]} px a {dpi} dpi")
