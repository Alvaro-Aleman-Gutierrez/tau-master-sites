"""Paso 6: FEL y FUBAR (HyPhy) sobre el alineamiento de codones de los 259 mamíferos.

FEL estima una tasa sinónima y otra no sinónima por codón y prueba si difieren. FUBAR coloca cada
codón en una rejilla de pares de tasas común a todos y devuelve la probabilidad posterior de que la
tasa no sinónima sea menor que la sinónima. Ambos usan la topología con restricción taxonómica del
paso 5 sin longitudes de rama, que HyPhy reestima.

El ejecutable se toma de la variable de entorno HYPHY (por defecto, "hyphy" en el PATH).
"""
import io
import os
import subprocess

from Bio import Phylo

HYPHY = os.environ.get("HYPHY", "hyphy")
CPU = os.environ.get("THREADS", "10")
OUT = "results/selection"
ALN = "results/mapt_2N4R_final.codon.fna"
os.makedirs(OUT, exist_ok=True)

tree = Phylo.read("results/phylo_protein/mapt_2N4R_constrained.treefile", "newick")
for c in tree.find_clades():
    c.branch_length = None
    if not c.is_terminal():
        c.name = None
        c.confidence = None
s = io.StringIO()
Phylo.write(tree, s, "newick", plain=True)
TREE = f"{OUT}/all.nwk"
open(TREE, "w").write(s.getvalue().strip() + "\n")

for method in ("FEL", "FUBAR"):
    js = f"{OUT}/all.{method}.json"
    if os.path.exists(js):
        print(f"{js} ya existe, no se vuelve a correr")
        continue
    cmd = [HYPHY, f"CPU={CPU}", method.lower(), "--alignment", ALN, "--tree", TREE, "--output", js]
    print(" ".join(cmd), flush=True)
    with open(f"{OUT}/all.{method}.log", "w") as log:
        subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)
    print(f"{method} terminado", flush=True)
