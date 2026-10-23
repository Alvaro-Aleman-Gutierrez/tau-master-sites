"""Paso 4: árbol de restricción a partir de los linajes de NCBI Taxonomy (todos los rangos con nombre).

Se usa con IQ-TREE (-g) para obtener un árbol ML consistente con los clados taxonómicos, que es la
topología de FEL, FUBAR y la reconstrucción ancestral. Los linajes se leen de la copia de NCBI
Taxonomy incluida en data/raw/ncbi_lineages.tsv (descargada el 21 de octubre de 2026), para que el
resultado no cambie si NCBI actualiza la taxonomía.
"""
import collections
import os

import pandas as pd

os.makedirs("results/phylo_protein", exist_ok=True)
sel = pd.read_csv("results/selected_2N4R.tsv", sep="\t")
keep = set(sel[sel.included].species.str.replace(" ", "_"))
lin_tab = pd.read_csv("data/raw/ncbi_lineages.tsv", sep="\t")
lineage = {r.id: r.lineage.split(";") for r in lin_tab.itertuples() if r.id in keep}
missing = keep - set(lineage)
assert not missing, f"especies sin linaje en la copia de NCBI Taxonomy: {sorted(missing)}"

# Árbol anidado: cada taxón con más de un descendiente se vuelve un nodo (politomías permitidas).
# Los nodos se identifican por la ruta completa porque hay nombres repetidos entre rangos
# (por ejemplo, género y subgénero Mustela).
children = collections.defaultdict(set)
for sp, lin in lineage.items():
    path = [tuple(["Mammalia"] + lin[:i]) for i in range(len(lin) + 1)]
    for a, b in zip(path, path[1:]):
        children[a].add(b)
    children[path[-1]].add(sp)


def newick(node):
    kids = [newick(k) if isinstance(k, tuple) else k for k in sorted(children[node], key=str)]
    return kids[0] if len(kids) == 1 else "(" + ",".join(kids) + ")"


open("results/phylo_protein/taxonomy_constraint.nwk", "w").write(newick(("Mammalia",)) + ";\n")
print(len(lineage), "especies en la restricción")
