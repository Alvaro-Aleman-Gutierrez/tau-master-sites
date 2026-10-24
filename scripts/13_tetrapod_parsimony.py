"""Paso 13: parsimonia de la posición 50 sobre la topología establecida de los tetrápodos.

El árbol de los pasos 5 y 8 solo tiene mamíferos, así que no ubica cambios fuera de ellos. Aquí se
usa la topología aceptada de los grandes grupos, (Amphibia, (Mammalia, Sauropsida)), con
Sauropsida = (Lepidosauria, (Testudines, (Crocodylia, Aves))), y el estado de cada grupo:

- anfibios y saurópsidos: los residuos observados en la posición 50 (regla por posición del paso 7);
- mamíferos: el estado reconstruido en su base, que no se puede estimar sin grupo externo, así que
  se toma como el conjunto de estados presentes en sus dos ramas basales (S en monotremas y T en
  el ancestro de Theria, posterior 0,86, paso 8).

Fitch sobre residuos asigna el estado de Amniota y cuenta los cambios. El resultado depende de los
pocos anfibios con exón 2 alineable, y así se informa.

Uso: python3 scripts/13_tetrapod_parsimony.py  → results/tetrapod_parsimony.txt
"""
import csv
import importlib.util
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("tm", ROOT / "scripts/07_site_table.py")
tm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tm)

SITE = 50
AMPHIBIA = {"Anura", "Caudata", "Gymnophiona"}


def fitch(node, states):
    """Primer paso de Fitch: devuelve (conjunto, cambios) para un árbol de tuplas anidadas."""
    if isinstance(node, str):
        return states[node], 0
    sets, cost = [], 0
    for child in node:
        s, c = fitch(child, states)
        sets.append(s)
        cost += c
    inter = set.intersection(*sets)
    return (inter, cost) if inter else (set.union(*sets), cost + 1)


def main():
    _, _, _, human = tm.mammal_data()
    tetra, groups = tm.tetrapod_residues(human)
    observed = {}
    for g in ("Lepidosauria", "Testudines", "Crocodylia", "Aves"):
        observed[g] = Counter(next(iter(d[SITE])) for sp, d in tetra.items() if groups[sp] == g and SITE in d)
    observed["Amphibia"] = Counter(next(iter(d[SITE])) for sp, d in tetra.items()
                                   if groups[sp] in AMPHIBIA and SITE in d)
    states = {g: set(c) for g, c in observed.items()}
    states["Mammalia"] = {"S", "T"}

    tree = ("Amphibia", ("Mammalia", ("Lepidosauria", ("Testudines", ("Crocodylia", "Aves")))))
    amniota = ("Mammalia", ("Lepidosauria", ("Testudines", ("Crocodylia", "Aves"))))
    sauropsida = ("Lepidosauria", ("Testudines", ("Crocodylia", "Aves")))
    root_set, cost = fitch(tree, states)
    lines = [f"Posición {SITE}, residuos observados por grupo:"]
    lines += [f"  {g}: {dict(c)} (n = {sum(c.values())})" for g, c in observed.items()]
    lines.append("  Mammalia (base): {S (monotremas), T (Theria, PP 0,86)}")
    lines.append(f"Fitch: Tetrapoda {sorted(root_set)}, Amniota {sorted(fitch(amniota, states)[0])}, "
                 f"Sauropsida {sorted(fitch(sauropsida, states)[0])}; cambios = {cost}")
    out = "\n".join(lines)
    (ROOT / "results/tetrapod_parsimony.txt").write_text(out + "\n")
    print(out)


if __name__ == "__main__":
    main()
