#!/usr/bin/env bash
# Reproduce todo el análisis desde los datos crudos de NCBI en data/raw/.
# Requiere Python 3 (Biopython, pandas, scipy, matplotlib, openpyxl), MAFFT, IQ-TREE 2 y HyPhy.
# Variables opcionales: THREADS (por defecto 10) y HYPHY (ruta al ejecutable de HyPhy).
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p results figures tables
FROM="${1:-1}"   # paso desde el que empezar (1 a 14)

step() { local n=$1; shift; if [ "$n" -ge "$FROM" ]; then echo "== paso $n: $*"; "$@"; fi; }

step 1  python3 scripts/01_human_exons.py
step 2  python3 scripts/02_select_2N4R.py
step 3  python3 scripts/03_alignment.py
step 4  python3 scripts/04_taxonomy_constraint.py
step 5  bash scripts/05_phylogeny.sh
step 6  python3 scripts/06_selection.py
step 7  python3 scripts/07_site_table.py
step 8  python3 scripts/08_ancestral_states.py
step 9  python3 scripts/09_figures.py
step 10 python3 scripts/10_tables.py
step 11 python3 scripts/11_graphical_abstract.py
step 12 python3 scripts/12_isoform_sensitivity.py
step 13 python3 scripts/13_tetrapod_parsimony.py
step 14 python3 scripts/14_figure3_schematic.py
echo "Listo. Resultados en results/, figuras en figures/ y tablas en tables/."
