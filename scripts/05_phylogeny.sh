#!/usr/bin/env bash
# Paso 5: filogenia ML de la proteína 2N4R.
# 5a. Árbol sin restricción: ModelFinder (BIC), 1000 réplicas de UFBoot (--bnni) y 1000 de SH-aLRT.
#     Solo se usa para informar el soporte de los clados discutidos.
# 5b. Árbol con la restricción taxonómica del paso 4 y el modelo elegido en 5a. Es la topología
#     de FEL, FUBAR y la reconstrucción ancestral.
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p results/phylo_protein
THREADS="${THREADS:-10}"

iqtree2 -s results/mapt_2N4R_final.aln.faa -st AA -m MFP -B 1000 --bnni --alrt 1000 \
        -T AUTO --threads-max "$THREADS" --seed 12345 --prefix results/phylo_protein/mapt_2N4R -redo

MODEL=$(grep -m1 "Best-fit model:" results/phylo_protein/mapt_2N4R.log | sed -E "s/Best-fit model: ([^ ]+) chosen.*/\1/")
echo "Modelo elegido por ModelFinder: $MODEL"

iqtree2 -s results/mapt_2N4R_final.aln.faa -st AA -m "$MODEL" -g results/phylo_protein/taxonomy_constraint.nwk \
        -T AUTO --threads-max "$THREADS" --seed 12345 --prefix results/phylo_protein/mapt_2N4R_constrained -redo
