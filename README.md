# Evolution of the master phosphorylation sites of tau

Data and code for the article

> Alemán-Gutiérrez, A. Threonine 50, a master phosphorylation site of tau, is widespread in placental mammals but absent from murid and cricetid rodents. Submitted to *Computational Biology and Chemistry*.

The repository rebuilds every result, table and figure of the article from public sequence data. It traces the 85 serines, threonines and tyrosines of human 2N4R tau (NP_005901.2) across 259 mammals and 240 non-mammalian tetrapods, estimates codon-level selection with FEL and FUBAR, reconstructs ancestral states on a species tree and annotates each site with published evidence of its function and of its phosphorylation in Alzheimer disease.

## Contents

| Folder | What it holds |
|---|---|
| `data/raw/` | *MAPT* orthologs, protein and coding sequences, gene product reports and taxonomy as downloaded from NCBI on 21 October 2026 (see `data/raw/README.md`) |
| `data/literature/` | Phosphorylation annotations transcribed from Hanger et al. (2007) and Wesseling et al. (2020) |
| `scripts/` | The analysis, numbered in the order it runs |
| `scripts/download/` | The scripts that produced `data/raw/` (not run by default) |
| `results/` | Intermediate and final results (created by the scripts) |
| `figures/` | Figures 1 and 2 and the graphical abstract (EPS, PDF, TIFF, PNG) |
| `tables/` | Table 1, Table S1 and Table S2 |

## Running the analysis

```bash
./run_all.sh          # every step, from data/raw/ to figures/ and tables/
./run_all.sh 7        # from step 7 onwards, reusing the trees and HyPhy output
```

Optional environment variables are `THREADS` (default 10) and `HYPHY`, the path to the HyPhy executable (default `hyphy` in the `PATH`). A full run takes about two hours on ten threads, most of it in steps 5 and 6.

| Step | Script | Output |
|---|---|---|
| 1 | `01_human_exons.py` | Human coding exons of *MAPT* and exon coverage of every mammalian isoform |
| 2 | `02_select_2N4R.py` | One 2N4R sequence per species (annotated, or edited in silico to the 2N4R exon set) |
| 3 | `03_alignment.py` | Final set of 259 species, MAFFT L-INS-i protein alignment and codon alignment |
| 4 | `04_taxonomy_constraint.py` | Constraint tree from NCBI Taxonomy |
| 5 | `05_phylogeny.sh` | Unconstrained ML tree with ModelFinder, UFBoot and SH-aLRT, and ML tree under the taxonomic constraint |
| 6 | `06_selection.py` | FEL and FUBAR on the codon alignment |
| 7 | `07_site_table.py` | `results/site_table.tsv`, one row per S/T/Y with annotations, conservation in mammals and tetrapods, and codon rates |
| 8 | `08_ancestral_states.py` | Marginal ancestral states (IQ-TREE) at selected nodes |
| 9 | `09_figures.py` | Figures 1 and 2 |
| 10 | `10_tables.py` | Table 1 and Table S1 |
| 11 | `11_graphical_abstract.py` | Graphical abstract |
| 12 | `12_isoform_sensitivity.py` | Table S2, the sensitivity of tetrapod counts to the isoform selection rule |

The column names of `results/site_table.tsv` are in Spanish. `tables/Table_S1.tsv` holds the same data with English column names.

## Software

Python 3.11 with Biopython 1.83, pandas 2.2.2, SciPy 1.13.1, matplotlib 3.8.4 and openpyxl 3.1.5; MAFFT 7.490; IQ-TREE 2.0.7; HyPhy 2.5.103. Figures use Arial and Courier New.

## Licence

Code is released under the MIT licence (`LICENSE`). The sequence data in `data/raw/` are public NCBI records and keep their original terms. The files in `data/literature/` are transcriptions of published supplementary tables and should be cited through the original articles.

## Citation

Please cite the article and the archived release of this repository (Zenodo DOI to be added on release).
