# Raw data from NCBI

Downloaded on 21 October 2026. The files are kept as downloaded so that the analysis gives the same result if NCBI later updates its annotations. The scripts in `scripts/download/` regenerate them from the current NCBI records, which may differ.

| File | Content | Source |
|---|---|---|
| `mapt_orthologs_ncbi.json` | Orthologs of human *MAPT* (NCBI Gene 4137) | NCBI Datasets v2, `gene/id/4137/orthologs` |
| `mammal_orthologs.tsv` | Mammalian orthologs with taxid, superorder, order and family (267 species) | Built by `d1_mammal_orthologs.py` from the JSON and NCBI Taxonomy |
| `mammal_protein.faa`, `mammal_cds.fna` | All annotated protein isoforms and their coding sequences (3 571 isoforms) | NCBI Datasets v2, gene download |
| `mammal_product_reports.json` | Transcript and exon coordinates of each isoform | NCBI Datasets v2, `product_report` |
| `human_MAPT_region.fasta` | GRCh38 NC_000017.11:45894000-46030000, used to translate each human exon | NCBI Nucleotide |
| `ncbi_lineages.tsv` | Full NCBI Taxonomy lineage below Mammalia for the 259 species analysed | NCBI Taxonomy (efetch) |
| `tetrapod_orthologs.tsv` | Non-mammalian tetrapod orthologs with their group (240 species) | Built by `d3_tetrapod_orthologs.py` |
| `tetrapod_protein.faa` | All annotated protein isoforms of those species (3 652 isoforms) | NCBI Datasets v2, gene download |

NCBI asks for a contact e-mail in Entrez queries. The download scripts carry a placeholder that must be replaced before running them.
