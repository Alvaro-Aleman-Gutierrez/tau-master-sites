"""Paso 1: ortólogos de MAPT (NCBI Gene 4137) restringidos a Mammalia, con su linaje."""
import json, sys, time, urllib.request
from Bio import Entrez

Entrez.email = "your.name@example.org"  # NCBI pide un correo de contacto; poner el propio
orth = json.load(open("data/raw/mapt_orthologs_ncbi.json"))["reports"]
genes = {str(r["gene"]["tax_id"]): r["gene"] for r in orth}

ranks = ("superorder", "order", "family")
rows = []
ids = list(genes)
for i in range(0, len(ids), 200):
    h = Entrez.efetch(db="taxonomy", id=",".join(ids[i:i + 200]), retmode="xml")
    for t in Entrez.read(h):
        lin = {x["Rank"]: x["ScientificName"] for x in t["LineageEx"]}
        if lin.get("class") != "Mammalia":
            continue
        g = genes[t["TaxId"]]
        rows.append([t["TaxId"], t["ScientificName"], lin.get("superorder", ""),
                     lin.get("order", ""), lin.get("family", ""), g["gene_id"],
                     str(g.get("protein_count", "")), str(g.get("transcript_count", ""))])
    time.sleep(0.4)

with open("data/raw/mammal_orthologs.tsv", "w") as f:
    f.write("taxid\tspecies\tsuperorder\torder\tfamily\tgene_id\tn_proteins\tn_transcripts\n")
    for r in sorted(rows, key=lambda r: (r[2], r[3], r[1])):
        f.write("\t".join(r) + "\n")
print(len(rows), "mamíferos")
