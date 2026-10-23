"""Paso 2: descarga proteínas, CDS y reporte de productos (exones) de todos los ortólogos de mamíferos."""
import csv, json, io, zipfile, urllib.request

API = "https://api.ncbi.nlm.nih.gov/datasets/v2/gene"
ids = [r["gene_id"] for r in csv.DictReader(open("data/raw/mammal_orthologs.tsv"), delimiter="\t")]

def post(url, body):
    req = urllib.request.Request(url, data=json.dumps(body).encode(),
                                 headers={"Content-Type": "application/json"})
    return urllib.request.urlopen(req, timeout=600).read()

import os
if not os.path.exists("data/raw/mammal_protein.faa"):
  z = zipfile.ZipFile(io.BytesIO(post(f"{API}/download",
        {"gene_ids": [int(i) for i in ids], "include_annotation_type": ["FASTA_PROTEIN", "FASTA_CDS"]})))
  for name in z.namelist():
    if name.endswith((".faa", ".fna")):
        out = "data/raw/mammal_" + name.split("/")[-1]
        open(out, "wb").write(z.read(name))
        print(out)

reports = []
for i in range(0, len(ids), 50):
    url = f"{API}/id/{','.join(ids[i:i + 50])}/product_report?page_size=1000"
    reports += json.loads(urllib.request.urlopen(url, timeout=600).read())["reports"]
json.dump(reports, open("data/raw/mammal_product_reports.json", "w"))
print(len(reports), "reportes de producto")
