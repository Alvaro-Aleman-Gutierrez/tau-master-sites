"""Paso 14: ortólogos de MAPT fuera de Mammalia (grupo externo para el exón 10).

El conjunto de ortólogos de NCBI no incluye peces con aleta radiada, porque la duplicación
teleósteo produjo mapta/maptb y el pipeline no los declara ortólogos. El grupo externo queda
entonces en tetrápodos más el pez pulmonado.
"""
import csv, io, json, os, time, urllib.request, zipfile
from Bio import Entrez

Entrez.email = "your.name@example.org"  # NCBI pide un correo de contacto; poner el propio
API = "https://api.ncbi.nlm.nih.gov/datasets/v2/gene"
orth = json.load(open("data/raw/mapt_orthologs_ncbi.json"))["reports"]
genes = {str(r["gene"]["tax_id"]): r["gene"] for r in orth}

# grupo por el primer clado de la lista que aparezca en el linaje
GROUPS = ["Aves", "Crocodylia", "Testudines", "Lepidosauria", "Anura", "Caudata",
          "Gymnophiona", "Dipnoi"]

rows, ids = [], list(genes)
for i in range(0, len(ids), 200):
    for t in Entrez.read(Entrez.efetch(db="taxonomy", id=",".join(ids[i:i + 200]), retmode="xml")):
        lin = {x["Rank"]: x["ScientificName"] for x in t["LineageEx"]}
        names = {x["ScientificName"] for x in t["LineageEx"]}
        if lin.get("class") == "Mammalia":
            continue
        grp = next((g for g in GROUPS if g in names), "")
        if not grp:
            continue
        g = genes[t["TaxId"]]
        rows.append([t["TaxId"], t["ScientificName"], grp, lin.get("order", ""),
                     lin.get("family", ""), g["gene_id"]])
    time.sleep(0.4)

rows.sort(key=lambda r: (GROUPS.index(r[2]), r[1]))
with open("data/raw/tetrapod_orthologs.tsv", "w") as f:
    f.write("taxid\tspecies\tgroup\torder\tfamily\tgene_id\n")
    for r in rows:
        f.write("\t".join(map(str, r)) + "\n")
print(len(rows), "ortólogos fuera de Mammalia")

def post(url, body):
    req = urllib.request.Request(url, data=json.dumps(body).encode(),
                                 headers={"Content-Type": "application/json"})
    return urllib.request.urlopen(req, timeout=900).read()

gids = [r[5] for r in rows]
if not os.path.exists("data/raw/tetrapod_protein.faa"):
    z = zipfile.ZipFile(io.BytesIO(post(f"{API}/download",
          {"gene_ids": [int(i) for i in gids],
           "include_annotation_type": ["FASTA_PROTEIN", "FASTA_CDS"]})))
    for name in z.namelist():
        if name.endswith((".faa", ".fna")):
            out = "data/raw/tetrapod_" + name.split("/")[-1]
            open(out, "wb").write(z.read(name))
            print(out)

reports = []
for i in range(0, len(gids), 50):
    url = f"{API}/id/{','.join(gids[i:i + 50])}/product_report?page_size=1000"
    reports += json.loads(urllib.request.urlopen(url, timeout=900).read())["reports"]
json.dump(reports, open("data/raw/tetrapod_product_reports.json", "w"))
print(len(reports), "reportes de producto")
