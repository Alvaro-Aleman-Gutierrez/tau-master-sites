"""Paso 3: filtro final del conjunto 2N4R, alineamiento proteico (MAFFT L-INS-i) y alineamiento de codones."""
import subprocess
import pandas as pd
from Bio import SeqIO

MIN_CORE, MIN_LEN, MAX_LEN = 0.85, 380, 480

sel = pd.read_csv("results/selected_2N4R.tsv", sep="\t")
def reason(r):
    if r.source == "excluida":
        return "sin isoforma con exones core completos editable en marco"
    if r.core_coverage < MIN_CORE:
        return f"modelo génico incompleto (cobertura core {r.core_coverage:.2f})"
    if not MIN_LEN <= r.length <= MAX_LEN:
        return f"longitud atípica ({int(r.length)} aa): extensión no homóloga en el modelo génico"
    return ""
sel["exclusion_reason"] = sel.apply(reason, axis=1)
sel["included"] = sel.exclusion_reason == ""
sel.to_csv("results/selected_2N4R.tsv", sep="\t", index=False)
keep = set(sel[sel.included].species.str.replace(" ", "_"))

prot = [r for r in SeqIO.parse("results/mapt_2N4R_protein.faa", "fasta") if r.id in keep]
nuc = {r.id: str(r.seq) for r in SeqIO.parse("results/mapt_2N4R_cds.fna", "fasta") if r.id in keep}
SeqIO.write(prot, "results/mapt_2N4R_final.faa", "fasta")
with open("results/mapt_2N4R_final.aln.faa", "w") as out:
    subprocess.run(["mafft", "--localpair", "--maxiterate", "1000", "--thread", "8", "--quiet",
                    "results/mapt_2N4R_final.faa"], stdout=out, check=True)

# alineamiento de codones guiado por la proteína (equivalente a PAL2NAL)
with open("results/mapt_2N4R_final.codon.fna", "w") as out:
    for r in SeqIO.parse("results/mapt_2N4R_final.aln.faa", "fasta"):
        nt, i, cod = nuc[r.id], 0, []
        for aa in str(r.seq):
            if aa == "-":
                cod.append("---")
            else:
                cod.append(nt[i:i + 3]); i += 3
        assert i == len(nt) or nt[i:] in ("TAA", "TAG", "TGA"), r.id
        out.write(f">{r.id}\n{''.join(cod)}\n")
print(len(prot), "especies incluidas;", (~sel.included).sum(), "excluidas")
print(sel[~sel.included][["species", "order", "exclusion_reason"]].to_string(index=False))
