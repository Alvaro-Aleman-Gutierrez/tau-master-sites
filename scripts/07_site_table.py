"""Paso 7: tabla maestra de sitios fosforilables de tau humana 2N4R.

Escribe results/site_table.tsv con una fila por S/T/Y de la tau humana 2N4R:
categoría bibliográfica, exón, región, motivo, conservación en 259 mamíferos,
tasas FEL/FUBAR por codón y conservación en tetrápodos no mamíferos.

Uso: python3 scripts/07_site_table.py  (después de los pasos 1 a 6)
"""
import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median

from Bio import Align
from Bio.Align import substitution_matrices

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)

ALN = OUT / "mapt_2N4R_final.aln.faa"
SPECIES = OUT / "selected_2N4R.tsv"
EXONS = OUT / "human_exons.tsv"
FEL = OUT / "selection/all.FEL.json"
FUBAR = OUT / "selection/all.FUBAR.json"
LINEAGES = ROOT / "data/raw/ncbi_lineages.tsv"
TETRA_PROT = ROOT / "data/raw/tetrapod_protein.faa"
TETRA_SP = ROOT / "data/raw/tetrapod_orthologs.tsv"
HANGER = ROOT / "data/literature/hanger2007_table_s10.tsv"
WESSELING = ROOT / "data/literature/wesseling2020_table_s4_frequencies.tsv"

EXONS_2N4R = ["E1", "E2", "E3", "E4", "E5", "E7", "E9", "E10", "E11", "E12", "E13"]
WINDOW = 5  # residuos a cada lado para el control de homología local
TETRA_MIN_WINDOW_ID = 0.5

# Regiones de 2N4R según Trushina et al. 2019 (Fig. 2): NTR 1–171, PRR 172–243, MBR 244–368, CTR 369–441.
REGIONS = [(1, 171, "NTR"), (172, 243, "PRR"), (244, 368, "MBR"), (369, 441, "CTR")]

# Categorías bibliográficas. Solo lo verificado en el texto de cada fuente.
CATEGORIES = {
    50: "master site (Stefanoska 2022)",
    69: "master site (Stefanoska 2022)",
    181: "master site (Stefanoska 2022); blood biomarker (Karikari 2020)",
    111: "weaker master site (Stefanoska 2022)",
    205: "weaker master site (Stefanoska 2022); protective (Trushina 2019)",
    199: "interdependence target (Stefanoska 2022)",
    212: "interdependence target (Stefanoska 2022)",
    231: "interdependence target (Stefanoska 2022); plasma biomarker (Ashton 2021)",
    235: "interdependence target (Stefanoska 2022)",
    396: "interdependence target (Stefanoska 2022)",
    404: "interdependence target (Stefanoska 2022)",
    422: "interdependence target (Stefanoska 2022)",
    113: "low conservation in vertebrates (Trushina 2019)",
    184: "low conservation in vertebrates (Trushina 2019)",
    217: "plasma biomarker (Palmqvist 2020)",
}


def hanger():
    """Hanger et al. 2007, Table S10 (data/literature): PHF-tau y cerebro control por posición."""
    phf, ctl = {}, {}
    for r in csv.DictReader(open(HANGER), delimiter="\t"):
        if r["PHF_tau"]:
            phf[int(r["position"])] = r["PHF_tau"]
        if r["control_brain"]:
            ctl[int(r["position"])] = r["control_brain"]
    return phf, ctl


def wesseling_frequencies():
    """Porcentaje de pacientes con Alzheimer (n = 49) y controles (n = 43) con cada fosfositio.

    Wesseling et al. 2020, Table S4. Solo entradas de un único residuo; si el sitio aparece en
    varias isoformas o péptidos (por ejemplo S113, S113_0N y S113_1N) se toma el máximo. Las
    entradas con dos residuos posibles ("T111orS113") o dos fosfatos se omiten.
    """
    out = {}
    for r in csv.DictReader(open(WESSELING), delimiter="\t"):
        h, a, c = r["entry"], float(r["AD_pct"]), float(r["control_pct"])
        if "Phospho" not in h or "or" in h:
            continue
        m = re.fullmatch(r"([STY])(\d+)_(?:\dN_)?Phospho(?:_\w+)?", h)
        if not m:
            continue
        pos = int(m.group(2))
        prev = out.get(pos, (0.0, 0.0))
        out[pos] = (max(prev[0], a), max(prev[1], c))
    return out


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:]
            seqs[name] = []
        elif name:
            seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def region(pos):
    return next(r for a, b, r in REGIONS if a <= pos <= b)


def exon_map(human):
    rows = {r["exon"]: r["peptide"] for r in csv.DictReader(open(EXONS), delimiter="\t")}
    pos2exon, seq = {}, ""
    for ex in EXONS_2N4R:
        for k in range(len(rows[ex])):
            pos2exon[len(seq) + k + 1] = ex
        seq += rows[ex]
    # Los péptidos se tradujeron en el marco de cada transcrito, así que en un codón
    # partido entre exones el residuo puede diferir; solo se exige el largo y como
    # mucho una diferencia por unión.
    assert len(seq) == len(human), "los exones no suman el largo de la 2N4R humana"
    assert sum(a != b for a, b in zip(seq, human)) < len(EXONS_2N4R)
    return pos2exon


def mammal_data():
    aln = {k.split()[0]: v for k, v in read_fasta(ALN).items()}
    order = {r["species"].replace(" ", "_"): r["order"]
             for r in csv.DictReader(open(SPECIES), delimiter="\t")}
    human_aln = aln["Homo_sapiens"]
    col, i = {}, 0
    for c, a in enumerate(human_aln):
        if a != "-":
            i += 1
            col[i] = c
    return aln, order, col, human_aln.replace("-", "")


def window_identity(seq, human_aln, cols):
    same = total = 0
    for c in cols:
        if seq[c] != "-" and human_aln[c] != "-":
            total += 1
            same += seq[c] == human_aln[c]
    return same / total if total else 0.0


def selection(col):
    fel = json.load(open(FEL))["MLE"]["content"]["0"]
    fub = json.load(open(FUBAR))["MLE"]["content"]["0"]
    out = {}
    for pos, c in col.items():
        alpha, beta, _, _, p, _ = fel[c][:6]
        out[pos] = {"FEL_alpha": alpha, "FEL_beta": beta, "FEL_p": p,
                    "FUBAR_beta_menos_alpha": fub[c][2], "FUBAR_P_neg": fub[c][3]}
    return out


def tetrapod_residues(human):
    """Residuo homólogo a cada posición humana en cada especie de tetrápodo.

    Cada isoforma se alinea globalmente con la 2N4R humana. Las isoformas con
    extremos N alternativos se alinean a la fuerza sobre los exones 2 y 3 con
    secuencia no homóloga, así que para cada posición se toma la isoforma con
    mayor identidad local (±WINDOW residuos, gaps como diferencia) y solo si esa
    identidad es al menos TETRA_MIN_WINDOW_ID.
    """
    groups = {r["species"]: r["group"] for r in csv.DictReader(open(TETRA_SP), delimiter="\t")}
    aligner = Align.PairwiseAligner(mode="global",
                                    substitution_matrix=substitution_matrices.load("BLOSUM62"),
                                    open_gap_score=-10, extend_gap_score=-0.5,
                                    end_gap_score=0)
    best = defaultdict(dict)  # especie -> posición -> (identidad local, residuo)
    for header, seq in read_fasta(TETRA_PROT).items():
        sp = header.split("[organism=")[1].split("]")[0]
        if sp not in groups:
            continue
        a = aligner.align(human, seq.rstrip("*"))[0]
        res = {}
        for (h0, h1), (s0, s1) in zip(*a.aligned):
            for k in range(h1 - h0):
                res[h0 + k + 1] = seq[s0 + k]
        for pos in res:
            win = [p for p in range(pos - WINDOW, pos + WINDOW + 1) if p != pos and 1 <= p <= len(human)]
            wid = sum(res.get(p) == human[p - 1] for p in win) / len(win)
            if wid >= TETRA_MIN_WINDOW_ID and wid > best[sp].get(pos, (-1, ""))[0]:
                best[sp][pos] = (wid, res[pos])
    per_species = {sp: {pos: {r} for pos, (_, r) in d.items()} for sp, d in best.items()}
    return per_species, groups


def main():
    aln, order, col, human = mammal_data()
    pos2exon = exon_map(human)
    sel = selection(col)
    tetra, tgroups = tetrapod_residues(human)
    wess = wesseling_frequencies()
    hanger_phf, hanger_ctl = hanger()
    n_tetra = len(tetra)

    rows = []
    for pos, aa in enumerate(human, 1):
        if aa not in "STY":
            continue
        c = col[pos]
        residues = {sp: s[c] for sp, s in aln.items()}
        cnt = Counter(residues.values())
        phospho = "ST" if aa in "ST" else "Y"
        by_order = defaultdict(Counter)
        for sp, r in residues.items():
            if r != aa:
                by_order[order.get(sp, "?")][r] += 1
        wcols = [col[p] for p in range(max(1, pos - WINDOW), min(len(human), pos + WINDOW) + 1) if p != pos]
        wid = [window_identity(s, aln["Homo_sapiens"], wcols) for s in aln.values()]

        t_same = t_phos = t_with = 0
        t_groups = defaultdict(lambda: [0, 0])
        for sp, sites in tetra.items():
            res = sites.get(pos, set())
            if not res:
                continue
            t_with += 1
            t_same += aa in res
            t_phos += bool(res & set(phospho))
            t_groups[tgroups[sp]][0] += 1
            t_groups[tgroups[sp]][1] += aa in res

        rows.append({
            "sitio": f"{aa}{pos}", "pos": pos, "residuo": aa,
            "categoria": CATEGORIES.get(pos, "background"),
            "exon": pos2exon[pos], "region": region(pos),
            "motivo_SP_TP": "yes" if pos < len(human) and human[pos] == "P" else "no",
            "mam_n": len(residues),
            "mam_identico": cnt[aa],
            "mam_fosforilable": sum(cnt[x] for x in phospho),
            "mam_gap": cnt["-"],
            "mam_otros": ";".join(f"{r}:{n}" for r, n in cnt.most_common() if r != aa),
            "mam_cambios_por_orden": ";".join(
                f"{o}({','.join(f'{r}{n}' for r, n in c2.most_common())})"
                for o, c2 in sorted(by_order.items(), key=lambda x: -sum(x[1].values()))),
            "ventana_id_mediana": round(median(wid), 3),
            "ventana_n_menor_0.5": sum(w < 0.5 for w in wid),
            **{k: round(v, 4) for k, v in sel[pos].items()},
            "PHF_Hanger2007": hanger_phf.get(pos, ""),
            "control_Hanger2007": hanger_ctl.get(pos, ""),
            "Wesseling_AD_pct": wess[pos][0] if pos in wess else "",
            "Wesseling_CTR_pct": wess[pos][1] if pos in wess else "",
            "tetra_n_con_posicion": t_with,
            "tetra_identico": t_same,
            "tetra_fosforilable": t_phos,
            "tetra_por_grupo": ";".join(f"{g}:{v[1]}/{v[0]}" for g, v in sorted(t_groups.items())),
        })

    with open(OUT / "site_table.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} sitios S/T/Y; {n_tetra} especies de tetrápodos -> {OUT / 'site_table.tsv'}")


if __name__ == "__main__":
    main()
