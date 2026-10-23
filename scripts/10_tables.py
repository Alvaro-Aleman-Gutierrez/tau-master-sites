"""Tabla 1 (sitios anotados) y Tabla S1 (los 85 S/T/Y) en inglés, desde results/site_table.tsv.

Uso: python3 scripts/10_tables.py  (después del paso 7)
"""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "tables"
OUT.mkdir(parents=True, exist_ok=True)

TABLE1_SITES = [50, 69, 181, 111, 205, 113, 184, 199, 212, 217, 231, 235, 262, 396, 403, 404, 422]

# Códigos de anotación de la Tabla 1, explicados en la nota al pie.
CODES = {
    50: "M", 69: "M", 181: "M, B", 111: "m", 205: "m", 113: "L", 184: "L",
    199: "I", 212: "I", 217: "B", 231: "I, B", 235: "I", 396: "I", 404: "I", 422: "I",
    262: "", 403: "",
}
FOOTNOTE = (
    "M, master site and m, weaker master site (Stefanoska et al., 2022). I, site whose phosphorylation "
    "depends on prior phosphorylation at other sites in the same study. B, blood or plasma biomarker "
    "(Karikari et al., 2020; Palmqvist et al., 2020; Ashton et al., 2021). L, low conservation across "
    "vertebrates in Trushina et al. (2019). PHF, detection in paired helical filament tau by mass "
    "spectrometry (MS) or only with a phosphospecific antibody (Ab) in Hanger et al. (2007). AD and "
    "control, percentage of 49 Alzheimer disease patients and 43 controls with the site phosphorylated "
    "(Wesseling et al., 2020); a dash means the site was not reported as a single residue. Tetrapods counts "
    "species whose best matching isoform aligns at the site. FEL p tests whether the nonsynonymous rate "
    "differs from the synonymous rate at the codon."
)


def fmt_p(p):
    p = float(p)
    return "<0.0001" if p < 1e-4 else f"{p:.2g}"


def phf(v):
    return {"MS": "MS", "antibody": "Ab", "MS ambiguous": "MS (ambiguous)",
            "MS ambiguous; antibody": "MS (ambiguous), Ab"}.get(v, "")


def pct(v):
    return f"{float(v):.0f}" if v else "–"


def main():
    rows = {int(r["pos"]): r for r in csv.DictReader(open(ROOT / "results/site_table.tsv"), delimiter="\t")}

    header = ["Site", "Exon", "Region", "SP/TP", "Annotation", "Mammals identical (of 259)",
              "Mammals S or T", "Tetrapods identical", "FEL p", "PHF", "AD (%)", "Control (%)"]
    t1 = []
    for pos in TABLE1_SITES:
        r = rows[pos]
        t1.append([r["sitio"], r["exon"].replace("E", ""), r["region"],
                   r["motivo_SP_TP"], CODES[pos],
                   r["mam_identico"], r["mam_fosforilable"],
                   f"{r['tetra_identico']}/{r['tetra_n_con_posicion']}", fmt_p(r["FEL_p"]),
                   phf(r["PHF_Hanger2007"]), pct(r["Wesseling_AD_pct"]), pct(r["Wesseling_CTR_pct"])])

    with open(OUT / "Table_1.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        w.writerows(t1)
    with open(OUT / "Table_1.md", "w") as fh:
        fh.write("**Table 1. Annotated phosphorylation sites of human 2N4R tau across tetrapods.**\n\n")
        fh.write("| " + " | ".join(header) + " |\n|" + "---|" * len(header) + "\n")
        for row in t1:
            fh.write("| " + " | ".join(row) + " |\n")
        fh.write("\n" + FOOTNOTE + "\n")

    s1_cols = {
        "sitio": "site", "pos": "position", "residuo": "residue", "exon": "exon", "region": "region",
        "motivo_SP_TP": "SP_TP_motif", "mam_n": "mammals_n", "mam_identico": "mammals_identical",
        "mam_fosforilable": "mammals_S_or_T", "mam_gap": "mammals_gap", "mam_otros": "mammals_other_residues",
        "mam_cambios_por_orden": "mammals_changes_by_order", "ventana_id_mediana": "local_identity_median",
        "ventana_n_menor_0.5": "local_identity_below_0.5", "FEL_alpha": "FEL_dS", "FEL_beta": "FEL_dN",
        "FEL_p": "FEL_p", "FUBAR_beta_menos_alpha": "FUBAR_dN_minus_dS", "FUBAR_P_neg": "FUBAR_P_negative",
        "PHF_Hanger2007": "PHF_Hanger2007", "control_Hanger2007": "control_Hanger2007",
        "Wesseling_AD_pct": "Wesseling_AD_pct", "Wesseling_CTR_pct": "Wesseling_control_pct",
        "tetra_n_con_posicion": "tetrapods_n", "tetra_identico": "tetrapods_identical",
        "tetra_fosforilable": "tetrapods_S_or_T", "tetra_por_grupo": "tetrapods_identical_by_group",
    }
    translate = {"antibody": "Ab", "MS ambiguous; antibody": "MS ambiguous; Ab"}
    with open(OUT / "Table_S1.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(list(s1_cols.values()))
        for pos in sorted(rows):
            r = rows[pos]
            w.writerow([translate.get(r[k], r[k]) for k in s1_cols])
    print(f"Table_1: {len(t1)} sitios; Table_S1: {len(rows)} sitios -> {OUT}")


if __name__ == "__main__":
    main()
