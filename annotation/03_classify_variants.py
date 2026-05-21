"""
03_classify_variants.py — Classify variants into functional categories

Parses VEP output (canonical, protein-coding transcripts), joins MPC scores,
and assigns each variant to one functional category:

  PTV   protein-truncating, LOFTEE high-confidence
  MisB  missense, MPC >= 2
  MisA  missense, 1 <= MPC < 2
  Mis0  missense, MPC < 1
  SYN   synonymous (calibration set)

Input:
  - vep_output.txt      (from step 02)
  - mpc.GRCh38.txt.gz   (MPC values; Samocha et al. 2017)

Output (variant ID + gene per line):
  - variants_PTV.txt
  - variants_MisB.txt
  - variants_MisA.txt
  - variants_Mis0.txt
  - variants_SYN.txt

NOTE: VEP output column layout varies by version/options. Adjust the header
handling and the MPC join key below to match your files.
"""

import pandas as pd

PTV_CONSEQUENCES = {
    "stop_gained", "frameshift_variant", "splice_donor_variant",
    "splice_acceptor_variant", "transcript_ablation",
}

def parse_extra(extra):
    out = {}
    for kv in str(extra).split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
    return out

# --- Load VEP output (skip '##' metadata lines) ---
vep = pd.read_csv("vep_output.txt", sep="\t", comment="#", header=None)

cols = ["Uploaded_variation", "Location", "Allele", "Gene", "Feature",
        "Feature_type", "Consequence", "cDNA_position", "CDS_position",
        "Protein_position", "Amino_acids", "Codons", "Existing_variation", "Extra"]

vep.columns = cols[:vep.shape[1]]

extra = vep["Extra"].apply(parse_extra)
vep["CANONICAL"] = extra.apply(lambda d: d.get("CANONICAL"))
vep["LoF"]       = extra.apply(lambda d: d.get("LoF"))
vep["SYMBOL"]    = extra.apply(lambda d: d.get("SYMBOL"))

# Canonical transcripts only
vep = vep[vep["CANONICAL"] == "YES"].copy()

# Primary consequence (VEP can list several, comma-separated)
vep["primary_csq"] = vep["Consequence"].str.split(",").str[0]

# --- Join MPC scores (gzip auto-detected by .gz extension) ---
mpc = pd.read_csv("mpc.GRCh38.txt.gz", sep="\t")  # adjust columns to match your file
vep = vep.merge(mpc, left_on="Uploaded_variation", right_on="variant", how="left")

is_ptv = vep["primary_csq"].isin(PTV_CONSEQUENCES) & (vep["LoF"] == "HC")
is_mis = vep["primary_csq"] == "missense_variant"
is_syn = vep["primary_csq"] == "synonymous_variant"

categories = {
    "PTV":  vep[is_ptv],
    "MisB": vep[is_mis & (vep["MPC"] >= 2)],
    "MisA": vep[is_mis & (vep["MPC"] >= 1) & (vep["MPC"] < 2)],
    "Mis0": vep[is_mis & (vep["MPC"] < 1)],
    "SYN":  vep[is_syn],
}

for name, df in categories.items():
    out = df[["Uploaded_variation", "SYMBOL", "Gene"]].drop_duplicates()
    out.to_csv(f"variants_{name}.txt", sep="\t", index=False)
    print(f"{name}: {len(out)} variants")
