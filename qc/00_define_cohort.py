"""
00_define_cohort.py — Define postpartum psychosis cases and controls

Reads pre-extracted AoU person and condition tables (built via the AoU
Researcher Workbench cohort builder), filters cases to psychiatric diagnoses
within 112 days (12 weeks) of delivery, and produces a unified phenotype file.

Input (place under ${WORKSPACE_BUCKET}/data/):
  - pps_cases_condition.csv : condition table for candidate PPS cases
  - pps_cases_person.csv    : person table for candidate PPS cases
  - controls_person.csv     : person table for parous controls
  - existing_pps_cases.csv  : previously curated PPS cases (optional)

Output:
  - cohort_phenotype.csv : person_id, phenotype (1=case, 0=control)

Case definition (paper Methods, Table S1):
  Mania or psychosis ICD codes within 0-3 months postpartum, including
  brief psychotic disorder, manic episodes, bipolar disorder, severe
  depressive episodes with psychotic features, and puerperal psychosis.
"""

import os
import pandas as pd

BUCKET = os.getenv("WORKSPACE_BUCKET")
DATA_DIR = f"{BUCKET}/data"

case_conditions = pd.read_csv(f"{DATA_DIR}/pps_cases_condition.csv")
controls        = pd.read_csv(f"{DATA_DIR}/controls_person.csv")

# Cases: psychiatric diagnosis within 112 days of delivery
deliveries = case_conditions[case_conditions["standard_concept_name"] == "Delivery normal"]
psych_dx   = case_conditions[case_conditions["standard_concept_name"] != "Delivery normal"]

paired = pd.merge(deliveries, psych_dx, on="person_id",
                  suffixes=("_delivery", "_dx"))

for col in ["condition_start_datetime_delivery", "condition_start_datetime_dx"]:
    paired[col] = pd.to_datetime(paired[col], format="ISO8601")

paired["days_from_delivery"] = (
    paired["condition_start_datetime_dx"] - paired["condition_start_datetime_delivery"]
).abs().dt.days

within_window = paired[paired["days_from_delivery"] <= 112]

# One row per case (earliest psychiatric event in the postpartum window)
idx = within_window.groupby("person_id")["days_from_delivery"].idxmin()
new_cases = within_window.loc[idx, ["person_id"]].drop_duplicates()

print(f"New cases within postpartum window: {new_cases['person_id'].nunique()}")

# Merge with previously curated PPS cases (if available)
existing_path = f"{DATA_DIR}/existing_pps_cases.csv"
if os.path.exists(existing_path):
    existing = pd.read_csv(existing_path)
    existing = existing[existing["standard_concept_name"] == "Postpartum psychosis"]
    existing_cases = existing[["person_id"]].drop_duplicates()
    all_cases = pd.concat([new_cases, existing_cases]).drop_duplicates(subset="person_id")
else:
    all_cases = new_cases

print(f"Total cases: {len(all_cases)}")

# Controls: drop any IDs overlapping with cases
case_ids = set(all_cases["person_id"])
controls_clean = (controls[~controls["person_id"].isin(case_ids)]
                  .drop_duplicates(subset="person_id"))

print(f"Total controls: {len(controls_clean)}")

# Combined phenotype file
all_cases = all_cases.assign(phenotype=1)
controls_clean = controls_clean.assign(phenotype=0)

cohort = pd.concat(
    [all_cases[["person_id", "phenotype"]],
     controls_clean[["person_id", "phenotype"]]],
    ignore_index=True,
)

cohort.to_csv("cohort_phenotype.csv", index=False)
os.system(f"gsutil -m cp cohort_phenotype.csv {DATA_DIR}/cohort_phenotype.csv")

print(f"Wrote cohort_phenotype.csv: "
      f"{(cohort.phenotype==1).sum()} cases, "
      f"{(cohort.phenotype==0).sum()} controls")
