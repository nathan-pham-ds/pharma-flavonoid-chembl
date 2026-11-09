"""fetch_chembl.py: re-pull xanthine oxidase IC50 data from the ChEMBL REST API.

Targets: CHEMBL1929 (human XDH/XO) and CHEMBL3649 (bovine XO; most natural-product papers use bovine milk XO).
Usage: python src/fetch_chembl.py   ->  data/raw/xo_ic50_activities.csv, xo_molecules.csv, chembl_funnel_counts.json
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import requests

API = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
TARGETS = ["CHEMBL1929", "CHEMBL3649"]
FIELDS = ("molecule_chembl_id,canonical_smiles,standard_value,assay_chembl_id,assay_type,target_chembl_id,"
          "document_year,data_validity_comment,potential_duplicate,document_chembl_id")


def total(params: dict) -> int:
    return requests.get(API, params={**params, "limit": 1}, timeout=60).json()["page_meta"]["total_count"]


def main() -> None:
    out = Path("data/raw"); out.mkdir(parents=True, exist_ok=True)
    funnel, rows = {}, []
    for t in TARGETS:
        base = {"target_chembl_id": t}
        funnel[t] = dict(all_activities=total(base), ic50=total({**base, "standard_type": "IC50"}),
                         ic50_nM=total({**base, "standard_type": "IC50", "standard_units": "nM"}),
                         ic50_nM_eq=total({**base, "standard_type": "IC50", "standard_units": "nM",
                                           "standard_relation__exact": "="}))
        url, params = API, {**base, "standard_type": "IC50", "standard_units": "nM",
                            "standard_relation__exact": "=", "limit": 1000, "only": FIELDS}
        while url:
            j = requests.get(url, params=params, timeout=120).json()
            rows += j["activities"]
            nxt = j["page_meta"]["next"]
            url, params = (f"https://www.ebi.ac.uk{nxt}", None) if nxt else (None, None)
    with open(out / "xo_ic50_activities.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["molecule_chembl_id", "target_chembl_id", "standard_value_nM", "assay_chembl_id", "assay_type",
                    "document_chembl_id", "document_year", "data_validity_flag", "potential_duplicate"])
        for r in rows:
            w.writerow([r["molecule_chembl_id"], r["target_chembl_id"], r["standard_value"], r["assay_chembl_id"],
                        r["assay_type"], r["document_chembl_id"], r["document_year"] or "",
                        "V" if r["data_validity_comment"] else "", "D" if r["potential_duplicate"] else ""])
    mols = {r["molecule_chembl_id"]: r["canonical_smiles"] for r in rows if r["canonical_smiles"]}
    with open(out / "xo_molecules.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["molecule_chembl_id", "canonical_smiles"]); w.writerows(mols.items())
    (out / "chembl_funnel_counts.json").write_text(json.dumps(funnel, indent=1))
    print(f"{len(rows)} activities, {len(mols)} molecules")


if __name__ == "__main__":
    main()
