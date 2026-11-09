"""run_pipeline.py: ChEMBL xanthine oxidase IC50 -> curated dataset -> flavonoid SAR -> QSAR model.

Inputs (pulled from the ChEMBL REST API; re-pull with src/fetch_chembl.py):
  data/raw/xo_ic50_activities.csv, data/raw/xo_molecules.csv, data/raw/chembl_funnel_counts.json
Outputs: data/curated/xo_curated.csv, outputs/tables/*, outputs/figures/*, outputs/key_numbers.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.ensemble import HistGradientBoostingRegressor, RandomForestRegressor
from sklearn.metrics import mean_squared_error, r2_score, roc_auc_score

sys.path.insert(0, str(Path(__file__).parent))
from chem import descriptors, flavonoid_class, morgan, sar_flags, scaffold, standardise  # noqa: E402
import pharmakit as pk  # noqa: E402

pk.mpl_theme()
RNG = 42
ACTIVE_PIC50 = 5.0  # IC50 <= 10 uM, the usual cut-off; allopurinol sits around 5-6 in these assays
for d in ["data/curated", "outputs/tables", "outputs/figures"]:
    Path(d).mkdir(parents=True, exist_ok=True)

acts = pd.read_csv("data/raw/xo_ic50_activities.csv", dtype=str)
mols = pd.read_csv("data/raw/xo_molecules.csv")
pre = json.loads(Path("data/raw/chembl_funnel_counts.json").read_text())

# ---------------------------------------------------------------- curation funnel
funnel = [("All activity records (human + bovine XO)", sum(v["all_activities"] for v in pre.values())),
          ("IC50 only", sum(v["ic50"] for v in pre.values())),
          ("Units nM", sum(v["ic50_nM"] for v in pre.values())),
          ("Exact values (relation '=')", sum(v["ic50_nM_eq"] for v in pre.values()))]
a = acts.copy()
a["value"] = pd.to_numeric(a["standard_value_nM"], errors="coerce")
a = a[a["value"] > 0]
funnel.append(("Positive numeric value", len(a)))
a = a[a["data_validity_flag"].isna()]
funnel.append(("No ChEMBL data-validity warning", len(a)))
a = a[a["potential_duplicate"].isna()]
funnel.append(("Not flagged as potential duplicate (citation)", len(a)))
a = a[a["assay_type"] == "B"]
funnel.append(("Binding/enzyme assays (type B)", len(a)))
a["pIC50"] = 9 - np.log10(a["value"])
a = a[(a["pIC50"] > 2) & (a["pIC50"] < 12)]
funnel.append(("Plausible range (pIC50 2–12)", len(a)))

# structure standardisation
smi = dict(zip(mols["molecule_chembl_id"], mols["canonical_smiles"]))
std = {}
for cid in a["molecule_chembl_id"].unique():
    m = standardise(smi.get(cid, "")) if cid in smi else None
    if m is not None and m.GetNumHeavyAtoms() >= 6:
        std[cid] = m
a = a[a["molecule_chembl_id"].isin(std)]
funnel.append(("Structure standardised (salts, charges)", len(a)))
from rdkit import Chem  # noqa: E402
a["inchikey"] = a["molecule_chembl_id"].map(lambda c: Chem.MolToInchiKey(std[c]))
a["std_smiles"] = a["molecule_chembl_id"].map(lambda c: Chem.MolToSmiles(std[c]))

# aggregate repeated measurements per compound (both species pooled; species mix recorded)
g = a.groupby("inchikey")
cur = g.agg(std_smiles=("std_smiles", "first"), chembl_ids=("molecule_chembl_id", lambda s: ";".join(sorted(set(s)))),
            pIC50=("pIC50", "median"), pIC50_sd=("pIC50", "std"), n_meas=("pIC50", "size"),
            n_docs=("document_chembl_id", "nunique"),
            species=("target_chembl_id", lambda s: "human" if set(s) == {"CHEMBL1929"} else ("bovine" if set(s) == {"CHEMBL3649"} else "both")),
            first_year=("document_year", "min")).reset_index()
cur["pIC50_sd"] = cur["pIC50_sd"].fillna(0)
cur["high_variance"] = cur["pIC50_sd"] > 1.0
funnel.append(("Unique compounds after aggregation (median pIC50)", len(cur)))
n_hv = int(cur["high_variance"].sum())
cur = cur[~cur["high_variance"]]
funnel.append(("Unique compounds, high-variance (SD>1) removed", len(cur)))

molmap = {r.inchikey: Chem.MolFromSmiles(r.std_smiles) for r in cur.itertuples()}
cur["flavonoid_class"] = cur["inchikey"].map(lambda k: flavonoid_class(molmap[k]))
cur["is_flavonoid"] = cur["flavonoid_class"] != "non-flavonoid"
sar = pd.DataFrame([sar_flags(molmap[k]) for k in cur["inchikey"]], index=cur.index)
cur = pd.concat([cur, sar], axis=1)
cur["scaffold"] = cur["inchikey"].map(lambda k: scaffold(molmap[k]))
cur["active"] = cur["pIC50"] >= ACTIVE_PIC50
cur.to_csv("data/curated/xo_curated.csv", index=False)
fun = pd.DataFrame(funnel, columns=["step", "n"])
fun.to_csv("outputs/tables/curation_funnel.csv", index=False)

# ---------------------------------------------------------------- flavonoid overview + SAR
cls = cur.groupby("flavonoid_class").agg(n=("pIC50", "size"), median_pIC50=("pIC50", "median"),
                                         pct_active=("active", "mean")).sort_values("n", ascending=False)
cls["pct_active"] = (100 * cls["pct_active"]).round(1)
cls.to_csv("outputs/tables/flavonoid_classes.csv")
fl = cur[cur["flavonoid_class"].isin(["flavone", "flavonol"])]
sar_rows = []
for feat, label in [("OH_5_7", "5,7-di-OH on A ring"), ("catechol_B", "3',4'-catechol on B ring"),
                    ("OH_3", "3-OH on C ring (flavonol)"), ("glycoside", "Sugar attached (glycoside)")]:
    w, wo = fl[fl[feat]], fl[~fl[feat]]
    if len(w) >= 3 and len(wo) >= 3:
        sar_rows.append(dict(feature=label, n_with=len(w), median_with=w["pIC50"].median(),
                             n_without=len(wo), median_without=wo["pIC50"].median(),
                             delta=w["pIC50"].median() - wo["pIC50"].median()))
sar_tbl = pd.DataFrame(sar_rows).round(2)
sar_tbl.to_csv("outputs/tables/flavone_sar.csv", index=False)
# planarity: C2=C3 double bond (flavone/flavonol) vs saturated (flavanone/flavan)
planar = cur[cur["flavonoid_class"].isin(["flavone", "flavonol"])]["pIC50"]
nonplanar = cur[cur["flavonoid_class"].isin(["flavanone", "flavan"])]["pIC50"]

# ---------------------------------------------------------------- features + scaffold split
X_fp = np.vstack([morgan(molmap[k]) for k in cur["inchikey"]])
D = pd.DataFrame([descriptors(molmap[k]) for k in cur["inchikey"]]).to_numpy()
X = np.hstack([X_fp, D])
y = cur["pIC50"].to_numpy()

# scaffold split: whole scaffold groups go to test until ~20% (largest groups stay in train)
rng = np.random.default_rng(RNG)
groups = cur["scaffold"].replace("", "acyclic").to_numpy()
uniq, counts = np.unique(groups, return_counts=True)
order = rng.permutation(len(uniq))
test_g, n_test = set(), 0
for i in order:
    if counts[i] > 0.05 * len(y):
        continue  # very common scaffolds stay in train (mirrors DeepChem's scaffold split)
    if n_test + counts[i] <= 0.2 * len(y):
        test_g.add(uniq[i]); n_test += counts[i]
test = np.array([g in test_g for g in groups])
rand_test = rng.random(len(y)) < 0.2


def fit_eval(train_mask, test_mask, model):
    model.fit(X[train_mask], y[train_mask])
    p = model.predict(X[test_mask])
    yt = y[test_mask]
    boots = []
    for _ in range(500):
        idx = rng.integers(0, len(yt), len(yt))
        boots.append(r2_score(yt[idx], p[idx]))
    auc = roc_auc_score(yt >= ACTIVE_PIC50, p) if 0 < (yt >= ACTIVE_PIC50).mean() < 1 else np.nan
    return dict(r2=r2_score(yt, p), r2_lcl=np.percentile(boots, 2.5), r2_ucl=np.percentile(boots, 97.5),
                rmse=mean_squared_error(yt, p) ** 0.5, auc=auc, n_train=int(train_mask.sum()), n_test=int(test_mask.sum())), p


res, preds = [], {}
for split_name, tmask in [("scaffold", test), ("random", rand_test)]:
    for mname, model in [("Random forest", RandomForestRegressor(n_estimators=500, min_samples_leaf=2, n_jobs=-1, random_state=RNG)),
                         ("Gradient boosting", HistGradientBoostingRegressor(max_iter=400, learning_rate=0.05, random_state=RNG))]:
        r, p = fit_eval(~tmask, tmask, model)
        res.append(dict(split=split_name, model=mname, **r))
        preds[(split_name, mname)] = p
# baseline: predict training mean
res.append(dict(split="scaffold", model="Mean baseline", r2=r2_score(y[test], np.full(test.sum(), y[~test].mean())),
                rmse=mean_squared_error(y[test], np.full(test.sum(), y[~test].mean())) ** 0.5,
                n_train=int((~test).sum()), n_test=int(test.sum())))
perf = pd.DataFrame(res).round(3)
perf.to_csv("outputs/tables/model_performance.csv", index=False)

# flavonoids in the scaffold test set
fl_test = cur["is_flavonoid"].to_numpy() & test
best = preds[("scaffold", "Random forest")]
fl_rmse = mean_squared_error(y[fl_test], best[fl_test[test]]) ** 0.5 if fl_test.sum() > 2 else np.nan

# ---------------------------------------------------------------- figures
fig, axes = plt.subplots(1, 2, figsize=(11, 4.3))
pca = PCA(n_components=2, random_state=RNG).fit_transform(X_fp.astype(float))
nf = ~cur["is_flavonoid"].to_numpy()
axes[0].scatter(pca[nf, 0], pca[nf, 1], s=9, c=pk.PALETTE["sage"], alpha=.6, label="Other chemotypes", edgecolor="none")
axes[0].scatter(pca[~nf, 0], pca[~nf, 1], s=16, c=pk.PALETTE["coral"], alpha=.9, label="Flavonoids", edgecolor="none")
axes[0].set_title("Chemical space of XO inhibitors (PCA, Morgan FP)")
axes[0].set_xlabel("PC1"); axes[0].set_ylabel("PC2"); axes[0].legend(loc="upper right")
yt = y[test]
axes[1].scatter(yt, best, s=14, c=np.where(fl_test[test], pk.PALETTE["coral"], pk.PALETTE["teal"]), alpha=.8, edgecolor="none")
lo, hi = min(yt.min(), best.min()) - .2, max(yt.max(), best.max()) + .2
axes[1].plot([lo, hi], [lo, hi], color=pk.PALETTE["grey"], lw=1, ls="--")
sc = perf.query("split=='scaffold' and model=='Random forest'").iloc[0]
axes[1].set_title(f"Scaffold-split test: R² {sc.r2:.2f}, RMSE {sc.rmse:.2f}")
axes[1].set_xlabel("Observed pIC50"); axes[1].set_ylabel("Predicted pIC50 (random forest)")
fig.tight_layout(); fig.savefig("outputs/figures/chemspace_and_model.png"); plt.close(fig)

fig, ax = plt.subplots(figsize=(6.6, 3.8))
order_cls = [c for c in ["flavonol", "flavone", "isoflavone", "flavanone", "flavan", "chalcone", "aurone", "non-flavonoid"] if c in cls.index]
data = [cur.loc[cur["flavonoid_class"] == c, "pIC50"] for c in order_cls]
bp = ax.boxplot(data, tick_labels=[f"{c}\n(n={len(d)})" for c, d in zip(order_cls, data)], patch_artist=True, widths=.6)
for patch, c in zip(bp["boxes"], order_cls):
    patch.set_facecolor(pk.PALETTE["panel"] if c == "non-flavonoid" else "#FBD9CD"); patch.set_edgecolor(pk.PALETTE["ink"])
ax.axhline(ACTIVE_PIC50, color=pk.PALETTE["coral"], ls="--", lw=1)
ax.set_ylabel("pIC50 (XO)"); ax.set_title("Potency by flavonoid class (dashed = 10 µM)")
ax.tick_params(axis="x", labelsize=8)
fig.tight_layout(); fig.savefig("outputs/figures/class_potency.png"); plt.close(fig)

key = dict(funnel=funnel, n_high_variance_removed=n_hv, n_compounds=len(cur), n_flavonoids=int(cur["is_flavonoid"].sum()),
           classes=cls.reset_index().to_dict("records"), sar=sar_tbl.to_dict("records"),
           planar_median=float(planar.median()), nonplanar_median=float(nonplanar.median()),
           n_planar=len(planar), n_nonplanar=len(nonplanar),
           performance=perf.to_dict("records"), flavonoid_test_n=int(fl_test.sum()), flavonoid_test_rmse=fl_rmse,
           species=cur["species"].value_counts().to_dict())
Path("outputs/key_numbers.json").write_text(json.dumps(key, indent=1, default=float))
print(fun.to_string(index=False)); print(cls); print(sar_tbl.to_string(index=False)); print(perf.to_string(index=False))
print("planar", planar.median(), len(planar), "non-planar", nonplanar.median(), len(nonplanar), "flav test", fl_test.sum(), fl_rmse)
