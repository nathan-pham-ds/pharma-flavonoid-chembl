# pharma-flavonoid-chembl

**Result:** I curated **4,443 ChEMBL activity records for xanthine oxidase (XO) down to 1,786 unique, standardised compounds**, and documented every filter. **170 of them are flavonoids.** Three structure–activity patterns known from bench work reappear in the data: 5,7-di-OH on the A ring (+0.5 log units), a planar C2=C3 bond (+0.3) and glycosylation (−0.6). A random-forest model reaches **R² 0.73 on a scaffold split**, compared with 0.83 on a random split, which shows how much random splits overstate performance.

> Path: R&D / lab data (data curation, cheminformatics) · Continues my BPharm thesis on flavonoids from *Elsholtzia ciliata* · Python + RDKit + scikit-learn

## Why xanthine oxidase

XO makes uric acid, and allopurinol and febuxostat are XO inhibitors used in gout. Flavonoids such as quercetin, luteolin and apigenin are classic natural XO inhibitors, so the target links the thesis chemistry (flavonoid isolation, HPLC) to a real drug target and a large public dataset.

## Data

- **ChEMBL REST API**, pulled 10 Nov 2026: human XO (CHEMBL1929) and bovine milk XO (CHEMBL3649). Most natural-product papers use bovine XO, so both species are pooled and the species mix is recorded per compound.
- Raw pulls are in `data/raw/`. The curated dataset is `data/curated/xo_curated.csv`, with one row per compound: InChIKey, SMILES, median pIC50, SD, n measurements, species, flavonoid class, SAR flags, scaffold.

## Curation funnel (`outputs/tables/curation_funnel.csv`)

| Step | Records |
|---|---|
| All activity records (human + bovine XO) | 4,443 |
| IC50 only | 2,552 |
| Units nM | 2,424 |
| Exact values (relation "=") | 2,174 |
| No ChEMBL data-validity warning | 2,094 |
| Not flagged as potential duplicate citation | 1,977 |
| Binding/enzyme assays, plausible range | 1,977 |
| Structure standardised (largest fragment, neutralised, cleaned) | 1,973 |
| Unique compounds (median pIC50 per InChIKey) | 1,797 |
| High-variance compounds removed (SD > 1 log unit) | **1,786** |

## Flavonoids: classes and SAR

Flavonoid classes are identified with bond-order-agnostic SMARTS (`src/chem.py`), which are unit-tested on quercetin, luteolin, apigenin, naringenin, genistein, catechin, isoliquiritigenin, allopurinol and febuxostat.

| Class | n | Median pIC50 | % active (≤ 10 µM) |
|---|---|---|---|
| Chalcone | 74 | 5.63 | 70.3 |
| Flavone | 37 | 5.05 | 54.1 |
| Flavonol | 28 | 5.19 | 64.3 |
| Aurone | 14 | 6.37 | 100 |
| Flavanone | 11 | 4.84 | 36.4 |
| Isoflavone / flavan | 4 / 2 | 5.03 / 5.59 | – |
| Non-flavonoid | 1,616 | 6.70 | 77.8 |

SAR within flavones and flavonols (median pIC50, with vs without the feature):

| Feature | With | Without | Δ | Literature |
|---|---|---|---|---|
| 5,7-di-OH on A ring | 5.48 (n=37) | 4.97 (n=28) | **+0.51** | agrees (Cos et al. 1998: C5/C7-OH important) |
| Planar C2=C3 (flavone/ol vs flavanone/flavan) | 5.12 (n=65) | 4.84 (n=13) | **+0.28** | agrees (planarity needed) |
| Sugar attached (glycoside) | 4.73 (n=12) | 5.28 (n=53) | **−0.56** | agrees (aglycones are more potent) |
| 3-OH on C ring | 5.43 (n=12) | 5.01 (n=53) | +0.42 | literature says neutral or slightly negative; confounded here by A-ring OH pattern |
| 3',4'-catechol on B ring | 5.12 | 5.12 | 0.00 | mixed in literature |

Flavonoids are moderate XO inhibitors (median about 10 µM), far weaker than optimised synthetic series such as the febuxostat-like compounds in the non-flavonoid group. That fits their role as dietary modulators rather than drugs.

## Model (`outputs/tables/model_performance.csv`)

Features: Morgan fingerprints (r = 2, 2,048 bits) plus 8 RDKit descriptors. Target: pIC50. The test set is about 20% of compounds.

| Split | Model | R² (95% bootstrap CI) | RMSE | ROC-AUC (active ≤ 10 µM) |
|---|---|---|---|---|
| **Scaffold** | Random forest | **0.73 (0.65–0.78)** | 0.73 | 0.91 |
| Scaffold | Gradient boosting | 0.72 (0.66–0.77) | 0.74 | 0.90 |
| Random | Random forest | 0.83 (0.77–0.87) | 0.59 | 0.94 |
| Scaffold | Mean baseline | −0.01 | 1.40 | – |

Under a random split, near-identical analogues end up in both train and test, so the model looks better than it is. The scaffold split is the honest estimate for new chemotypes. On the 27 flavonoids in the scaffold test set, RMSE is 0.52.

## How to run

```bash
pip install rdkit scikit-learn pandas numpy matplotlib requests pytest
python -m pytest -q              # SMARTS / standardisation tests
python src/fetch_chembl.py       # optional: re-pull from ChEMBL (raw pull included)
python src/run_pipeline.py       # curation -> SAR -> model -> figures
```

## Limitations

1. Human and bovine XO are pooled. Potency can differ by species and by assay (uric acid formation vs superoxide).
2. Only exact IC50 values are used. Censored values (">100 µM", the true inactives) are dropped, which biases the dataset towards actives (78% active).
3. IC50 values from different labs are not strictly comparable (substrate concentration, enzyme source). Median aggregation plus the SD > 1 filter only partly fixes this.
4. The SAR comparisons are univariate and have small n (12–65). They are descriptive and haven't been through a matched-molecular-pair analysis.
5. There is no applicability-domain filter yet. Predictions for scaffolds far from the training data should be flagged.

## Interview talking points

Why curation matters more than the model · InChIKey vs SMILES for deduplication · random vs scaffold split · censored data in bioactivity sets · how bench knowledge (hydroxylation pattern, glycosides from extraction work) helps judge whether a model's SAR is plausible.

## CV bullets

- EN: Curated 4,443 ChEMBL xanthine-oxidase bioactivity records into 1,786 unique compounds (RDKit standardisation, InChIKey deduplication, documented 10-step curation funnel), flagging 170 flavonoids with unit-tested SMARTS.
- EN: Trained a random-forest QSAR model on Morgan fingerprints with a scaffold split (R² 0.73, ROC-AUC 0.91) and recovered known flavonoid SAR (5,7-di-OH +0.5 log units; glycosides −0.6).
- VI: Tuyển chọn 4.443 dữ liệu hoạt tính trên xanthine oxidase từ ChEMBL thành 1.786 hợp chất duy nhất (chuẩn hóa cấu trúc bằng RDKit, loại trùng theo InChIKey, ghi nhận quy trình lọc 10 bước), xác định 170 flavonoid bằng SMARTS đã kiểm thử.
- VI: Xây dựng mô hình QSAR rừng ngẫu nhiên với fingerprint Morgan, chia dữ liệu theo khung cấu trúc (R² 0,73; ROC-AUC 0,91), tái hiện quan hệ cấu trúc–tác dụng đã biết của flavonoid (5,7-di-OH +0,5 log; glycosid −0,6).
