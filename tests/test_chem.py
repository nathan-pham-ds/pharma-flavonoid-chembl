"""Known structures must classify correctly.  Run: python -m pytest -q"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from chem import flavonoid_class, sar_flags, standardise  # noqa: E402

KNOWN = {
    "quercetin": ("O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12", "flavonol"),
    "quercetin_kekule": ("O=C1C(O)=C(OC2=CC(O)=CC(O)=C12)C1=CC(O)=C(O)C=C1", "flavonol"),
    "luteolin": ("O=c1cc(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12", "flavone"),
    "apigenin": ("O=c1cc(-c2ccc(O)cc2)oc2cc(O)cc(O)c12", "flavone"),
    "naringenin": ("O=C1CC(c2ccc(O)cc2)Oc2cc(O)cc(O)c21", "flavanone"),
    "genistein": ("O=c1c(-c2ccc(O)cc2)coc2cc(O)cc(O)c12", "isoflavone"),
    "catechin": ("Oc1ccc(cc1O)[C@H]1Oc2cc(O)cc(O)c2C[C@@H]1O", "flavan"),
    "isoliquiritigenin": ("O=C(/C=C/c1ccc(O)cc1)c1ccc(O)cc1O", "chalcone"),
    "allopurinol": ("O=c1[nH]cnc2[nH]ncc12", "non-flavonoid"),
    "febuxostat": ("Cc1nc(-c2ccc(OCC(C)C)c(C#N)c2)sc1C(=O)O", "non-flavonoid"),
}


def test_classes():
    for name, (smi, cls) in KNOWN.items():
        assert flavonoid_class(standardise(smi)) == cls, name


def test_sar_flags():
    q = sar_flags(standardise(KNOWN["quercetin"][0]))
    assert q["OH_5_7"] and q["catechol_B"] and q["OH_3"]
    a = sar_flags(standardise(KNOWN["apigenin"][0]))
    assert a["OH_5_7"] and not a["catechol_B"] and not a["OH_3"]


def test_salt_stripping():
    m = standardise("[Na+].[O-]C(=O)c1sc(-c2ccc(OCC(C)C)c(C#N)c2)nc1C")
    assert m.GetNumAtoms() == 22 and "Na" not in [a.GetSymbol() for a in m.GetAtoms()]
