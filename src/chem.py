"""chem.py: structure standardisation, flavonoid classification, features."""
from __future__ import annotations

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

# Bond-order-agnostic SMARTS (=,: and -,:) so Kekulé and aromatic inputs both match.
_B = "[c,C]"
FLAVONOID_SMARTS = {
    # 2-phenylchromen-4-one (C2=C3 double bond): flavones + flavonols
    "flavone": "[#8]=[#6]1[#6]=,:[#6](-c2ccccc2)[#8][#6]2:[#6]:[#6]:[#6]:[#6]:[#6]12",
    # 3-hydroxy flavone
    "flavonol": "[#8]=[#6]1[#6](-[OX2H1,OX2])=,:[#6](-c2ccccc2)[#8][#6]2:[#6]:[#6]:[#6]:[#6]:[#6]12",
    # 2,3-dihydro (no C2=C3): flavanones / flavanonols
    "flavanone": "[#8]=[#6]1[#6;X4][#6;X4](-c2ccccc2)[#8][#6]2:[#6]:[#6]:[#6]:[#6]:[#6]12",
    # 3-phenylchromen-4-one
    "isoflavone": "[#8]=[#6]1[#6](-c2ccccc2)=,:[#6][#8][#6]2:[#6]:[#6]:[#6]:[#6]:[#6]12",
    # flavan-3-ols (catechins) and other flavans
    "flavan": "[#6;X4]1(-c2ccccc2)[#8][#6]2:[#6]:[#6]:[#6]:[#6]:[#6]2[#6;X4][#6;X4]1",
    # open-chain chalcones (biosynthetic precursors, usually grouped with flavonoids)
    "chalcone": "[#8]=[#6](-[#6]=[#6]-c1ccccc1)-c1ccccc1",
    # aurones
    "aurone": "[#8]=[#6]1[#6](=[#6]-c2ccccc2)[#8][#6]2:[#6]:[#6]:[#6]:[#6]:[#6]12",
}
_PATTERNS = {k: Chem.MolFromSmarts(v) for k, v in FLAVONOID_SMARTS.items()}
# SAR features on the flavone/flavonol core (A ring 5,7-OH; B ring catechol; C ring 3-OH)
SAR_SMARTS = {
    "OH_5_7": "[OX2H1]c1cc([OX2H1])c2c(c1)[#8][#6](-c)=,:[#6][#6]2=[#8]",
    "catechol_B": "[#6]1(=,:[#6][#6](=[#8])c2ccccc2[#8]1)-c1ccc([OX2H1])c([OX2H1])c1",
    "OH_3": "[OX2H1][#6]1=,:[#6](-c)[#8]c2ccccc2[#6]1=[#8]",
    "glycoside": "[#6;R1]1([#8][#6;R1])[#8][#6;R1][#6;R1][#6;R1][#6;R1]1",
}
_SAR = {k: Chem.MolFromSmarts(v) for k, v in SAR_SMARTS.items()}

_uncharger = rdMolStandardize.Uncharger()
_chooser = rdMolStandardize.LargestFragmentChooser(preferOrganic=True)
_te = rdMolStandardize.TautomerEnumerator()


def standardise(smiles: str, tautomer: bool = False) -> Chem.Mol | None:
    """Strip salts/solvents (largest organic fragment), neutralise, cleanup, optional canonical tautomer."""
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    m = rdMolStandardize.Cleanup(m)
    m = _chooser.choose(m)
    m = _uncharger.uncharge(m)
    if tautomer:
        m = _te.Canonicalize(m)
    Chem.SanitizeMol(m)
    return m


def flavonoid_class(m: Chem.Mol) -> str:
    order = ["flavonol", "flavone", "isoflavone", "flavanone", "flavan", "aurone", "chalcone"]
    for k in order:
        if m.HasSubstructMatch(_PATTERNS[k]):
            return k
    return "non-flavonoid"


def sar_flags(m: Chem.Mol) -> dict:
    return {k: m.HasSubstructMatch(p) for k, p in _SAR.items()}


def scaffold(m: Chem.Mol) -> str:
    try:
        return MurckoScaffold.MurckoScaffoldSmiles(mol=m)
    except Exception:
        return ""


def morgan(m: Chem.Mol, radius: int = 2, n_bits: int = 2048) -> np.ndarray:
    gen = AllChem.GetMorganGenerator(radius=radius, fpSize=n_bits)
    return gen.GetFingerprintAsNumPy(m).astype(np.uint8)


def descriptors(m: Chem.Mol) -> dict:
    return dict(MW=Descriptors.MolWt(m), cLogP=Descriptors.MolLogP(m), TPSA=rdMolDescriptors.CalcTPSA(m),
                HBD=rdMolDescriptors.CalcNumHBD(m), HBA=rdMolDescriptors.CalcNumHBA(m),
                RotB=rdMolDescriptors.CalcNumRotatableBonds(m), AromRings=rdMolDescriptors.CalcNumAromaticRings(m),
                FracSP3=rdMolDescriptors.CalcFractionCSP3(m))
