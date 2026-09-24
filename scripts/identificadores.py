"""SMILES/InChI/InChIKey canônicos de um candidato (data/candidatos/*.json), a partir do grafo exato.

Moléculas com modos zero não têm estrutura de Kekulé completa: usamos um emparelhamento máximo
(ligações duplas) e marcamos os carbonos que sobram como radicais.
Uso: python3 scripts/identificadores.py data/candidatos/C51H25_eta3.json
"""
import json, sys
import networkx as nx
from rdkit import Chem
from rdkit.Chem import Descriptors


def molecule(struct):
    g = nx.Graph()
    g.add_nodes_from(a["id"] for a in struct["carbonos"])
    g.add_edges_from(tuple(b) for b in struct["ligacoes"])
    matching = nx.max_weight_matching(g, maxcardinality=True)
    doubles = {frozenset(e) for e in matching}
    m = Chem.RWMol()
    for a in struct["carbonos"]:
        m.AddAtom(Chem.Atom(6))
    for u, v in g.edges():
        m.AddBond(u, v, Chem.BondType.DOUBLE if frozenset((u, v)) in doubles else Chem.BondType.SINGLE)
    matched = {x for e in matching for x in e}
    radicals = 0
    for a in m.GetAtoms():
        deg = a.GetDegree()
        a.SetNoImplicit(True)
        a.SetNumExplicitHs(1 if deg == 2 else 0)
        if a.GetIdx() not in matched:
            a.SetNumRadicalElectrons(1)
            radicals += 1
    mol = m.GetMol()
    Chem.SanitizeMol(mol)
    return mol, len(matching), radicals


if __name__ == "__main__":
    rec = json.load(open(sys.argv[1]))
    mol, doubles, radicals = molecule(rec["estrutura"])
    print("fórmula      ", Chem.rdMolDescriptors.CalcMolFormula(mol))
    print("duplas/radicais (emparelhamento máximo):", doubles, radicals)
    print("SMILES       ", Chem.MolToSmiles(mol))
    print("InChI        ", Chem.MolToInchi(mol))
    print("InChIKey     ", Chem.MolToInchiKey(mol))
    print("massa (g/mol)", round(Descriptors.MolWt(mol), 3))
