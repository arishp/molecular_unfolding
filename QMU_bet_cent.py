import networkx as nx
from rdkit import Chem

mol_folder = './Molecules/'
mol_fname = '50A12Rot_3.mol2'
mol_filepath = mol_folder + mol_fname
mol = Chem.rdmolfiles.MolFromMol2File(mol_filepath)
n_atoms_mol = mol.GetNumAtoms()
atoms_list = [x for x in range(n_atoms_mol)]

G = nx.Graph()
G.add_nodes_from(atoms_list)

n_bonds_mol = mol.GetNumBonds()
for i in range(n_bonds_mol):
    G.add_edge(mol.GetBondWithIdx(i).GetBeginAtomIdx(), mol.GetBondWithIdx(i).GetEndAtomIdx())
bet_cent = nx.algorithms.betweenness_centrality(G)
for item in bet_cent.items():
    print('atom:', item[0], '\tb_c value: ', item[1])
print(G)
