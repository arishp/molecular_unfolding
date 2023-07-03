import time
import sympy as sp
import copy
from rdkit import Chem
from sympy.matrices import ones, eye


n_angles = 8  # no. of discrete angles


def rotate_coordinates(rotation_matrix, old_coords):
    coord_vector = ones(4, 1)
    coord_vector[0, 0] = old_coords[0]
    coord_vector[1, 0] = old_coords[1]
    coord_vector[2, 0] = old_coords[2]
    coord_rot_vector = sp.expand(rotation_matrix * coord_vector)
    return [coord_rot_vector[0, 0], coord_rot_vector[1, 0], coord_rot_vector[2, 0]]


def rotation_matrix_new_coords(first_coords, second_coords, bond_no, soln_theta):
    x_dash, y_dash, z_dash = first_coords[0], first_coords[1], first_coords[2]
    x_ddash, y_ddash, z_ddash = second_coords[0], second_coords[1], second_coords[2]
    dx = x_ddash - x_dash
    dy = y_ddash - y_dash
    dz = z_ddash - z_dash
    l_sq = dx ** 2 + dy ** 2 + dz ** 2
    l = sp.sqrt(l_sq)
    c_theta = sp.cos(soln_theta)
    s_theta = sp.sin(soln_theta)
    # index = n_angles * bond_no
    rotation_matrix = eye(4)
    rotation_matrix[0, 0] = ((dx ** 2 + (dy ** 2 + dz ** 2) * c_theta) / l_sq).evalf()
    rotation_matrix[0, 1] = ((dx * dy * (1 - c_theta) - dz * l * s_theta) / l_sq).evalf()
    rotation_matrix[0, 2] = ((dx * dz * (1 - c_theta) + dy * l * s_theta) / l_sq).evalf()
    rotation_matrix[0, 3] = (((x_dash * (dy ** 2 + dz ** 2) - dx * (y_dash * dy + z_dash * dz)) * (1 - c_theta) + (
            y_dash * dz - z_dash * dy) * l * s_theta) / l_sq).evalf()
    rotation_matrix[1, 0] = ((dx * dy * (1 - c_theta) + dz * l * s_theta) / l_sq).evalf()
    rotation_matrix[1, 1] = ((dy ** 2 + (dx ** 2 + dz ** 2) * c_theta) / l_sq).evalf()
    rotation_matrix[1, 2] = ((dy * dz * (1 - c_theta) - dx * l * s_theta) / l_sq).evalf()
    rotation_matrix[1, 3] = (((y_dash * (dx ** 2 + dz ** 2) - dy * (x_dash * dx + z_dash * dz)) * (1 - c_theta) + (
            z_dash * dx - x_dash * dz) * l * s_theta) / l_sq).evalf()
    rotation_matrix[2, 0] = ((dx * dz * (1 - c_theta) - dy * l * s_theta) / l_sq).evalf()
    rotation_matrix[2, 1] = ((dy * dz * (1 - c_theta) + dx * l * s_theta) / l_sq).evalf()
    rotation_matrix[2, 2] = ((dz ** 2 + (dx ** 2 + dy ** 2) * c_theta) / l_sq).evalf()
    rotation_matrix[2, 3] = (((z_dash * (dx ** 2 + dy ** 2) - dz * (x_dash * dx + y_dash * dy)) * (1 - c_theta) + (
            x_dash * dy - y_dash * dx) * l * s_theta) / l_sq).evalf()
    return rotation_matrix


def generate_thetas():
    angle_incr = 2 * sp.pi / n_angles
    thetas = [i*angle_incr for i in range(n_angles)]
    return thetas


def distance_squared(first_coords, second_coords):
    dis_sq = (second_coords[0] - first_coords[0])**2 + (second_coords[1] - first_coords[1])**2 + (second_coords[2] - first_coords[2])**2
    return dis_sq


start = time.time()
thetas = generate_thetas()
mol_folder = './Molecules/'
mol_fname = '25A12Rot_37.mol2'
mol_filepath = mol_folder + mol_fname
mol = Chem.rdmolfiles.MolFromMol2File(mol_filepath)
n_atoms_mol = mol.GetNumAtoms()  # no. of atoms in the molecule (excluding hydrogen atoms)
# print('No. of atoms: ', n_atoms_mol)
conformers = mol.GetConformers()
conf = conformers[0]
coords_mol = {}  # Coordinates of the atoms
for i in range(n_atoms_mol):
    coords_mol[i] = list(conf.GetAtomPosition(i))
input_folder = './Input_Rot_2/'
input_fname = mol_fname[:-5] + '_input.txt'
input_filepath = input_folder + input_fname
input_lines = open(input_filepath, 'r').readlines()
torsional_bonds_mol = eval(input_lines[0])
coords_rotation_mol = eval(input_lines[2])
torsional_config = {}
for bond in torsional_bonds_mol.keys():
    torsional_config[bond] = 0
old_volume = 0
for i in range(n_atoms_mol - 1):
    for j in range(i + 1, n_atoms_mol):
        old_volume += distance_squared(coords_mol[i], coords_mol[j])
# print('old volume: ', old_volume)
best_volume = old_volume
best_soln = copy.deepcopy(torsional_config)

old_best_volume = best_volume
iterations = 0
while True:
    iterations = iterations + 1
    for bond_no in torsional_bonds_mol:
        temp_soln = copy.deepcopy(best_soln)
        # print('\nbond no: ', bond_no)
        for angle_index in range(8):
            # print('angle: ', thetas[angle_index])
            temp_soln[bond_no] = thetas[angle_index]
            # print('temp solution: ', temp_soln)
            new_coords = copy.deepcopy(coords_mol)
            for atom, bonds in coords_rotation_mol.items():
                rot_mat = eye(4, 4)
                for bond in bonds:
                    temp_rot_mat = rotation_matrix_new_coords(coords_mol[torsional_bonds_mol[bond][0]],
                                                              coords_mol[torsional_bonds_mol[bond][1]],
                                                              bond, temp_soln[bond])
                    rot_mat = rot_mat * temp_rot_mat
                new_coords[atom] = rotate_coordinates(rot_mat, coords_mol[atom])
            new_volume = 0
            for u in range(n_atoms_mol - 1):
                for v in range(u+1, n_atoms_mol):
                    new_volume += distance_squared(new_coords[u], new_coords[v])
            # print('new volume: ', new_volume)
            if new_volume > best_volume:
                # print('BEST VOLUME!!!')
                best_volume = new_volume
                best_soln = copy.deepcopy(temp_soln)
    if best_volume - old_best_volume < 10:
        break
    else:
        old_best_volume = best_volume
torsional_config = best_soln
end = time.time()
print('no. of iterations: ', iterations)
print('best solution: ', torsional_config)
print('volume change: ', best_volume - old_volume)
print("The time of execution: ", (end-start), "s")

