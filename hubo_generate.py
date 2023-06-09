import copy
import sympy as sp
from sympy.matrices import ones, eye
from sympy import Point3D
import dimod
import neal
from dwave.system import DWaveSampler, EmbeddingComposite

#########################
# INPUTS FOR AUTOMATION #
#########################

# dictionary of coordinates involved
coords_original = {
    0: [1.0, -0.5, 0.0],
    1: [0.0, 0.0, 0.0],
    2: [0.0, 1.0, 0.0],
    3: [1.0, 1.5, 0.0],
    4: [2.0, 1.0, 0.0]
}

# dictionary of bond numbers
torsional_bonds = {
    0: (1, 2),
    1: (2, 3)
}

# dictionary of bonds affecting coordinates
coords_rotation_dict = {
    0: [],
    1: [],
    2: [],
    3: [0],
    4: [0, 1]
}

# pair of coordinates to find distances
distance_pairs = [
    (0, 3),
    (0, 4),
    (1, 3),
    (1, 4),
    (2, 4)]

coords_dict = copy.deepcopy(coords_original)
final_coords = copy.deepcopy(coords_original)

n_bonds = len(torsional_bonds)  # no. of bonds
n_angles = 4  # no. of discrete angles

x = []  # global variables for hubo variables


#################
# END OF INPUTS #
#################


def generate_hard_constraint(include_a=False, a_value=1.0):
    for i in range(n_bonds):
        x.append(sp.symbols(f'x_{str(i)}_(0:{n_angles})'))
    hard_constraint = 0
    for i in range(n_bonds):
        summation = 0
        for j in range(n_angles):
            summation += x[i][j]
        hard_constraint += (summation - 1) ** 2
    if include_a:
        a_const = sp.Symbol('A_const')
        hard_constraint *= a_const
    else:
        hard_constraint *= a_value
    return hard_constraint.expand()


def generate_thetas():
    angle_incr = 2 * sp.pi / n_angles
    angle = 0.0
    thetas = [angle]
    for i in range(1, n_angles):
        angle += angle_incr
        thetas.append(angle)
    return thetas


def distance_squared(first_coords, second_coords):
    dis_sq = sp.expand((second_coords[0] - first_coords[0])**2 + (second_coords[1] - first_coords[1])**2 +(second_coords[2] - first_coords[2])**2)
    return dis_sq


def generate_distance_hubo():
    distance_sq = 0
    for pair in distance_pairs:
        distance_sq += distance_squared(coords_dict[pair[0]], coords_dict[pair[1]])
    return distance_sq.expand()


def generate_rotation_matrix(first_coords, second_coords, bond_no):
    x_dash, y_dash, z_dash = first_coords[0], first_coords[1], first_coords[2]
    x_ddash, y_ddash, z_ddash = second_coords[0], second_coords[1], second_coords[2]
    dx = x_ddash - x_dash
    dy = y_ddash - y_dash
    dz = z_ddash - z_dash
    l_sq = sp.expand(dx ** 2 + dy ** 2 + dz ** 2)
    l = sp.expand(sp.sqrt(l_sq))
    thetas = generate_thetas()
    c_theta = 0.0
    s_theta = 0.0
    for i in range(n_angles):
        c_theta += sp.expand(sp.cos(thetas[i]) * x[bond_no][i])
        s_theta += sp.expand(sp.sin(thetas[i]) * x[bond_no][i])
    rotation_matrix = eye(4)
    rotation_matrix[0, 0] = sp.expand((dx ** 2 + (dy ** 2 + dz ** 2) * c_theta) / l_sq)
    rotation_matrix[0, 1] = sp.expand((dx * dy * (1 - c_theta) - dz * l * s_theta) / l_sq)
    rotation_matrix[0, 2] = sp.expand((dx * dz * (1 - c_theta) + dy * l * s_theta) / l_sq)
    rotation_matrix[0, 3] = sp.expand(((x_dash * (dy ** 2 + dz ** 2) - dx * (y_dash * dy + z_dash * dz)) * (1 - c_theta) + (
            y_dash * dz - z_dash * dy) * l * s_theta) / l_sq)
    rotation_matrix[1, 0] = sp.expand((dx * dy * (1 - c_theta) + dz * l * s_theta) / l_sq)
    rotation_matrix[1, 1] = sp.expand((dy ** 2 + (dx ** 2 + dz ** 2) * c_theta) / l_sq)
    rotation_matrix[1, 2] = sp.expand((dy * dz * (1 - c_theta) - dx * l * s_theta) / l_sq)
    rotation_matrix[1, 3] = sp.expand(((y_dash * (dx ** 2 + dz ** 2) - dy * (x_dash * dx + z_dash * dz)) * (1 - c_theta) + (
            z_dash * dx - x_dash * dz) * l * s_theta) / l_sq)
    rotation_matrix[2, 0] = sp.expand((dx * dz * (1 - c_theta) - dy * l * s_theta) / l_sq)
    rotation_matrix[2, 1] = sp.expand((dy * dz * (1 - c_theta) + dx * l * s_theta) / l_sq)
    rotation_matrix[2, 2] = sp.expand((dz ** 2 + (dx ** 2 + dy ** 2) * c_theta) / l_sq)
    rotation_matrix[2, 3] = sp.expand(((z_dash * (dx ** 2 + dy ** 2) - dz * (x_dash * dx + y_dash * dy)) * (1 - c_theta) + (
            x_dash * dy - y_dash * dx) * l * s_theta) / l_sq)
    return sp.expand(rotation_matrix)


def rotate_coordinates(rotation_matrix, old_coords):
    coord_vector = ones(4, 1)
    coord_vector[0, 0] = old_coords[0]
    coord_vector[1, 0] = old_coords[1]
    coord_vector[2, 0] = old_coords[2]
    coord_rot_vector = sp.expand(rotation_matrix * coord_vector)
    return [coord_rot_vector[0, 0], coord_rot_vector[1, 0], coord_rot_vector[2, 0]]


def rotate_all_coordinates():
    for i in coords_original.keys():
        if len(coords_rotation_dict[i]) > 0:
            rot_mat = eye(4, 4)
            for bond_no in coords_rotation_dict[i]:
                temp_rot_mat = generate_rotation_matrix(coords_original[torsional_bonds[bond_no][0]],
                                                        coords_original[torsional_bonds[bond_no][1]], bond_no)
                rot_mat = rot_mat * temp_rot_mat
            coords_dict[i] = rotate_coordinates(rot_mat, coords_original[i])


def extract_hubo_dict(hubo_expr_str):
    hubo_expr_str = hubo_expr_str.replace('+ ', '+')
    hubo_expr_str = hubo_expr_str.replace('- ', '-')
    hubo_expr_list = hubo_expr_str.split()
    hubo_dict = {}
    for mono in hubo_expr_list:
        mono_list = mono.split("*")
        if len(mono_list) > 1:
            dict_index = (mono_list[1],)
            if len(mono_list) > 2:
                if mono_list[2] != '':
                    temp_index_list = []
                    for item in mono_list[2:]:
                        temp_index_list.append(item)
                    if len(temp_index_list) > 0:
                        dict_index += tuple(temp_index_list)
                else:
                    repeat_tuple = dict_index
                    for i in range(int(mono_list[3])-1):
                        dict_index += repeat_tuple

            hubo_dict[dict_index] = float(mono_list[0])
        else:
            hubo_dict[()] = float(mono_list[0])
    return hubo_dict


def generate_final_rotation_matrix(first_coords, second_coords, bond_no, torsional_config):
    x_dash, y_dash, z_dash = first_coords[0], first_coords[1], first_coords[2]
    x_ddash, y_ddash, z_ddash = second_coords[0], second_coords[1], second_coords[2]
    dx = x_ddash - x_dash
    dy = y_ddash - y_dash
    dz = z_ddash - z_dash
    l_sq = dx ** 2 + dy ** 2 + dz ** 2
    l = sp.sqrt(l_sq)
    c_theta = sp.cos(torsional_config[bond_no])
    s_theta = sp.sin(torsional_config[bond_no])
    rotation_matrix = eye(4)
    rotation_matrix[0, 0] = (dx ** 2 + (dy ** 2 + dz ** 2) * c_theta) / l_sq
    rotation_matrix[0, 1] = (dx * dy * (1 - c_theta) - dz * l * s_theta) / l_sq
    rotation_matrix[0, 2] = (dx * dz * (1 - c_theta) + dy * l * s_theta) / l_sq
    rotation_matrix[0, 3] = ((x_dash * (dy ** 2 + dz ** 2) - dx * (y_dash * dy + z_dash * dz)) * (1 - c_theta) + (
            y_dash * dz - z_dash * dy) * l * s_theta) / l_sq
    rotation_matrix[1, 0] = (dx * dy * (1 - c_theta) + dz * l * s_theta) / l_sq
    rotation_matrix[1, 1] = (dy ** 2 + (dx ** 2 + dz ** 2) * c_theta) / l_sq
    rotation_matrix[1, 2] = (dy * dz * (1 - c_theta) - dx * l * s_theta) / l_sq
    rotation_matrix[1, 3] = ((y_dash * (dx ** 2 + dz ** 2) - dy * (x_dash * dx + z_dash * dz)) * (1 - c_theta) + (
            z_dash * dx - x_dash * dz) * l * s_theta) / l_sq
    rotation_matrix[2, 0] = (dx * dz * (1 - c_theta) - dy * l * s_theta) / l_sq
    rotation_matrix[2, 1] = (dy * dz * (1 - c_theta) + dx * l * s_theta) / l_sq
    rotation_matrix[2, 2] = (dz ** 2 + (dx ** 2 + dy ** 2) * c_theta) / l_sq
    rotation_matrix[2, 3] = ((z_dash * (dx ** 2 + dy ** 2) - dz * (x_dash * dx + y_dash * dy)) * (1 - c_theta) + (
            x_dash * dy - y_dash * dx) * l * s_theta) / l_sq
    return rotation_matrix


def print_new_coords(solution):
    torsional_config = {}
    thetas = generate_thetas()
    for key, value in solution.items():
        if value == 1:
            key_list = key.split('_')
            torsional_config[int(key_list[1])] = thetas[int(key_list[2])]
    for i in coords_dict.keys():
        if len(coords_rotation_dict[i]) > 0:
            rot_mat = eye(4, 4)
            for bond_no in coords_rotation_dict[i]:
                temp_rot_mat = generate_final_rotation_matrix(coords_original[torsional_bonds[bond_no][0]],
                                                              coords_original[torsional_bonds[bond_no][1]], bond_no,
                                                              torsional_config)
                rot_mat = rot_mat * temp_rot_mat
            final_coords[i] = rotate_coordinates(rot_mat, final_coords[i])
    print(final_coords)


def main():
    hubo_expr = generate_hard_constraint(include_a=False, a_value=20.0)
    # print("\nHARD CONSTRAINT")
    # print("---- ----------")
    # sp.pprint(hubo_expr)
    #
    # rotate_all_coordinates()
    #
    # print("\n\nOPTIMIZATION CONSTRAINT: ")
    # print("------------ ----------")
    # distance_hubo = generate_distance_hubo()
    # sp.pprint(distance_hubo)
    #
    # print('\nFINAL HUBO: ')
    # print('----- ----')
    # hubo_expr -= distance_hubo
    # sp.pprint(hubo_expr)
    #
    # hubo_expr_str = str(hubo_expr.expand())
    #
    # print("\n\nFile hubo_expr.txt is getting ready...")
    # f = open("hubo_expr.txt", "w")
    # f.write(hubo_expr_str)
    # f.close()
    # print('\n\nFile hubo_expr.txt created!')

    # print('\nHUBO DICTIONARY')
    # print('---- ----------')
    # hubo_expr_str = open("hubo_expr.txt", "r").read()
    # hubo_dict = extract_hubo_dict(hubo_expr_str)
    # print(hubo_dict)

    # read hubo_dict from a file
    # bqm = dimod.make_quadratic(hubo_dict, 12.0, dimod.BINARY)
    # sampler = neal.SimulatedAnnealingSampler()
    # sample_size = 10
    # sampleset = sampler.sample(bqm, num_reads=sample_size)
    # print("\nSA RESULTS:\n-- -------\n", sampleset)
    # sa_solution = sampleset.first.sample
    # print("\nBEST SA RESULT:\n---- -- ------\n", sa_solution)
    
    # sa_solution = {'x_0_2': 1, 'x_1_0': 1}
    # print_new_coords(sa_solution)


if __name__ == "__main__":
    main()
