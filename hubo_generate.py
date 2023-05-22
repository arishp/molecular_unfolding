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
coords_dict = {
    1: [1.0, -0.5, 0.0],
    2: [0.0, 0.0, 0.0],
    3: [0.0, 1.0, 0.0],
    4: [1.0, 1.5, 0.0]
}

# dictionary of bond numbers
torsional_bonds = {
    0: (2, 3)
}

# dictionary of bonds affecting coordinates
coords_rotation_dict = {
    1: [],
    2: [],
    3: [],
    4: [0]
}

# pair of coordinates to find distances
distance_pairs = [
    (1, 4),
    (2, 4)]

n_bonds = len(torsional_bonds)  # no. of bonds
n_angles = 4  # no. of discrete angles

x = []  # global variables for hubo variables

#################
# END OF INPUTS #
#################


def generate_hard_constraint(include_a=False, a_value=1.0):
    for i in range(n_bonds):
        x.append(sp.symbols(f'x{str(i)}(0:{n_angles})'))
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
    return hard_constraint


def distance_squared(first_coords, second_coords):
    p1 = Point3D(first_coords[0], first_coords[1], first_coords[2])
    p2 = Point3D(second_coords[0], second_coords[1], second_coords[2])
    return p2.distance(p1) ** 2


def generate_distance_hubo():
    distance_sq = 0
    for pair in distance_pairs:
        distance_sq += distance_squared(coords_dict[pair[0]], coords_dict[pair[1]])
    return distance_sq


def generate_rotation_matrix(first_coords, second_coords, bond_no):
    x_dash, y_dash, z_dash = first_coords[0], first_coords[1], first_coords[2]
    x_ddash, y_ddash, z_ddash = second_coords[0], second_coords[1], second_coords[2]
    dx = x_ddash - x_dash
    dy = y_ddash - y_dash
    dz = z_ddash - z_dash
    l_sq = dx ** 2 + dy ** 2 + dz ** 2
    l = sp.sqrt(l_sq)
    angle_incr = 2 * sp.pi / n_angles
    angle = 0.0
    thetas = [angle]
    for i in range(1, n_angles):
        angle += angle_incr
        thetas.append(angle)
    c_theta = 0.0
    s_theta = 0.0
    for i in range(n_angles):
        c_theta += (sp.cos(thetas[i]) * x[bond_no][i])
        s_theta += (sp.sin(thetas[i]) * x[bond_no][i])
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


def rotate_coordinates(rotation_matrix, old_coords):
    coord_vector = ones(4, 1)
    coord_vector[0, 0] = old_coords[0]
    coord_vector[1, 0] = old_coords[1]
    coord_vector[2, 0] = old_coords[2]
    coord_rot_vector = rotation_matrix * coord_vector
    return [coord_rot_vector[0, 0], coord_rot_vector[1, 0], coord_rot_vector[2, 0]]


def rotate_all_coordinates():
    for i in coords_dict.keys():
        if len(coords_rotation_dict[i]):
            rot_mat = eye(4, 4)
            for bond_no in coords_rotation_dict[i]:
                temp_rot_mat = generate_rotation_matrix(coords_dict[torsional_bonds[bond_no][0]],
                                                        coords_dict[torsional_bonds[bond_no][1]], bond_no)
                rot_mat = rot_mat * temp_rot_mat
            coords_dict[i] = rotate_coordinates(rot_mat, coords_dict[i])


def extract_hubo_dict(hubo_expr):
    hubo_expr_str = str(hubo_expr.expand())
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


def main():
    hubo_expr = generate_hard_constraint(include_a=False, a_value=20.0)
    print("\nHARD CONSTRAINT")
    print("---- ----------")
    sp.pprint(hubo_expr)

    rotate_all_coordinates()

    print("\n\nOPTIMIZATION CONSTRAINT: ")
    print("------------ ----------")
    distance_hubo = generate_distance_hubo()
    sp.pprint(distance_hubo)

    print('\nFINAL HUBO: ')
    print('----- ----')
    hubo_expr -= distance_hubo
    sp.pprint(hubo_expr)

    print("\nHUBO EXPANDED")
    print("---- --------")
    print(hubo_expr.expand())

    print('\nHUBO DICTIONARY')
    print('---- ----------')
    hubo_dict = extract_hubo_dict(hubo_expr)
    print(hubo_dict)

    bqm = dimod.make_quadratic(hubo_dict, 12.0, dimod.BINARY)

    sampler = neal.SimulatedAnnealingSampler()
    sample_size=10
    sampleset = sampler.sample(bqm, num_reads=sample_size)
    print("\nSA RESULTS:\n",sampleset)

    sampler = EmbeddingComposite(DWaveSampler())
    sampleset = sampler.sample(bqm, num_reads=1000)
    print("\nQA RESULTS:\n",sampleset)


if __name__ == "__main__":
    main()
