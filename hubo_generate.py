import sympy as sp
from sympy.matrices import zeros, ones, eye
from sympy import Point3D

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


def generate_hard_constraint():
    A_const = sp.Symbol('A_const')
    for i in range(n_bonds):
        x.append(sp.symbols(f'x{str(i)}(0:{n_angles})'))
    hard_constraint = 0
    for i in range(n_bonds):
        summation = 0
        for j in range(n_angles):
            summation += x[i][j]
        hard_constraint += (summation - 1) ** 2
    hard_constraint *= A_const
    return hard_constraint


def distance_squared(first_coords, second_coords):
    p1 = Point3D(first_coords[0], first_coords[1], first_coords[2])
    p2 = Point3D(second_coords[0], second_coords[1], second_coords[2])
    return p1.distance(p2) ** 2


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
    rotation_matrix = zeros(4, 4)
    rotation_matrix[0, 0] = (dx ** 2 + (dy ** 2 + dz ** 2) * c_theta) / l_sq
    rotation_matrix[0, 1] = (dx * dy * (1 - c_theta) - dz * l * s_theta) / l_sq
    rotation_matrix[0, 2] = (dx * dz * (1 - c_theta) - dy * l * s_theta) / l_sq
    rotation_matrix[0, 3] = ((x_dash * (dy ** 2 + dz ** 2) - dx * (y_dash * dy + z_dash * dz)) * (1 - c_theta) + (
            y_dash * dz - z_dash * dy) * l * s_theta) / l_sq
    rotation_matrix[1, 0] = (dx * dy * (1 - c_theta) + dz * l * s_theta) / l_sq
    rotation_matrix[1, 1] = (dy * dy + (dx * dx + dz * dz) * c_theta) / l_sq
    rotation_matrix[1, 2] = (dy * dz * (1 - c_theta) - dx * l * s_theta) / l_sq
    rotation_matrix[1, 3] = ((y_dash * (dx * dx + dz * dz) - dy * (x_dash * dx + z_dash * dz)) * (1 - c_theta) + (
            z_dash * dx - x_dash * dz) * l * s_theta) / l_sq
    rotation_matrix[2, 0] = (dx * dz * (1 - c_theta) - dy * l * s_theta) / l_sq
    rotation_matrix[2, 1] = (dy * dz * (1 - c_theta) + dx * l * s_theta) / l_sq
    rotation_matrix[2, 2] = (dz * dz + (dx * dx + dy * dy) * c_theta) / l_sq
    rotation_matrix[2, 3] = ((z_dash * (dx * dx + dy * dy) - dz * (x_dash * dx + y_dash * dy)) * (1 - c_theta) + (
            x_dash * dy - y_dash * dx) * l * s_theta) / l_sq
    rotation_matrix[3, 0] = 0.0
    rotation_matrix[3, 1] = 0.0
    rotation_matrix[3, 2] = 0.0
    rotation_matrix[3, 3] = 1.0
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


def main():
    hubo_expr = generate_hard_constraint()
    print("\nHARD CONSTRAINT")
    print("---- ----------")
    sp.pprint(hubo_expr)

    print("\n\nCO-ORDINATES BEFORE ROTATION")
    print("------------ ------ --------")
    sp.pprint(coords_dict)

    rotate_all_coordinates()

    print("\n\nCO-ORDINATES AFTER ROTATION")
    print("------------ ----- --------")
    sp.pprint(coords_dict)

    print("\n\nOPTIMIZATION CONSTRAINT: ")
    print("------------ ----------")
    distance_hubo = generate_distance_hubo()
    sp.pprint(distance_hubo)

    print('\nFINAL HUBO: ')
    print('----- ----')
    hubo_expr += distance_hubo
    sp.pprint(hubo_expr)
    print("\nHUBO EXPANDED")
    print("---- --------")
    print(hubo_expr.expand())


if __name__ == "__main__":
    main()

