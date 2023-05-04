from sympy import *
import numpy as np


d = 4
angle = 0.0
angle_incr = 2*np.pi/d
d_angles = [angle]
for i in range(1, d):
    angle += angle_incr
    d_angles.append(angle)
#print(d_angles)



x11, x12, x13, x14 = symbols('x11, x12, x13, x14')
init_printing(use_unicode=False, wrap_line=False)
e = (x11 - x13)*(x11 - x13)
print(e.as_poly())

hubo = {}

A_const= 7.0
### Hard constraint
hubo[('x11', 'x11')] = 1 * A_const
hubo[('x11', 'x12')] = 2.0 * A_const
hubo[('x11', 'x13')] = 2.0 * A_const
hubo[('x11', 'x14')] = 2.0 * A_const
hubo[('x11', )] = -2.0 * A_const
hubo[('x12', 'x12')] = 1 * A_const
hubo[('x12', 'x13')] = 2.0 * A_const
hubo[('x12', 'x14')] = 2.0 * A_const
hubo[('x12', )] = -2.0 * A_const
hubo[('x13', 'x13')] = 1 * A_const
hubo[('x13', 'x14')] = 2.0 * A_const
hubo[('x13', )] = -2.0 * A_const
hubo[('x14', 'x14')] = 1 * A_const
hubo[('x14', )] = -2.0 * A_const
hubo[()] = 1.0 * A_const

### optimizing constraint
#C1, C4 - X
hubo[('x11', 'x11')] -= 1
hubo[('x11', 'x13')] -= 2.0
hubo[('x11', )] -= -2.0
hubo[('x13', 'x13')] -= 1
hubo[('x13', )] -= 2.0
hubo[()] -= 1.0
#C1, C4 - Y
hubo[('x11', 'x11')] -= 2.25
hubo[('x11', 'x12')] -= 4.5
hubo[('x11', 'x13')] -= 4.5
hubo[('x11', 'x14')] -= 4.5
hubo[('x11', )] -= -1.5
hubo[('x12', 'x12')] -= 2.25
hubo[('x12', 'x13')] -= 4.5
hubo[('x12', 'x14')] -= 4.5
hubo[('x12', )] -= 1.5
hubo[('x13', 'x13')] -= 2.25
hubo[('x13', 'x14')] -= 4.5
hubo[('x13', )] -= 1.5
hubo[('x14', 'x14')] -= 2.25
hubo[('x14', )] -= 1.5
hubo[()] -= 0.25
#C1, C4 - Z
hubo[('x12', 'x12')] -= 1
hubo[('x12', 'x14')] -= -2.0
hubo[('x14', 'x14')] -= 1


#C2, C4 - X
hubo[('x11', 'x11')] -= 1
hubo[('x11', 'x13')] -= -2.0
hubo[('x13', 'x13')] -= 1
#C2, C4 - Y
hubo[('x11', 'x11')] -= 2.25
hubo[('x11', 'x12')] -= 4.5
hubo[('x11', 'x13')] -= 4.5
hubo[('x11', 'x14')] -= 4.5
hubo[('x12', 'x12')] -= 2.25
hubo[('x12', 'x13')] -= 4.5
hubo[('x12', 'x14')] -= 4.5
hubo[('x13', 'x13')] -= 2.25
hubo[('x13', 'x14')] -= 4.5
hubo[('x14', 'x14')] -= 2.25
#C2, C4 - Z
hubo[('x12', 'x12')] -= 1
hubo[('x12', 'x14')] -= -2.0
hubo[('x14', 'x14')] -= 1

print(hubo)
