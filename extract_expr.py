import sympy as sp

f = open("full_hubo_expr.txt", "r")
hubo_str = f.read()
hubo_str = hubo_str.replace('+ ', '+')
hubo_str = hubo_str.replace('- ', '-')
hubo_list = hubo_str.split()
print(hubo_list,'\n\n')
for monom in hubo_list:
    if monom != '+' and monom != '-':
        print(monom.split(['*','**']))
        