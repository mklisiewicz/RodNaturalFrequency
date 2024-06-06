import sympy as sp
from sympy import latex

n, l, E, J, m = sp.symbols('n l E J m')
k = (2 * n - 1) * sp.pi / (2 * l)
omega = sp.sqrt(l * E * J / m)/sp.pi * k**2
# omega= sp.expand(omega)
omega_latex = latex(omega)

with open('omega.tex', 'w') as f:
    f.write(omega_latex)