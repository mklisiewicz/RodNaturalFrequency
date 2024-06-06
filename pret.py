# import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib.widgets import Slider

# # Experimental values
# exp_values_b = np.array([130, 366, 586])
# exp_values_gf = np.array([240, 700, 1660])
# n_exp = np.array([1, 2, 3])
# n_lin = np.linspace(1, 3, 100)
# # Default parameters
# E_b_init = 130*10**9
# d_b_init = 0.009
# l_b_init = 2.89
# m_b_init = 1.55
# J_b_init = np.pi * d_b_init**4 / 64

# E_gf_init = 69.5*10**9
# d_gf_init = 0.01
# l_gf_init = 2.6
# m_gf_init = 0.390
# J_gf_init = np.pi * d_gf_init**4 / 64

# # Function to calculate omega values
# def get_omega(n, l, E, J, m):
#     def get_k(n, l):
#         return (2 * n - 1) * np.pi / (2 * l)
#     k = get_k(n, l)
#     return np.sqrt(l * E * J / m) * k**2

# # Plot initial data
# fig, ax = plt.subplots()
# plt.subplots_adjust(left=0.1, bottom=0.25)
# line1, = plt.plot(n_exp, exp_values_b/2, 'o', label='Eksperymentalne wartości $\\omega(n)$ dla pręta z brązu')
# line2, = plt.plot(n_exp, exp_values_gf/2, 'o', label='Eksperymentalne wartości $\\omega(n)$ dla pręta z włókna szklanego')
# line3, = plt.plot([], [], 'o', label='Teoretyczne wartości $\\omega(n)$ dla pręta z brązu')
# line4, = plt.plot([], [], 'o', label='Teoretyczne wartości $\\omega(n)$ dla pręta z włókna szklanego')
# line5, = plt.plot([], [], '--', label='Linia trendu teoretycznych wartości $\\omega(n)$ dla pręta z brązu')
# line6, = plt.plot([], [], '--', label='Linia trendu teoretycznych wartości $\\omega(n)$ dla pręta z włókna szklanego')
# plt.xlabel('$n$')
# plt.ylabel('$\\omega$')
# plt.title('$\\omega(n)$')
# plt.grid(True)
# plt.legend()

# # Create sliders
# axcolor = 'lightgoldenrodyellow'
# ax_E_b = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
# ax_E_gf = plt.axes([0.25, 0.3, 0.65, 0.03], facecolor=axcolor)
# ax_J_gf = plt.axes([0.25, 0.4, 0.65, 0.03], facecolor=axcolor)
# ax_J_b = plt.axes([0.25, 0.5, 0.65, 0.03], facecolor=axcolor)
# ax_r_b = plt.axes([0.25, 0.6, 0.65, 0.03], facecolor=axcolor)
# ax_r_gf = plt.axes([0.25, 0.7, 0.65, 0.03], facecolor=axcolor)
# s_r_b = Slider(ax_r_b, 'r_b', 0, 10, valinit=1)
# s_r_gf = Slider(ax_r_gf, 'r_gf', 0, 10, valinit=1)
# s_E_b = Slider(ax_E_b, 'E_b', 0, 300*10**9, valinit=E_b_init)
# s_E_gf = Slider(ax_E_gf, 'E_gf', 0, 100*10**9, valinit=E_gf_init)
# s_J_gf = Slider(ax_J_gf, 'J_gf', 0, 100*10**9, valinit=J_gf_init)
# s_J_b = Slider(ax_J_b, 'J_b', 0, 100*10**9, valinit=J_b_init)

# # Update function
# def update(val):
#     E_b = s_E_b.val
#     E_gf = s_E_gf.val
#     J_gf = s_J_gf.val
#     J_b = s_J_b.val
#     r_b=s_r_b.val
#     r_gf=s_r_gf.val
    
#     line3.set_data(n_exp, r_b*get_omega(n_exp, l_b_init, E_b, J_b, m_b_init))
#     line4.set_data(n_exp, r_gf*get_omega(n_exp, l_gf_init, E_gf, J_gf, m_gf_init))
#     line5.set_data(n_lin, r_b*get_omega(n_lin, l_b_init, E_b, J_b, m_b_init))
#     line6.set_data(n_lin, r_gf*get_omega(n_lin, l_gf_init, E_gf, J_gf, m_gf_init))
#     fig.canvas.draw_idle()

# # Link sliders to update function
# s_E_b.on_changed(update)
# s_J_b.on_changed(update)
# s_J_gf.on_changed(update)
# s_E_gf.on_changed(update)
# s_r_b.on_changed(update)
# s_r_gf.on_changed(update)

# plt.show()

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage[T1]{fontenc}'

exp_values_b = np.array([65, 183, 293])
exp_values_gf = np.array([120, 350, 830])
exp_values_b = exp_values_b*2*np.pi/60
exp_values_gf = exp_values_gf*2*np.pi/60


n = np.array([1, 2, 3])
n_lin = np.linspace(1, 3, 100)

E_b = 50*10**9
d_b = 0.009
l_b = 2.89
m_b = 1.55
J_b = np.pi * d_b**4 / 64

E_gf = 69.5*10**9
d_gf = 0.01
l_gf = 2.6
m_gf = 0.390
J_gf = np.pi * d_gf**4 / 64
print(J_gf*E_gf)
print(J_b*E_b)


def get_omega(n, l, E, J, m):
    
    return np.pi**2*np.sqrt(l * E * J / m)/(4*l**2)*(2 * n - 1)**2

exp_interp_b = interp1d(n, exp_values_b, kind='quadratic')
exp_interp_gf = interp1d(n, exp_values_gf, kind='quadratic')
print(get_omega(n, l_b, E_b, J_b, m_b)/np.pi)
print(get_omega(n, l_gf, E_gf, J_gf, m_gf))
plt.figure(figsize=(10, 6))
plt.scatter(n, get_omega(n, l_b, E_b, J_b, m_b)/np.pi, label='Teoretyczne wartości $\\omega(n)$ dla pręta z brązu')
plt.plot(n_lin, get_omega(n_lin, l_b, E_b, J_b, m_b)/np.pi, '--', label='Linia trendu teoretycznych wartości $\\omega(n)$ dla pręta z brązu')
plt.scatter(n, exp_values_b/np.pi, label='Eksperymentalne wartości $\\omega(n)$ dla pręta z brązu')
plt.plot(n_lin, exp_interp_b(n_lin)/np.pi, '--', label='Interpolacja trendu eksperymentalnych wartości $\\omega(n)$ dla pręta z brązu')
# plt.scatter(n, get_omega(n, l_gf, E_gf, J_gf, m_gf), label='Teoretyczne wartości $\\omega(n)$ dla pręta z włókna szklanego')
# plt.plot(n_lin, get_omega(n_lin, l_gf, E_gf, J_gf, m_gf), '--', label='Linia trendu teoretycznych wartości $\\omega(n)$ dla pręta z włókna szklanego')
# plt.scatter(n, exp_values_gf, label='Eksperymentalne wartości $\\omega(n)$ dla pręta z włókna szklanego')
# plt.plot(n_lin, exp_interp_gf(n_lin), '--', label='Interpolacja trendu eksperymentalnych wartości $\\omega(n)$ dla pręta z włókna szklanego')
plt.title('$\\omega(n)$')
plt.xlabel('$n$')
plt.ylabel('$\\omega [\\pi\\frac{rad}{s}]$')
plt.grid(True)
plt.legend()
plt.show()