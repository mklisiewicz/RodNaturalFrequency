import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp
from scipy.optimize import fsolve

# Given parameters
E_b_init = 213.5 * 10**9  # Young's modulus (Pa)
d_b_init = 0.009  # Diameter (m)
l_b_init = 2.89  # Length (m)
m_b_init = 1.55  # Mass (kg)
J_b_init = np.pi * d_b_init**4 / 64  # Moment of inertia (m^4)
A_cross = np.pi * (d_b_init / 2)**2  # Cross-sectional area (m^2)
rho_b_init = m_b_init / (A_cross * l_b_init)  # Density (kg/m^3)
N = 1000  # Number of spatial points
x = np.linspace(0, l_b_init, N)

# Natural frequencies for a cantilever beam (beta values)
beta_guess = [1.875, 4.694, 7.855, 10.996]  # Initial guesses for beta L values
beta = []

def char_eq(beta):
    return np.cos(beta) * np.cosh(beta) + 1

for guess in beta_guess:
    beta.append(fsolve(char_eq, guess)[0])

# Natural frequencies (omega_n)
omega_n = [((beta[i]/l_b_init) ** 2) * np.sqrt(E_b_init * J_b_init * l_b_init / m_b_init) for i in range(len(beta))]

# Displacement function for given omega and time t
def displacement(x, omega, t, mode):
    return A * np.sin(mode * np.pi * x / l_b_init) * np.cos(omega * t)

# Plot settings
A = 0.1  # Amplitude of the forced oscillation (m)

# Set up the figure and axis
fig, axes = plt.subplots(len(omega_n), 1, figsize=(6, 15))

# Plotting the motion at natural frequencies
for i, omega in enumerate(omega_n):
    T = 2 * np.pi / omega  # Period of the current natural frequency (s)
    t = np.linspace(0, T, 200)  # Time array for the animation
    ax = axes[i]
    
    for ti in t:
        u = displacement(x, omega, ti, i+1)
        ax.plot(u, x, 'b-')
        ax.set_xlabel('Displacement (m)')
        ax.set_ylabel('Length (m)')
        ax.set_ylim(0, l_b_init)
        ax.set_xlim(-A, A)
        ax.set_title(f'Mode {i+1}, Natural Frequency = {omega/np.pi:.2f}$\pi$ rad/s')
        plt.pause(0.05)
        ax.cla()  # Clear the axis for the next frame

plt.tight_layout()
plt.show()
