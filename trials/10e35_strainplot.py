import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Constants
G = 6.67430e-11  # Gravitational constant
c = 299792458
m = 1.989e30  # solar mass
m1 = 1.4 * m
m2 = 1.4 * m
M = m1 + m2
Mc = (m1 * m2) ** (0.6) / M ** 0.2
µ = m1 * m2 / M
AU = 150 * 10 ** 9
a = AU / 10 ** 4.5  # in AU
e = 0  # eccentricity
T = 2 * np.pi * np.sqrt(a ** 3 / (G * M))
tm = 50  # merger time
fusion_distance = a / 27  # distance below which stars will merge
R = 500 * 3.086 * 10 ** (19)  # distance of observation 500Mpc
α = 4.1
β = 4
γ = 1


B = 10e11
magnetar_radius = 20
mu_naught = 4e-7 * np.pi
magnetic_moment = (B * magnetar_radius ** 3) / (2 * mu_naught)
reduced_mass = (m1 ** 2) / (m1 + m2)

# Dipole interaction function
def dipole_interaction(dist):
    return ((6 * mu_naught) / (4 * np.pi * dist ** 4)) * (magnetic_moment ** 2)


def v(a):
    return 2 * np.pi * a / T


def relative_velocity(v1, v2):
    return np.sqrt((v1[0] - v2[0]) ** 2 + (v1[1] - v2[1]) ** 2)


def distance(x1, y1, x2, y2):
    return np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)


def friction_coefficient(r, t):
    return β / ((t + γ * tm) ** (3 / 5) * np.exp(α * r / a))  # Avoid division by zero


def equation(t, r):
    x1, y1, vx1, vy1, x2, y2, vx2, vy2 = r

    R1 = np.sqrt(x1 ** 2 + y1 ** 2)
    R2 = np.sqrt(x2 ** 2 + y2 ** 2)

    # Gravitational forces
    dvx1dt = -G * m2 * x1 / (R1 ** 3 + 1e-10)  # Avoid division by zero
    dvy1dt = -G * m2 * y1 / (R1 ** 3 + 1e-10)

    dvx2dt = -G * m1 * x2 / (R2 ** 3 + 1e-10)
    dvy2dt = -G * m1 * y2 / (R2 ** 3 + 1e-10)

    # Friction forces
    v_rel = np.array([vx1 - vx2, vy1 - vy2])
    dist = distance(x1, y1, x2, y2)
    k = friction_coefficient(dist, tm)

    dvx1dt -= (k * v_rel[0] + (dipole_interaction(dist) * 1.26) / (10 ** 4))
    dvy1dt -= (k * v_rel[1] + (dipole_interaction(dist) * 1.26) / (10 ** 4))
    dvx2dt += (k * v_rel[0] + (dipole_interaction(dist) * 1.26) / (10 ** 4))
    dvy2dt += (k * v_rel[1] + (dipole_interaction(dist) * 1.26) / (10 ** 4))


    return [vx1, vy1, dvx1dt, dvy1dt, vx2, vy2, dvx2dt, dvy2dt]


def h(t, distances):
    f_GW = [np.sqrt(G * M / (distances[i] + 1e-10) ** 3) / (np.pi) for i, dist in enumerate(distances)]
    h0 = [
        4 * G * Mc / (c ** 2 * R) * (G / c ** 3 * np.pi * f_GW[i] * Mc) ** (2 / 3) for i, dist in enumerate(distances)
    ]
    dfdt = [
        96 / 5 * c ** 3 / G * f_GW[i] / Mc * (G / c ** 3 * np.pi * f_GW[i] * Mc) ** (8 / 3)
        for i, dist in enumerate(distances)
    ]
    ϕ = [2 * np.pi * (f_GW[i] * t[i] + 0.5 * dfdt[i] * t[i] ** 2) for i, dist in enumerate(distances)]

    h = [0 if dist < fusion_distance else h0[i] * np.cos(ϕ[i]) for i, dist in enumerate(distances)]
    return h


initial_conditions = [a * (1 - e) * m1 / M, 0, 0, v(a), -a * (1 - e) * m2 / M, 0, 0, -v(a)]

# Integrate using solve_ivp with adaptive step sizing
sol = solve_ivp(equation, [0, tm / 2], initial_conditions, method='RK45', t_eval=np.arange(0.0, tm / 2, 0.01))
solution = sol.y.T

# Use sol.t (the time points actually returned by solve_ivp) instead of t_values
distances = [distance(solution[i, 0], solution[i, 1], solution[i, 4], solution[i, 5]) for i in range(len(sol.t))]

# Create the figure and axes
fig = plt.figure(figsize=(10, 12))
fig.set_facecolor("black")
gs = fig.add_gridspec(2, hspace=0.20)

# Plot the strain (gravitational wave signal)
ax2 = fig.add_subplot(gs[0])
ax2.set_xlim(0, tm / 2)

ax2.set_ylim(min(h(sol.t, distances)), max(h(sol.t, distances)))
ax2.set_facecolor("black")

# Plotting the gravitational wave strain over time
ax2.plot(sol.t, h(sol.t, distances), color="white")
ax2.set_title("Strain over time", color="white")
ax2.set_xlabel("Time", color="white")
ax2.set_ylabel("Strain", color="white")

# Customize ticks and spines for the signal plot
ax2.tick_params(axis="x", colors="white")
ax2.tick_params(axis="y", colors="white")
ax2.spines["bottom"].set_color("white")
ax2.spines["left"].set_color("white")

# Show plot
plt.show()
