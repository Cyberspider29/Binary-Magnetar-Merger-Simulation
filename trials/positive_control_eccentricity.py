import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Constants
G = 6.67430e-11  # Gravitational constant
c = 299792458
m = 1.989e30  # solar mass
m1 = 1.4 * m
m2 = 1.4 * m
M = m1 + m2
Mc = (m1 * m2) ** 0.6 / M ** 0.2
µ = m1 * m2 / M
AU = 150 * 10 ** 9
a = AU / 10 ** 4.5  # in AU
e = 0  # initial eccentricity
T = 2 * np.pi * np.sqrt(a ** 3 / (G * M))
tm = 50  # merger time
fusion_distance = a / 27  # distance below which stars will merge
R = 500 * 3.086 * 10 ** 19  # distance of observation 500 Mpc
α = 4.1
β = 4
γ = 1

# EM Interaction Constants
B = 10e11
magnetar_radius = 20
mu_naught = 4e-7 * np.pi
magnetic_moment = (B * magnetar_radius ** 3) / (2 * mu_naught)
reduced_mass = (m1 ** 2) / (m1 + m2)

# Dipole interaction function
def dipole_interaction(dist):
    return ((6 * mu_naught) / (4 * np.pi * dist ** 4)) * (magnetic_moment ** 2)

# Velocity function
def v(a):
    return 2 * np.pi * a / T

# Distance calculation function
def distance(x1, y1, x2, y2):
    return np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)

# Friction coefficient function
def friction_coefficient(r, t):
    return β / ((t + γ * tm) ** (3/5) * np.exp(α * r / a))

# Equations of motion
def equation(t, r):
    x1, y1, vx1, vy1, x2, y2, vx2, vy2 = r

    R1 = np.sqrt(x1 ** 2 + y1 ** 2)
    R2 = np.sqrt(x2 ** 2 + y2 ** 2)

    # Gravitational forces
    dvx1dt = -G * m2 * x1 / (R1 ** 3 + 1e-10)
    dvy1dt = -G * m2 * y1 / (R1 ** 3 + 1e-10)
    dvx2dt = -G * m1 * x2 / (R2 ** 3 + 1e-10)
    dvy2dt = -G * m1 * y2 / (R2 ** 3 + 1e-10)

    # Friction forces and magnetic interaction
    v_rel = np.array([vx1 - vx2, vy1 - vy2])
    dist = distance(x1, y1, x2, y2)
    k = friction_coefficient(dist, tm)

    dvx1dt -= (k * v_rel[0] + dipole_interaction(dist) / 10e10)
    dvy1dt -= (k * v_rel[1] + dipole_interaction(dist) / 10e10)
    dvx2dt += (k * v_rel[0] + dipole_interaction(dist) / 10e10)
    dvy2dt += (k * v_rel[1] + dipole_interaction(dist) / 10e10)

    return [vx1, vy1, dvx1dt, dvy1dt, vx2, vy2, dvx2dt, dvy2dt]

# Initial conditions
initial_conditions = [a * (1 - e) * m1 / M, 0, 0, v(a), -a * (1 - e) * m2 / M, 0, 0, -v(a)]
t_span = (0.0, tm / 2)
t_values = np.linspace(0.0, tm / 2, 1000)

# Solve the differential equations using solve_ivp
solution = solve_ivp(equation, t_span, initial_conditions, t_eval=t_values, method='RK45', rtol=1e-8, atol=1e-10)
x1, y1, x2, y2 = solution.y[0], solution.y[1], solution.y[4], solution.y[5]

# Calculate separation distances for each frame
separation_distances = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)

# Animation function
def update(frame):
    dist = separation_distances[frame]
    
    r_a = np.max(separation_distances[:frame+1])  # current maximum separation distance 
    r_p = np.min(separation_distances[:frame+1])  # current minimum separation distance

    # Calculate eccentricity
    if (r_a + r_p) != 0:  # Avoid division by zero
        eccentricity = (r_a - r_p) / (r_a + r_p)
    else:
        eccentricity = 0  # Default to 0 if both distances are zero
    print("Eccentricity:", eccentricity)
    # Check if merger should happen
    if dist < fusion_distance:
        star1.set_data(0, 0)
        star2.set_data(0, 0)
        star_merged.set_data(0, 0)
        star_merged.set_markersize(33)  # Adjust size
        star_merged.set_color("gray")  # Set color to gray
    else:
        star1.set_data(x1[frame], y1[frame])
        star2.set_data(x2[frame], y2[frame])
        star_merged.set_data([], [])  # Clear merged star data

    return star1, star2, star_merged

# Create the figure and axes
fig = plt.figure(figsize=(10, 12))
fig.set_facecolor("black")
gs = fig.add_gridspec(2, hspace=0.20)

# Upper plot for animation
ax1 = fig.add_subplot(gs[0])
ax1.set_xlim(-1.1 * a / 2, 1.1 * a / 2)
ax1.set_ylim(-1.1 * a / 2, 1.1 * a / 2)
ax1.set_facecolor("black")

star1, = ax1.plot([], [], "o", markersize=20, color="lightblue")
star2, = ax1.plot([], [], "o", markersize=20, color="white")
star_merged, = ax1.plot([], [], "o", markersize=33, color="gray")

# Adjust the size of the plots
ax1.set_position([0.1, 0.35, 0.8, 0.63])  # Increase the size of the upper plot
ani = FuncAnimation(fig, update, frames=len(t_values), interval=1, blit=True)

# Animation
plt.show()
