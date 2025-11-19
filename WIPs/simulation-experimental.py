# fixed_merge_simulation.py
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# ================================
# USER SETTINGS
# ================================
USE_DRAG = True
USE_MAG  = False

# Set the INITIAL SEPARATION between the magnetars (meters)
INITIAL_SEPARATION = 0.4e8   # try 4e8, 4e7, 4e6 to experiment

# Choose fusion threshold as fraction of initial separation (or set absolute)
FUSION_FRACTION = 2.0e4 / INITIAL_SEPARATION  # fusion_distance = INITIAL_SEPARATION * FUSION_FRACTION

# ================================
# CONSTANTS
# ================================
G = 6.67430e-11
m_sun = 1.989e30
m1 = 1.4 * m_sun
m2 = 1.4 * m_sun
M  = m1 + m2

tm = 200.0                 # used in drag model if enabled
mu0 = 4e-7 * np.pi

# Magnetar magnetic moment (used only if USE_MAG=True)
B = 1e12
R_mag = 20.0
magnetic_moment = (B * R_mag**3) / (2 * mu0)

# drag parameters
alpha = 4.1
beta  = 1e6
gamma = 1.0

# ================================
# HELPER FUNCTIONS
# ================================
def friction_coefficient(r, t):
    # keep dependence on current separation r and current time t
    # note use of INITIAL_SEPARATION to normalize exponent scale
    return beta / ((t + gamma * tm)**(3/5) * np.exp(alpha * r / INITIAL_SEPARATION))


# ================================
# EQUATIONS OF MOTION
# ================================
def eom(t, S):
    """
    State vector S = [x1, y1, vx1, vy1, x2, y2, vx2, vy2]
    Returns derivatives: [vx1, vy1, ax1, ay1, vx2, vy2, ax2, ay2]
    """
    x1, y1, vx1, vy1, x2, y2, vx2, vy2 = S
    
    r1 = np.array([x1, y1], dtype=float)
    r2 = np.array([x2, y2], dtype=float)
    r_vec = r1 - r2                       # vector from star2 --> star1
    dist = np.linalg.norm(r_vec) + 1e-15  # avoid div-by-zero
    r_hat = r_vec / dist

    v1 = np.array([vx1, vy1], dtype=float)
    v2 = np.array([vx2, vy2], dtype=float)
    v_rel = v1 - v2

    # gravitational acceleration (two-body Newtonian)
    a1 = -G * m2 * r_vec / dist**3
    a2 =  G * m1 * r_vec / dist**3   # this equals -G*m1*(r2-r1)/dist**3

    # optional drag (dissipative)
    if USE_DRAG:
        k = friction_coefficient(dist, t)
        a1 += (-k * v_rel) / m1
        a2 += (+k * v_rel) / m2

    # optional magnetic dipole repulsion (axial approx)
    if USE_MAG:
        mu1 = mu2 = magnetic_moment
        # axial approximation scalar (repulsive for parallel aligned dipoles)
        Fmag = (3.0 * mu0 * mu1 * mu2) / (2.0 * np.pi * (dist**4 + 1e-15))
        a1 += (Fmag / m1) * r_hat
        a2 += -(Fmag / m2) * r_hat

    return [vx1, vy1, a1[0], a1[1], vx2, vy2, a2[0], a2[1]]


# ================================
# EVENT: stop integration when separation < fusion_distance
# ================================
def fusion_event(t, S):
    x1, y1, vx1, vy1, x2, y2, vx2, vy2 = S
    dist = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
    return dist - fusion_distance

fusion_event.terminal = True   # stop integration
fusion_event.direction = -1    # trigger when dist is decreasing through threshold


# ================================
# INITIAL CONDITIONS
# ================================
# Put the stars at +/- separation/2 on x-axis
r1x = +INITIAL_SEPARATION / 2.0
r2x = -INITIAL_SEPARATION / 2.0

# compute circular orbital speed for separation = INITIAL_SEPARATION (approx)
v_rel = np.sqrt(G * (m1 + m2) / INITIAL_SEPARATION)

# split velocities about center of mass
v1 =  v_rel * (m2 / M)     # +y
v2 = -v_rel * (m1 / M)     # -y

state0 = [r1x, 0.0, 0.0, v1,
          r2x, 0.0, 0.0, v2]

fusion_distance = R_mag * 5000000

print("Starting ODE integration with:")
print("  INITIAL_SEPARATION =", INITIAL_SEPARATION)
print("  fusion_distance    =", fusion_distance)
print("  USE_DRAG =", USE_DRAG, "USE_MAG =", USE_MAG)


# ================================
# INTEGRATE (with event stopping if fusion occurs)
# ================================
# choose a long enough timespan but integration will stop early if fusion_event triggers
t_final = 5e5   # seconds (you can change this)
t_eval = np.linspace(0.0, t_final, 5000)

sol = solve_ivp(eom, [0.0, t_final], state0, t_eval=t_eval,
                rtol=1e-9, atol=1e-12, events=fusion_event)

print("solve_ivp finished. nsteps:", sol.t.size)
if sol.t_events and sol.t_events[0].size > 0:
    t_fuse = sol.t_events[0][0]
    print("Fusion event detected at t =", t_fuse, "s")
else:
    t_fuse = None
    print("No fusion event detected during integration.")

x1 = sol.y[0]
y1 = sol.y[1]
x2 = sol.y[4]
y2 = sol.y[5]

sep = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
min_sep = np.min(sep)
min_idx = np.argmin(sep)
print("min separation during integration: ", min_sep, " at t =", sol.t[min_idx])


# ================================
# ANIMATION (no-blit to keep drawing reliable)
# ================================
fig, ax = plt.subplots(figsize=(6,6))
ax.set_facecolor('black')
lim = max(1.1 * INITIAL_SEPARATION/2.0, np.max(sep)*1.1)
ax.set_xlim(-lim, lim)
ax.set_ylim(-lim, lim)

star1, = ax.plot([],[], 'o', color='cyan',  markersize=10)
star2, = ax.plot([],[], 'o', color='white', markersize=10)
star_merged, = ax.plot([],[], 'o', color='gray', markersize=18)  # merged marker (initially hidden)
star_merged.set_visible(False)

def update(i):
    # called every frame; i indexes into t_eval/sol.t
    if sep[i] < fusion_distance:
        # fusion case: hide individuals and show merged marker
        star1.set_data([], [])
        star2.set_data([], [])
        star_merged.set_data(0.0, 0.0)   # put merged object at COM (choose 0,0)
        star_merged.set_visible(True)
    else:
        star_merged.set_visible(False)
        star1.set_data(x1[i], y1[i])
        star2.set_data(x2[i], y2[i])
    return star1, star2, star_merged

ani = FuncAnimation(fig, update, frames=len(sol.t), interval=20, blit=False)
plt.show()
