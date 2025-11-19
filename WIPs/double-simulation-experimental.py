import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# ================================
# USER SETTINGS
# ================================
USE_DRAG = True
USE_MAG  = False

# Set the INITIAL SEPARATION between the magnetars
INITIAL_SEPARATION = 0.4e8   # meters (try 4e8, 4e7, etc)

# ================================
# CONSTANTS
# ================================
G = 6.67430e-11
m_sun = 1.989e30
m1 = 1.4 * m_sun
m2 = 1.4 * m_sun
M  = m1 + m2

tm = 200.0                 # used only in drag model
mu0 = 4e-7 * np.pi

# Magnetar magnetic moment (used only if USE_MAG=True)
B = 1e12
R_mag = 20.0
magnetic_moment = (B * R_mag**3) / (2 * mu0)

# drag parameters
alpha = 4.1
beta  = 4.0
gamma = 1.0


# ================================
# HELPER FUNCTIONS
# ================================
def friction_coefficient(r, t):
    return beta / ((t + gamma * tm)**(3/5) * np.exp(alpha * r / INITIAL_SEPARATION))


# ================================
# EQUATIONS OF MOTION
# ================================
def eom(t, S):
    ##print("ODE called at t =", t)
    x1, y1, vx1, vy1, x2, y2, vx2, vy2 = S
    
    r1 = np.array([x1, y1])
    r2 = np.array([x2, y2])
    r_vec = r1 - r2
    dist = np.linalg.norm(r_vec) + 1e-12
    r_hat = r_vec / dist

    v1 = np.array([vx1, vy1])
    v2 = np.array([vx2, vy2])
    v_rel = v1 - v2

    # gravitational
    a1 = -G * m2 * r_vec / dist**3
    a2 =  G * m1 * r_vec / dist**3

    # optional drag
    if USE_DRAG:
        k = friction_coefficient(dist, t)
        a1 += (-k * v_rel) / m1
        a2 += (+k * v_rel) / m2

    # optional magnetic repulsion
    if USE_MAG:
        mu1 = mu2 = magnetic_moment
        Fmag = (3*mu0*mu1*mu2)/(2*np.pi*(dist**4))
        a1 += (Fmag/m1)*r_hat
        a2 -= (Fmag/m2)*r_hat

    return [vx1, vy1, *a1, vx2, vy2, *a2]


# ================================
# INITIAL CONDITIONS
# ================================

# Put the stars at ± separation/2
r1x = +INITIAL_SEPARATION/2
r2x = -INITIAL_SEPARATION/2

# Correct circular velocity for this new separation
v_rel = np.sqrt(G * M / INITIAL_SEPARATION)

# velocity split about COM
v1 =  v_rel * (m2/M)     # +y
v2 = -v_rel * (m1/M)     # -y

state0 = [r1x, 0, 0, v1,
          r2x, 0, 0, v2]

fusion_distance = INITIAL_SEPARATION / 30


# ================================
# SOLVE LONG ENOUGH FOR MERGER
# ================================
t_final = 5e4             # long integration allowed
t_eval = np.linspace(0, t_final, 5000)
print("Starting IVP...")
sol = solve_ivp(eom, [0, t_final], state0,
                t_eval=t_eval, rtol=1e-9, atol=1e-12)
print("Finished SOLVE_IVP")
print("Processing ODE Solution")
x1, y1, x2, y2 = sol.y[0], sol.y[1], sol.y[4], sol.y[5]
print("Finished processing ODE Solution")
sep = np.sqrt((x2-x1)**2 + (y2-y1)**2)
print("Length of t_eval:", len(t_eval))

# ================================
# ANIMATION
# ================================
fig, ax = plt.subplots(figsize=(6,6))
ax.set_facecolor("black")

lim = np.max(sep)*1.2
ax.set_xlim(-lim, lim)
ax.set_ylim(-lim, lim)

star1, = ax.plot([],[], 'o', color='cyan',  markersize=10)
star2, = ax.plot([],[], 'o', color='white', markersize=10)

def update(i):
    print("update called with i =", i)
    if sep[i] < fusion_distance:
        print(f"Fusion at frame {i}, sep={sep[i]}")
        # merged object
        star1.set_data(0,0)
        star2.set_data(0,0)
    else:
        star1.set_data(x1[i], y1[i])
        star2.set_data(x2[i], y2[i])
    return star1, star2

print("Creating animation now...")
ani = FuncAnimation(fig, update, frames=len(t_eval), interval=10, blit=True)
print("Created animation")
print("Calling plt.show()")
plt.show()
