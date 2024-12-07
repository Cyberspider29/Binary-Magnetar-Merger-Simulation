#2D-Binary-Magnetar-Simulation#

This code simulates the 2D merger of two hypercharged binary neutron stars (magnetars) and the resulting gravitational wave signal:
The code and derivation for the gravitational interaction can be found via HugoGW's repository (https://github.com/HugoGW/Star-merging-and-gravitational-wave-signal).
However, I have made some major changes regarding the code since I needed to calculate the orbital eccentricity for each frame. Thus, in this code I have replaced the `odeint` differential equation solver with `solve_ivp`.

The eccentricity plots return the orbital eccentricity and progression of the binary magnetar system.
The strain plots return the resulting gravitational strain as a result of the interactions of the magnetars during each point in time.

To calculate the magnetic force interaction I will do the following:

We first start off with the formula for the dipole-dipole interaction: $U= \frac{\mu_0 \cdot \overrightarrow[\mu_1] \cdot \overrightarrow[\mu_2] -3(\overrightarrow[\mu_1] \cdot \hat[r])(\overrightarrow[\mu_2] \cdot \hat[r])}{4\pi \cdot r^3} $

To simplify calculations I will assume that both magnetic dipoles have the same magnetic moment. Thus the equation simplifies to: 
$U_r = \frac{3\mu_0 \cdot 2\mu_1\mu_2}{4\pi \cdot r^3}$

Next, since the program needs to calculate the repulsive force between the two magnetars, we can take the derivative of the above repulsive energy with respect to the distance between the magnetars: $F_r = -\frac{dU}{dr}$

We obtain the following repulsive force equation to be used in our calculations: $F_r = \frac{3\mu_0 \cdot 2\mu_1\mu_2}{4\pi \cdot r^4}$

To calculate the orbital eccentricity we use Kepler's 3rd Law to obtain the following: $e = \frac{r_a - r_p}{r_a + r_p}$
The `solve_ivp` integrator will allow me to get the maximum and minimum distances for each magnetar at each frame, allowing me to return the eccentricity for each frame as well.

