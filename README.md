This code simulates the 2D merger of two hypercharged binary neutron stars and the resulting gravitational wave signal:
The code and derivation for the gravitational interaction can be found via HugoGW's repository (https://github.com/HugoGW/Star-merging-and-gravitational-wave-signal).

To calculate the magnetic force interaction I will do the following:

We first start off with the formula for the dipole-dipole interaction: 
$U= \frac{\mu_0 \cdot \overrightarrow[\mu_1] \cdot \overrightarrow[\mu_2] -3(\overrightarrow[\mu_1] \cdot \hat[r])(\overrightarrow[\mu_2] \cdot \hat[r])}{4\pi \cdot r^3} $

To simplify calculations I will assume that both magnetic dipoles have the same magnetic moment. Thus the equation simplifies to:
$U_repulsive = \frac{3\mu_0 \cdot 2\mu_1\mu_2}{4\pi \cdot r^3}$

Next, since the program needs to calculate the repulsive force between the two magnetars, we can take the derivative of the above repulsive energy with respect to the distance between the magnetars:
$F_repulsive = \-frac{dU}{dr}$

We obtain the following repulsive force equation to be used in our calculations:![image](https://github.com/user-attachments/assets/0cc39ff4-d240-4748-9ce9-af9199267b50)


