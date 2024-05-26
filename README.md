# MATLAB Code Source for Course Paper

This repository contains the MATLAB code source files for our course paper. The following is a list of the files and their respective functions:

1. `RK_FK1D_Classical.m`: Numerical simulation of a one-dimensional classical Frenkel-Kontorova (FK) model, run to generate particle position images and relative position change plots.

2. `RK_FK2D_Classical.m`: Using the Runge-Kutta (RK) method to simulate the classical two-dimensional FK model. The initial condition `m0` needs to be carefully adjusted.

3. `RK_FK2D_alpha_F_relationship.m`: Generates an image of `Fs` (static force) versus `Fc` (critical force) with respect to `alpha`.

4. `RK_FK2D_Ek_Ep_t_relationship.m`: Simulation to verify the energy homogenization theorem.

5. `RK_FK2D_T_F_relationship.m`: An interesting code that generates images of `Fs` with respect to temperature `T` for different `kappa`, producing bizarre image results.

6. `RK_FK2D_Temper_distribution_t.m`: A simulation of microscopic heat transfer, providing a dynamic image of the temperature map.

7. `RK_FK2D_fig1.m`: Generates plots of `v-F` and `beta-F`.

8. `RK_FK2D_strange.m`: Investigates the relationship between `alpha` and `Fs` and `Fc`, including the presence of singular peaks at certain temperatures.

9. `RK_FK2D_v_FcFs_parfor.m`: Generates "force hysteresis loops".

Feel free to explore and utilize these code files for your own research and projects. If you have any questions or need further assistance, please don't hesitate to contact us.
