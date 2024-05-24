# code
This is our course paper MATLAB code source file, we list the name of the file with specific functions below:

RK_FK1D_Classical.m             ----            Numerical simulation of a one-dimensional classical FK model, 
                                                run to generate particle position images and relative 
                                                position change Plot

RK_FK2D_Classical.m             ----            Using the RK method to simulate the classical two-dimensional
                                                FK model, the initial condition m0 needs to be carefully
                                                adjusted.

RK_FK2D_alpha_F_relationship.m  ----            An image of Fs versus Fc with respect to alpha is generated 
                                                after running the program.

RK_FK2D_Ek_Ep_t_relationship.m  ----            Simulation to verify the energy homogenization theorem.

RK_FK2D_T_F_relationship.m      ----            An interesting code. From the results of this code run we get
                                                images of Fs with respect to temperature T for different kappa,
                                                producing bizarre image results.

RK_FK2D_Temper_distribution_t.m ----            A nice simulations of microscopic heat transfer. We will get
                                                a dynamic image of the temperature map.

RK_FK2D_fig1.m                  ----            v-F and beta-F

RK_FK2D_strange.m               ----            Alpha and Fs, Fc singular peaks at presence temperature

RK_FK2D_v_FcFs_parfor.m         ----            Generation of “force hysteresis loops”
