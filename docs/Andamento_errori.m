%% Fixed Fraction

clear; close all; clc;

global_error = [0.00257948 0.00066308 0.000161263 3.75914e-05 1.39086e-05 4.10352e-06 7.66208e-07 1.45803e-07 8.33394e-08];

local_error = [0.00257953 0.00066309 0.000161264 3.75917e-05 1.39088e-05 4.10365e-06 7.66198e-07 1.45792e-07 8.3385e-08];

n_dof= [8741 8951 9451 10916 14267 21200 35017 73580 154548];

n_refinement_cycles = 0:8;

plot(n_refinement_cycles,global_error,".-b")
hold on
plot(n_refinement_cycles,local_error,"or")
legend("Global error", "Global error as sum of local contributions")
xlabel("Number of refinements")
ylabel("errors  [Nm^2/C]")
title("Error estimates")

%%
close all

semilogy(n_refinement_cycles,global_error./global_error(1),".-b")
hold on
semilogy(n_refinement_cycles,local_error./local_error(1),"or")

semilogy(n_refinement_cycles,exp(-1.3*n_refinement_cycles),':k')

legend("Global error", "Global error as sum of local contributions","e^{-1.3 n}")
xlabel("Number of refinements")
ylabel("errors  [Nm^2/C]")
title("Error estimates in logaritmic scale")


%%
close all;

loglog(n_dof,global_error,".-b")
hold on
loglog(n_dof,local_error,"or")
legend("Global error", "Global error as sum of local contributions")
xlabel("Number of DoF")
ylabel("errors  [Nm^2/C]")
title("Error estimates")


%% Comparison between two grid refinment algorithms
close all; clear; clc;

n_dof_fraction = [8741 8951 9451 10916 14267 21200 35017 73580 154548];
fraction = [0.00257948 0.00066308 0.000161263 3.75914e-05 1.39086e-05 4.10352e-06 7.66208e-07 1.45803e-07 8.33394e-08];

n_dof_optimize = [8741 9343 11183 17500 25910 35106 86317 135037 371782  ];
optimize = [0.00257948  0.00065952 0.000184072 5.30751e-05 1.76e-05 6.87238e-06 3.91509e-06 2.10714e-06 1.27066e-06 ];

loglog(n_dof_fraction,fraction,".-b")
hold on
loglog(n_dof_optimize,optimize,".-r")
legend("Fixed fraction refinement algorithm", "Optimize refinement algorithm")
xlabel("Number of DoF")
ylabel("errors  [Nm^2/C]")
title("Error estimates comparison between the two refinement algorithms")
