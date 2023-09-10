clear; close all; clc;

global_error = [0.000195323 8.85697e-05 5.59325e-05 3.05931e-05 1.86367e-05 9.4473e-06 4.85154e-06 2.54274e-06];

local_error = [0.000195439 8.82656e-05 5.60238e-05 3.05629e-05 1.86367e-05 9.44257e-06 4.85023e-06 2.54266e-06];

n_refinement_cycles = 0:6;


semilogy(n_refinement_cycles,global_error,"-ob")
hold on
semilogy(n_refinement_cycles,local_error,"-or")