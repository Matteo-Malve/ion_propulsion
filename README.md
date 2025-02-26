# Ion propulsion in atmosphere: Goal-oriented adaptive mesh refinements for corona discharge simulations

*A Master Thesis by Matteo Malvestiti*

---
This is a further development of the homonymous project for the Master course "Numerical Analysis of PDEs", NAPDE in short.

That code can be found in branch "NAPDE_project"

---

### Instructions to build the Makefiles

Move to the "build" directory of the project folder
> cd ion-propulsion/build

Generate your Makefiles, and store them orderly inside the build directory, with the following command:
> cmake ..

You can now choose to build the code with deal.II in debug or release mode by typing either:
> make debug

or

> make release

---
### Instructions for QUICK RUN of the code
To reproduce the results of the thesis, inside the `build` folder run with python 
> python3 ./make_plots.py -n 6

where after -n you specify the number of processors you want to run the MPI code with.
All the necessary simulations will be run and all the plots will be saved in the `build` folder named "Figure_X.png", 
where the name corresponds to the same figure label in the manuscript of the thesis.

You're done 🎉

---

### Instructions to run the code with user defined setup

Open the .yaml configuration file in the root folder and set up the parameters you prefer.
Choose one of the many predefined setups, for examples the ones used inside the manuscript of the thesis are:
- 2: [deal.II's step-14](https://www.dealii.org/current/doxygen/deal.II/step_14.html) recreation (only available in serial code)
- 5: Annulus ratio 1:100, homogeneous forcing term, 20kV on the emitter, 0kV on the collector
- 16: Annulus ratio 1:2, f and BCs as previous
- 18: Real case mesh for corona discharge. Here the voltages can be set up in the configuration file.

Many more options can be chosen in the setup file, like the mapping degree, the refinement algorithm, the dual functional 
to be used and the exact point and flux values that will be used for the convergence analysis.

NOTE: Due to time constraints, not all features have been made available in the development, most notably the manual lifting and the point-value dual functional.

Finally run it with:
> mpirun -np N ./ion_propulsion

### Further developing the code for unlimited new scenarios 

Geometry comes with multilayered information it isn't possible to reduce it to few parameters.
If you want to use a new mesh, with new boundary condition, or a new manifold description or a new forcing term, you can follow these very easy steps!

You will need to change very little in code, namely only `data.h` and `data.cpp` in the include and src directories.

In `data.h`, take setup-0 and change it to your liking, every setup works indeed by overloading a DataSetup class.
Define your rhs function, your exact solution if you need convergence results and a function describing your boundary conditions.

In `data.cpp` change setup-0's create_coarse_grid() function and create/import your coarse mesh. You can also define a manifold description to your liking.

Finally, select setup 0 in the config file and you're ready to go!
