# Ion Propulsion:
### Simulation of the electrical field with Goal Oriented adaptive refinement

*A project by Claudia Mallimaci e Matteo Malvestiti*

---

### Instructions to build the Makefiles the first time.

Move to the root directory of the project folder
> cd ion-propulsion 

Generate your Makefiles, and store them orderly inside the build directory, with the following command:
> cmake -S . -B build

---

### Everyday use of the code

To run the code, **move to subdirectory /build**
> cd ion_propulsion/build

Here just type
> make run

Before running this code you may want to open the **data_setup** file that you find in the root directory of the project.
Here you can manipulate all the main parameters of the code, like:
- Max number of refinements
- FE degree
- Mesh geometrical information if you don't use a mesh of your own instead of the standard one
- grid_option:            1. input_mesh                2. doughnut grid for validation
- refinement_algorithm:   1. (default) fixed_fraction, 2. optimize


---
### Understanding the tree of folders

There are many directories in the root folder:

- **build**: Where you go to run the code. Also contains Makefiles and cache files.
- **Results** Here you'll find all the results of the simulation, at each refinement loop, in .vtu format.
You can use software like Paraview to inspect them.
Files named *solution-#.vtu*, with # number of refinement loop, contain the Primal or geometrical solution.
Files named *dual_solution-#.vtu* have no geometrical interpretation and were mostly a tool for us developers. 
Still they can be of interest, and they don't burden a high computational cost to generate.
- **gmsh_grids**: Here you're supposed to place all the meshes you generate with Gmsh and want to be imported in the code.
Our mesh is already present, under the name *horizontal_strip.msh*.
NOTE: Put here only meshes in *.msh* format.
- **src**: All the scripts, subdivided in root, **Foundamentals** and **Mesh**
- **src/Foundamentals** contains most of the scripts
- **src/Mesh** actually contains just one class, named GridForge, which is in charge of the loading of grids

They were divided in the beginning because they worked independently. 
This division lost much of its meaning, after the current inheritance structure was introduced.


---

To use another custom mesh:
- It must have the same macroscopic shape: the emitter must be circular
- Name it input_mesh.msh and save it in directory gmsh_grid
- Or, if in .vtu, add it directly to directory mesh_storage
- Select grid_option = 1

Note: We didn't test such portability.