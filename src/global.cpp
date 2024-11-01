#include "globals.h"

const double eps_0 = 8.854; // 10^-12 [F/m]= [C^2 s^2 / kg / m^3]
const double eps_r = 1.0006;
const double Ve = 2.e+4; // [V]
const double E_ON = 3.31e+6; // [V/m] corona inception threshold

const double nn = 2;
//const double l = nn*R; // [m]
const double l = 0.0004;
const double L = 0.004;
const double R = std::sqrt(2.)*l;
const double X = 0.;//-l/2.; // [m]
const double g = 0.2; // [m]
const double mesh_height = 0.1;; // [m]

const double Rc = 3.0*R;
const double dR2 = Rc*Rc - R*R;
double AC = - Ve / (dR2*dR2*dR2);
double AD = 3. * Ve / (dR2*dR2);
double AE = -3. * Ve / dR2;
double AF = Ve;

//std::string PATH_TO_MESH = "../mesh/input_mesh.msh";
//std::string PATH_TO_MESH = "../mesh/TestSquare.msh";
std::string PATH_TO_MESH = "../mesh/FullTestSquare.msh";
//std::string PATH_TO_MESH = "../mesh/cerchi_concentrici.msh";

const unsigned int NUM_PRELIMINARY_REF = 4; 
const unsigned int NUM_PRELIMINARY_GLOBAL_REF = 0; 

int NUM_REFINEMENT_CYCLES = 12;

const std::string TEST_NAME = "4_0_12"; 

