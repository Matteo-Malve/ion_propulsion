#include "globals.h"

const double eps_0 = 8.854; // 10^-12 [F/m]= [C^2 s^2 / kg / m^3]
const double eps_r = 1.0006;
const double Ve = 2.e+4; // [V]
const double E_ON = 3.31e+6; // [V/m] corona inception threshold

const double nn = 2;
//const double L = nn*R; // [m]
const double L = 0.0004;
const double R = std::sqrt(2.)*L;
const double X = 0.;//-L/2.; // [m]
const double g = 0.2; // [m]
const double mesh_height = 0.1;; // [m]

const double Rc = 3.0*R;
const double dR2 = Rc*Rc - R*R;
double factor1 = - Ve / (dR2*dR2*dR2);
double factor2 = 3. * Ve / (dR2*dR2);
double factor3 = -3. * Ve / dR2;


std::string PATH_TO_MESH = "../mesh/TestSquare.msh";
const unsigned int NUM_PRELIMINARY_REF = 4; // Meglio 7 ma poi troppo lento
int NUM_REFINEMENT_CYCLES = 4;
