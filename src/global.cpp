#include "globals.h"

//const double eps_0 = 8.854; // 10^-12 [F/m]= [C^2 s^2 / kg / m^3]
//const double eps_r = 1.0006;
// For step-14:
const double eps_0 = 1.; // 10^-12 [F/m]= [C^2 s^2 / kg / m^3]
const double eps_r = 1.;

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

const double Ri = l;
const double Rc = 5.0*Ri;
const double dR2 = Rc*Rc - Ri*Ri;
double AC = - Ve / (dR2*dR2*dR2);
double AD = 3. * Ve / (dR2*dR2);
double AE = -3. * Ve / dR2;
double AF = Ve;


// ######################################################################
//    MESH
// ######################################################################

const bool READ_FROM_MESH_FILE = true;
//std::string PATH_TO_MESH = "../mesh/input_mesh.msh";
//std::string PATH_TO_MESH = "../mesh/TestSquare.msh";
//std::string PATH_TO_MESH = "../mesh/FullTestSquare.msh";
//std::string PATH_TO_MESH = "../mesh/cerchi_concentrici.msh";
//std::string PATH_TO_MESH = "../mesh/TinyStep14.msh";
std::string PATH_TO_MESH = "../mesh/TinyStep14_1_99.msh";


const unsigned int NUM_PRELIMINARY_REF = 0; 
const unsigned int NUM_PRELIMINARY_GLOBAL_REF = 0; 

// ######################################################################
//    RUN 
// ######################################################################

int NUM_REFINEMENT_CYCLES = 20;

//const std::string REFINEMENT_STRATEGY = "GlobRef";
const std::string REFINEMENT_STRATEGY = "GO";

const std::string GOAL_FUNCTIONAL = "PointValue";
//const std::string GOAL_FUNCTIONAL = "PointYDerivative";
//const std::string GOAL_FUNCTIONAL = "AreaEvaluation";
//const std::string GOAL_FUNCTIONAL = "BoundaryFluxEvaluation";

// ######################################################################
//    CONVERGENCE TESTS
// ######################################################################

const bool ENABLE_CONVERGENCE_ANALYSIS = true;
const bool ENABLE_FLUX_EVALUATION = false;

//const dealii::Point<2> EVALUATION_POINT(0.00025, 0.0005); 
//const dealii::Point<2> EVALUATION_POINT(0.0, 0.001); 
//const dealii::Point<2> EVALUATION_POINT(0.00025, 0.0004125); 
//const dealii::Point<2> EVALUATION_POINT(0.75, 0.75);  // step-14
//const dealii::Point<2> EVALUATION_POINT(0.004, 0.004);  // step-14 ratio 1:100
const dealii::Point<2> EVALUATION_POINT(0.0039, 0.0039);  // step-14 ratio 1:99

//const double EXACT_VALUE = 0.0334473;
//const double EXACT_VALUE = 1.767446e-05;    // step-14 ratio 1:100
const double EXACT_VALUE = 1.742630e-05; // step-14 ratio 1:99


const double EVALUATION_RADIUS = 0.0002;

// ######################################################################
//    OUTPUTS
// ######################################################################

const std::string TEST_NAME = REFINEMENT_STRATEGY + "_" + "replicate-step14_1:99"; 
const unsigned int DUAL_OUTPUT_PATCHES = 1;