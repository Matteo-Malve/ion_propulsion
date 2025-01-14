#include "Constants.h"

#include <deal.II/base/exceptions.h>

GlobalConstants* GlobalConstants::instance = nullptr;

// Declare global variables as uninitialized

std::string OUTPUT_PATH;

double eps_r;
double eps_0;

double Ve;
double Vc;
std::string RHS_EXPRESSION;

//unsigned int NUM_PRELIMINARY_GLOBAL_REF;
std::string PATH_TO_MESH;
unsigned int NUM_CONCENTRIC_REF;
unsigned int LOAD_FROM_SETUP;

bool MANUAL_LIFTING_ON;
unsigned int REFINEMENT_CRITERION;
unsigned int DUAL_FUNCTIONAL;

double EVALUATION_POINT_X;
double EVALUATION_POINT_Y;

double EXACT_POINT_VALUE;
double EXACT_FLUX;

// Initialize global constants after `GlobalConstants::initialize`
void useGlobalConstants() {
  LOAD_FROM_SETUP = static_cast<unsigned int>(GlobalConstants::getInstance().get("LOAD_FROM_SETUP"));

  eps_r = GlobalConstants::getInstance().get("eps_r");
  eps_0 = GlobalConstants::getInstance().get("eps_0");

  if (LOAD_FROM_SETUP==0) {
    Ve = GlobalConstants::getInstance().get("Ve");
    Vc = GlobalConstants::getInstance().get("Vc");
    RHS_EXPRESSION = GlobalConstants::getInstance().getString("RHS_EXPRESSION");
  }

  PATH_TO_MESH = GlobalConstants::getInstance().getString("PATH_TO_MESH");
  NUM_CONCENTRIC_REF = GlobalConstants::getInstance().get("NUM_CONCENTRIC_REF",true);

  MANUAL_LIFTING_ON = static_cast<bool>(GlobalConstants::getInstance().get("MANUAL_LIFTING_ON"));
  REFINEMENT_CRITERION = static_cast<unsigned int>(GlobalConstants::getInstance().get("REFINEMENT_CRITERION"));
  if (REFINEMENT_CRITERION>1)
    DUAL_FUNCTIONAL = static_cast<unsigned int>(GlobalConstants::getInstance().get("DUAL_FUNCTIONAL"));

  EVALUATION_POINT_X = GlobalConstants::getInstance().get("EVALUATION_POINT_X");
  EVALUATION_POINT_Y = GlobalConstants::getInstance().get("EVALUATION_POINT_Y");

  EXACT_POINT_VALUE = GlobalConstants::getInstance().get("EXACT_POINT_VALUE");
  EXACT_FLUX = GlobalConstants::getInstance().get("EXACT_FLUX");
}

//unsigned int NUM_PRELIMINARY_REF = GlobalConstants::getInstance().get("NUM_PRELIMINARY_REF");;


void printParsedConstants() {
  std::string setup_name;
  if(LOAD_FROM_SETUP == 0)
    setup_name="none";
  else if(LOAD_FROM_SETUP == 1)
    setup_name="CurvedRidges";
  else if(LOAD_FROM_SETUP == 2)
    setup_name="Exercise_2_3";
  else if (LOAD_FROM_SETUP == 3)
    setup_name="Rectangle_1_99";
  else if (LOAD_FROM_SETUP == 4)
    setup_name="LogCircular_1_10";
  else if (LOAD_FROM_SETUP == 5)
    setup_name="LogCircular_1_100";
  else
    DEAL_II_NOT_IMPLEMENTED();

  std::string functional_name;
  if(DUAL_FUNCTIONAL == 1)
    setup_name="Point evaluation";
  else if(DUAL_FUNCTIONAL == 2)
    setup_name="std. Flux Evaluation";
  else if (DUAL_FUNCTIONAL == 3)
    setup_name="cons. Flux Evaluation";
  else
    DEAL_II_NOT_IMPLEMENTED();

  std::cout << "\n===================== Parsed Constants =====================\n";
  std::cout << std::left << std::setw(30) << "Variable Name"
            << std::setw(20) << "Value" << "\n";
  std::cout << "------------------------------------------------------------\n";

  // Print the actual variables
  std::cout << std::left << std::setw(30) << "eps_r"
            << std::setw(20) << eps_r << "\n";
  std::cout << std::left << std::setw(30) << "eps_0"
            << std::setw(20) << eps_0 << "\n";

  std::cout << "------------------------------------------------------------\n";
  std::cout << std::left << std::setw(30) << "LOAD_FROM_SETUP"
              << std::setw(20) << setup_name << "\n";
  //std::cout << std::left << std::setw(30) << "NUM_PRELIMINARY_GLOBAL_REF"
  //          << std::setw(20) << NUM_PRELIMINARY_GLOBAL_REF << "\n";
  std::cout << std::left << std::setw(30) << "PATH_TO_MESH"
           << std::setw(20) << PATH_TO_MESH << "\n";
  std::cout << std::left << std::setw(30) << "NUM_CONCENTRIC_REF"
            << std::setw(20) << NUM_CONCENTRIC_REF << "\n";

  if (LOAD_FROM_SETUP==0){
    std::cout << "------------------------------------------------------------\n";
    std::cout << std::left << std::setw(30) << "Ve"
              << std::setw(20) << Ve << "\n";
    std::cout << std::left << std::setw(30) << "Vc"
              << std::setw(20) << Vc << "\n";
    std::cout << std::left << std::setw(30) << "RHS_EXPRESSION"
              << std::setw(20) << RHS_EXPRESSION << "\n";
    //std::cout << std::left << std::setw(30) << "UEX_EXPRESSION"
    //          << std::setw(20) << UEX_EXPRESSION << "\n";
    }

  std::cout << "------------------------------------------------------------\n";
  std::cout << std::left << std::setw(30) << "MANUAL_LIFTING_ON"
            << std::setw(20) << (MANUAL_LIFTING_ON ? "true" : "false") << "\n";
  std::cout << std::left << std::setw(30) << "REFINEMENT_CRITERION"
            << std::setw(20) << ((REFINEMENT_CRITERION==1) ? "global_refinement" : "dual_weighted_error_estimator") << "\n";
  if (REFINEMENT_CRITERION>1) {
    std::cout << std::left << std::setw(30) << "DUAL_FUNCTIONAL" << std::setw(20) << functional_name << "\n";
  }
  std::cout << "------------------------------------------------------------\n";
  std::string eval_point_alltogether = "( " + std::to_string(EVALUATION_POINT_X) + " , " + std::to_string(EVALUATION_POINT_Y) + " )";
  std::cout << std::left << std::setw(30) << "EVALUATION_POINT"
             << std::setw(20) << eval_point_alltogether << "\n";
  std::cout << "------------------------------------------------------------\n";
  if (std::fabs(EXACT_POINT_VALUE)>1.e-10) {}
    std::cout << std::left << std::setw(30) << "EXACT_POINT_VALUE"
              << std::setw(20) << EXACT_POINT_VALUE << "\n";
  if (std::fabs(EXACT_FLUX)>1.e-10)
    std::cout << std::left << std::setw(30) << "EXACT_FLUX"
              << std::setw(20) << EXACT_FLUX << "\n";

  std::cout << "============================================================\n\n";
}
/*
double eps_r = 1.0;     // 8.854
double eps_0 = 1.0; // 1.0006

//double EXACT_POINT_VALUE = 0.0334473; // step14 classic (0.75, 0.75)
//double EXACT_POINT_VALUE = 1.742630e-05;       // 1:99 (0.0039, 0.0039)
//double EXACT_POINT_VALUE = 9030.9;        // LogCircular_1_10, (0.001,0.001)
//double EXACT_POINT_VALUE = 6466.13;      // LogCircular_1_10, (0.0019,0)
double EXACT_POINT_VALUE = 3148.182801 ;       // LogCircular_1_100, (0.019375,0.)

//double EXACT_FLUX = -3.9252746598790566e-05;  // Rectangle_1_99 by extrapolation
//double EXACT_FLUX = -2.193245422464e+00;   // CircularZeroDirichlet exact
double EXACT_FLUX = 54575.1;             // LogCircular 1:10 exact
double EXACT_FLUX = 27287.5;             // LogCircular 1:100 exact
//double EXACT_FLUX = 13962.6;               // Circular

unsigned int NUM_PRELIMINARY_GLOBAL_REF = 0;
unsigned int NUM_PRELIMINARY_REF = 0;

bool MANUAL_LIFTING_ON = 1;    // 1: ON, 0: OFF

*/