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

unsigned int MAPPING_DEGREE;
bool MANUAL_LIFTING_ON;
unsigned int REFINEMENT_CRITERION;
unsigned int DUAL_FUNCTIONAL;

double EVALUATION_POINT_X;
double EVALUATION_POINT_Y;

double EXACT_POINT_VALUE;
double EXACT_FLUX;

unsigned int REFINEMENT_STRATEGY;
double TOP_FRACTION;
double BOTTOM_FRACTION;
unsigned int OPTIMIZE_ORDER;

// Special ones:   (only to generate specific setups that were wrong, on purpose)
unsigned int MANIFOLD_IS_APPLIED;        // 0: no, 1: yes, 2: only on boundary

// Initialize global constants after `GlobalConstants::initialize`
void useGlobalConstants() {
  LOAD_FROM_SETUP = static_cast<unsigned int>(GlobalConstants::getInstance().get("LOAD_FROM_SETUP"));

  eps_r = GlobalConstants::getInstance().get("eps_r");
  eps_0 = GlobalConstants::getInstance().get("eps_0");

  if (LOAD_FROM_SETUP==0) {
    Ve = GlobalConstants::getInstance().get("Ve");
    Vc = GlobalConstants::getInstance().get("Vc");
    //RHS_EXPRESSION = GlobalConstants::getInstance().getString("RHS_EXPRESSION");
  }

  PATH_TO_MESH = GlobalConstants::getInstance().getString("PATH_TO_MESH");
  NUM_CONCENTRIC_REF = static_cast<unsigned int>(GlobalConstants::getInstance().get("NUM_CONCENTRIC_REF",0.0));

  MAPPING_DEGREE = static_cast<unsigned int>(GlobalConstants::getInstance().get("MAPPING_DEGREE",1));
  MANUAL_LIFTING_ON = static_cast<bool>(GlobalConstants::getInstance().get("MANUAL_LIFTING_ON"));
  REFINEMENT_CRITERION = static_cast<unsigned int>(GlobalConstants::getInstance().get("REFINEMENT_CRITERION"));
  DUAL_FUNCTIONAL = static_cast<unsigned int>(GlobalConstants::getInstance().get("DUAL_FUNCTIONAL",1));

  if (LOAD_FROM_SETUP > 0) {
    EVALUATION_POINT_X = GlobalConstants::getInstance().get("EVALUATION_POINT_X");
    EVALUATION_POINT_Y = GlobalConstants::getInstance().get("EVALUATION_POINT_Y");

    EXACT_POINT_VALUE = GlobalConstants::getInstance().get("EXACT_POINT_VALUE");
    EXACT_FLUX = GlobalConstants::getInstance().get("EXACT_FLUX");
  }

  REFINEMENT_STRATEGY = static_cast<unsigned int>(GlobalConstants::getInstance().get("REFINEMENT_STRATEGY", 1));
  if (REFINEMENT_STRATEGY == 1 || REFINEMENT_STRATEGY==2) {
    TOP_FRACTION = GlobalConstants::getInstance().get("TOP_FRACTION", 0.8);
    BOTTOM_FRACTION = GlobalConstants::getInstance().get("BOTTOM_FRACTION", 0.02);
  }
  if (REFINEMENT_STRATEGY==3) {
    OPTIMIZE_ORDER= static_cast<unsigned int>(GlobalConstants::getInstance().get("OPTIMIZE_ORDER", 2));
  }

  MANIFOLD_IS_APPLIED = static_cast<unsigned int>(GlobalConstants::getInstance().get("MANIFOLD_IS_APPLIED",2));
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
  else if (LOAD_FROM_SETUP == 6)
    setup_name="Rectangle_1_99_manifold";
  else if (LOAD_FROM_SETUP == 7)
    setup_name="angle_step14_forced";
  else if (LOAD_FROM_SETUP == 8)
    setup_name="angle_Rectangle_1_100_forced";
  else if (LOAD_FROM_SETUP == 9)
    setup_name="CircularStep14";
  else if (LOAD_FROM_SETUP == 10)
    setup_name="LogCircular_1_2";
  else
    DEAL_II_NOT_IMPLEMENTED();

  std::string refinement_criterion;
  if(REFINEMENT_CRITERION == 1)
    refinement_criterion="global_refinement";
  else if(REFINEMENT_CRITERION == 2)
    refinement_criterion="dual_weighted_error_estimator";
  else if (REFINEMENT_CRITERION == 3)
    refinement_criterion="global* (also evaluate dual-weighted residual)";
  else
    DEAL_II_NOT_IMPLEMENTED();

  std::string functional_name;
  if(DUAL_FUNCTIONAL == 1)
    functional_name="Point evaluation";
  else if(DUAL_FUNCTIONAL == 2)
    functional_name="std. Flux Evaluation";
  else if (DUAL_FUNCTIONAL == 3)
    functional_name="cons. Flux Evaluation";
  else
    DEAL_II_NOT_IMPLEMENTED();

  std::string refinement_strategy;
  if(REFINEMENT_STRATEGY == 1)
    refinement_strategy="fixed_fraction (default)";
  else if(REFINEMENT_STRATEGY == 2)
    refinement_strategy="fixed_number";
  else if (REFINEMENT_STRATEGY == 3)
    refinement_strategy="optimize";
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
  std::cout << std::left << std::setw(30) << "MAPPING_DEGREE"
            << std::setw(20) << MAPPING_DEGREE << "\n";
  std::cout << std::left << std::setw(30) << "MANUAL_LIFTING_ON"
            << std::setw(20) << (MANUAL_LIFTING_ON ? "true" : "false") << "\n";
  std::cout << std::left << std::setw(30) << "REFINEMENT_CRITERION"
            << std::setw(20) << refinement_criterion << "\n";
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
  std::cout << "------------------------------------------------------------\n";
  if (REFINEMENT_CRITERION==2) {
    std::cout << std::left << std::setw(30) << "REFINEMENT_STRATEGY"
              << std::setw(20) << refinement_strategy << "\n";
    if (REFINEMENT_STRATEGY == 1 || REFINEMENT_STRATEGY==2) {
      std::cout << std::left << std::setw(30) << "TOP_FRACTION"
               << std::setw(20) << TOP_FRACTION << "\n";
      std::cout << std::left << std::setw(30) << "BOTTOM_FRACTION"
               << std::setw(20) << BOTTOM_FRACTION << "\n";
    }
    if (REFINEMENT_STRATEGY==3) {
      std::cout << std::left << std::setw(30) << "OPTIMIZE_ORDER"
              << std::setw(20) << OPTIMIZE_ORDER << "\n";
    }
  }
  if (MANIFOLD_IS_APPLIED!=1) {
    std::cout << std::left << std::setw(30) << "MANIFOLD"
              << std::setw(20) << ((MANIFOLD_IS_APPLIED==0) ? "not applied" : "partially applied (only boundary)") << "\n";
  }
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