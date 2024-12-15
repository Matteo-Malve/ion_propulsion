#include "Constants.h"
void useGlobalConstants();


double eps_r = GlobalConstants::getInstance().get("eps_r");
double eps_0 = GlobalConstants::getInstance().get("eps_0");

double EXACT_POINT_VALUE = GlobalConstants::getInstance().get("EXACT_POINT_VALUE");
double EXACT_FLUX = GlobalConstants::getInstance().get("EXACT_FLUX");

unsigned int NUM_PRELIMINARY_GLOBAL_REF = GlobalConstants::getInstance().get("NUM_PRELIMINARY_GLOBAL_REF");;
unsigned int NUM_PRELIMINARY_REF = GlobalConstants::getInstance().get("NUM_PRELIMINARY_REF");;

std::string GOAL_FUNCTIONAL = GlobalConstants::getInstance().getString("GOAL_FUNCTIONAL");