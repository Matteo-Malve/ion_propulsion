#include "Constants.h"
void useGlobalConstants();


double eps_r = GlobalConstants::getInstance().get("eps_r");
double eps_0 = GlobalConstants::getInstance().get("eps_0");

double EXACT_POINT_VALUE = GlobalConstants::getInstance().get("EXACT_POINT_VALUE");
double EXACT_FLUX = GlobalConstants::getInstance().get("EXACT_FLUX");

std::string GOAL_FUNCTIONAL = GlobalConstants::getInstance().getString("GOAL_FUNCTIONAL");