#include "Constants.h"
void useGlobalConstants();
double eps_r = GlobalConstants::getInstance().get("eps_r");
double eps_0 = GlobalConstants::getInstance().get("eps_0");