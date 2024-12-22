#include "Constants.h"
void useGlobalConstants();


double eps_r = GlobalConstants::getInstance().get("eps_r");
double eps_0 = GlobalConstants::getInstance().get("eps_0");

double EXACT_POINT_VALUE = 0.0334473;
 // GlobalConstants::getInstance().get("EXACT_POINT_VALUE");
double EXACT_FLUX = GlobalConstants::getInstance().get("EXACT_FLUX");

unsigned int NUM_PRELIMINARY_GLOBAL_REF = GlobalConstants::getInstance().get("NUM_PRELIMINARY_GLOBAL_REF");;
unsigned int NUM_PRELIMINARY_REF = GlobalConstants::getInstance().get("NUM_PRELIMINARY_REF");;

bool MANUAL_LIFTING_ON = GlobalConstants::getInstance().get("MANUAL_LIFTING_ON");;