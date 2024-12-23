#include "Constants.h"
void useGlobalConstants();


double eps_r = GlobalConstants::getInstance().get("eps_r");
double eps_0 = GlobalConstants::getInstance().get("eps_0");

double EXACT_POINT_VALUE = 0.0334473;
 // GlobalConstants::getInstance().get("EXACT_POINT_VALUE");
//double EXACT_FLUX = -3.9252746598790566e-05;  // Rectangle_1_99 by extrapolation
double EXACT_FLUX = -2.193245422464e+00;   // CircularZeroDirichlet exact


unsigned int NUM_PRELIMINARY_GLOBAL_REF = GlobalConstants::getInstance().get("NUM_PRELIMINARY_GLOBAL_REF");;
unsigned int NUM_PRELIMINARY_REF = GlobalConstants::getInstance().get("NUM_PRELIMINARY_REF");;

bool MANUAL_LIFTING_ON = GlobalConstants::getInstance().get("MANUAL_LIFTING_ON");;