#include "Constants.h"
void useGlobalConstants();


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
//double EXACT_FLUX = 13962.6;               // Circular

unsigned int NUM_PRELIMINARY_GLOBAL_REF = 0;
unsigned int NUM_PRELIMINARY_REF = 0;

bool MANUAL_LIFTING_ON = 1;    // 1: ON, 0: OFF