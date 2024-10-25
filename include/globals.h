#ifndef GLOBALS_H
#define GLOBALS_H

#include <string>
#include <cmath>
#include <iostream>    

using std::cout;
using std::endl;

// Declare global variables using `extern`
extern const double eps_0 ; 
extern const double eps_r;
extern const double Vmax; 
extern const double E_ON; 
//const double R = 4.0e-4; // [m]
//const double Rc = 10.0*R;
//const double dR2 = Rc*Rc - R*R;
//double factor1 = - Vmax / (dR2*dR2*dR2);
//double factor2 = 3. * Vmax / (dR2*dR2);
//double factor3 = -3. * Vmax / dR2;
extern const double nn;
//const double L = nn*R; // [m]
extern const double L ;
extern const double R ;
extern const double X ;
extern const double g; 
extern const double mesh_height ; 

extern int NUM_REFINEMENT_CYCLES;
extern std::string PATH_TO_MESH;
extern const unsigned int NUM_PRELIMINARY_REF; 


#endif
