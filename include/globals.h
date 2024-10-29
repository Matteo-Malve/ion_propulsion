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
extern const double Ve; 
extern const double E_ON; 

extern const double nn;
//const double L = nn*R; // [m]
extern const double L ;
extern const double R ;
extern const double X ;
extern const double g; 
extern const double mesh_height ; 

extern const double R; 
extern const double Rc;
extern const double dR2;
extern double AC;
extern double AD;
extern double AE;
extern double AF;

extern int NUM_REFINEMENT_CYCLES;
extern std::string PATH_TO_MESH;
extern const unsigned int NUM_PRELIMINARY_REF; 


#endif
