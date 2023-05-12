Function Naca

   /*
    * X, Y, Z, lc
    * point   dX    dY
    * 1       0     0
    * 2       dX2   dY2
    * 3       dX3   dY3
    * 4       dX4   dY4
    *
    * Transfinite:
    * Horizontal: N_H
    * Vertical: N_V
    */
   p1 = newp; Point(p1) = {X,Y,Z,lc}; 
   p2 = newp; Point(p2) = {X+dX2,Y+dY2,Z,lc};
   p3 = newp; Point(p3) = {X+dX3,Y+dY3,Z,lc};
   p4 = newp; Point(p4) = {X+dX4,Y+dY4,Z,lc};

   l1= newl; Line(l1) = {p1, p2} ;
   l2= newl; Line(l2) = {p2, p3} ;
   l3= newl; Line(l3) = {p3, p4} ;
   l4= newl; Line(l4) = {p4, p1} ;

   Transfinite Line {l1, l3} = N_H;
   Transfinite Line {l2, l4} = N_V;

   ll1 = newll; Line Loop(ll1) = {l1, l2, l3, l4};

   pl1 = news;

   Plane Surface (pl1) = {ll1};

   Transfinite Surface {pl1};

   Recombine Surface {pl1};

Return







double a0 = 0.2969;
		double a1 = -0.126;
		double a2 = -0.3516;
		double a3 = 0.2843;
		double a4 = -0.1036;
		double t = 0.5; // Last 2 digits of the NACA divided by 20

		y = t*( a0 * std::sqrt(x) + a1 * x + a2 * pow(x,2.0) + a3 * pow(x,3.0) + a4 * pow(x,4.0) );