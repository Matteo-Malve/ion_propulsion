// Gmsh project created on Tue Sep 17 11:36:43 2024
//+
Point(1) = {-2, 0, 0, 1.0};
//+
Point(2) = {-2, 0.1, 0, 1.0};
//+
Point(3) = {-0.004, 0.004, 0, 1.0};
//+
Point(4) = {0.004, 0.004, 0, 1.0};
//+
Point(5) = {0.2, 0.1, 0, 1.0};
//+
Point(6) = {0.2, 0, 0, 1.0};
//+
Point(7) = {0.004, 0, 0, 1.0};
//+
Point(8) = {0.00040, 0, 0, 1.0};
//+
Point(9) = {0.00040, 0.00040, 0, 1.0};
//+
Point(10) = {-0.00040, 0.00040, 0, 1.0};
//+
Point(11) = {-0.00040, 0, 0, 1.0};
//+
Point(12) = {-0.004, 0, 0, 1.0};

//+
Line(1) = {1, 2};
//+
Line(2) = {2, 5};
//+
Line(3) = {3, 4};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 7};
//+
Line(7) = {7, 8};
//+
Line(8) = {8, 9};
//+
Line(9) = {9, 10};
//+
Line(10) = {10, 11};
//+
Line(11) = {11, 12};
//+
Line(12) = {12, 1};
//+
Physical Curve("Emitter", 1) = {8, 9, 10};
//+
Physical Curve("Collector", 2) = {5};

//+
Line(13) = {3, 12};
//+
Line(14) = {4, 7};

//+
Curve Loop(1) = {12, 1, 2, 5, 6, -14, -3, 13};
//+
Plane Surface(1) = {1};
//Physical Surface("Surface exrternal") = {1};

//+
Curve Loop(2) = {11, -13, 3, 14, 7, 8, 9, 10};
//+
Plane Surface(2) = {2};
//Physical Surface("Surface center") = {2};

Characteristic Length {1,2,5,6,7,4,3,11} = 0.01;   // External coarse mesh
Characteristic Length {12, 3, 4, 7, 8, 9, 10, 11} = 0.0001;   // Center finer mesh

Mesh.Algorithm = 8;  // Use the frontal-Delaunay algorithm for quadrilaterals

Recombine Surface {1};
Recombine Surface {2};

Mesh 2;