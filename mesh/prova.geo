// Gmsh project created on Thu Sep 19 14:57:21 2024
//+
Point(1) = {-1, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {1, 0.5, 0, 1.0};
//+
Point(4) = {1, 1, 0, 1.0};
//+
Point(5) = {-1, 0.5, 0, 1.0};
//+
Point(6) = {-1, 1, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 6};
//+
Line(5) = {6, 5};
//+
Line(6) = {5, 1};
//+
Line(7) = {5, 3};
//+
Curve Loop(1) = {5, 7, 3, 4};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {7, -2, -1, -6};
//+
Plane Surface(2) = {2};
//+
Physical Surface(8) = {2};
//+
Physical Surface(9) = {1};
