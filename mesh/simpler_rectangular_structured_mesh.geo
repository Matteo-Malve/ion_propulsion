
Point(6) = {0.00040, 0, 0, 1.0};
Point(7) = {0.00040, 0.00040, 0, 1.0};
Point(8) = {-0.00040, 0.00040, 0, 1.0};
Point(9) = {-0.00040, 0, 0, 1.0};
Point(10) = {0.004, 0, 0, 1.0};
Point(11) = {0.004, 0.004, 0, 1.0};
Point(12) = {-0.004, 0.004, 0, 1.0};
Point(13) = {-0.004, 0, 0, 1.0};
Point(14) = {-0.004, 0.0004, 0, 1.0};
Point(15) = {0.004, 0.0004, 0, 1.0};
Point(19) = {0.0004, 0.004, 0, 1.0};
Point(20) = {-0.0004, 0.004, 0, 1.0};
//+
Line(1) = {13, 14};
//+
Line(2) = {14, 12};
//+
Line(3) = {12, 20};
//+
Line(4) = {20, 19};
//+
Line(5) = {19, 11};
//+
Line(6) = {11, 15};
//+
Line(7) = {15, 10};
//+
Line(8) = {10, 6};
//+
Line(9) = {6, 7};
//+
Line(10) = {7, 8};
//+
Line(11) = {8, 9};
//+
Line(12) = {9, 13};
//+
Line(13) = {14, 8};
//+
Line(14) = {7, 15};
//+
Line(15) = {19, 7};
//+
Line(16) = {20, 8};
//+
Curve Loop(1) = {2, 3, 16, -13};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {1, 13, 11, 12};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {16, -10, -15, -4};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {15, 14, -6, -5};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {9, 14, 7, 8};
//+
Plane Surface(5) = {5};
//+
Transfinite Surface {1};
//+
Transfinite Surface {3};
//+
Transfinite Surface {4};
//+
Transfinite Surface {5};
//+
Transfinite Surface {2};
//+
Transfinite Curve {2, 3, 16, 15, 5, 6, 14, 13, 12, 8} = 37 Using Progression 1;
//+
Transfinite Curve {1, 11, 9, 7} = 5 Using Progression 1;
//+
Transfinite Curve {10, 4} = 9 Using Progression 1;
//+
Physical Curve("Collector", 2) = {6, 7};
//+
Physical Curve("Emitter", 1) = {11, 10, 9};


Mesh.Algorithm = 6;
Mesh.RecombineAll = 1;
//Mesh.CharacteristicLengthFactor = 2.5;
Mesh.SubdivisionAlgorithm = 1;
Mesh.Smoothing = 20;
Show "*";
Mesh 2;