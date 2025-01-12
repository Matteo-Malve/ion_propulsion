Point(1) = {0, 0, 0, 1.0};
Point(2) = {-0.04, 0, 0, 1.0};
Point(3) = {0.04, 0, 0, 1.0};
Point(4) = {-0.0004, 0, 0, 1.0};
Point(5) = {0.0004, 0, 0, 1.0};
//+
Circle(1) = {2, 1, 3};
//+
Circle(2) = {3, 1, 2};
//+
Circle(3) = {4, 1, 5};
//+
Circle(4) = {5, 1, 4};
//+
Line(5) = {2, 4};
//+
Line(6) = {5, 3};
//+
Curve Loop(1) = {2, 5, -4, 6};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {5, 3, 6, -1};
//+
Plane Surface(2) = {2};
//+
Transfinite Surface {1};
//+
Transfinite Surface {2};
//+
Transfinite Curve {4, 2, 1, 3} =31 Using Progression 1;
//+
Transfinite Curve {-5, 6} = 97 Using Progression 1;
//+
Physical Surface(7) = {1, 2};
//+
Physical Curve("emitter", 1) = {4, 3};
//+
Physical Curve("outer_shell", 9) = {2, 1};
