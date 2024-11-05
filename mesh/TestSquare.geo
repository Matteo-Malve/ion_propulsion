
Point(8) = {0.00040, 0, 0, 1.0};
Point(9) = {0.00040, 0.00040, 0, 1.0};
Point(10) = {-0.00040, 0.00040, 0, 1.0};
Point(11) = {-0.00040, 0, 0, 1.0};

Point(13) = {0.004, 0, 0, 1.0};
Point(14) = {0.004, 0.004, 0, 1.0};
Point(15) = {-0.004, 0.004, 0, 1.0};
Point(16) = {-0.004, 0, 0, 1.0};

Point(17) = {0.004, 0.0004, 0, 1.0};
Point(18) = {-0.004, 0.0004, 0, 1.0};
Point(19) = {0.0004, 0.004, 0, 1.0};
Point(20) = {-0.0004, 0.004, 0, 1.0};
//+
Line(1) = {16, 18};
//+
Line(2) = {18, 15};
//+
Line(3) = {15, 20};
//+
Line(4) = {20, 19};
//+
Line(5) = {19, 14};
//+
Line(6) = {14, 17};
//+
Line(7) = {17, 13};
//+
Line(8) = {13, 8};
//+
Line(9) = {8, 9};
//+
Line(10) = {9, 10};
//+
Line(11) = {10, 11};
//+
Line(12) = {11, 16};
//+
Line(13) = {10, 18};
//+
Line(14) = {10, 20};
//+
Line(15) = {9, 19};
//+
Line(16) = {9, 17};
//+
Curve Loop(1) = {13, 2, 3, -14};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {10, 14, 4, -15};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {15, 5, 6, -16};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {16, 7, 8, 9};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {13, -1, -12, -11};
//+
Plane Surface(5) = {5};
//+
Physical Curve("emitter", 1) = {9, 10, 11};
Physical Curve("floor", 3) = {12,8};
Physical Curve("others", 9) = {1, 2, 3, 4, 5, 6, 7};
//+
Physical Surface("all", 9) = {5, 1, 2, 3, 4};
//+
Transfinite Surface {5} = {16, 18, 10, 11};
//+
Transfinite Surface {1};
//+
Transfinite Surface {2};
//+
Transfinite Surface {3};
//+
Transfinite Surface {4};
//+
Transfinite Curve {11, 9, 7, 1} = 5 Using Progression 1;
//+
Transfinite Curve {2, 3, 13, 12, 14, 15, 5, 6, 16, 8} = 37 Using Progression 1;
//+
Transfinite Curve {10, 4} = 9 Using Progression 1;
