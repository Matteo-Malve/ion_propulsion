//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {20, 0, 0, 1.0};
//+
Point(3) = {20, 10, 0, 1.0};
//+
Point(4) = {-200, 10, 0, 1.0};
//+
Point(5) = {-200, 0, 0, 1.0};
//+
Point(6) = {0.01, 0, 0, 1.0};
//+
Point(7) = {-0.01, 0, 0, 1.0};
//+
Point(10) = {-10, 10, 0, 1.0};
//+
Point(11) = {10, 10, 0, 1.0};
//+
Point(14) = {-10, 0, 0, 1.0};
//+
Point(15) = {10, 0, 0, 1.0};
//+
Point(16) = {0.01*Cos(Pi/4.), 0.01*Sin(Pi/4.), 0, 1.0};
//+
Point(17) = {-0.01*Cos(Pi/4.), 0.01*Sin(Pi/4.), 0, 1.0};

//+
Line(6) = {2, 3};
//+
Line(8) = {4, 5};
//+
Line(10) = {3, 11};
//+
Line(11) = {11, 10};
//+
Line(12) = {10, 4};
//+
Line(15) = {11, 15};
//+
Line(16) = {10, 14};
//+
Circle(17) = {6, 1, 16};
//+
Circle(18) = {16, 1, 17};
//+
Circle(19) = {17, 1, 7};
//+
Line(25) = {2, 15};
//+
Line(28) = {14, 5};
//+
Line(29) = {6, 15};
//+
Line(30) = {16, 11};
//+
Line(31) = {17, 10};
//+
Line(32) = {7, 14};

//+
Curve Loop(1) = {29, -15, -30, -17};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {30, 11, -31, -18};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {19, 32, -16, -31};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {28, -8, -12, 16};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {25, -15, -10, -6};
//+
Plane Surface(5) = {5};

Transfinite Surface{2};
Transfinite Line {30} = 32 Using Progression 1.2;
Transfinite Line {31} = 32 Using Progression 1.2;
Transfinite Line {18} = 12 Using Progression 1;
Transfinite Line {11} = 12 Using Progression 1;
Recombine Surface{2};

Transfinite Surface{1};
Transfinite Line {29} = 32 Using Progression 1.2;
Transfinite Line {17} = 6 Using Progression 1;
Transfinite Line {15} = 6 Using Progression 1;
Recombine Surface{1};

Transfinite Surface{3};
Transfinite Line {32} = 32 Using Progression 1.2;
Transfinite Line {19} = 6 Using Progression 1;
Transfinite Line {16} = 6 Using Progression 1;
Recombine Surface{3};

Transfinite Surface{4};
Transfinite Line {12} = 120 Using Progression 1;
Transfinite Line {28} = 120 Using Progression 1;
Transfinite Line {8} = 6 Using Progression 1;
Recombine Surface{4};

Transfinite Surface{5};
Transfinite Line {25} = 6 Using Progression 1;
Transfinite Line {10} = 6 Using Progression 1;
Transfinite Line {6} = 6 Using Progression 1;
Recombine Surface{5};

Mesh.Algorithm = 6;
Mesh.RecombineAll = 1;
Mesh.CharacteristicLengthFactor = 2.5;
Mesh.SubdivisionAlgorithm = 1;
Mesh.Smoothing = 20;
Show "*";
Mesh 2;

//+
Physical Curve("emettitor", 1) = {17, 18, 19};
//+
Physical Curve("collector", 2) = {6};

//+
Physical Surface("right", 33) = {5};
//+
Physical Surface("middle", 34) = {1, 2, 3};
//+
Physical Surface("left", 35) = {4};
