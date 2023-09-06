
//+
Point(2) = {20, 0, 0, 1.0};
//+
Point(3) = {20, 10, 0, 1.0};
//+
Point(4) = {-200, 10, 0, 1.0};
//+
Point(5) = {-200, 0, 0, 1.0};
//+
Point(6) = {-0.1, 0, 0, 1.0};
//+
Point(7) = {-0.1, 0.02, 0, 1.0};
//+
Point(8) = {0.05, 0.02, 0, 1.0};
//+
Point(9) = {0.1, 0, 0, 1.0};

//+
Point(10) = {0.05, 0.04, 0, 1.0};
//+
Point(11) = {0.1, 0.04, 0, 1.0};
//+
Point(12) = {-1.25, 1, 0, 1.0};
//+
Point(13) = {1.25, 1, 0, 1.0};
//+
Point(14) = {1.25, 0, 0, 1.0};
//+
Point(15) = {-1.25, 0, 0, 1.0};

//+
Point(16) = {-0.1, 0.04, 0, 1.0};
//+
Point(17) = {-1.25, 0.04, 0, 1.0};
//+
Point(18) = {1.25, 0.04, 0, 1.0};

//+
Point(21) = {-200, 1, 0, 1.0};
//+
Point(22) = {20, 1, 0, 1.0};
//+
Point(23) = {20, 0.04, 0, 1.0};
//+
Point(24) = {-200, 0.04, 0, 1.0};

//+
Line(1) = {6, 7};
//+
Line(2) = {7, 8};
//+
Line(3) = {8, 9};
//+
Line(4) = {9, 11};
//+
Line(5) = {11, 10};
//+
Line(6) = {10, 8};
//+
Line(7) = {10, 16};
//+
Line(8) = {16, 7};
//+
Line(9) = {16, 17};
//+
Line(10) = {17, 15};
//+
Line(11) = {17, 12};
//+
Line(12) = {12, 13};
//+
Line(14) = {11, 18};
//+
Line(15) = {18, 13};
//+
Line(17) = {9, 14};
//+
Line(18) = {14, 18};
//+
Line(19) = {14, 2};
//+
Line(21) = {3, 4};
//+
Line(23) = {5, 15};
//+
Line(24) = {15, 6};

//+
Line(26) = {22, 3};
//+
Line(27) = {5, 21};
//+
Line(28) = {21, 4};
//+
Line(29) = {12, 21};
//+
Line(30) = {13, 22};

//+
Line(31) = {16, 6};
//+
Line(32) = {17, 18};
//+
Line(33) = {22, 21};
//+
Line(34) = {15, 12};
//+
Line(35) = {14, 13};
//+
Line(36) = {23, 22};
//+
Line(37) = {2, 23};
//+
Line(38) = {18, 23};
//+
Line(39) = {24, 5};
//+
Line(40) = {24, 21};
//+
Line(41) = {24, 17};


//+
Curve Loop(1) = {2, -6, 7, 8};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {3, 4, 5, 6};
//+
Plane Surface(2) = {2};
//+
Curve Loop(4) = {11, 12, -15, -32};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {21, -28, -33, 26};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {14, -18, -17, 4};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {9, 10, 24, -31};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {10, -23, -39, 41};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {36, -30, -15, 38};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {-18, 19, 37, -38};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {40, -29, -11, -41};
//+
Plane Surface(11) = {11};

Transfinite Surface{:};
Recombine Surface{:};

//Transfinite Line {2,7} = 16 Using Progression 1;
//Transfinite Line {8,6} = 3 Using Progression 1;
//Transfinite Line {3,5} = 6 Using Progression 1;
//Transfinite Line {31,10,39} = 5 Using Progression 1;
//Transfinite Line {14,17} = 116 Using Progression 1;
//Transfinite Line {9,24} = 116 Using Progression 1;
//Transfinite Line {32,12} = 251 Using Progression 1;
//Transfinite Line {11,15,36,40} = 9 Using Progression 1.5;
//Transfinite Line {30,38,19} = 16 Using Progression 1;
//Transfinite Line {29,23,41} = 160 Using Progression 1;

Transfinite Line {2,7} = 10 Using Progression 1;
Transfinite Line {8,6,4,18,37} = 2 Using Progression 1;
Transfinite Line {3,5} = 4 Using Progression 1;
Transfinite Line {31,10,39} = 3 Using Progression 1;
Transfinite Line {14,17} = 70 Using Progression 1;
Transfinite Line {9,24} = 70 Using Progression 1;
Transfinite Line {32,12} = 151 Using Progression 1;
Transfinite Line {11,15,36,40} = 7 Using Progression 1.6;
Transfinite Line {30,38,19} = 16 Using Progression 1;
Transfinite Line {29,23,41} = 160 Using Progression 1;

//+
Physical Line("emitter", 1) = {2, 1, 3};
//+
Physical Line("collector", 2) = {26, 36, 37};

//+
Physical Surface("left", 33) = {11, 8};
//+
Physical Surface("middle", 34) = {4, 7, 1, 2, 6};
//+
Physical Surface("right", 35) = {9, 10};
//+
Physical Surface("big", 36) = {5};

Mesh.Algorithm = 4;
Mesh.RecombineAll = 1;
Mesh.CharacteristicLengthFactor = 2.5;
Mesh.SubdivisionAlgorithm = 1;
Mesh.Smoothing = 0;
Show "*";
Mesh 2;




