//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {20, 0, 0, 1.0};
//+
Point(3) = {20, 100, 0, 1.0};
//+
Point(4) = {-200, 100, 0, 1.0};
//+
Point(5) = {-200, 0, 0, 1.0};
//+
Point(6) = {0.025, 0, 0, 1.0};
//+
Point(7) = {-0.025, 0, 0, 1.0};
//+
Circle(1) = {6, 1, 7};

//+
Point(8) = {1, 0, 0, 1.0};
//+
Point(9) = {-1, 0, 0, 1.0};
//+
Point(20) = {0, 1, 0, 1.0};
//+
Circle(2) = {8, 1, 20};
//+
Circle(3) = {20, 1, 9};
//+
Circle(26) = {8, 1, 9};


//+
Point(10) = {-4, 0, 0, 1.0};
//+
Point(11) = {0, 4, 0, 1.0};
//+
Point(12) = {20, 4, 0, 1.0};


//+
Point(13) = {-7, 0, 0, 1.0};
//+
Point(14) = {0, 7, 0, 1.0};
//+
Point(15) = {20, 7, 0, 1.0};

//+
Point(16) = {7, 4, 0, 1.0};
//+
Point(17) = {13, 4, 0, 1.0};
//+
Point(18) = {7, 0, 0, 1.0};
//+
Point(19) = {13, 0, 0, 1.0};
//+
Point(21) = {20, 1.2, 0, 1.0};


//+
Line(4) = {21, 12};
//+
Line(5) = {12, 15};
//+
Line(6) = {15, 3};
//+
Line(7) = {3, 4};
//+
Line(8) = {4, 5};
//+
Line(9) = {5, 13};
//+
Line(10) = {13, 10};
//+
Line(11) = {10, 9};
//+
Line(12) = {9, 7};
//+
Line(13) = {6, 8};
//+
Circle(14) = {10, 1, 11};
//+
Circle(15) = {13, 1, 14};


//+
Line(16) = {11, 16};
//+
Line(17) = {16, 17};
//+
Line(18) = {17, 12};
//+
Line(19) = {14, 15};
//+
Line(20) = {16, 18};
//+
Line(21) = {17, 19};
//+
Line(22) = {8, 18};
//+
Line(23) = {18, 19};
//+
Line(24) = {19, 2};

//+
Line(25) = {20, 11};


//+
Line(27) = {2, 21};


//+
Curve Loop(1) = {12, -1, 13, 2,3};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {11, -3, 25, -14};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {25, 16, 20, -22, 2};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {20, 23, -21, -17};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {21, 24, 27, 4, -18};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {15, 19, -5, -18, -17, -16, -14, -10};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {9, 15, 19, 6, 7, 8};
//+
Plane Surface(7) = {7};


Transfinite Line {13} = 25 Using Progression 1.3;
Transfinite Line {-12} = 25 Using Progression 1.3;
Transfinite Line {1} = 17 Using Progression 1;
Transfinite Line {26} = 17 Using Progression 1;
Recombine Surface {1};

Transfinite Line {-11} = 13 Using Progression 1.1;
Transfinite Line {25} = 13 Using Progression 1.1;
Transfinite Line {14} = 9 Using Progression 1;
Transfinite Line {3} = 9 Using Progression 1;
Recombine Surface {2};

Transfinite Line {2} = 9 Using Progression 1;
Transfinite Line {22} = 16 Using Progression 1.1;
Transfinite Line {16} = 18 Using Progression 1;
Transfinite Line {20} = 8 Using Progression 1;
Recombine Surface {3};

Transfinite Line {21} = 8 Using Progression 1;
Transfinite Line {23} = 10 Using Progression 1;
Transfinite Line {17} = 10 Using Progression 1;
Recombine Surface {4};

Transfinite Line {-18} = 18 Using Progression 1;
Transfinite Line {-24} = 30 Using Progression 1.1;
Transfinite Line {4} = 10 Using Progression 1.0U;
Transfinite Line {27} = 20 Using Progression 1.12;
Recombine Surface {5};

Transfinite Line {15} = 9 Using Progression 1;
Transfinite Line {-10} = 6 Using Progression 1.1;
Transfinite Line {5} = 7 Using Progression 1.1;
Transfinite Line {19} = 25 Using Progression 1;
Recombine Surface {6};

Recombine Surface {7};


//+
Physical Curve("emettitor", 1) = {1};
//+
Physical Curve("collector", 2) = {27};

//+
Physical Surface("unique", 3) = {1,2,3,4,5,6, 7};



Mesh.Algorithm = 6;
Mesh.RecombineAll = 1;
Mesh.CharacteristicLengthFactor = 10;
Mesh.SubdivisionAlgorithm = 1;
Mesh.Smoothing = 20;
Show "*";
Mesh 2; 


