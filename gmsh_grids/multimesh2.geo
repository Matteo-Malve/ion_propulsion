//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {20, 0, 0, 1.0};
//+
Point(3) = {20, 100, 0, 1.0};
//+
Point(401) = {1, 0, 0, 1.0};
//+
Point(402) = {0.025, 0, 0, 1.0};

//+
Point(5) = {-200, 0, 0, 1.0};
//+
Point(6) = {-200, 100, 0, 1.0};
//+
Point(701) = {-1, 0, 0, 1.0};
//+
Point(702) = {-0.025, 0, 0, 1.0};
//+
Point(8) = {20, 1.2, 0, 1.0};

//+
Point(703) = {-4, 0, 0, 1.0};
//+
Point(9) = {0, 4, 0, 1.0};

//+
Point(10) = {20, 4, 0, 1.0};

//+
Circle(1) = {402, 1, 702};
//+
Line(201) = {402, 401};
//+
Line(202) = {401, 2};
//+
Line(3) = {2, 8};
//+
Line(401) = {8, 10};
//+
Line(402) = {10, 3};
//+
Line(5) = {3, 6};
//+
Line(6) = {6, 5};
//+
Line(701) = {5, 703};
//+
Line(702) = {703, 701};
//+
Line(703) = {701, 702};
//+
Circle(8) = {401, 1, 701};

//+
Circle(9) = {703, 1, 9};
//+
Line(10) = {9, 10};

//+
Curve Loop(1) = {701, 9, 10, 402, 5, 6};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {8, 703, -1, 201};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {702, -8, 202, 3, 401, -10, -9};
//+
Plane Surface(3) = {3};

Recombine Surface {1};

Transfinite Surface{2};
Transfinite Line {201} = 15 Using Progression 1.5;
Transfinite Line {-703} = 15 Using Progression 1.5;
Transfinite Line {8} = 16 Using Progression 1;
Transfinite Line {1} = 16 Using Progression 1;
Recombine Surface {2};

Transfinite Line {9} = 8 Using Progression 1;
Transfinite Line {-10} = 30 Using Progression 1.1;
Transfinite Line {-202} = 26 Using Progression 1.1;
Transfinite Line {401} = 12 Using Progression 1.2;
Transfinite Line {3} = 20 Using Progression 1;
Recombine Surface {3};

Mesh.Algorithm = 6;
Mesh.RecombineAll = 1;
Mesh.CharacteristicLengthFactor = 10;
Mesh.SubdivisionAlgorithm = 1;
Mesh.Smoothing = 20;
Show "*";
Mesh 2;

//+
Physical Curve("emettitor", 1) = {1};
//+
Physical Curve("collector", 2) = {3};



//+
Physical Surface("outer", 704) = {1};
//+
Physical Surface("middle", 705) = {3};
//+
Physical Surface("inner", 706) = {2};
