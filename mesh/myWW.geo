//+
SetFactory("Built-in");
re = 110e-6;
Ly  = 12e-2;
Lx =8e-2;
rc = 3.175e-3;

//Set a spacing///////
d=60e-3;

//SET A GAP
gap=50e-3;

//////DOMAIN/////////

Point(11) = {-Lx,    -Ly/2,  0};
Point(12) = {-Lx,     Ly/2,  0};
Point(41) = {Lx+gap,  Ly/2,  0};
Point(43) = {Lx+gap, -Ly/2,  0};


////EMITTERS//////
h  = d/2;
k = 100;
kk = 5;
//First Emitter

Point(21) = {0,   h,    0};
Point(22) = {re,  h,    0};
Point(23) = {0,   re+h, 0};
Point(24) = {-re, h,    0};
Point(25) = {0,   h-re, 0};
Circle(21) = {22, 21, 23};
Circle(22) = {23, 21, 24};
Circle(23) = {24, 21, 25};
Circle(24) = {25, 21, 22};

Point(122) = {k*re,  h,    0};
Point(123) = {0,   k*re+h, 0};
Point(124) = {-k*re, h,    0};
Point(125) = {0,   h-k*re, 0};
Circle(121) = {122, 21, 123};
Circle(122) = {123, 21, 124};
Circle(123) = {124, 21, 125};
Circle(124) = {125, 21, 122};

//Second Emitter
Point(31) = {0,   -h,     0};
Point(32) = {re,  -h,     0};
Point(33) = {0,   -h+re,  0};
Point(34) = {-re, -h,     0};
Point(35) = {0,  -(re+h), 0};
Circle(31) = {32, 31, 33};
Circle(32) = {33, 31, 34};
Circle(33) = {34, 31, 35};
Circle(34) = {35, 31, 32};

Point(132) = {k*re,  -h,     0};
Point(133) = {0,   -h+k*re,  0};
Point(134) = {-k*re, -h,     0};
Point(135) = {0,  -(k*re+h), 0};
Circle(131) = {132, 31, 133};
Circle(132) = {133, 31, 134};
Circle(133) = {134, 31, 135};
Circle(134) = {135, 31, 132};

////Collectors//////
//First Collector

Point(6) = {gap+rc,  h,    0};
Point(7) = {gap+2*rc,h,    0};
Point(8) = {gap+rc,  h+rc, 0};
Point(9) = {gap,     h,    0};
Point(10) = {gap+rc, h-rc, 0};
Circle(5) = {7, 6, 8};
Circle(6) = {8, 6, 9};
Circle(7) = {9, 6, 10};
Circle(8) = {10, 6, 7};

Point(107) = {gap+rc+k*re,h,    0};
Point(108) = {gap+rc,  k*re+h, 0};
Point(109) = {gap+rc-k*re,     h,    0};
Point(110) = {gap+rc, h-k*re, 0};
Circle(105) = {107, 6, 108};
Circle(106) = {108, 6, 109};
Circle(107) = {109, 6, 110};
Circle(108) = {110, 6, 107};

//Second Collector
Point(50) = {gap+rc,  -h,     0};
Point(51) = {gap+2*rc,-h,     0};
Point(52) = {gap+rc, -(h+rc), 0};
Point(53) = {gap,     -h,     0};
Point(54) = {gap+rc, -(h-rc), 0};
Circle(10) = {54, 50, 53};
Circle(11) = {53, 50, 52};
Circle(12) = {52, 50, 51};
Circle(13) = {51, 50, 54};

Point(151) = {gap+rc+k*re,-h,     0};
Point(152) = {gap+rc, -h+k*re, 0};
Point(153) = {gap+rc-k*re,     -h,     0};
Point(154) = {gap+rc, -(k*re+h), 0};
Circle(110) = {154, 50, 153};
Circle(111) = {153, 50, 152};
Circle(112) = {152, 50, 151};
Circle(113) = {151, 50, 154};


Line(135) = {23, 123};
Line(136) = {22, 122};
Line(137) = {25, 125};
Line(138) = {24, 124};
Line(139) = {33, 133};
Line(140) = {32, 132};
Line(141) = {35, 135};
Line(142) = {34, 134};
Line(143) = {8, 108};
Line(144) = {9, 109};
Line(145) = {10, 110};
Line(146) = {152, 54};
Line(147) = {51, 151};
Line(148) = {52, 154};
Line(149) = {53, 153};
Point(155) = {-0.1, -0, 0, 1.0};
Point(156) = {-0.08, -0.03, 0, 1.0};
Point(157) = {-0.08, 0, 0, 1.0};
Point(158) = {-0.08, 0.03, 0, 1.0};
Point(159) = {0.13, -0.03, 0, 1.0};
Point(160) = {0.13, 0.03, 0, 1.0};
Point(161) = {0.13, 0., 0, 1.0};
Point(162) = {0, -0.06, 0, 1.0};
Point(163) = {0.0265875, -0.06, 0, 1.0};
Point(164) = {0.0265875, 0.06, 0, 1.0};
Point(165) = {0.0265875, -0.03, 0, 1.0};
Point(166) = {0.0265875, 0.03, 0, 1.0};
Point(167) = {0.0265875, 0, 0, 1.0};
Point(168) = {0, 0, 0, 1.0};
Point(169) = {0.053175, 0, 0, 1.0};
Point(170) = {0.053175, -0.06, 0, 1.0};
Point(171) = {0.053175, 0.06, 0, 1.0};
Point(172) = {0, 0.06, 0, 1.0};

Line(150) = {158, 124};
Line(151) = {123, 172};
Line(152) = {122, 166};
Line(153) = {166, 109};
Line(154) = {166, 164};
Line(155) = {108, 171};
Line(156) = {7, 107};
Line(157) = {107, 160};
Line(158) = {110, 169};
Line(159) = {169, 152};
Line(160) = {169, 161};
Line(161) = {151, 159};
Line(162) = {154, 170};
Line(163) = {153, 165};
Line(164) = {165, 163};
Line(165) = {165, 132};
Line(166) = {135, 162};
Line(167) = {167, 169};
Line(168) = {165, 167};
Line(169) = {167, 166};
Line(170) = {125, 168};
Line(171) = {168, 133};
Line(172) = {167, 168};
Line(173) = {168, 157};
Line(174) = {134, 156};
Line(175) = {12, 158};
Line(176) = {158, 157};
Line(177) = {157, 156};
Line(178) = {156, 11};
Line(179) = {11, 162};
Line(180) = {162, 163};
Line(181) = {163, 170};
Line(182) = {170, 43};
Line(183) = {43, 159};
Line(184) = {159, 161};
Line(185) = {161, 160};
Line(186) = {160, 41};
Line(187) = {41, 171};
Line(188) = {171, 164};
Line(189) = {164, 172};
Line(190) = {172, 12};
Curve Loop(2) = {190, 175, 150, -122, 151};
Plane Surface(1) = {2};
Curve Loop(3) = {150, 123, 170, 173, -176};
Plane Surface(2) = {3};
Curve Loop(4) = {173, 177, -174, -132, -171};
Plane Surface(3) = {4};
Curve Loop(5) = {174, 178, 179, -166, -133};
Plane Surface(4) = {5};
Curve Loop(6) = {166, 180, -164, 165, -134};
Plane Surface(5) = {6};
Curve Loop(7) = {164, 181, -162, 110, 163};
Plane Surface(6) = {7};
Curve Loop(8) = {162, 182, 183, -161, 113};
Plane Surface(7) = {8};
Curve Loop(9) = {159, 112, 161, 184, -160};
Plane Surface(8) = {9};
Curve Loop(10) = {158, 160, 185, -157, -108};
Plane Surface(9) = {10};
Curve Loop(11) = {105, 155, -187, -186, -157};
Plane Surface(10) = {11};
Curve Loop(12) = {188, -154, 153, -106, 155};
Plane Surface(11) = {12};
Curve Loop(13) = {189, -151, -121, 152, 154};
Plane Surface(12) = {13};
Curve Loop(14) = {170, -172, 169, -152, -124};
Plane Surface(13) = {14};
Curve Loop(15) = {169, 153, 107, 158, -167};
Plane Surface(14) = {15};
Curve Loop(16) = {172, 171, -131, -165, 168};
Plane Surface(15) = {16};
Curve Loop(17) = {167, 159, -111, 163, 168};
Plane Surface(16) = {17};
Curve Loop(18) = {122, -138, -22, 135};
Plane Surface(17) = {18};
Curve Loop(19) = {138, 123, -137, -23};
Plane Surface(18) = {19};
Curve Loop(20) = {24, 136, -124, -137};
Plane Surface(19) = {20};
Curve Loop(21) = {135, -121, -136, 21};
Plane Surface(20) = {21};
Curve Loop(22) = {132, -142, -32, 139};
Plane Surface(21) = {22};
Curve Loop(23) = {142, 133, -141, -33};
Plane Surface(22) = {23};
Curve Loop(24) = {141, 134, -140, -34};
Plane Surface(23) = {24};
Curve Loop(25) = {139, -131, -140, 31};
Plane Surface(24) = {25};
Curve Loop(26) = {111, 146, 10, 149};
Plane Surface(25) = {26};
Curve Loop(27) = {149, -110, -148, -11};
Plane Surface(26) = {27};
Curve Loop(28) = {12, 147, 113, -148};
Plane Surface(27) = {28};
Curve Loop(29) = {146, -13, 147, -112};
Plane Surface(28) = {29};
Curve Loop(30) = {106, -144, -6, 143};
Plane Surface(29) = {30};
Curve Loop(31) = {144, 107, -145, -7};
Plane Surface(30) = {31};
Curve Loop(32) = {8, 156, -108, -145};
Plane Surface(31) = {32};
Curve Loop(33) = {143, -105, -156, 5};
Plane Surface(32) = {33};


Physical Surface("Domain", 1) = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32};
Physical Curve("inlet", 11) = {175,176,177,178};
Physical Curve("outlet", 12) = {183,184,185,186};
Physical Curve("walls", 13) = {179,180,181,182,187,188,189,190};
Physical Curve("EmitterUp", 1) = {22, 21, 24, 23};
Physical Curve("EmitterDown", 2) = {32, 31, 33, 34};
Physical Curve("CollectorUp", 3) = {6, 7, 8, 5};
Physical Curve("CollectorDown", 4) = {10, 11, 12, 13};

Transfinite Curve {175, 176, 177, 178, 183, 184, 185, 186} = 2 Using Progression 1;
Transfinite Curve {190, 150, 187, 182, 179, 173, 160} = 4 Using Progression 1;
Transfinite Curve {150, 174, 157, 161} = 3 Using Progression 1;
Transfinite Curve {164, 168, 169, 154} = 5 Using Progression 1;
Transfinite Curve {151, 152, 153, 155, 158, 170, 171, 159, 163, 165, 166, 162} = 3 Using Progression 1;
Transfinite Curve {122, 121, 124, 123, 21, 24, 23, 22, 131, 134, 133, 132, 31, 34, 33, 32} = 3 Using Progression 1;
Transfinite Curve {105, 108, 107, 106, 5, 8, 7, 6, 112, 113, 110, 111, 13, 12, 11, 10} = 3 Using Progression 1;
Transfinite Curve {135, 136, 137, 138, 139, 140, 141, 142} = 5 Using Progression 2.9;
Transfinite Curve {143, 156, 145, 144, 146, 147, 148, 149} = 3 Using Progression 1;
Transfinite Curve {189, 188, 172, 167, 180, 181} = 3 Using Progression 1;

Transfinite Surface {20};
Transfinite Surface {19};
Transfinite Surface {18};
Transfinite Surface {17};
Transfinite Surface {21};
Transfinite Surface {24};
Transfinite Surface {23};
Transfinite Surface {22};
Transfinite Surface {29};
Transfinite Surface {32};
Transfinite Surface {31};
Transfinite Surface {30};
Transfinite Surface {28};
Transfinite Surface {27};
Transfinite Surface {26};
Transfinite Surface {25};

Mesh.Algorithm=1;
Mesh.RecombineAll=1;

Mesh 2;
