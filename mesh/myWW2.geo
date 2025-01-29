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
Point(173) = {-0.08, 0.019, 0, 1.0};
Point(174) = {-0.08, -0.019, 0, 1.0};
Point(175) = {0.13, 0.019, 0, 1.0};
Point(176) = {0.13, -0.019, 0, 1.0};
Point(177) = {-0.08, 0.041, 0, 1.0};
Point(178) = {-0.08, -0.041, 0, 1.0};
Point(179) = {0.13, 0.041, 0, 1.0};
Point(180) = {0.13, -0.041, 0, 1.0};

Line(150) = {7, 107};
Line(151) = {122, 109};
Line(152) = {123, 108};
Line(153) = {125, 110};
Line(154) = {110, 152};
Line(155) = {152, 133};
Line(156) = {133, 125};
Line(157) = {154, 135};
Line(158) = {132, 153};
Line(159) = {151, 159};
Line(160) = {134, 156};
Line(161) = {124, 158};
Line(162) = {123, 172};
Line(163) = {108, 171};
Line(164) = {172, 171};
Line(165) = {135, 162};
Line(166) = {162, 170};
Line(167) = {170, 154};
Line(168) = {125, 173};
Line(169) = {173, 158};
Line(170) = {158, 177};
Line(171) = {177, 123};
Line(172) = {110, 175};
Line(173) = {175, 160};
Line(174) = {160, 107};
Line(175) = {160, 179};
Line(176) = {179, 108};
Line(177) = {179, 41};
Line(178) = {41, 171};
Line(179) = {172, 12};
Line(180) = {12, 177};
Line(181) = {173, 174};
Line(182) = {174, 133};
Line(183) = {174, 156};
Line(184) = {156, 178};
Line(185) = {178, 11};
Line(186) = {11, 162};
Line(187) = {135, 178};
Line(188) = {152, 176};
Line(189) = {176, 159};
Line(190) = {159, 180};
Line(191) = {180, 154};
Line(192) = {180, 43};
Line(193) = {43, 170};
Line(194) = {176, 175};


Curve Loop(1) = {122, -138, -22, 135};
Plane Surface(1) = {1};
Curve Loop(2) = {138, 123, -137, -23};
Plane Surface(2) = {2};
Curve Loop(3) = {24, 136, -124, -137};
Plane Surface(3) = {3};
Curve Loop(4) = {135, -121, -136, 21};
Plane Surface(4) = {4};
Curve Loop(5) = {132, -142, -32, 139};
Plane Surface(5) = {5};
Curve Loop(6) = {142, 133, -141, -33};
Plane Surface(6) = {6};
Curve Loop(7) = {141, 134, -140, -34};
Plane Surface(7) = {7};
Curve Loop(8) = {139, -131, -140, 31};
Plane Surface(8) = {8};
Curve Loop(9) = {111, 146, 10, 149};
Plane Surface(9) = {9};
Curve Loop(10) = {149, -110, -148, -11};
Plane Surface(10) = {10};
Curve Loop(11) = {148, -113, -147, -12};
Plane Surface(11) = {11};
Curve Loop(12) = {13, -146, 112, -147};
Plane Surface(12) = {12};
Curve Loop(13) = {106, -144, -6, 143};
Plane Surface(13) = {13};
Curve Loop(14) = {144, 107, -145, -7};
Plane Surface(14) = {14};
Curve Loop(15) = {145, 108, -150, -8};
Plane Surface(15) = {15};
Curve Loop(16) = {143, -105, -150, 5};
Plane Surface(16) = {16};
Curve Loop(17) = {179, 180, 171, 162};
Plane Surface(17) = {17};
Curve Loop(18) = {171, 122, 161, 170};
Plane Surface(18) = {18};
Curve Loop(19) = {161, -169, -168, -123};
Plane Surface(19) = {19};
Curve Loop(20) = {168, 181, 182, 156};
Plane Surface(20) = {20};
Curve Loop(21) = {182, 132, 160, -183};
Plane Surface(21) = {21};
Curve Loop(22) = {160, 184, -187, -133};
Plane Surface(22) = {22};
Curve Loop(23) = {187, 185, 186, -165};
Plane Surface(23) = {23};
Curve Loop(24) = {165, 166, 167, 157};
Plane Surface(24) = {24};
Curve Loop(25) = {134, 158, -110, 157};
Plane Surface(25) = {25};
Curve Loop(26) = {131, -155, -111, -158};
Plane Surface(26) = {26};
Curve Loop(27) = {156, 153, 154, 155};
Plane Surface(27) = {27};
Curve Loop(28) = {151, 107, -153, 124};
Plane Surface(28) = {28};
Curve Loop(29) = {152, 106, -151, 121};
Plane Surface(29) = {29};
Curve Loop(30) = {164, -163, -152, 162};
Plane Surface(30) = {30};
Curve Loop(31) = {178, -163, -176, 177};
Plane Surface(31) = {31};
Curve Loop(32) = {176, -105, -174, 175};
Plane Surface(32) = {32};
Curve Loop(33) = {174, -108, 172, 173};
Plane Surface(33) = {33};
Curve Loop(34) = {172, -194, -188, -154};
Plane Surface(34) = {34};
Curve Loop(35) = {188, 189, -159, -112};
Plane Surface(35) = {35};
Curve Loop(36) = {159, 190, 191, -113};
Plane Surface(36) = {36};
Curve Loop(37) = {191, -167, -193, -192};
Plane Surface(37) = {37};

Physical Surface("Domain", 1) = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37};
Physical Curve("inlet", 11) = {180,170,169,181,183,184,185};
Physical Curve("outlet", 12) = {177,175,173,194,189,190,192};
Physical Curve("walls", 13) = {179,164,178,186,166,193};
Physical Curve("EmitterUp", 1) = {22, 21, 24, 23};
Physical Curve("EmitterDown", 2) = {32, 31, 33, 34};
Physical Curve("CollectorUp", 3) = {6, 7, 8, 5};
Physical Curve("CollectorDown", 4) = {10, 11, 12, 13};


Transfinite Curve {122, 121, 124, 123, 21, 24, 23, 22, 131, 134, 133, 132, 31, 34, 33, 32,  170,169,183,184,175,173,189,190} = 2 Using Progression 1;
Transfinite Curve {105, 108, 107, 106, 5, 8, 7, 6, 112, 113, 110, 111, 13, 12, 11, 10} = 2 Using Progression 1;
Transfinite Curve {135, 136, 137, 138, 139, 140, 141, 142} = 5 Using Progression 2.9;
Transfinite Curve {143, 156, 145, 144, 146, 147, 148, 149} = 3 Using Progression 1;
Transfinite Surface {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37};

Transfinite Curve {180, 185, 165, 167, 192, 177, 163, 162} = 2 Using Progression 1;
Transfinite Curve {181, 156, 154, 194} = 3 Using Progression 1;
Transfinite Curve {164, 152, 153, 155, 157, 166} = 3 Using Progression 1;
Transfinite Curve {151, 158} = 3 Using Progression 1;
Transfinite Curve {179, 171, 168, 182, 187, 186, 178, 176, 172, 188, 191, 193} = 2 Using Progression 1;
Transfinite Curve {161, 160, 159, 174} = 2 Using Progression 1;


Mesh.Algorithm=1;
Mesh.RecombineAll=1;

Mesh 2;
