
Point(1) = {-2, 0, 0, 1.0};
Point(2) = {-2, 0.1, 0, 1.0};
Point(3) = {-0.04, 0.04, 0, 1.0};
Point(4) = {0.04, 0.04, 0, 1.0};
Point(5) = {0.2, 0.1, 0, 1.0};
Point(6) = {0.2, 0, 0, 1.0};
Point(7) = {0.04, 0, 0, 1.0};
Point(8) = {0.00040, 0, 0, 1.0};
Point(9) = {0.00040, 0.00040, 0, 1.0};
Point(10) = {-0.00040, 0.00040, 0, 1.0};
Point(11) = {-0.00040, 0, 0, 1.0};
Point(12) = {-0.04, 0, 0, 1.0};

Point(13) = {0.004, 0, 0, 1.0};
Point(14) = {0.004, 0.004, 0, 1.0};
Point(15) = {-0.004, 0.004, 0, 1.0};
Point(16) = {-0.004, 0, 0, 1.0};

Line(2) = {1, 2};
Line(3) = {2, 5};
Line(4) = {5, 6};
Line(5) = {6, 7};
Line(6) = {7, 4};
Line(7) = {4, 3};
Line(8) = {3, 12};
Line(9) = {12, 1};
Line(10) = {7, 13};
Point(17) = {0.004, 0.0004, 0, 1.0};
Point(18) = {-0.004, 0.0004, 0, 1.0};
Point(19) = {0.0004, 0.004, 0, 1.0};
Point(20) = {-0.0004, 0.004, 0, 1.0};
Line(11) = {13, 17};
Line(12) = {17, 14};
Line(13) = {14, 19};
Line(14) = {19, 20};
Line(15) = {20, 15};
Line(16) = {15, 18};
Line(17) = {18, 16};
Line(18) = {16, 12};
Line(19) = {13, 8};
Line(20) = {8, 9};
Line(21) = {9, 10};
Line(22) = {10, 11};
Line(23) = {11, 16};
Line(24) = {17, 9};
Line(25) = {19, 9};
Line(26) = {20, 10};
Line(27) = {18, 10};


Physical Curve("Emitter", 1) = {20,21,22};
Physical Curve("Collector", 2) = {4};
Physical Curve("Others", 9) =  {};


Curve Loop(1) = {9, 2, 3, 4, 5, 6, 7, 8};
Plane Surface(1) = {1};
Curve Loop(2) = {18, -8, -7, -6, 10, 11, 12, 13, 14, 15, 16, 17};
Plane Surface(2) = {2};
Curve Loop(3) = {17, -23, -22, -27};
Plane Surface(3) = {3};
Curve Loop(4) = {27, -26, 15, 16};
Plane Surface(4) = {4};
Curve Loop(5) = {26, -21, -25, 14};
Plane Surface(5) = {5};
Curve Loop(6) = {25, -24, 12, 13};
Plane Surface(6) = {6};
Curve Loop(7) = {20, -24, -11, 19};
Plane Surface(7) = {7};

Transfinite Curve {17} = 5 Using Progression 1;
Transfinite Curve {22} = 5 Using Progression 1;
Transfinite Curve {20} = 5 Using Progression 1;
Transfinite Curve {11} = 5 Using Progression 1;
Transfinite Curve {21} = 9 Using Progression 1;
Transfinite Curve {14} = 9 Using Progression 1;
Transfinite Curve {16} = 37 Using Progression 1;
Transfinite Curve {15} = 37 Using Progression 1;
Transfinite Curve {13} = 37 Using Progression 1;
Transfinite Curve {12} = 37 Using Progression 1;
Transfinite Curve {26} = 37 Using Progression 1;
Transfinite Curve {25} = 37 Using Progression 1;
Transfinite Curve {27} = 37 Using Progression 1;
Transfinite Curve {24} = 37 Using Progression 1;
Transfinite Curve {19} = 37 Using Progression 1;
Transfinite Curve {23} = 37 Using Progression 1;
Transfinite Surface {4};
Transfinite Surface {5};
Transfinite Surface {6};
Transfinite Surface {3};
Transfinite Surface {7};

Characteristic Length {1, 2, 5, 6, 7, 4, 3, 12} = 0.01;
Characteristic Length {12, 3, 4, 7, 13, 17, 14, 19, 20, 15, 18, 16} = 0.0005;

Recombine Surface {1};
Recombine Surface {2};
Recombine Surface {3};
Recombine Surface {4};
Recombine Surface {5};
Recombine Surface {6};
Recombine Surface {7};

