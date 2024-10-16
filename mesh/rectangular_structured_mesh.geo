
Point(1) = {-2, 0, 0, 1.0};
Point(2) = {-2, 0.1, 0, 1.0};
Point(3) = {-0.1, 0.1, 0, 1.0};
Point(4) = {0.1, 0.1, 0, 1.0};
Point(5) = {0.2, 0.1, 0, 1.0};
Point(6) = {0.2, 0, 0, 1.0};
Point(7) = {0.1, 0, 0, 1.0};
Point(8) = {0.00040, 0, 0, 1.0};
Point(9) = {0.00040, 0.00040, 0, 1.0};
Point(10) = {-0.00040, 0.00040, 0, 1.0};
Point(11) = {-0.00040, 0, 0, 1.0};
Point(12) = {-0.1, 0, 0, 1.0};

Point(13) = {0.004, 0, 0, 1.0};
Point(14) = {0.004, 0.004, 0, 1.0};
Point(15) = {-0.004, 0.004, 0, 1.0};
Point(16) = {-0.004, 0, 0, 1.0};

Line(2) = {1, 2};
Line(301) = {2, 3};
Line(303) = {4, 5};
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
Line(28) = {15, 3};
Line(29) = {14, 4};

Curve Loop(1) = {9, 2, 301, 8};
Plane Surface(1) = {1};
Curve Loop(2) = {303,4,5,6};
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
Curve Loop(8) = {8, -18, -17, -16, 28};
Plane Surface(8) = {8};
Curve Loop(9) = {28, -7, -29, 13, 14, 15};
Plane Surface(9) = {9};
Curve Loop(10) = {29, -6, 10, 11, 12};
Plane Surface(10) = {10};

Physical Curve("Emitter", 1) = {20,21,22};
Physical Curve("Collector", 2) = {4};
Physical Curve("Others", 9) =  {9,2,301,7,303,5,10,19,23,18};

Physical Surface(100) = {1,2,3,4,5,6,7,8,9,10};

Transfinite Curve {11,17,20,22} = 3 Using Progression 1;
Transfinite Curve {14,21} = 5 Using Progression 1;
Transfinite Curve {12,13,15,16,19,23,24,25,26,27} = 19 Using Progression 1;
Transfinite Surface {3,4,5,6,7};

Transfinite Curve {6,8} = 21 Using Progression 1;
Transfinite Curve {7} = 41 Using Progression 1;
Transfinite Curve {18,28,29,-10} = 41 Using Progression 1.1;
Transfinite Surface {8} = {16,12,3,15};
Transfinite Surface {9} = {15,3,4,14};
Transfinite Surface {10} = {13,14,4,7};

Transfinite Curve {2,4} = 21 Using Progression 1;
Transfinite Curve {9,-301} = 61 Using Progression 1.01;
Transfinite Curve {5,303} = 16 Using Progression 1;
Transfinite Surface {1};
Transfinite Surface {2};

Characteristic Length {1, 2, 5, 6, 7, 4, 3, 12} = 0.01;
Characteristic Length {12, 3, 4, 7, 13, 17, 14, 19, 20, 15, 18, 16} = 0.0005;

Recombine Surface {1};
Recombine Surface {2};
Recombine Surface {3};
Recombine Surface {4};
Recombine Surface {5};
Recombine Surface {6};
Recombine Surface {7};
Recombine Surface {8};
Recombine Surface {9};
Recombine Surface {10};

