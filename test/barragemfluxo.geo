// Gmsh project created on Thu Aug 27 09:18:33 2015
Point(1) = {-1.9, 0.4, -0, 1.0};
Point(2) = {-1.9, -1, 0, 1.0};
Point(3) = {2.1, -1.1, 0, 1.0};
Point(4) = {2, 0.4, -0, 1.0};
Point(5) = {1.5, 0.4, -0, 1.0};
Point(6) = {0.7, 1.1, -0, 1.0};
Point(7) = {-0.3, 1.1, -0, 1.0};
Point(8) = {-0.9, 0.5, -0, 1.0};
Point(9) = {-1, 0.3, -0, 1.0};
Point(10) = {0.4, 1.1, -0, 1.0};
Point(11) = {0.3, 1.1, -0, 1.0};
Point(12) = {0.3, 0.3, -0, 1.0};
Point(13) = {1.3, 0.4, -0, 1.0};
Point(14) = {1.3, 0.5, -0, 1.0};
Point(15) = {0.4, 0.5, -0, 1.0};
Line(1) = {1, 2};
Line(2) = {3, 3};
Line(3) = {2, 3};
Line(4) = {3, 4};
Line(5) = {4, 5};
Line(6) = {5, 14};
Line(7) = {14, 6};
Line(8) = {6, 10};
Line(9) = {10, 11};
Line(10) = {11, 7};
Line(11) = {8, 8};
Line(12) = {7, 8};
Line(13) = {8, 9};
Line(14) = {9, 1};
Line(15) = {10, 15};
Line(16) = {15, 14};
Line(17) = {14, 13};
Line(18) = {13, 12};
Line(19) = {12, 11};
Line(20) = {5, 13};
Line(21) = {13, 12};
Line(22) = {12, 9};
Line Loop(23) = {14, 1, 3, 4, 5, 20, 18, 22};
Plane Surface(24) = {23};
Line Loop(25) = {7, 8, 15, 16};
Plane Surface(26) = {25};
Line Loop(27) = {6, 17, -20};
Plane Surface(28) = {27};
Line Loop(29) = {22, -13, -12, -10, -19};
Plane Surface(30) = {29};
Line Loop(31) = {16, 17, 18, 19, -9, 15};
Plane Surface(32) = {31};
Physical Line(33) = {1};
Physical Line(34) = {3};
Physical Line(35) = {4};
Physical Line(36) = {5};
Physical Line(37) = {6};
Physical Line(38) = {7};
Physical Line(39) = {8};
Physical Line(40) = {9};
Physical Line(41) = {10};
Physical Line(42) = {12};
Physical Line(43) = {13};
Physical Line(44) = {14};
Physical Surface(45) = {24};
Physical Surface(46) = {30};
Physical Surface(47) = {26};
Physical Surface(48) = {28};
Physical Surface(49) = {32};
