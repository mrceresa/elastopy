cl__1 = 1;
Point(1) = {0, 0, 0, 1};
Point(2) = {200, 0, 0, 1};
Point(3) = {200, 70, 0, 1};
Point(4) = {147.9, 70, 0, 1};
Point(5) = {77.28700000000001, 169.06, 0, 1};
Point(6) = {73.477, 176.68, 0, 5};
Point(7) = {73.477, 191.92, 0, 1};
Point(8) = {63.727, 191.92, 0, 1};
Point(9) = {63.727, 169.06, 0, 1};
Point(10) = {52.1, 70, 0, 0.2};
Point(11) = {0, 70, 0, 1};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 9};
Line(9) = {9, 10};
Line(10) = {10, 11};
Line(11) = {11, 1};
Line(12) = {4, 10};
Line Loop(14) = {4, 5, 6, 7, 8, 9, -12};
Plane Surface(14) = {14};
Line Loop(16) = {10, 11, 1, 2, 3, 12};
Plane Surface(16) = {16};
Physical Line(17) = {1};
Physical Line(18) = {2};
Physical Line(19) = {3};
Physical Line(20) = {4};
Physical Line(21) = {5};
Physical Line(22) = {6};
Physical Line(23) = {7};
Physical Line(24) = {8};
Physical Line(25) = {9};
Physical Line(26) = {10};
Physical Line(27) = {11};
Physical Surface(28) = {16};
Physical Surface(29) = {14};
