res1 = 0.01;
res2 = 0.005;
// Tube domain
Point(1) = {0, 0, 0, 0.04};
Point(2) = {0, 0.41, 0, 0.04};
Point(3) = {2.5, 0.41, 0, 0.1};
Point(4) = {2.5, 0, 0, 0.1};

// Circle
Point(5) = {0.2, 0.2, 0, 1};
Point(7) = {0.2, 0.15, 0, 0.004};
Point(8) = {0.2, 0.25, 0, 0.004};
Point(9) = {0.15, 0.2, 0, 0.004};
Point(10) = {0.2489897949, 0.21, 0.0, 0.01};
Point(11) = {0.2489897949, 0.19, 0.0, 0.01};

Point(31) = {0.175, 0.243301, 0.0, 0.002};
Point(30) = {0.155, 0.221794, 0.0, 0.002};
/*
Point(22) = {0.175, 0.243301, 0.0, 0.002};
Point(24) = {0.155, 0.221794, 0.0, 0.002};
Point(23) = {0.16, 0.23, 0.0, 0.002};
Point(25) = {0.165, 0.235707, 0.0, 0.002};
Point(26) = {0.152, 0.214, 0.0, 0.002};
*/
//Flag end points
Point(12) = {0.6, 0.19, 0, 0.008};
Point(13) = {0.6, 0.21, 0, 0.008};

//Extra point Flag
Point(22) = {0.32, 0.21, 0.0, 0.008};
Point(23) = {0.32, 0.19, 0.0, 0.008};

Point(14) = {0.43, 0.19, 0, 0.008};
Point(15) = {0.43, 0.21, 0, 0.008};

//Extra points boundary
Point(16) = {0.43, 0.0, 0, 0.05};
Point(17) = {0.43, 0.41, 0, 0.05};

Point(18) = {1.5, 0.0, 0, 0.1};
Point(19) = {1.5, 0.41, 0, 0.1};

Point(20) = {0.9, 0.0, 0, 0.1};
Point(21) = {0.9, 0.41, 0, 0.1};

//+
Line(1) = {2, 17};
//+
Line(2) = {17, 21};
//+
Line(3) = {21, 19};
//+
Line(4) = {19, 3};
//+
Line(5) = {3, 4};
//+
Line(6) = {4, 18};
//+
Line(7) = {18, 20};
//+
Line(8) = {20, 16};
//+
Line(9) = {16, 1};
//+
Line(10) = {1, 2};
//+
Line(11) = {11, 23};
//+
Line(12) = {23, 14};
//+
Line(13) = {14, 12};
//+
Line(14) = {12, 13};
//+
Line(15) = {13, 15};
//+
Line(16) = {15, 22};
//+
Line(17) = {22, 10};
//+
Circle(18) = {11, 5, 7};
//+
Circle(19) = {7, 5, 9};
//+
Circle(20) = {9, 5, 30};
//+
Circle(21) = {30, 5, 31};
//+
Circle(22) = {31, 5, 8};
//+
Circle(23) = {8, 5, 10};
//+
Line Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
//+
Line Loop(2) = {23, -17, -16, -15, -14, -13, -12, -11, 18, 19, 20, 21, 22};
//+
Circle(24) = {10, 5, 11};
//+
Plane Surface(1) = {1, 2};
//+
Line Loop(3) = {17, 24, 11, 12, 13, 14, 15, 16};
//+
Plane Surface(2) = {3};


Field[1] = BoundaryLayer;
Field[1].EdgesList = {18,19,20,21,22};
Field[1].hwall_n = 0.0005;
Field[1].ratio = 3.1;
Field[1].thickness = 0.003;
Field[1].Quads = 0;
BoundaryLayer Field = 1;
