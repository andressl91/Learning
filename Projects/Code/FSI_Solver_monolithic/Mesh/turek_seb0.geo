// Gmsh project created on Tue May 16 10:57:48 2017
Point(1) = {0.0, 0, 0, 0.5};
Point(2) = {0.0, 0.41, 0, 0.5};
Point(3) = {2.5, 0.41, 0, 0.5};
Point(4) = {2.5, 0, 0, 0.5};
Point(5) = {0.2, 0.2, 0, 0.5};
Point(6) = {0.2, 0.05, 0, 0.5};
Delete {
  Point{6};
}
Point(6) = {0.2, 0.25, 0, 0.5};
Point(7) = {0.2, 0.15, 0, 0.5};
Point(8) = {0.25, 0.2, 0, 0.5};
Point(9) = {0.15, 0.2, 0, 0.5};