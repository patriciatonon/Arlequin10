// Gmsh project created on Fri Apr 28 14:41:38 2023
SetFactory("OpenCASCADE");
//+
Point(1) = {1, 0, 0, 1.0};
//+
Point(2) = {1, 0.5, 0, 1.0};
//+
Point(3) = {0, 0.5, 0, 1.0};
//+
Point(4) = {0, 0, 0, 1.0};
//+
Line(1) = {4, 1};
//+
Line(2) = {1, 2};
//+
Line(3) = {2, 3};
//+
Line(4) = {3, 4};
//+
Line Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Physical Line("down") = {1};
//+
Physical Line("leftandright") = {4, 2};
//+
Physical Line("up") = {3};
//+
Transfinite Line {1} = 11 Using Progression 1;
//+
Transfinite Line {2} = 6 Using Progression 1;
//+
Transfinite Line {4} = 6 Using Progression 1;
//+
Transfinite Line {3} = 11 Using Progression 1;
//+
Transfinite Surface {1} Right;
//+
Physical Surface("surf") = {1};
