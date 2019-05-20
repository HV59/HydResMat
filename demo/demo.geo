/*****************************************************************************
 *
 *  Gmsh demo case for a cylinder with a concave and a convex spherical end
 *
 *****************************************************************************/

/* Copyright (C) 2018-2019 Johannes Voss, Julian Jeggle, Raphael Wittkowski

    This file is part of HydResMat.

    HydResMat is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    HydResMat is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with HydResMat. If not, see <http://www.gnu.org/licenses/>.*/


// Defining the mesh size lcar1 at the outer surface of the cuboid-shaped simulation domain, the mesh size lcar2 at the particle surface and the edge length (2 times box_size) of the cuboid-shaped domain
//lcar1 = 20;
//lcar2 = 0.5;
lcar1 = 80;
lcar2 = 2;
box_size = 500;
// Resolution for good results
//lcar1 = 15;
//lcar2 = 0.1;

// Defining the corner points of the cuboid-shaped domain with edge length = 2 times box_size
Point(8) = {-box_size,box_size,-box_size,lcar1};
Point(9) = {box_size,box_size,-box_size,lcar1};
Point(10) = {-box_size,-box_size,box_size,lcar1};
Point(11) = {-box_size,box_size,box_size,lcar1};
Point(12) = {box_size,box_size,box_size,lcar1};
Point(13) = {box_size,-box_size,box_size,lcar1};
Point(14) = {box_size,-box_size,-box_size,lcar1};
Point(15) = {-box_size,-box_size,-box_size,lcar1};

// Defining the edges of the cuboid-shaped domain
Line(1) = {8,9};     Line(2) = {9,12};    Line(3) = {12,11};
Line(4) = {11,8};    Line(5) = {9,14};    Line(6) = {14,13};
Line(7) = {13,12};   Line(8) = {11,10};   Line(9) = {10,13};
Line(10) = {10, 15}; Line(11) = {14, 15}; Line(12) = {8, 15};

// Defining the faces of the cuboid-shaped domain
Line Loop(30) = {7, 3, 8, 9};	    Plane Surface(31) = {30};
Line Loop(26) = {-4, -3, -2, -1};   Plane Surface(27) = {26};
Line Loop(28) = {7, -2, 5, 6};      Plane Surface(29) = {28};
Line Loop(20) = {9, -6, 11, -10};   Plane Surface(21) = {20};
Line Loop(22) = {-4, 8, 10, -12};   Plane Surface(23) = {22};
Line Loop(24) = {-5, -1,  12, -11}; Plane Surface(25) = {24};

// If you want to have periodic faces of the cuboid-shaped domain, you need to uncomment the following three lines
// Periodic Surface 23{-8, 4, 12, -10} = 29{7, -2, 5, 6};
// Periodic Surface 21{-10, 9, -6, 11} = 27{-4, -3, -2, -1};
// Periodic Surface 25{-5, -1, 12, -11} = 31{7, 3, 8, 9};

// Declaring the outer surface of the cuboid-shaped domain
theloops[0] = 102;
Surface Loop(theloops[0]) = {21, 23, 25, 27, 29, 31};

// Declaring some parameters related to the particle
// x, y, z define a shift of the particle position
// sigma is the diameter of the particle
// L is the length of the cylindrical part of the particle
// h is the height of the spherical caps of the particle
// r is the radius of curvature of the caps of the particle
x = 0; y = 0; z = 0; sigma = 0.1; L = 5 * sigma; h = 0.048; r = (sigma * sigma / 4 + h * h) / (2 * h);

// Declaring points that are used to define the particle shape
p1 = newp; Point(p1) = {x-L/2, y, z, lcar2};
p2 = newp; Point(p2) = {x-L/2, y-sigma/2, z, lcar2};
p3 = newp; Point(p3) = {x-L/2, y+sigma/2, z, lcar2};
p4 = newp; Point(p4) = {x-L/2, y, z-sigma/2, lcar2};
p5 = newp; Point(p5) = {x-L/2, y, z+sigma/2, lcar2};
p6 = newp; Point(p6) = {x+L/2, y, z, lcar2};
p7 = newp; Point(p7) = {x+L/2, y-sigma/2, z, lcar2};
p8 = newp; Point(p8) = {x+L/2, y+sigma/2, z, lcar2};
p9 = newp; Point(p9) = {x+L/2, y, z-sigma/2, lcar2};
p10 = newp; Point(p10) = {x+L/2, y, z+sigma/2, lcar2};
p11 = newp; Point(p11) = {x-L/2+h, y, z, lcar2};
p12 = newp; Point(p12) = {x+L/2+h, y, z, lcar2};
p13 = newp; Point(p13) = {x-L/2+h-r, y, z, lcar2};
p14 = newp; Point(p14) = {x+L/2+h-r, y, z, lcar2};

// Defining lines and circles of the particle shape
l1 = newreg; Line(l1) = {p2, p7};
l2 = newreg; Line(l2) = {p3, p8};
l3 = newreg; Line(l3) = {p4, p9};
l4 = newreg; Line(l4) = {p5, p10};

c1 = newreg; Circle(c1) = {p2, p1, p4};
c2 = newreg; Circle(c2) = {p4, p1, p3};
c3 = newreg; Circle(c3) = {p3, p1, p5};
c4 = newreg; Circle(c4) = {p5, p1, p2};
c5 = newreg; Circle(c5) = {p7, p6, p9};
c6 = newreg; Circle(c6) = {p9, p6, p8};
c7 = newreg; Circle(c7) = {p8, p6, p10};
c8 = newreg; Circle(c8) = {p10, p6, p7};
c9 = newreg; Circle(c9) = {p2, p13, p11};
c10 = newreg; Circle(c10) = {p11, p13, p3};
c11 = newreg; Circle(c11) = {p4, p13, p11};
c12 = newreg; Circle(c12) = {p11, p13, p5};
c13 = newreg; Circle(c13) = {p7, p14, p12};
c14 = newreg; Circle(c14) = {p12, p14, p8};
c15 = newreg; Circle(c15) = {p9, p14, p12};
c16 = newreg; Circle(c16) = {p12, p14, p10};

// Defining parts of the surface of the particle (the order of the lines is important, because they have to create a closed curve)
ll1 = newreg; Line Loop(ll1) = {-c1, l1, c5, -l3};	Surface(newreg) = {ll1};
ll2 = newreg; Line Loop(ll2) = {-c2, l3, c6, -l2};	Surface(newreg) = {ll2};
ll3 = newreg; Line Loop(ll3) = {-c3, l2, c7, -l4};	Surface(newreg) = {ll3};
ll4 = newreg; Line Loop(ll4) = {-c4, l4, c8, -l1};	Surface(newreg) = {ll4};
ll5 = newreg; Line Loop(ll5) = {c8, c13, c16};		Surface(newreg) = {ll5};
ll6 = newreg; Line Loop(ll6) = {c7, -c16, c14};		Surface(newreg) = {ll6};
ll7 = newreg; Line Loop(ll7) = {c6, -c14, -c15};	Surface(newreg) = {ll7};
ll8 = newreg; Line Loop(ll8) = {c5, c15, -c13};		Surface(newreg) = {ll8};
ll9 = newreg; Line Loop(ll9) = {c9, c12, c4};		Surface(newreg) = {ll9};
ll10 = newreg; Line Loop(ll10) = {c1, c11, -c9};	Surface(newreg) = {ll10};
ll11 = newreg; Line Loop(ll11) = {c2, -c10, -c11};	Surface(newreg) = {ll11};
ll12 = newreg; Line Loop(ll12) = {c3, -c12, c10};	Surface(newreg) = {ll12};

// Declaring the whole surface of the particle as sum of the surface parts
theloops[1] = 100;
Surface Loop(theloops[1]) = {ll1+1, ll2+1, ll3+1, ll4+1, ll5+1, ll6+1, ll7+1, ll8+1, ll9+1, ll10+1, ll11+1, ll12+1};

// Defining the volume that has to be meshed, i.e., the domain outside of the particle
Volume(186) = {theloops[]};

// Defining some numbers including those that are later used in FEniCS to describe the surface parts in the boundary conditions or surface integrals
Physical Volume (10) = 186;
// The surface of the particle gets the number 32
Physical Surface (32) = {ll1+1, ll2+1, ll3+1, ll4+1, ll5+1, ll6+1, ll7+1, ll8+1, ll9+1, ll10+1, ll11+1, ll12+1};
// The faces of the cuboid-shaped domain get the numbers 21, 23, 25, 27, 29 and 31
Physical Surface(21) = {21};
Physical Surface(23) = {23};
Physical Surface(25) = {25};
Physical Surface(27) = {27};
Physical Surface(29) = {29};
Physical Surface(31) = {31};
