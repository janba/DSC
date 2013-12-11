//
//  Deformabel Simplicial Complex (DSC) method
//  Copyright (C) 2013  Technical University of Denmark
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  See licence.txt for a copy of the GNU General Public License.

#pragma once

#include "util.h"


inline void test_distance_triangle_triangle()
{
    std::cout << "Testing utility functions:";
    real d = distance_triangle_triangle<real>(vec3(0.), vec3(0., 1., 0.), vec3(1., 0., 0.), vec3(1., 1., -1.), vec3(1.,1.,2.), vec3(4., 2., 4.));
    assert(abs(d - sqrt(2.)/2) < EPSILON);
    d = distance_triangle_triangle<real>(vec3(0.), vec3(0., 1., 0.), vec3(1., 0., 0.), vec3(1., 1., 0.), vec3(4.,1.,2.), vec3(4., 2., 4.));
    assert(abs(d - sqrt(2.)/2) < EPSILON);
    std::cout << " PASSED" << std::endl;
}

inline void simplex_set_test()
{
    std::cout << "Testing simplex set class: ";
    SimplexSet<int> A = {1,3,9,4};
    SimplexSet<int> B = {1,7,5,3,10};
    
    SimplexSet<int> U = {1,3,9,4,7,5,10};
    assert((A+B) == U);
    
    SimplexSet<int> C = {9,4};
    assert((A-B) == C);
    
    SimplexSet<int> I = {1,3};
    assert((A&B) == I);
    
    A -= 3;
    A += 9;
    A += 11;
    SimplexSet<int> E = {1,9,4,11};
    assert(A == E);
    
    std::cout << "PASSED" << std::endl;
}