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

class ObjectGenerator {    
    
public:
    
    static void create_sphere(const std::vector<vec3>& points, const std::vector<int>& tets, const vec3& center, const real& radius, int label, std::vector<int>& tet_labels);
    
    static void create_cube(const std::vector<vec3>& points, const std::vector<int>& tets, const vec3& origin, const vec3& size, int label, std::vector<int>& tet_labels);
};