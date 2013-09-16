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

namespace DSC {
    
    class Tetralizer {
        
    private:        
        static void tetralize_cube1(int i, int j, int k, int Ni, int Nj, int Nk, std::vector<int>& tets);
        
        static void tetralize_cube2(int i, int j, int k, int Ni, int Nj, int Nk, std::vector<int>& tets);
        
        static void create_tets(int Ni, int Nj, int Nk, std::vector<int>& tets);
        
        static void create_points(const vec3& size, real avg_edge_length, int Ni, int Nj, int Nk, std::vector<real>& points);
        
    public:
        static void tetralize(const vec3& size, real avg_edge_length, std::vector<real>& points, std::vector<int>& tets)
        {
            int Ni = std::ceil(size[0]/avg_edge_length) + 1;
            int Nj = std::ceil(size[1]/avg_edge_length) + 1;
            int Nk = std::ceil(size[2]/avg_edge_length) + 1;
            
            create_points(size, avg_edge_length, Ni, Nj, Nk, points);
            create_tets(Ni, Nj, Nk, tets);
        }
    };
    
}
