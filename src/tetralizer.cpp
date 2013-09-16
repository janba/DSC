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

#include "tetralizer.h"

namespace DSC {
    
    int get_index(int i, int j, int k, int Ni, int Nj, int Nk)
    {
        return i + j*Ni + k*Ni*Nj;
    }
    
    void Tetralizer::tetralize_cube1(int i, int j, int k, int Ni, int Nj, int Nk, std::vector<int>& tets)
    {
        // First tetrahedron:
        tets.push_back(get_index(i, j, k, Ni, Nj, Nk));
        tets.push_back(get_index(i, j, k+1, Ni, Nj, Nk));
        tets.push_back(get_index(i+1, j, k+1, Ni, Nj, Nk));
        tets.push_back(get_index(i, j+1, k+1, Ni, Nj, Nk));
        
        // Second tetrahedron:
        tets.push_back(get_index(i, j, k, Ni, Nj, Nk));
        tets.push_back(get_index(i, j+1, k+1, Ni, Nj, Nk));
        tets.push_back(get_index(i+1, j+1, k, Ni, Nj, Nk));
        tets.push_back(get_index(i, j+1, k, Ni, Nj, Nk));
        
        // Third tetrahedron:
        tets.push_back(get_index(i, j, k, Ni, Nj, Nk));
        tets.push_back(get_index(i+1, j, k, Ni, Nj, Nk));
        tets.push_back(get_index(i+1, j+1, k, Ni, Nj, Nk));
        tets.push_back(get_index(i+1, j, k+1, Ni, Nj, Nk));
        
        // Fourth tetrahedron:
        tets.push_back(get_index(i+1, j+1, k+1, Ni, Nj, Nk));
        tets.push_back(get_index(i+1, j, k+1, Ni, Nj, Nk));
        tets.push_back(get_index(i, j+1, k+1, Ni, Nj, Nk));
        tets.push_back(get_index(i+1, j+1, k, Ni, Nj, Nk));
        
        // Fifth tetrahedron:
        tets.push_back(get_index(i, j, k, Ni, Nj, Nk));
        tets.push_back(get_index(i+1, j+1, k, Ni, Nj, Nk));
        tets.push_back(get_index(i, j+1, k+1, Ni, Nj, Nk));
        tets.push_back(get_index(i+1, j, k+1, Ni, Nj, Nk));
    }
    
    void Tetralizer::tetralize_cube2(int i, int j, int k, int Ni, int Nj, int Nk, std::vector<int>& tets)
    {
        // First tetrahedron:
        tets.push_back(get_index(i, j, k, Ni, Nj, Nk));
        tets.push_back(get_index(i, j, k+1, Ni, Nj, Nk));
        tets.push_back(get_index(i+1, j, k, Ni, Nj, Nk));
        tets.push_back(get_index(i, j+1, k, Ni, Nj, Nk));
        
        // Second tetrahedron:
        tets.push_back(get_index(i, j+1, k, Ni, Nj, Nk));
        tets.push_back(get_index(i+1, j+1, k, Ni, Nj, Nk));
        tets.push_back(get_index(i+1, j+1, k+1, Ni, Nj, Nk));
        tets.push_back(get_index(i+1, j, k, Ni, Nj, Nk));
        
        // Third tetrahedron:
        tets.push_back(get_index(i, j+1, k+1, Ni, Nj, Nk));
        tets.push_back(get_index(i, j+1, k, Ni, Nj, Nk));
        tets.push_back(get_index(i+1, j+1, k+1, Ni, Nj, Nk));
        tets.push_back(get_index(i, j, k+1, Ni, Nj, Nk));
        
        // Fourth tetrahedron:
        tets.push_back(get_index(i+1, j, k+1, Ni, Nj, Nk));
        tets.push_back(get_index(i+1, j, k, Ni, Nj, Nk));
        tets.push_back(get_index(i, j, k+1, Ni, Nj, Nk));
        tets.push_back(get_index(i+1, j+1, k+1, Ni, Nj, Nk));
        
        // Fifth tetrahedron:
        tets.push_back(get_index(i, j, k+1, Ni, Nj, Nk));
        tets.push_back(get_index(i+1, j, k, Ni, Nj, Nk));
        tets.push_back(get_index(i, j+1, k, Ni, Nj, Nk));
        tets.push_back(get_index(i+1, j+1, k+1, Ni, Nj, Nk));
    }
    
    void Tetralizer::create_tets(int Ni, int Nj, int Nk, std::vector<int>& tets)
    {
        for (int k = 0; k < Nk-1; k++) {
            for (int j = 0; j < Nj-1; j++) {
                for (int i = 0; i < Ni-1; i++)
                {
                    if((i + j + k)%2 == 0)
                    {
                        tetralize_cube1(i, j, k, Ni, Nj, Nk, tets);
                    }
                    else {
                        tetralize_cube2(i, j, k, Ni, Nj, Nk, tets);
                    }
                }
            }
        }
    }
    
    void Tetralizer::create_points(const vec3& size, real avg_edge_length, int Ni, int Nj, int Nk, std::vector<real>& points)
    {
        for (int k = 0; k < Nk; k++) {
            for (int j = 0; j < Nj; j++) {
                for (int i = 0; i < Ni; i++)
                {
                    points.push_back(std::min(i*avg_edge_length, size[0]) - size[0]/2.);
                    points.push_back(std::min(j*avg_edge_length, size[1]) - size[1]/2.);
                    points.push_back(std::min(k*avg_edge_length, size[2]) - size[2]/2.);
                }
            }
        }
    }

}
