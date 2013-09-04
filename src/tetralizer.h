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

#ifndef __DSC__tetralizer__
#define __DSC__tetralizer__

#include "util.h"

template<typename MT>
class Tetralizer {
    
    typedef typename MT::real_type      T;
    
    T DEPTH;
    T WIDTH;
    T HEIGHT;
    T AVG_EDGE_LENGTH;
    
    int Ni, Nj, Nk;
    
public:
    Tetralizer(T width, T height, T depth, T avg_edge_length): WIDTH(width), HEIGHT(height), DEPTH(depth), AVG_EDGE_LENGTH(avg_edge_length)
    {
        Ni = std::ceil(WIDTH/AVG_EDGE_LENGTH) + 1;
        Nj = std::ceil(HEIGHT/AVG_EDGE_LENGTH) + 1;
        Nk = std::ceil(DEPTH/AVG_EDGE_LENGTH) + 1;
    }
    
    int get_index(int i, int j, int k)
    {
        return i + j*Ni + k*Ni*Nj;
    }
    
    void tetralize_cube(int i, int j, int k, std::vector<int>& tets)
    {        
        // First tetrahedron:
        tets.push_back(get_index(i, j, k));
        tets.push_back(get_index(i, j, k+1));
        tets.push_back(get_index(i+1, j, k+1));
        tets.push_back(get_index(i, j+1, k+1));
        
        // Second tetrahedron:
        tets.push_back(get_index(i, j, k));
        tets.push_back(get_index(i, j+1, k+1));
        tets.push_back(get_index(i+1, j+1, k));
        tets.push_back(get_index(i, j+1, k));
        
        // Third tetrahedron:
        tets.push_back(get_index(i, j, k));
        tets.push_back(get_index(i+1, j, k));
        tets.push_back(get_index(i+1, j+1, k));
        tets.push_back(get_index(i+1, j, k+1));
        
        // Fourth tetrahedron:
        tets.push_back(get_index(i+1, j+1, k+1));
        tets.push_back(get_index(i+1, j, k+1));
        tets.push_back(get_index(i, j+1, k+1));
        tets.push_back(get_index(i+1, j+1, k));
        
        // Fifth tetrahedron:
        tets.push_back(get_index(i, j, k));
        tets.push_back(get_index(i+1, j+1, k));
        tets.push_back(get_index(i, j+1, k+1));
        tets.push_back(get_index(i+1, j, k+1));
    }
    
    void create_tets(std::vector<int>& tets)
    {
        for (int k = 0; k < Nk-1; k++) {
            for (int j = 0; j < Nj-1; j++) {
                for (int i = 0; i < Ni-1; i++)
                {
                    tetralize_cube(i, j, k, tets);
                }
            }
        }
    }
    
    void create_points(std::vector<T>& points)
    {        
        for (int k = 0; k < Nk; k++) {
            for (int j = 0; j < Nj; j++) {
                for (int i = 0; i < Ni; i++)
                {
                    points.push_back(std::min(i*AVG_EDGE_LENGTH, WIDTH) - WIDTH/2.);
                    points.push_back(std::min(j*AVG_EDGE_LENGTH, HEIGHT) - HEIGHT/2.);
                    points.push_back(std::min(k*AVG_EDGE_LENGTH, DEPTH) - DEPTH/2.);
                }
            }
        }
    }
    
    void label_tets(std::vector<int>& tet_labels)
    {
        for (int k = 0; k < Nk-1; k++) {
            for (int j = 0; j < Nj-1; j++) {
                for (int i = 0; i < Ni-1; i++)
                {
//                    if(i > Ni/4. && i < Ni*3./4. && j > Nj/4. && j < Nj*3./4. && k > Nk/4. && k < Nk*3./4.)
//                    {
//                        for(int t = 0; t < 5; t++)
//                        {
//                            tet_labels.push_back(1);
//                        }
//                    }
//                    else {
//                        for(int t = 0; t < 5; t++)
//                        {
//                            tet_labels.push_back(0);
//                        }
//                    }
                    for(int t = 0; t < 5; t++)
                    {
                        tet_labels.push_back(1);
                    }
                }
            }
        }
    }
    
    void tetralize(std::vector<T>& points, std::vector<int>& tets, std::vector<int>& tet_labels)
    {
        create_points(points);
        create_tets(tets);
        label_tets(tet_labels);
    }
};

#endif
