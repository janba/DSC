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
    
    // TODO:
    template <typename DeformableSimplicialComplex>
    static void fit_mesh_to_object(DeformableSimplicialComplex& dsc)
    {
        
    }
    
    
public:
    
    template <typename DeformableSimplicialComplex>
    static void create_sphere(DeformableSimplicialComplex& dsc, const vec3& center, const real& radius, int label)
    {
        for (auto tit = dsc.tetrahedra_begin(); tit != dsc.tetrahedra_end(); tit++)
        {
            bool inside = true;
            auto verts = dsc.get_pos(dsc.get_nodes(tit.key()));
            for (auto &v : verts)
            {
                if((center - v).length() > radius)
                {
                    inside = false;
                    break;
                }
            }
            if(inside)
            {
                dsc.set_label(tit.key(), label);
            }
        }
    }
    
    template <typename DeformableSimplicialComplex>
    static void create_cube(DeformableSimplicialComplex& dsc, const vec3& origin, const vec3& size, int label)
    {
        vec3 max_pos = origin + size;
        for (auto tit = dsc.tetrahedra_begin(); tit != dsc.tetrahedra_end(); tit++)
        {
            bool inside = true;
            auto verts = dsc.get_pos(dsc.get_nodes(tit.key()));
            for (auto &v : verts)
            {
                for(int i = 0; i < 3; i++)
                {
                    if(v[i] < origin[i] || v[i] > max_pos[i])
                    {
                        inside = false;
                        break;
                    }
                }
            }
            if(inside)
            {
                dsc.set_label(tit.key(), label);
            }
        }
    }
};