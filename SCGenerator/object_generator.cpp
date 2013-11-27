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

#include "object_generator.h"

void ObjectGenerator::create_sphere(const std::vector<vec3>& points, const std::vector<int>& tets, const vec3& center, const real& radius, int label, std::vector<int>& tet_labels)
{
    for (unsigned int k = 0; k < tets.size()/4; k++)
    {
        bool inside = true;
        for (int j = 0; j < 4; j++)
        {
            vec3 v = points[tets[4*k + j]];
            for(int i = 0; i < 3; i++)
            {
                if((center - v).length() > radius)
                {
                    inside = false;
                    break;
                }
            }
        }
        if(inside)
        {
            tet_labels[k] = label;
        }
    }
}

void ObjectGenerator::create_cube(const std::vector<vec3>& points, const std::vector<int>& tets, const vec3& origin, const vec3& size, int label, std::vector<int>& tet_labels)
{
    vec3 max_pos = origin + size;
    for (unsigned int k = 0; k < tets.size()/4; k++)
    {
        bool inside = true;
        for (int j = 0; j < 4; j++)
        {
            vec3 v = points[tets[4*k + j]];
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
            tet_labels[k] = label;
        }
    }
}