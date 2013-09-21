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

#include "design_domain.h"

namespace DSC {
        
    void DesignDomain::clamp_vector(const vec3& p, vec3& v) const
    {
        if(!is_inside(p+v))
        {
            for (auto &plane : planes) {
                
                real t = Util::intersection_ray_plane<real>(p, v, plane.p, plane.n);
                if(t >= 0. && t < 1.)
                {
                    v = t*v;
                }
            }
        }
    }
    
    bool DesignDomain::is_inside(const vec3& p) const
    {
        for (auto &plane : planes) {
            if (!Util::is_inside(p, plane.p, plane.n)) {
                return false;
            }
        }
        return true;
    }
    
    
    
    bool DesignDomain::is_inside(const std::vector<vec3>& verts) const
    {
        for (auto &p : verts) {
            if (!is_inside(p)) {
                return false;
            }
        }
        return true;
    }
    
}