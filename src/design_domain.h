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
    
    /**
     A domain class for specifying the design domain.
     */
    class DesignDomain
    {
        struct Plane {
            vec3 p;
            vec3 n;
        };
        
        std::vector<Plane> planes;
        
    public:
        enum DESIGN_DOMAIN_TYPE {CUBE};
        
        /**
         Creates a design domain defined by the design domain type and size. It is possible to specify a boundary gap which translates the entire domain by the amount specified by the input parameter.
         */
        DesignDomain(DESIGN_DOMAIN_TYPE design, const vec3& size)
        {
            vec3 center(0.);
            switch (design) {
                case CUBE:
                    for (int i = 0; i < 3; i++) {
                        vec3 point1 = center;
                        vec3 normal1(0.);
                        point1[i] += size[i]/2.;
                        normal1[i] = 1.;
                        planes.push_back({point1, normal1});
                        
                        vec3 point2 = center;
                        vec3 normal2(0.);
                        point2[i] -= size[i]/2.;
                        normal2[i] = -1.;
                        planes.push_back({point2, normal2});
                    }
                    break;
            }
        }
        
        /**
         Clamps the position pos to be within the domain.
         */
        void clamp_position(vec3& p) const;
        
        /**
         Clamps p + v to be within the domain by scaling the vector v if p is inside the domain. The position p is therefore not garanteed to be within the domain.
         */
        void clamp_vector(const vec3& p, vec3& v) const;
        
        /**
         Returns whether the position p is inside the domain.
         */
        bool is_inside(const vec3& p) const;
    };
    
}