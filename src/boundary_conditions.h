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
        std::vector<vec3> corners; // Specified in a clockwise order
        real volume = -1.;
        
    public:
        enum DESIGN_DOMAIN_TYPE {RECTANGLE, L, ESO};
        
        /**
         Creates a design domain defined by the design domain type and size. It is possible to specify a boundary gap which translates the entire domain by the amount specified by the input parameter.
         */
        DesignDomain(DESIGN_DOMAIN_TYPE design, int SIZE_X, int SIZE_Y, real boundary)
        {
            switch (design) {
                case RECTANGLE:
//                    corners.push_back(vec2(0.,0.));
//                    corners.push_back(vec2(0., SIZE_Y));
//                    corners.push_back(vec2(SIZE_X, SIZE_Y));
//                    corners.push_back(vec2(SIZE_X, 0.));
                    break;
                case L:
//                    corners.push_back(vec2(0.,0.));
//                    corners.push_back(vec2(0., SIZE_Y));
//                    corners.push_back(vec2(SIZE_X/2., SIZE_Y));
//                    corners.push_back(vec2(SIZE_X/2., SIZE_Y/2.));
//                    corners.push_back(vec2(SIZE_X, SIZE_Y/2.));
//                    corners.push_back(vec2(SIZE_X, 0.));
                    break;
                case ESO:
//                    corners.push_back(vec2(0.,0.));
//                    corners.push_back(vec2(0., 3.*SIZE_Y/7.));
//                    corners.push_back(vec2(30.*SIZE_X/32., 3.*SIZE_Y/7.));
//                    corners.push_back(vec2(30.*SIZE_X/32., SIZE_Y));
//                    corners.push_back(vec2(31.*SIZE_X/32., SIZE_Y));
//                    corners.push_back(vec2(31.*SIZE_X/32., 3.*SIZE_Y/7.));
//                    corners.push_back(vec2(SIZE_X, 3.*SIZE_Y/7.));
//                    corners.push_back(vec2(SIZE_X, 0.));
                    break;
            }
            
            for(auto &c : corners)
            {
                c[0] += boundary;
                c[1] += boundary;
                c[2] += boundary;
            }
        }
        
        /**
         Returns the corners of the design domain.
         */
        std::vector<vec3> get_corners() const;
        
        /**
         Returns an approximate center of the design domain.
         */
        vec3 get_center();
        
        /**
         Returns the total volume of the domain.
         */
        real get_volume();
        
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