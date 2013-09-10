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

#include <CGLA/Vec3d.h>
#include <CGLA/Vec4d.h>
#include <CGLA/Mat4x4d.h>
#include <CGLA/Mat3x3d.h>
#include <CGLA/Mat2x3d.h>
#include <CGLA/eigensolution.h>

class GELTypes
{
public:
    
    typedef double            real_type;
    typedef CGLA::Vec3d       vector3_type;
    typedef CGLA::Vec4d       vector4_type;
    typedef CGLA::Mat4x4d     matrix4x4_type;
    typedef CGLA::Mat3x3d     matrix3x3_type;
    typedef CGLA::Mat2x3d     matrix2x3_type;
    
    typedef CGLA::Axis       AXIS;
    
public:
    
    static real_type dot(vector3_type const & v1, vector3_type const & v2)
    {
        return CGLA::dot(v1, v2);
    }
    
    static vector3_type cross(vector3_type const & v1, vector3_type const & v2)
    {
        return CGLA::cross(v1, v2);
    }
    
    static real_type length(vector3_type const & v)
    {
        return CGLA::length(v);
    }
    
    static bool is_nan(real_type const & t)
    {
        return CGLA::isnan(t);
    }
    
    static real_type determinant(matrix3x3_type const & m)
    {
        return CGLA::determinant(m);
    }
    
    static real_type determinant(matrix4x4_type const & m)
    {
        return CGLA::determinant(m);
    }
    
    
    static matrix4x4_type transpose(matrix4x4_type const & m)
    {
        return CGLA::transpose(m);
    }
    
    static matrix4x4_type invert(matrix4x4_type const & m)
    {
        return CGLA::invert(m);
    }
    
    
    static real_type sqr_length(vector3_type const & v)
    {
        return CGLA::sqr_length(v);
    }
    
    static vector3_type normalize(vector3_type const & v)
    {
        return CGLA::normalize(v);
    }
    
    static void orthogonal(vector3_type const & v1, vector3_type & v2, vector3_type & v3)
    {
        CGLA::orthogonal(v1, v2, v3);
    }
    
    static void eigen(matrix3x3_type AA, matrix3x3_type Q, matrix3x3_type L)
    {
        CGLA::power_eigensolution(AA, Q, L);
    }
    
    static matrix3x3_type get_rotation_matrix(const CGLA::Axis& axis, const real_type& angle)
    {
        return CGLA::rotation_Mat3x3d(axis, angle);
    }
};
