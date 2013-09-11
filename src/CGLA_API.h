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

#include <CGLA/Vec2d.h>
#include <CGLA/Vec3d.h>
#include <CGLA/Vec4d.h>
#include <CGLA/Mat3x3d.h>
#include <CGLA/Mat4x4d.h>

typedef double            real;
typedef CGLA::Vec2d       vec2;
typedef CGLA::Vec3d       vec3;
typedef CGLA::Vec4d       vec4;
typedef CGLA::Mat3x3d     mat3;
typedef CGLA::Mat4x4d     mat4;

typedef CGLA::Axis       AXIS;


inline bool is_nan(real const & t)
{
    return CGLA::isnan(t);
}

inline real dot(vec2 const & v1, vec2 const & v2)
{
    return CGLA::dot(v1, v2);
}

inline real length(vec2 const & v)
{
    return CGLA::length(v);
}

inline real sqr_length(vec2 const & v)
{
    return CGLA::sqr_length(v);
}

inline vec2 normalize(vec2 const & v)
{
    return CGLA::normalize(v);
}

inline real dot(vec3 const & v1, vec3 const & v2)
{
    return CGLA::dot(v1, v2);
}

inline vec3 cross(vec3 const & v1, vec3 const & v2)
{
    return CGLA::cross(v1, v2);
}

inline real length(vec3 const & v)
{
    return CGLA::length(v);
}

inline real sqr_length(vec3 const & v)
{
    return CGLA::sqr_length(v);
}

inline vec3 normalize(vec3 const & v)
{
    return CGLA::normalize(v);
}

inline real determinant(mat3 const & m)
{
    return CGLA::determinant(m);
}

inline real determinant(mat4 const & m)
{
    return CGLA::determinant(m);
}

inline mat3 transpose(mat3 const & m)
{
    return CGLA::transpose(m);
}

inline mat4 transpose(mat4 const & m)
{
    return CGLA::transpose(m);
}

inline mat3 invert(mat3 const & m)
{
    return CGLA::invert(m);
}

inline mat4 invert(mat4 const & m)
{
    return CGLA::invert(m);
}

inline mat3 rotation_mat3(const CGLA::Axis& axis, const real& angle)
{
    return CGLA::rotation_Mat3x3d(axis, angle);
}

inline mat4 rotation_mat4(const CGLA::Axis& axis, const real& angle)
{
    return CGLA::rotation_Mat4x4d(axis, angle);
}
