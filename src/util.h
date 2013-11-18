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

#include <vector>
#include <list>
#include <sstream>
#include <cmath>
#include <cassert>

#include <CGLA/Vec3d.h>
#include <CGLA/Vec4d.h>
#include <CGLA/Mat3x3d.h>
#include <CGLA/Mat4x4d.h>

namespace DSC {
    
    typedef double            real;
    typedef CGLA::Vec3d       vec3;
    typedef CGLA::Vec4d       vec4;
    typedef CGLA::Mat3x3d     mat3;
    typedef CGLA::Mat4x4d     mat4;
    
    constexpr real EPSILON = 1e-8;
    
#ifndef INFINITY
    constexpr real INFINITY = 1e32;
#endif
    
    namespace Util
    {
        using CGLA::isnan;
        using CGLA::dot;
        using CGLA::cross;
        using CGLA::length;
        using CGLA::sqr_length;
        using CGLA::normalize;
        
        struct Plane {
            vec3 p;
            vec3 n;
        };
        
        inline real min(real x, real y)
        {
            return std::min(x, y);
        }
        
        inline real max(real x, real y)
        {
            return std::max(x, y);
        }
        
        template <typename vec3>
        inline vec3 normal_direction(const vec3& a, const vec3& b, const vec3& c);
        
        
        template <typename real>
        inline int sign(real val)
        {
            return (0. < val) - (val < 0.);
        }
        
        /**
         * Computes the signed area of the triangle spanned by vertices with positions v0, v1 and v2.
         */
        template <typename real, typename vec3>
        inline real signed_area(const vec3& v0, const vec3& v1, const vec3& v2)
        {
            vec3 n = cross(v1-v0, v2-v0);
            return 0.5 * n.length();
        }
        
        /**
         * Computes the area of the triangle spanned by vertices with positions v0, v1 and v2.
         */
        template <typename real, typename vec3>
        inline real area(const vec3& v0, const vec3& v1, const vec3& v2)
        {
            return std::abs(signed_area<real>(v0, v1, v2));
        }
        
        template <typename real, typename vec3>
        inline real signed_volume(const vec3& a, const vec3& b, const vec3& c, const vec3& d)
        {
            return dot(a-d, cross(b-d, c-d))/6.;
        }
        
        template <typename real, typename vec3>
        inline real volume(const vec3& a, const vec3& b, const vec3& c, const vec3& d)
        {
            return std::abs(signed_volume<real>(a, b, c, d));
        }
        
        /**
         * Calculates the cosine of the angle between the line segments |ab| and |ac|.
         */
        template <typename real, typename vec3>
        inline real cos_angle(const vec3& a, const vec3& b, const vec3& c)
        {
            vec3 ab = normalize(b - a);
            vec3 ac = normalize(c - a);
            return dot(ab, ac);
        }
        
        /**
         * Calculates the angle between the line segments |ab| and |ac|.
         */
        template <typename real, typename vec3>
        inline real angle(const vec3& a, const vec3& b, const vec3& c)
        {
            return acos(cos_angle<real>(a, b, c));
        }
        
        /**
         * Calculate the cosine of angles in the triangle defined by the vertices a, b and c.
         */
        template <typename real, typename vec3>
        inline std::vector<real> cos_angles(const vec3& a, const vec3& b, const vec3& c)
        {
            std::vector<real> cosines(3);
            cosines[0] = cos_angle<real>(a, b, c);
            cosines[1] = cos_angle<real>(b, c, a);
            cosines[2] = cos_angle<real>(c, a, b);
            return cosines;
        }
        
        template <typename real, typename vec3>
        inline real cos_min_angle(const vec3& a, const vec3& b, const vec3& c)
        {
            std::vector<real> cosines = cos_angles<real>(a, b, c);
            real max_cos = -1.;
            for(auto cos : cosines)
            {
                max_cos = std::max(cos, max_cos);
            }
            return max_cos;
        }
        
        template <typename real, typename vec3>
        inline real min_angle(const vec3& a, const vec3& b, const vec3& c)
        {
            return acos(cos_min_angle<real>(a, b, c));
        }
        
        template <typename real, typename vec3>
        inline real cos_max_angle(const vec3& a, const vec3& b, const vec3& c)
        {
            std::vector<real> cosines = cos_angles<real>(a, b, c);
            real min_cos = 1.;
            for(auto cos : cosines)
            {
                min_cos = std::min(cos, min_cos);
            }
            return min_cos;
        }
        
        template <typename real, typename vec3>
        inline real max_angle(const vec3& a, const vec3& b, const vec3& c)
        {
            return std::acos(cos_max_angle<real>(a, b, c));
        }
        
        /**
         * Returns the cosine to the dihedral angle between face |abc| and face |abd|.
         */
        template<typename real, typename vec3>
        inline real cos_dihedral_angle(const vec3& a, const vec3& b, const vec3& c, const vec3& d)
        {
            vec3 n0 = normal_direction(a, b, c);
            vec3 n1 = normal_direction(b, a, d);
            real angle = dot(n0, n1);
#ifdef DEBUG
            assert(angle < 1. + EPSILON);
            assert(angle > -1. - EPSILON);
#endif
            return angle;
        }
        
        /**
         * Returns the dihedral angle between face |abc| and face |abd|.
         */
        template<typename real, typename vec3>
        inline real dihedral_angle(const vec3& a, const vec3& b, const vec3& c, const vec3& d)
        {
            return std::acos(cos_dihedral_angle<real>(a, b, c, d));
        }
        
        template <typename vec3>
        inline vec3 barycenter(const vec3& a, const vec3& b)
        {
            return (a + b)*0.5;
        }
        
        template <typename vec3>
        inline vec3 barycenter(const vec3& a, const vec3& b, const vec3& c)
        {
            return (a + b + c)/3.;
        }
        
        template <typename vec3>
        inline vec3 barycenter(const vec3& a, const vec3& b, const vec3& c, const vec3& d)
        {
            return (a + b + c + d)*0.25;
        }
        
        /**
         * Finds the barycentric coordinates of point v in a triangle spanned by the vertices a, b and c.
         */
        template <typename real, typename vec3>
        inline std::vector<real> barycentric_coords(const vec3& p, const vec3& a, const vec3& b, const vec3& c)
        {
            std::vector<real> coords(3);
            
            vec3 v0 = b - a;
            vec3 v1 = c - a;
            vec3 v2 = p - a;
            real d00 = dot(v0, v0);
            real d01 = dot(v0, v1);
            real d11 = dot(v1, v1);
            real d20 = dot(v2, v0);
            real d21 = dot(v2, v1);
            real denom = d00 * d11 - d01 * d01;
#ifdef DEBUG
            assert(denom != 0.);
#endif
            coords[0] = (d11 * d20 - d01 * d21) / denom;
            coords[1] = (d00 * d21 - d01 * d20) / denom;
            coords[2] = 1. - coords[0] - coords[1];
            
            return coords;
        }
        
        /**
         * Calculates the barycentric coordinates of a point v in a tetrahedron spanned by the four vertices a, b, c and d.
         */
        template <typename real, typename vec3>
        inline std::vector<real> barycentric_coords(const vec3& p, const vec3& a, const vec3& b, const vec3& c, const vec3& d)
        {
            std::vector<real> coords(4);
            coords[0] = signed_volume<real>(p, b, c, d);
            coords[1] = signed_volume<real>(a, p, c, d);
            coords[2] = signed_volume<real>(a, b, p, d);
            coords[3] = signed_volume<real>(a, b, c, p);
            
            real s = coords[0] + coords[1] + coords[2] + coords[3];
            for (unsigned int i = 0; i < 4; ++i)
            {
                coords[i] /= s;
            }
            return coords;
        }
        
        template <typename vec3>
        inline vec3 normal_direction(const vec3& a, const vec3& b, const vec3& c)
        {
            vec3 ab = b - a;
            vec3 ac = c - a;
            vec3 n = cross(ab, ac);
#ifdef DEBUG
            assert(!isnan(n[0]) && !isnan(n[1]) && !isnan(n[2]));
#endif
            return normalize(n);
        }
        
        template <typename real, typename vec3>
        inline vec3 normal_direction(const vec3& a, const vec3& b, const vec3& c, const vec3& d)
        {
            vec3 n = normal_direction<real>(a, b, c);
            vec3 bf = barycenter<real>(a, b, c);
            vec3 bt = barycenter<real>(a, b, c, d);
            vec3 v_out = bf - bt;
            if (dot(v_out, n) > 0)
                return n;
            else
                return -n;
        }
        
        /**
         Returns p projected onto the line spanned by the two points a and b.
         */
        template <typename vec3>
        inline vec3 project(const vec3& p, const vec3& a, const vec3& b)
        {
            vec3 v1 = p - a;
            vec3 v2 = b - a;
            return a + v2 * dot(v1,v2)/dot(v2, v2);
        }
        
        /**
         * Project the point p onto the plane spanned by the three points a, b and c.
         */
        template<typename vec3>
        inline vec3 project(const vec3& p, const vec3& a, const vec3& b, const vec3& c)
        {
            vec3 normal = normal_direction(a, b, c);
            return p - normal * dot(p - a, normal);
        }
        
        template<typename real, typename vec3>
        inline real ms_length(const vec3& a, const vec3& b, const vec3& c, const vec3& d)
        {
            real result = 0.;
            result += sqr_length(a - b);
            result += sqr_length(a - c);
            result += sqr_length(a - d);
            result += sqr_length(b - c);
            result += sqr_length(b - d);
            result += sqr_length(c - d);
            return result / 6.;
        }
        
        template<typename real, typename vec3>
        inline real rms_length(const vec3& a, const vec3& b, const vec3& c, const vec3& d)
        {
            return sqrt(ms_length<real>(a, b, c, d));
        }
        
        template<typename real, typename vec3>
        inline real quality(const vec3& a, const vec3& b, const vec3& c, const vec3& d)
        {
            real v = signed_volume<real>(a, b, c, d);
            real lrms = rms_length<real>(a, b, c, d);
            
            real q = 8.48528 * v / (lrms * lrms * lrms);
#ifdef DEBUG
            assert(!isnan(q));
#endif
            return q;
        }
        
        /**
         * Returns whether the point p is on the inside (in the direction away from the normal) from the plane defined by the point a and the normal.
         */
        template<typename vec3>
        inline bool is_inside(const vec3& p, const vec3& a, const vec3& normal)
        {
            if(sqr_length(p-a) == 0.)
            {
                return true;
            }
            return dot(p-a, normal) <= 0.;
        }
        
        /**
         Returns whether the point p is between the points in the vector corners.
         */
        template<typename real, typename vec3>
        inline bool is_inside(const vec3& p, std::vector<vec3> corners)
        {
            while(corners.size() > 2)
            {
                auto bc = barycentric_coords<real>(p, corners[0], corners[1], corners[2]);
                if(bc[0] > -EPSILON && bc[1] > -EPSILON && bc[2] > -EPSILON)
                {
                    return true;
                }
                corners.erase(corners.begin() + 1);
            }
            return false;
        }
        
        /**
         * Returns whether points p1 and p2 lies on the same side of the plane spanned by points a, b and c. If p1 or p2 lies on the plane, the method returns false.
         */
        template<typename real, typename vec3>
        inline bool is_on_same_side(const vec3& p1, const vec3& p2, const vec3& a, const vec3& b, const vec3& c)
        {
            auto normal = normal_direction(a, b, c);
            auto d1 = dot(p1 - a, normal);
            auto d2 = dot(p2 - a, normal);
            if(std::abs(d1) > EPSILON && std::abs(d2) > EPSILON && sign(d1) == sign(d2))
            {
                return true;
            }
            return false;
        }
        
        /**
         * Returns the shortest distance from the point p to the plane defined by the point a and the normal.
         */
        template<typename real, typename vec3>
        inline real distance_plane(const vec3& p, const vec3& a, const vec3& normal)
        {
            vec3 v = p - a;
            return std::abs(dot(v, normal));
        }
        
        /**
         * Returns the shortest distance from the point p to the plane spanned by the points a, b and c.
         */
        template<typename real, typename vec3>
        inline real distance_plane(const vec3& p, const vec3& a, const vec3& b, const vec3& c)
        {
            vec3 n = normal_direction(a, b, c);
            return distance_plane<real>(p, a, n);
        }
        
        /**
         * Calculates the intersection between the line segment defined by p + t*r where 0 <= t <= 1 and the plane defined by the point a and the normal. The intersection point is defined by p + t*r and the function returns t. Returns infinity if it does not intersect.
         */
        template<typename real, typename vec3>
        inline real intersection_ray_plane(const vec3& p, const vec3& r, const vec3& a, const vec3& normal)
        {
            real n = dot(normal, a - p);
            real d = dot(normal, r);
            
            if (std::abs(d) < EPSILON) // Plane and line are parallel if true.
            {
                if (std::abs(n) < EPSILON)
                {
                    return 0.; // Line intersection
                }
                return INFINITY; // No intersection.
            }
            
            // Compute the t value for the directed line ray intersecting the plane.
            return n / d;
        }
        
        /**
         * Calculates the intersection between the line segment defined by p + t*r where 0 <= t <= 1 and the plane spanned by the vertices a, b and c. The intersection point is defined by p + t*r and the function returns t. Returns infinity if it does not intersect.
         */
        template<typename real, typename vec3>
        inline real intersection_ray_plane(const vec3& p, const vec3& r, const vec3& a, const vec3& b, const vec3& c)
        {
            vec3 normal = normal_direction(a, b, c);
            return intersection_ray_plane<real>(p, r, a, normal);
        }
        
        /**
         * Calculates the intersection between the line segment defined by p + t*r where 0 <= t <= 1 and the triangle |a b c|. The intersection point is defined by p + t*r and the function returns t. Returns infinity if it does not intersect.
         */
        template<typename real, typename vec3>
        inline real intersection_ray_triangle(const vec3& p, const vec3& r, const vec3& a, const vec3& b, const vec3& c)
        {
            real t = intersection_ray_plane<real>(p, r, a, b, c);
            if(t < 0.) // The ray goes away from the triangle
            {
                return t;
            }
            
            std::vector<real> coords = barycentric_coords<real>(p + t*r, a, b, c);
            if(coords[0] > EPSILON && coords[1] > EPSILON && coords[2] > EPSILON) // The intersection happens inside the triangle.
            {
                return t;
            }
            return INFINITY; // The intersection happens outside the triangle.
        }
        
        /**
         Concatonates the integer number to the string name.
         */
        inline std::string concat4digits(std::string name, int number)
        {
            std::ostringstream s;
            if (number < 10)
                s << name << "000" << number;
            else if (number < 100)
                s << name << "00" << number;
            else if (number < 1000)
                s << name << "0" << number;
            else
                s << name << number;
            return s.str();
        }
        
        /**
         Returns the maximum difference between the values x[i] and y[i] for all i less than the size of x and y.
         */
        template<typename real>
        inline real max_diff(const std::vector<real>& x, const std::vector<real>& y)
        {
            real max_diff = -INFINITY;
            for (int i = 0; i < x.size(); i++)
            {
                if (i < y.size()) {
                    max_diff = max(std::abs(x[i] - y[i]), max_diff);
                }
            }
            return max_diff;
        }
    }
    
}
