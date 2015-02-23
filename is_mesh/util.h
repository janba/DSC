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
#include <map>
#include <sstream>
#include <cmath>
#include <cassert>
#include <limits>

#include <CGLA/Vec2d.h>
#include <CGLA/Vec3d.h>
#include <CGLA/Vec4d.h>
#include <CGLA/Mat3x3d.h>
#include <CGLA/Mat4x4d.h>

using vec2 = CGLA::Vec2d;
using vec3 = CGLA::Vec3d;
using vec4 = CGLA::Vec4d;
using mat3 = CGLA::Mat3x3d;
using mat4 = CGLA::Mat4x4d;

static const double EPSILON = 1e-8;

#undef INFINITY
static const double INFINITY = std::numeric_limits<double>::max();;


#ifdef __GNUC__
#define DEPRECATED __attribute__((deprecated))
#elif defined(_MSC_VER)
#define DEPRECATED __declspec(deprecated)
#else
#pragma message("WARNING: You need to implement DEPRECATED for this compiler")
#define DEPRECATED
#endif

namespace Util
{
    using CGLA::isnan;
    using CGLA::dot;
    using CGLA::cross;
    using CGLA::length;
    using CGLA::sqr_length;
    using CGLA::normalize;
    
    inline double min(double x, double y)
    {
        return std::min(x, y);
    }
    
    inline double max(double x, double y)
    {
        return std::max(x, y);
    }

    inline vec3 normal_direction(const vec3& a, const vec3& b, const vec3& c){
        vec3 ab = b - a;
        vec3 ac = c - a;
        vec3 n = cross(ab, ac);
#ifdef DEBUG
        assert(!isnan(n[0]) && !isnan(n[1]) && !isnan(n[2]));
#endif
        return normalize(n);
    }
    
    
    inline int sign(double val)
    {
        return (0. < val) - (val < 0.);
    }
    
    /**
     * Computes the area of the triangle spanned by vertices with positions v0, v1 and v2.
     */
    inline double area(const vec3& v0, const vec3& v1, const vec3& v2)
    {
        return 0.5 * length(cross(v1-v0, v2-v0));
    }
    
    inline double signed_volume(const vec3& a, const vec3& b, const vec3& c, const vec3& d)
    {
        return dot(a-d, cross(b-d, c-d))/6.;
    }
    
    inline double volume(const vec3& a, const vec3& b, const vec3& c, const vec3& d)
    {
        return std::abs(signed_volume(a, b, c, d));
    }
    
    /**
     * Calculates the cosine of the angle between the line segments |ab| and |ac|.
     */
    inline double cos_angle(const vec3& a, const vec3& b, const vec3& c)
    {
        vec3 ab = normalize(b - a);
        vec3 ac = normalize(c - a);
        return dot(ab, ac);
    }
    
    /**
     * Calculates the angle between the line segments |ab| and |ac|.
     */
    inline double angle(const vec3& a, const vec3& b, const vec3& c)
    {
        return acos(cos_angle(a, b, c));
    }
    
    /**
     * Returns a rotation matrix.
     */
    inline mat3 rotation_matrix(int dimension, double angle)
    {
        mat3 m(0.);
        
        switch(dimension)
        {
            case 0:
                m[0][0] = 1.0;
                m[1][1] = cos(angle);
                m[1][2] = sin(angle);
                m[2][1] = -sin(angle);
                m[2][2] = cos(angle);
                break;
            case 1:
                m[0][0] = cos(angle);
                m[0][2] = -sin(angle);
                m[2][0] = sin(angle);
                m[2][2] = cos(angle);
                m[1][1] = 1.0;
                break;
            case 2:
                m[0][0] = cos(angle);
                m[0][1] = sin(angle);
                m[1][0] = -sin(angle);
                m[1][1] = cos(angle);
                m[2][2] = 1.0;
                break;
        }
        
        return m;
    }
    
    /**
     * Returns the derivative of a rotation matrix.
     */
    inline mat3 d_rotation_matrix(int dimension, double angle)
    {
        mat3 m(0.);
        
        switch(dimension)
        {
            case 0:
                m[1][1] = -sin(angle);
                m[1][2] = cos(angle);
                m[2][1] = -cos(angle);
                m[2][2] = -sin(angle);
                break;
            case 1:
                m[0][0] = -sin(angle);
                m[0][2] = -cos(angle);
                m[2][0] = cos(angle);
                m[2][2] = -sin(angle);
                break;
            case 2:
                m[0][0] = -sin(angle);
                m[0][1] = cos(angle);
                m[1][0] = -cos(angle);
                m[1][1] = -sin(angle);
                break;
        }
        
        return m;
    }
    
    /**
     * Calculate the cosine of angles in the triangle defined by the vertices a, b and c.
     */
    inline std::vector<double> cos_angles(const vec3& a, const vec3& b, const vec3& c)
    {
        std::vector<double> cosines(3);
        cosines[0] = cos_angle(a, b, c);
        cosines[1] = cos_angle(b, c, a);
        cosines[2] = cos_angle(c, a, b);
        return cosines;
    }
    
    inline double cos_min_angle(const vec3& a, const vec3& b, const vec3& c)
    {
        std::vector<double> cosines = cos_angles(a, b, c);
        double max_cos = -1.;
        for(auto cos : cosines)
        {
            max_cos = std::max(cos, max_cos);
        }
        return max_cos;
    }
    
    inline double min_angle(const vec3& a, const vec3& b, const vec3& c)
    {
        return acos(cos_min_angle(a, b, c));
    }
    
    inline double cos_max_angle(const vec3& a, const vec3& b, const vec3& c)
    {
        std::vector<double> cosines = cos_angles(a, b, c);
        double min_cos = 1.;
        for(auto cos : cosines)
        {
            min_cos = std::min(cos, min_cos);
        }
        return min_cos;
    }
    
    inline double max_angle(const vec3& a, const vec3& b, const vec3& c)
    {
        return std::acos(cos_max_angle(a, b, c));
    }
    
    /**
     * Returns the cosine to the dihedral angle between face |abc| and face |abd|.
     */
    inline double cos_dihedral_angle(const vec3& a, const vec3& b, const vec3& c, const vec3& d)
    {
        vec3 n0 = normal_direction(a, b, c);
        vec3 n1 = normal_direction(b, a, d);
        double angle = dot(n0, n1);
#ifdef DEBUG
        assert(angle < 1. + EPSILON);
        assert(angle > -1. - EPSILON);
#endif
        return angle;
    }
    
    /**
     * Returns the dihedral angle between face |abc| and face |abd|.
     */
    inline double dihedral_angle(const vec3& a, const vec3& b, const vec3& c, const vec3& d)
    {
        return std::acos(cos_dihedral_angle(a, b, c, d));
    }
    
    inline vec3 barycenter(const vec3& a, const vec3& b)
    {
        return (a + b)*0.5;
    }
    
    inline vec3 barycenter(const vec3& a, const vec3& b, const vec3& c)
    {
        return (a + b + c)/3.;
    }
    
    inline vec3 barycenter(const vec3& a, const vec3& b, const vec3& c, const vec3& d)
    {
        return (a + b + c + d)*0.25;
    }
    
    /**
     * Finds the barycentric coordinates of point v in a triangle spanned by the vertices a, b and c.
     */
    inline std::vector<double> barycentric_coords(const vec3& p, const vec3& a, const vec3& b, const vec3& c)
    {
        vec3 v0 = b - a;
        vec3 v1 = c - a;
        vec3 v2 = p - a;
        double d00 = dot(v0, v0);
        double d01 = dot(v0, v1);
        double d11 = dot(v1, v1);
        double d20 = dot(v2, v0);
        double d21 = dot(v2, v1);
        double denom = d00 * d11 - d01 * d01;
#ifdef DEBUG
        assert(std::abs(denom) > EPSILON);
#endif
        std::vector<double> coords(3);
        coords[1] = (d11 * d20 - d01 * d21) / denom;
        coords[2] = (d00 * d21 - d01 * d20) / denom;
        coords[0] = 1. - coords[1] - coords[2];
        
        return coords;
    }
    
    /**
     * Calculates the barycentric coordinates of a point v in a tetrahedron spanned by the four vertices a, b, c and d.
     */
    inline std::vector<double> barycentric_coords(const vec3& p, const vec3& a, const vec3& b, const vec3& c, const vec3& d)
    {
        std::vector<double> coords(4);
        coords[0] = signed_volume(p, b, c, d);
        coords[1] = signed_volume(a, p, c, d);
        coords[2] = signed_volume(a, b, p, d);
        coords[3] = signed_volume(a, b, c, p);
        
        double s = coords[0] + coords[1] + coords[2] + coords[3];
        for (unsigned int i = 0; i < 4; ++i)
        {
            coords[i] /= s;
        }
        return coords;
    }

    
    inline vec3 normal_direction(const vec3& a, const vec3& b, const vec3& c, const vec3& d)
    {
        vec3 n = normal_direction(a, b, c);
        vec3 bf = barycenter(a, b, c);
        vec3 bt = barycenter(a, b, c, d);
        vec3 v_out = bf - bt;
        if (dot(v_out, n) > 0)
            return n;
        else
            return -n;
    }
    
    /**
     Returns p projected onto the line spanned by the point a and the vector v.
     */
    inline vec3 project_point_line(const vec3& p, const vec3& a, const vec3& v)
    {
        vec3 v1 = p - a;
        return a + v * dot(v1,v)/dot(v, v);
    }
    
    /**
     Returns p projected onto the linesegment spanned by the points a and b.
     */
    inline vec3 project_point_linesegment(const vec3& p, const vec3& a, const vec3& b)
    {
        vec3 v1 = p - a;
        vec3 v2 = b - a;
        double d = dot(v1,v2)/dot(v2,v2);
        if(d < 0.)
        {
            return a;
        }
        if(d > 1.)
        {
            return b;
        }
        return a + v2 * d;
    }
    
    /**
     * Project the point p onto the plane defined by the point a and the normal n.
     */
    inline vec3 project_point_plane(const vec3& p, const vec3& a, const vec3& n)
    {
        return p - n * dot(p - a, n);
    }
    
    /**
     * Project the point p onto the plane spanned by the three points a, b and c.
     */
    inline vec3 project_point_plane(const vec3& p, const vec3& a, const vec3& b, const vec3& c)
    {
        vec3 normal = normal_direction(a, b, c);
        return project_point_plane(p, a, normal);
    }
    
    /**
     * Returns the closest point on the triangle |abc| from the point p.
     */
    inline vec3 closest_point_on_triangle(const vec3& p, const vec3& a, const vec3& b, const vec3& c)
    {
        auto bc = Util::barycentric_coords(p, a, b, c);
        bool b1 = bc[0] > -EPSILON, b2 = bc[1] > -EPSILON, b3 = bc[2] > -EPSILON;
        if(b1 && b2 && b3)
        {
            return project_point_plane(p, a, b, c);
        }
        if(b2 && b3)
        {
            return project_point_linesegment(p, b, c);
        }
        if(b1 && b3)
        {
            return project_point_linesegment(p, a, c);
        }
        if(b1 && b2)
        {
            return project_point_linesegment(p, a, b);
        }
        if(b1)
        {
            return a;
        }
        if(b2)
        {
            return b;
        }
        if(b3)
        {
            return c;
        }
        assert(false);
        return c;
    }
    
    inline double ms_length(const vec3& a, const vec3& b, const vec3& c, const vec3& d)
    {
        double result = 0.;
        result += sqr_length(a - b);
        result += sqr_length(a - c);
        result += sqr_length(a - d);
        result += sqr_length(b - c);
        result += sqr_length(b - d);
        result += sqr_length(c - d);
        return result / 6.;
    }
    
    inline double rms_length(const vec3& a, const vec3& b, const vec3& c, const vec3& d)
    {
        return sqrt(ms_length(a, b, c, d));
    }

    inline double quality(const vec3& a, const vec3& b, const vec3& c, const vec3& d)
    {
        double v = signed_volume(a, b, c, d);
        double lrms = rms_length(a, b, c, d);
        
        double q = 8.48528 * v / (lrms * lrms * lrms);
#ifdef DEBUG
        assert(!isnan(q));
#endif
        return q;
    }
    
    /**
     * Returns whether the point p is on the inside (in the direction away from the normal) from the plane defined by the point a and the normal.
     */
    inline bool is_inside(const vec3& p, const vec3& a, const vec3& normal)
    {
        if(sqr_length(p-a) == 0.)
        {
            return true;
        }
        return dot(p-a, normal) <= 0.;
    }
    
    /**
     Returns whether the point p is inside the triangle with corners a, b and c. The point p and the triangle must lie in the same plane.
     */
    inline bool is_inside(const vec3& p, const vec3& a, const vec3& b, const vec3& c)
    {
        auto bc = barycentric_coords(p, a, b, c);
        return bc[0] > -EPSILON && bc[1] > -EPSILON && bc[2] > -EPSILON;
    }
    
    /**
     Returns whether the point p is between the points in the vector corners.
     */
    inline bool is_inside(const vec3& p, std::vector<vec3> corners)
    {
        while(corners.size() > 2)
        {
            if(is_inside(p, corners[0], corners[1], corners[2]))
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
    inline bool is_on_same_side(const vec3& p1, const vec3& p2, const vec3& a, const vec3& normal)
    {
        auto d1 = dot(p1 - a, normal);
        auto d2 = dot(p2 - a, normal);
        if(std::abs(d1) > EPSILON && std::abs(d2) > EPSILON && sign(d1) == sign(d2))
        {
            return true;
        }
        return false;
    }
    
    /**
     * Returns the shortest distance from the point p to the line segment defined by the two end points a and b.
     */
    inline double distance_point_linesegment(const vec3& p, const vec3& a, const vec3& b)
    {
        vec3 v1 = p - a;
        vec3 v2 = b - a;
        double d = dot(v1, v2)/sqr_length(v2);
        if(d <= 0.)
        {
            return length(p - a);
        }
        if(d >= 1.)
        {
            return length(p - b);
        }
        return length(v1 - v2 * d);
    }
    
    /**
     * Returns the shortest distance from the line segment |ab| to the line segment |cd|.
     */
    inline double distance_linesegment_linesegment(const vec3& a, const vec3& b, const vec3& c, const vec3& d)
    {
        vec3 u = b - a;
        vec3 v = d - c;
        vec3 w = a - c;
        double duu = dot(u,u);         // always >= 0
        double duv = dot(u,v);
        double dvv = dot(v,v);         // always >= 0
        double duw = dot(u,w);
        double dvw = dot(v,w);
        double D = duu*dvv - duv*duv;        // always >= 0
        double sc, sN, sD = D;       // sc = sN / sD, default sD = D >= 0
        double tc, tN, tD = D;       // tc = tN / tD, default tD = D >= 0
        
        // compute the line parameters of the two closest points
        if (D < EPSILON) { // the lines are almost parallel
            sN = 0.0;         // force using point P0 on segment S1
            sD = 1.0;         // to prevent possible division by 0.0 later
            tN = dvw;
            tD = dvv;
        }
        else {                 // get the closest points on the infinite lines
            sN = (duv*dvw - dvv*duw);
            tN = (duu*dvw - duv*duw);
            if (sN < 0.0) {        // sc < 0 => the s=0 edge is visible
                sN = 0.0;
                tN = dvw;
                tD = dvv;
            }
            else if (sN > sD) {  // sc > 1  => the s=1 edge is visible
                sN = sD;
                tN = dvw + duv;
                tD = dvv;
            }
        }
        
        if (tN < 0.0) {            // tc < 0 => the t=0 edge is visible
            tN = 0.0;
            // recompute sc for this edge
            if (-duw < 0.0)
                sN = 0.0;
            else if (-duw > duu)
                sN = sD;
            else {
                sN = -duw;
                sD = duu;
            }
        }
        else if (tN > tD) {      // tc > 1  => the t=1 edge is visible
            tN = tD;
            // recompute sc for this edge
            if ((-duw + duv) < 0.0)
                sN = 0;
            else if ((-duw + duv) > duu)
                sN = sD;
            else {
                sN = (-duw +  duv);
                sD = duu;
            }
        }
        // finally do the division to get sc and tc
        sc = (abs(sN) < EPSILON ? 0.0 : sN / sD);
        tc = (abs(tN) < EPSILON ? 0.0 : tN / tD);
        
        // get the difference of the two closest points
        vec3 dP = w + (sc * u) - (tc * v);  // =  S1(sc) - S2(tc)
        
        return length(dP);   // return the closest distance
    }
    
    /**
     * Returns the shortest distance from the point p to the line defined by the point a and the vector v.
     */
    inline double distance_point_line(const vec3& p, const vec3& a, const vec3& v)
    {
        vec3 v1 = p - a;
        return length(v1 - v * dot(v1,v));
    }
    
    /**
     * Returns the shortest distance from the point p to the plane defined by the point a and the normal n.
     */
    inline double distance_point_plane(const vec3& p, const vec3& a, const vec3& n)
    {
        return std::abs(dot(p - a, n));
    }
    
    /**
     * Returns the shortest distance from the point p to the triangle spanned by the points a, b and c.
     */
    inline double distance_point_triangle(const vec3& p, const vec3& a, const vec3& b, const vec3& c)
    {
        vec3 p_proj = project_point_plane(p, a, b, c);
        if(is_inside(p_proj, a, b, c))
        {
            return length(p_proj - p);
        }
        
        double d_line_ab = distance_point_linesegment(p, a, b);
        double d_line_bc = distance_point_linesegment(p, b, c);
        double d_line_ca = distance_point_linesegment(p, c, a);
        
        return std::min(std::min(d_line_ab, d_line_bc), d_line_ca);
    }
    
    /**
     * Returns the shortest distance from the triangle |abc| to the triangle |def|.
     */
    inline double distance_triangle_triangle(const vec3& a, const vec3& b, const vec3& c, const vec3& d, const vec3& e, const vec3& f)
    {
        // Calculate 6 point-triangle distances
        double dist = distance_point_triangle(a, d, e, f);
        dist = std::min(dist, distance_point_triangle(b, d, e, f));
        dist = std::min(dist, distance_point_triangle(c, d, e, f));
        dist = std::min(dist, distance_point_triangle(d, a, b, c));
        dist = std::min(dist, distance_point_triangle(e, a, b, c));
        dist = std::min(dist, distance_point_triangle(f, a, b, c));
        
        // Calculate 9 line segment-line segment distances
        dist = std::min(dist, distance_linesegment_linesegment(a, b, d, e));
        dist = std::min(dist, distance_linesegment_linesegment(a, b, d, f));
        dist = std::min(dist, distance_linesegment_linesegment(a, b, e, f));
        dist = std::min(dist, distance_linesegment_linesegment(a, c, d, e));
        dist = std::min(dist, distance_linesegment_linesegment(a, c, d, f));
        dist = std::min(dist, distance_linesegment_linesegment(a, c, e, f));
        dist = std::min(dist, distance_linesegment_linesegment(b, c, d, e));
        dist = std::min(dist, distance_linesegment_linesegment(b, c, d, f));
        dist = std::min(dist, distance_linesegment_linesegment(b, c, e, f));
        
        return dist;
    }
    
    /**
     * Calculates the intersection between the line segment defined by p + t*r where 0 <= t <= 1 and the plane defined by the point a and the normal. The intersection point is defined by p + t*r and the function returns t. Returns infinity if it does not intersect.
     */
    inline double intersection_ray_plane(const vec3& p, const vec3& r, const vec3& a, const vec3& normal)
    {
        double n = dot(normal, a - p);
        double d = dot(normal, r);
        
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
    inline double intersection_ray_plane(const vec3& p, const vec3& r, const vec3& a, const vec3& b, const vec3& c)
    {
        vec3 normal = normal_direction(a, b, c);
        return intersection_ray_plane(p, r, a, normal);
    }
    
    /**
     * Calculates the intersection between the line segment defined by p + t*r where 0 <= t <= 1 and the triangle |a b c|. The intersection point is defined by p + t*r and the function returns t. Returns infinity if it does not intersect.
     */
    inline double intersection_ray_triangle(const vec3& p, const vec3& r, const vec3& a, const vec3& b, const vec3& c)
    {
        double t = intersection_ray_plane(p, r, a, b, c);
        if(t < 0.) // The ray goes away from the triangle
        {
            return t;
        }
        
        std::vector<double> coords = barycentric_coords(p + t*r, a, b, c);
        if(coords[0] > EPSILON && coords[1] > EPSILON && coords[2] > EPSILON) // The intersection happens inside the triangle.
        {
            return t;
        }
        return INFINITY; // The intersection happens outside the triangle.
    }
    
    /**
     Concatenates the integer number to the string name.
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
    inline double max_diff(const std::vector<double>& x, const std::vector<double>& y)
    {
        double max_diff = -INFINITY;
        for (unsigned int i = 0; i < x.size(); i++)
        {
            if (i < y.size()) {
                max_diff = max(std::abs(x[i] - y[i]), max_diff);
            }
        }
        return max_diff;
    }
}
