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

#include "CGLA_API.h"

constexpr double EPSILON = 1e-8;

#ifndef INFINITY
constexpr double INFINITY = 1e32;
#endif

namespace Util
{
    template <typename MT>
    inline typename MT::vector3_type normal_direction(typename MT::vector3_type const & a, typename MT::vector3_type const & b, typename MT::vector3_type const & c);
    
    
    template <typename MT>
    inline int sign(typename MT::real_type val)
    {
        return (0. < val) - (val < 0.);
    }
    
    template <typename MT>
    inline typename MT::real_type sqr_length(std::vector<typename MT::vector3_type> const & vv, int i, int j)
    {
        return MT::sqr_length(vv[i]-vv[j]);
    }
    
    template <typename MT>
    inline typename MT::real_type length(std::vector<typename MT::vector3_type> const & vv, int i = 0, int j = 1)
    {
        return (vv[i]-vv[j]).length();
    }
    
    template <typename MT>
    inline typename MT::vector3_type grad_length(std::vector<typename MT::vector3_type> const & vv, int i)
    {
        typedef typename MT::vector3_type V;
        
        V v = vv[i]-vv[1-i];
        V l = v / vv.length();
#ifdef DEBUG
        assert(!MT::is_nan(l.length()));
#endif
        return l;
    }
    
    template <typename MT>
    inline typename MT::vector3_type grad_length(std::vector<typename MT::vector3_type> const & vv, int i, typename MT::real_type lng)
    {
        typedef typename MT::vector3_type V;
        
        V v = vv[i]-vv[1-i];
        V l = v / lng;
#ifdef DEBUG
        assert(!MT::is_nan(l.length()));
#endif
        return l;
    }
    
    /**
     * Computes the signed area of the triangle spanned by vertices with positions v0, v1 and v2.
     */
    template <typename MT>
    inline typename MT::real_type signed_area(const typename MT::vector3_type& v0, const typename MT::vector3_type& v1, const typename MT::vector3_type& v2)
    {
        typedef typename MT::vector3_type V;
        V n = MT::cross(v1-v0, v2-v0);
        return 0.5 * n.length();
    }
    
    /**
     * Computes the area of the triangle spanned by vertices with positions v0, v1 and v2.
     */
    template <typename MT>
    inline typename MT::real_type area(const typename MT::vector3_type& v0, const typename MT::vector3_type& v1, const typename MT::vector3_type& v2)
    {
        return std::abs(signed_area<MT>(v0, v1, v2));
    }
    
    /**
     * Computes the signed area of the triangle spanned by vertices with positions at i, j, k in vv.
     */
    template <typename MT>
    inline typename MT::real_type signed_area(std::vector<typename MT::vector3_type> const & vv, int i = 0, int j = 1, int k = 2)
    {
        return signed_area<MT>(vv[i], vv[j], vv[k]);
    }
    
    /**
     * Computes the area of the triangle spanned by vertices with positions at i, j, k in vv.
     */
    template <typename MT>
    inline typename MT::real_type area(std::vector<typename MT::vector3_type> const & vv, int i = 0, int j = 1, int k = 2)
    {
        return std::abs(signed_area<MT>(vv, i, j, k));
    }
    
    template <typename MT>
    inline typename MT::vector3_type grad_area(std::vector<typename MT::vector3_type> const & vv, int i)
    {
        typedef typename MT::vector3_type V;
        
        int j = (i+1)%3;
        int k = (i+2)%3;
        V n = MT::cross(vv[0]-vv[2], vv[1]-vv[2]);
        V ar = -MT::cross(n, vv[j]-vv[k])/(4.0 * area(vv));
#ifdef DEBUG
        assert(!MT::is_nan(ar.length()));
#endif
        return ar;
    }
    
    template <typename MT>
    inline typename MT::vector3_type grad_area(std::vector<typename MT::vector3_type> const & vv, int i, double a)
    {
        typedef typename MT::vector3_type V;
        
        int j = (i+1)%3;
        int k = (i+2)%3;
        V n = MT::cross(vv[0]-vv[2], vv[1]-vv[2]);
        V ar = -MT::cross(n, vv[j]-vv[k]) / (4.0 * a);
#ifdef DEBUG
        assert(!MT::is_nan(ar.length()));
#endif
        return ar;
    }
    
    template <typename MT>
    inline typename MT::real_type signed_volume(typename MT::vector3_type const & a,
                                                  typename MT::vector3_type const & b,
                                                  typename MT::vector3_type const & c,
                                                  typename MT::vector3_type const & d)
    {
        typename MT::matrix3x3_type m(a-d,c-d,b-d);
        typename MT::real_type v = MT::determinant(m);
        return v/6.0;
    }
    
    template <typename MT>
    inline typename MT::real_type volume(typename MT::vector3_type const & a,
                                         typename MT::vector3_type const & b,
                                         typename MT::vector3_type const & c,
                                         typename MT::vector3_type const & d)
    {
        return std::abs(signed_volume<MT>(a, b, c, d));
    }
    
//    template <typename MT>
//    inline typename MT::vector3_type grad_volume(std::vector<typename MT::vector3_type> const & vv, int i)
//    {
//        typedef typename MT::real_type T;
//        typedef typename MT::vector3_type V;
//        
//        int j = (i+1)%4;
//        int k = (i+2)%4;
//        int l = (i+3)%4;
//        T sgn = (i%2==1)?(1):(-1);
//        V g = sgn*cross(vv[j]-vv[l], vv[k]-vv[l]);
//        return g / 6.0;
//    }
    
    /**
     * Calculates the cosine of the angle between the line segments |ab| and |ac|.
     */
    template <typename MT>
    inline typename MT::real_type cos_angle(typename MT::vector3_type const & a,
                                            typename MT::vector3_type const & b,
                                            typename MT::vector3_type const & c)
    {
        typedef typename MT::vector3_type V;
        V ab = MT::normalize(b - a);
        V ac = MT::normalize(c - a);
        return MT::dot(ab, ac);
    }
    
    /**
     * Calculates the angle between the line segments |ab| and |ac|.
     */
    template <typename MT>
    inline typename MT::real_type angle(typename MT::vector3_type const & a,
                                        typename MT::vector3_type const & b,
                                        typename MT::vector3_type const & c)
    {
        return acos(cos_angle<MT>(a, b, c));
    }
    
    /**
     * Calculate the cosine of angles in the triangle defined by the vertices a, b and c.
     */
    template <class MT>
    inline std::vector<typename MT::real_type> cos_angles(typename MT::vector3_type const & a,
                                                          typename MT::vector3_type const & b,
                                                          typename MT::vector3_type const & c)
    {
        typedef typename MT::real_type T;
        
        std::vector<T> cosines(3);
        cosines[0] = cos_angle<MT>(a, b, c);
        cosines[1] = cos_angle<MT>(b, c, a);
        cosines[2] = cos_angle<MT>(c, a, b);
        return cosines;
    }
    
    template<typename MT>
    inline typename MT::real_type min_angle(typename MT::vector3_type const & a,
                                            typename MT::vector3_type const & b,
                                            typename MT::vector3_type const & c)
    {
        typedef typename MT::real_type T;
        std::vector<T> cosines = cos_angles<MT>(a, b, c);
        T max_cos = -1.;
        for(auto cos : cosines)
        {
            max_cos = std::max(cos, max_cos);
        }
        return acos(max_cos);
    }
    
    template<typename MT>
    inline typename MT::real_type max_angle(typename MT::vector3_type const & a,
                                            typename MT::vector3_type const & b,
                                            typename MT::vector3_type const & c)
    {
        typedef typename MT::real_type T;
        std::vector<T> cosines = cos_angles<MT>(a, b, c);
        double min_cos = 1.;
        for(auto cos : cosines)
        {
            min_cos = std::min(cos, min_cos);
        }
        return acos(min_cos);
    }
    
    /**
     * Find cosines of angles between edges connecting v with triangle's vertices and triangle's edges.
     * Used for determining whether v lies very close to the edge -- in that case we do not want to use cap or sliver removal
     * (it would introduce a very low quality tetrahedron).
     */
    template <typename MT>
    inline void cosines(typename MT::vector3_type & v,
                            std::vector<typename MT::vector3_type> & verts,
                            std::vector<typename MT::real_type> & cosines)
    {
        typedef typename MT::vector3_type V;
        
        for (int i = 0; i < 3; ++i)
        {
            int k = i;
            int l = (i+1)%3;
            V p, q;
            
            p = v - verts[k];
            p.normalize();
            q = verts[l] - verts[k];
            q.normalize();
            cosines[2*k] = MT::dot(p,q);
            
            p = v - verts[l];
            p.normalize();
            q = verts[k] - verts[l];
            q.normalize();
            cosines[2*k+1] = MT::dot(p,q);
        }
    }
    
    /**
     * Returns the cosine to the dihedral angle between face |abc| and face |abd|.
     */
    template<typename MT>
    inline typename MT::real_type cos_dihedral_angle(const typename MT::vector3_type& a, const typename MT::vector3_type& b, const typename MT::vector3_type& c, const typename MT::vector3_type& d)
    {
        typedef typename MT::real_type      T;
        typedef typename MT::vector3_type   V;
        
        V n0 = normal_direction<MT>(a, b, c);
        V n1 = normal_direction<MT>(b, a, d);
        T angle = MT::dot(n0, n1);
#ifdef DEBUG
        assert(angle < 1. + EPSILON);
        assert(angle > -1. - EPSILON);
#endif
        return angle;
    }
    
    /**
     * Returns the dihedral angle between face |abc| and face |abd|.
     */
    template<typename MT>
    inline typename MT::real_type dihedral_angle(const typename MT::vector3_type& a, const typename MT::vector3_type& b, const typename MT::vector3_type& c, const typename MT::vector3_type& d)
    {
        return MT::acos(cos_dihedral_angle<MT>(a, b, c, d));
    }
    
    template <typename MT>
    inline typename MT::vector3_type barycenter(typename MT::vector3_type const & a, typename MT::vector3_type const & b)
    {
        typedef typename MT::vector3_type V;
        V result(0.);
        result += a + b;
        result /= 2.;
        return result;
    }
    
    template <typename MT>
    inline typename MT::vector3_type barycenter(typename MT::vector3_type const & a,
                                                typename MT::vector3_type const & b,
                                                typename MT::vector3_type const & c)
    {
        typedef typename MT::vector3_type V;
        V result(0.);
        result += a + b + c;
        result /= 3.;
        return result;
    }
    
    template <typename MT>
    inline typename MT::vector3_type barycenter(typename MT::vector3_type const & a,
                                                typename MT::vector3_type const & b,
                                                typename MT::vector3_type const & c,
                                                typename MT::vector3_type const & d)
    {
        typedef typename MT::vector3_type V;
        V result(0.);
        result += a + b + c + d;
        result /= 4.;
        return result;
    }
    
    /**
     * Finds the barycentric coordinates of point v in a triangle spanned by points v0, v1, v2.
     */
    template <typename MT>
    inline std::vector<typename MT::real_type> barycentric_coords(typename MT::vector3_type const& p, typename MT::vector3_type const& a, typename MT::vector3_type const& b, typename MT::vector3_type const& c)
    {
        typedef typename MT::real_type      T;
        typedef typename MT::vector3_type   V;
        
        std::vector<T> coords(3);
        
        V v0 = b - a;
        V v1 = c - a;
        V v2 = p - a;
        T d00 = MT::dot(v0, v0);
        T d01 = MT::dot(v0, v1);
        T d11 = MT::dot(v1, v1);
        T d20 = MT::dot(v2, v0);
        T d21 = MT::dot(v2, v1);
        T denom = d00 * d11 - d01 * d01;        
#ifdef DEBUG
        assert(denom != 0.);
#endif
        coords[0] = (d11 * d20 - d01 * d21) / denom;
        coords[1] = (d00 * d21 - d01 * d20) / denom;
        coords[2] = 1. - coords[0] - coords[1];
           
        return coords;
    }
    
    /**
     * Calculates the barycentric coordinates of a point v in a tetrahedron spanned by the four vertices in verts.
     */
    template <typename MT>
    inline void barycentric_coords(const typename MT::vector3_type& v, const std::vector<typename MT::vector3_type>& verts, std::vector<typename MT::real_type> & coords)
    {        
        coords[0] = signed_volume<MT>(v       , verts[1], verts[2], verts[3]);
        coords[1] = signed_volume<MT>(verts[0], v       , verts[2], verts[3]);
        coords[2] = signed_volume<MT>(verts[0], verts[1], v       , verts[3]);
        coords[3] = signed_volume<MT>(verts[0], verts[1], verts[2], v       );
        
        typename MT::real_type s = coords[0] + coords[1] + coords[2] + coords[3];
        for (unsigned int i = 0; i < 4; ++i)
        {
            coords[i] /= s;
        }
    }
    
    template <typename MT>
    inline typename MT::vector3_type normal_direction(typename MT::vector3_type const & a,
                                                      typename MT::vector3_type const & b,
                                                      typename MT::vector3_type const & c)
    {
        typedef typename MT::vector3_type V;
        V ab = b - a;
        V ac = c - a;
        V n = MT::cross(ab, ac);
#ifdef DEBUG
        assert(!MT::is_nan(n[0]) && !MT::is_nan(n[1]) && !MT::is_nan(n[2]));
#endif
        return MT::normalize(n);
    }
    
    template <typename MT>
    inline typename MT::vector3_type normal_direction(typename MT::vector3_type const & a,
                                                      typename MT::vector3_type const & b,
                                                      typename MT::vector3_type const & c,
                                                      typename MT::vector3_type const & d)
    {
        typedef typename MT::vector3_type V;
        V n = normal_direction<MT>(a, b, c);
        V bf = barycenter<MT>(a, b, c);
        V bt = barycenter<MT>(a, b, c, d);
        V v_out = bf - bt;
        if (dot(v_out, n) > 0)
            return n;
        else
            return -n;
    }
    
    /**
     Returns v projected onto the line spanned by the two points v1 and v2.
     */
    template <typename MT>
    inline typename MT::vector3_type project(typename MT::vector3_type const & v, typename MT::vector3_type const & v1,typename MT::vector3_type const & v2)
    {
        typedef typename MT::vector3_type V;
        V a = v - v1;
        V b = v2 - v1;
        return v1 + b * MT::dot(a,b)/MT::dot(b, b);
    }
    
    /**
     * Project the point v onto the plane spanned by the three points in verts.
     */
    template<typename MT>
    inline typename MT::vector3_type project(typename MT::vector3_type const & v, const std::vector<typename MT::vector3_type>& verts)
    {
        typename MT::vector3_type normal = Util::normal_direction<MT>(verts[0], verts[1], verts[2]);
        return v - normal * MT::dot(v - verts[0], normal);
    }
    
    template <typename MT>
    inline typename MT::real_type calc_flatness(typename MT::vector3_type const & a,
                                                typename MT::vector3_type const & b,
                                                typename MT::vector3_type const & c,
                                                typename MT::vector3_type const & d)
    {
        typedef typename MT::vector3_type V;
        V normal0 = normal_direction<MT>(d, a, b);
        V normal1 = normal_direction<MT>(c, b, a);
        
        return MT::dot(normal0, normal1);
    }
    
    
    
    template <typename MT>
    inline typename MT::real_type ms_length(std::vector<typename MT::vector3_type> const & vv)
    {
        typedef typename MT::real_type T;
        
        T result = 0.0;
        result += sqr_length<MT>(vv, 0, 1);
        result += sqr_length<MT>(vv, 0, 2);
        result += sqr_length<MT>(vv, 0, 3);
        result += sqr_length<MT>(vv, 1, 2);
        result += sqr_length<MT>(vv, 1, 3);
        result += sqr_length<MT>(vv, 2, 3);
        return result / 6.0;
    }
    
    template <typename MT>
    inline typename MT::real_type rms_length(typename MT::vector3_type const a, typename MT::vector3_type const b, typename MT::vector3_type const c, typename MT::vector3_type const d)
    {        
        std::vector<typename MT::vector3_type> verts(4);
        verts[0] = a;
        verts[1] = b;
        verts[2] = c;
        verts[3] = d;
        return sqrt(ms_length<MT>(verts));
    }
    
    template <typename MT>
    inline typename MT::real_type rms_length(const std::vector<typename MT::vector3_type>& verts)
    {
        return sqrt(ms_length<MT>(verts));
    }
    
    template <typename MT>
    inline typename MT::vector3_type grad_rms_length(std::vector<typename MT::vector3_type> const & vv, int i)
    {
        typedef typename MT::real_type    T;
        typedef typename MT::vector3_type V;
        
        std::vector<V> ve(4);
        T const rmsl = rms_length<MT>(vv);
        V result(0.0);
        ve[0] = vv[i];
        
        ve[1] = vv[(i+1)%4];
        T l = length<MT>(ve);
        result += l * grad_length<MT>(vv, 0, l);
        
        ve[1] = vv[(i+2)%4];
        l = length<MT>(ve);
        result += l * grad_length<MT>(vv, 0, l);
        
        ve[1] = vv[(i+3)%4];
        l = length<MT>(ve);
        result += l * grad_length<MT>(vv, 0, l);
        
        result /= rmsl;
#ifdef DEBUG
        assert(!MT::is_nan(result.length()) || !"????");
#endif
        return result;
    }
    
    template<typename MT>
    inline typename MT::real_type quality(typename MT::vector3_type const a, typename MT::vector3_type const b, typename MT::vector3_type const c, typename MT::vector3_type const d)
    {
        typedef typename MT::real_type      T;
        
        T v = Util::signed_volume<MT>(a, b, c, d);
        T lrms = rms_length<MT>(a, b, c, d);
        
        T q = 8.48528 * v / (lrms * lrms * lrms);
#ifdef DEBUG
        assert(!MT::is_nan(q));
#endif
        return q;
    }
    
    
    /**
     * Finds the center of a smallest circle containing the triangle specified by vertices a, b, c.
     * For an acute or right triangle, this is the circumcircle. For an obtuse triangle this is the midpoint of the longest edge.
     */
//    template <typename MT>
//    inline typename MT::vector3_type min_circle_center(typename MT::vector3_type & a,
//                                                       typename MT::vector3_type & b,
//                                                       typename MT::vector3_type & c)
//    {
//        typedef typename MT::real_type    T;
//        typedef typename MT::vector3_type V;
//        
//        V eba = b-a,
//        eca = c-a,
//        ecb = c-b;
//        
//        T c2 = MT::sqr_length(eba),
//        b2 = MT::sqr_length(eca),
//        a2 = MT::sqr_length(ecb);
//        
//        T alpha = a2 * (b2 + c2 - a2);
//        T beta  = b2 * (a2 + c2 - b2);
//        T gamma = c2 * (a2 + b2 - c2);
//        
//        T sum = alpha + beta + gamma;
//        alpha /= sum;	beta /= sum;	gamma /= sum;
//        
//        if (alpha <= 0)
//            return (b+c)/2.0;
//        if (beta <= 0)
//            return (a+c)/2.0;
//        if (gamma <= 0)
//            return (a+b)/2.0;
//        
//        return alpha * a + beta * b + gamma * c;
//    }
    
    /**
     * Computes the determinant of a 4-by-4 matrix specified by four 4D vectors a, b, c, d
     */
    template <typename MT>
    inline typename MT::real_type determinant(typename MT::vector4_type & a,
                                              typename MT::vector4_type & b,
                                              typename MT::vector4_type & c,
                                              typename MT::vector4_type & d)
    {
        typedef typename MT::real_type T;
        
        T a1, a2, a3, a4, b1, b2, b3, b4, c1, c2, c3, c4, d1, d2, d3, d4;
        
        a1 = a[0];	b1 = b[0];	c1 = c[0];	d1 = d[0];
        a2 = a[1];	b2 = b[1];	c2 = c[1];	d2 = d[1];
        a3 = a[2];	b3 = b[2];	c3 = c[2];	d3 = d[2];
        a4 = a[3];	b4 = b[3];	c4 = c[3];	d4 = d[3];
        
        return   a1 * (b2*(c3*d4-d3*c4)-c2*(b3*d4-d3*b4)+d2*(b3*c4-c3*b4))
        - b1 * (a2*(c3*d4-d3*c4)-c2*(a3*d4-d3*a4)+d2*(a3*c4-c3*a4))
        + c1 * (a2*(b3*d4-d3*b4)-b2*(a3*d4-d3*a4)+d2*(a3*b4-b3*a4))
        - d1 * (a2*(b3*c4-c3*b4)-b2*(a3*c4-c3*a4)+c2*(a3*b4-b3*a4));
    }
    
    template <typename MT>
    inline std::vector<typename MT::vector3_type> find_basis(std::vector<typename MT::vector3_type> & basis, std::vector<typename MT::vector3_type> & points)
    {
        typedef typename MT::real_type      T;
        typedef typename MT::vector3_type   V;
        V vp = points[0];
        std::vector<V> b,m;
        for (int i = 1; i < points.size(); ++i)
        {
            m.push_back(points[i]);
        }
        
        if (points.size() == 1)
        {
            if (basis.size() == 0)
            {
                return points;
            }
            for (unsigned int i = 0; i < basis.size(); ++i)
            {
                b.push_back(basis[i]);
            }
        }
        else
            b = find_basis<MT>(basis, m);
        
        if (b.size() == 1)
        {
            V vq = b[0];
            if (MT::dot(vq,vp-vq) >= 0)
            {
                return b;
            }
        }
        else if (b.size() == 2)
        {
            V vq = b[0];
            V vr = b[1];
            V vs = vp - vr;
            V vt = vq - vr;
            if (MT::dot(MT::cross(vs,vt),MT::cross(vr,vt)) >= 0)
            {
                return b;
            }
        }
        else if (b.size() == 3)
        {
            V vq = b[0];
            V vr = b[1];
            V vs = b[2];
            if (signed_volume<MT>(vp, vq, vr, vs) * signed_volume<MT>(V(0.0), vq, vr, vs) <= 0)
            {
                return b;
            }
        }
        else
        {
            return b;
        }
        
        basis.push_back(vp);
        if (points.size() == 1 || basis.size() == 3)
        {
            return basis;
        }
        else
        {
            return find_basis<MT>(basis, m);
        }
    }
    
    /**
     * Finds the minimum convex hull point.
     */
    template <typename MT>
    inline typename MT::vector3_type min_convex_hull_point(std::vector<typename MT::vector3_type> & points)
    {
        typedef typename MT::real_type      T;
        typedef typename MT::vector3_type   V;
        
        std::vector<V> basis;
        std::vector<V> b = find_basis<MT>(basis, points);
        if (b.size() == 1)
        {
            return b[0];
        }
        else if (b.size() == 2)
        {
            V vp = b[0];
            V vq = b[1];
            return vq - (vp-vq)*(MT::dot(vq,vp-vq)/MT::sqr_length(vp-vq));
        }
        else if (b.size() == 3)
        {
            V vp = b[0];
            V vq = b[1];
            V vr = b[2];
            V vs = vp-vr;
            V vt = vq-vr;
            return vr - (MT::dot(MT::cross(vs,vt), MT::cross(vr,vt))/MT::dot(MT::cross(vs,vt),MT::cross(vs,vt)))*vs - (MT::dot(MT::cross(vs,vt),MT::cross(vs,vr))/MT::dot(MT::cross(vs,vt),MT::cross(vs,vt)))*vt;
        }
        else
            return V(0.0);
    }
    
    /**
     * Returns the shortest distance from the point p to the plane spanned by the points a, b and c.
     */
    template<typename MT>
    inline typename MT::real_type distance(const typename MT::vector3_type& p, const typename MT::vector3_type& a, const typename MT::vector3_type& b, const typename MT::vector3_type& c)
    {
        typedef typename MT::vector3_type   V;
        
        V v = p - a;
        V n = normal_direction<MT>(a, b, c);
        
        return std::abs(MT::dot(v, n));
    }
    
    /**
     Returns whether you have to turn left when going from a to b to c.
     */
    template<typename MT>
    inline bool is_left_of(const typename MT::vector3_type& a, const typename MT::vector3_type& b, const typename MT::vector3_type& c)
    {
        if(signed_area<MT>(a, b, c) > 0.)
        {
            return true;
        }
        return false;
    }
    
    template<typename MT>
    inline bool is_between(const typename MT::vector3_type& p, const std::vector<typename MT::vector3_type>& verts)
    {
        bool is_l1 = is_left_of<MT>(verts[0], verts[1], p);
        bool is_l2 = is_left_of<MT>(verts[1], verts[2], p);
        bool is_l3 = is_left_of<MT>(verts[2], verts[0], p);
        return (is_l1 && is_l2 && is_l3) | (!is_l1 && !is_l2 && !is_l3);
    }
    
    /**
     * Calculates the intersection between the line segment |p0 p1| and the plane spanned by the vertices v0, v1 and v2. The intersection point is defined by p0 + t*(p1 - p0) and the function returns t. Returns infinity if it does not intersect.
     */
    template<typename MT>
    typename MT::real_type intersection_ray_plane(const typename MT::vector3_type& p0, const typename MT::vector3_type& p1, const typename MT::vector3_type& v0, const typename MT::vector3_type& v1, const typename MT::vector3_type& v2)
    {
        typedef typename MT::real_type      T;
        typedef typename MT::vector3_type   V;
        
        V normal = normal_direction<MT>(v0, v1, v2);
        T n = MT::dot(normal, v0 - p0);
        T d = MT::dot(normal, p1 - p0);
        
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
     * Calculates the intersection between the line segment |p0 p1| and the triangle |v0 v1 v2|. The intersection point is defined by p0 + t*(p1 - p0) and the function returns t. Returns infinity if it does not intersect.
     */
    template<typename MT>
    typename MT::real_type intersection_ray_triangle(const typename MT::vector3_type& p0, const typename MT::vector3_type& p1, const typename MT::vector3_type& v0, const typename MT::vector3_type& v1, const typename MT::vector3_type& v2)
    {
        typedef typename MT::real_type      T;
        typedef typename MT::vector3_type   V;
        
        T t = intersection_ray_plane<MT>(p0, p1, v0, v1, v2);
        if(t < 0.) // The ray goes away from the triangle
        {
            return t;
        }
        V p = p0 + t*(p1 - p0);
        
        std::vector<T> coords = barycentric_coords<MT>(p, v0, v1, v2);
        if(coords[0] > EPSILON && coords[1] > EPSILON && coords[2] > EPSILON) // The intersection happens inside the triangle.
        {
            return t;
        }
        return INFINITY; // The intersection happens outside the triangle.
    }
    
    /**
     * Implies ordering in the space of binary vectors of given size.
     *
     * @param n     Size of binary vectors.
     * @param v1    First vector.
     * @param v2    Second vector.
     *
     * @return      True if the first vector is smaller than the second, false otherwise.
     */
    inline bool compare(int n, const std::vector<bool> & v1, const std::vector<bool> & v2)
    {
        for (int i = 0; i < n; ++i)
        {
            if (v1[i] && !v2[i]) return false;
            if (!v1[i] && v2[i]) return true;
        }
        return false;
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
    
}
