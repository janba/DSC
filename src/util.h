#ifndef UTIL_H
#define UTIL_H

#include <vector>
#include <list>
#include <sstream>
#include <cmath>
#include <cassert>


#define EPSILON 1e-12

namespace Util
{
    
    template <typename MT>
    inline int sgn(typename MT::real_type val)
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
        assert(!MT::is_nan(l.length()));
        return l;
    }
    
    template <typename MT>
    inline typename MT::vector3_type grad_length(std::vector<typename MT::vector3_type> const & vv, int i, typename MT::real_type lng)
    {
        typedef typename MT::vector3_type V;
        
        V v = vv[i]-vv[1-i];
        V l = v / lng;
        assert(!MT::is_nan(l.length()));
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
        assert(!MT::is_nan(ar.length()));
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
        assert(!MT::is_nan(ar.length()));
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
//    inline typename MT::real_type volume(std::vector<typename MT::vector3_type> const & vv)
//    {
//        return std::abs(signed_volume<MT>(vv[0], vv[1], vv[2], vv[3]));
//    }
    
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
        double max_cos = -1.;
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
     * Calculate the maximum cosine of angles in the triangle defined by the a, b and c. The index parameter tells which angle (a, b or c) is maximal.
     */
//    template <class MT>
//    inline typename MT::real_type max_cos_angle(typename MT::vector3_type const & a,
//                                                typename MT::vector3_type const & b,
//                                                typename MT::vector3_type const & c, int& index)
//    {
//        typedef typename MT::real_type T;
//        T max_cos = -1.;
//        std::vector<T> cosines = cos_angles<MT>(a, b, c);
//        for (int i = 0; i < cosines.size(); i++)
//        {
//            if(cosines[i] > max_cos)
//            {
//                max_cos = cosines[i];
//                index = i;
//            }
//        }
//        return max_cos;
//        
//    }
    
    /**
     * Find cosines of angles between edges connecting v with triangle's vertices and triangle's edges.
     * Used for determining whether v lies very close to the edge -- in that case we do not want to use cap or sliver removal
     * (it would introduce a very low quality tetrahedron).
     */
    template <typename MT>
    inline void get_cosines(typename MT::vector3_type & v,
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
    
//    template <typename MT>
//    inline typename MT::real_type cos_dihedral_angle(typename MT::vector3_type const & a,
//                                                     typename MT::vector3_type const & b,
//                                                     typename MT::vector3_type const & c,
//                                                     typename MT::vector3_type const & d)
//    {
//        typedef typename MT::vector3_type V;
//        
//        V n0 = normal_direction(a, b, c, d);
//        V n1 = normal_direction(a, b, d, c);
//        
//        return -MT::dot(n0, n1);
//    }
//    
//    template <typename MT>
//    inline typename MT::real_type sin_dihedral_angle(typename MT::vector3_type const & a,
//                                                     typename MT::vector3_type const & b,
//                                                     typename MT::vector3_type const & c,
//                                                     typename MT::vector3_type const & d)
//    {
//        typedef typename MT::real_type T;
//        
//        T cda = cos_dihedral_angle(a, b, c, d);
//        return sqrt(1 - cda * cda);
//    }
//    
//    template <typename MT>
//    inline std::list<typename MT::real_type> sin_dihedral_angles(typename MT::vector3_type const & a,
//                                                                 typename MT::vector3_type const & b,
//                                                                 typename MT::vector3_type const & c,
//                                                                 typename MT::vector3_type const & d)
//    {
//        typedef typename MT::real_type T;
//        
//        std::list<T> result;
//        
//        result.push_back(sin_dihedral_angle(a, b, c, d));
//        result.push_back(sin_dihedral_angle(a, c, b, d));
//        result.push_back(sin_dihedral_angle(a, d, b, c));
//        result.push_back(sin_dihedral_angle(b, c, a, d));
//        result.push_back(sin_dihedral_angle(b, d, a, c));
//        result.push_back(sin_dihedral_angle(c, d, a, b));
//        
//        result.sort();
//        return result;
//    }
    
    //	template <class VT>
    //	double min_sin_dihedral_angle(VT const & a, VT const & b, VT const & c, VT const & d)
    //	{
    //		double result, t;
    //
    //		t = sin_dihedral_angle(a, b, c, d);
    //		result = t;
    //		t = sin_dihedral_angle(a, c, b, d);
    //		if (t < result) result = t;
    //		t = sin_dihedral_angle(a, d, b, c);
    //		if (t < result) result = t;
    //		t = sin_dihedral_angle(b, c, a, d);
    //		if (t < result) result = t;
    //		t = sin_dihedral_angle(b, d, a, c);
    //		if (t < result) result = t;
    //		t = sin_dihedral_angle(c, d, a, b);
    //		if (t < result) result = t;
    //
    //		return result;
    //	}
    
//    
//    template <typename MT>
//    inline typename MT::real_type sin_dihedral_angle(std::vector<typename MT::vector3_type> const & vv, int i, int j, int k, int l)
//    {
//        typedef typename MT::real_type T;
//        typedef typename MT::vector3_type V;
//        
//        std::vector<V> edge(2);
//        edge[0] = vv[i];
//        edge[1] = vv[j];
//        T lng = length(edge);
//        
//        std::vector<V> face_k(3);
//        face_k[0] = vv[i];
//        face_k[1] = vv[j];
//        face_k[2] = vv[l];
//        T a_k = area(face_k);
//        
//        std::vector<V> face_l(3);
//        face_l[0] = vv[i];
//        face_l[1] = vv[j];
//        face_l[2] = vv[k];
//        T a_l = area(face_l);
//        
//        T v = volume(vv);
//        
//        T a = 0.f;
//        if (a_k != 0 && a_l != 0)
//            a = 1.5 * v * lng / (a_k * a_l);
//        
//        assert(!MT::is_nan(a));
//        return a;
//    }
//    
//    template <typename MT>
//    inline typename MT::vector3_type grad_sin_dihedral_angle(std::vector<typename MT::vector3_type> const & vv, int i, int j, int alpha)
//    {
//        typedef typename MT::real_type T;
//        typedef typename MT::vector3_type V;
//        
//        int k,l,n;
//        for (n=0; n<4; ++n)
//            if ((n!=i) && (n!=j))
//            {
//                k = n;
//                ++n;
//                break;
//            }
//        for ( ; n<4; ++n)
//            if ((n!=i) && (n!=j))
//            {
//                l = n;
//                break;
//            }
//        
//        T s = sin_dihedral_angle(vv, i, j, k, l);
//        
//        T v = volume(vv);
//        V g1 = grad_volume(vv,alpha);
//        
//        V g2(0.0);
//        int beta = -1;
//        if (i == alpha)
//            beta = 0;
//        else if (j == alpha)
//            beta = 1;
//        
//        std::vector<V> edge(2);
//        edge[0] = vv[i];
//        edge[1] = vv[j];
//        T lng = length(edge);
//        
//        if (beta > -1)
//        {
//            g2 = grad_length(edge,beta);
//        }
//        
//        V g3(0.0);
//        beta = -1;
//        if (i == alpha)
//            beta = 0;
//        else if (j == alpha)
//            beta = 1;
//        else if (k == alpha)
//            beta = 2;
//        
//        std::vector<V> face_l(3);
//        face_l[0] = vv[i];
//        face_l[1] = vv[j];
//        face_l[2] = vv[k];
//        T a_l = area(face_l);
//        
//        if (beta > -1)
//        {
//            g3 = grad_area(face_l,beta);
//        }
//        
//        V g4(0.0);
//        beta = -1;
//        if (i == alpha)
//            beta = 0;
//        else if (j == alpha)
//            beta = 1;
//        else if (l == alpha)
//            beta = 2;
//        
//        std::vector<V> face_k(3);
//        face_k[0] = vv[i];
//        face_k[1] = vv[j];
//        face_k[2] = vv[l];
//        T a_k = area(face_k);
//        
//        if (beta > -1)
//        {
//            g4 = grad_area(face_k, beta);
//        }
//        
//        return s * (g1/v + g2/lng - g3/a_l - g4/a_k);
//    }
//    
//    template <typename MT>
//    inline typename MT::vector3_type grad_sin_dihedral_angle(std::vector<typename MT::vector3_type> const & vv,
//                                                             int i, int j, int alpha, typename MT::real_type vol,
//                                                             typename MT::real_type ak,
//                                                             typename MT::real_type al,
//                                                             typename MT::real_type lij)
//    {
//        typedef typename MT::real_type    T;
//        typedef typename MT::vector3_type V;
//        
//        int k = -1, l = -1, n = -1;
//        for (n=0; n<4; ++n)
//            if ((n!=i) && (n!=j))
//            {
//                k = n;
//                ++n;
//                break;
//            }
//        for ( ; n<4; ++n)
//            if ((n!=i) && (n!=j))
//            {
//                l = n;
//                break;
//            }
//        
//        T s = 1.5 * vol * lij / (ak * al);
//        assert(!MT::is_nan(s));
//        
//        V g1 = grad_volume(vv,alpha);
//        
//        V g2(0.0);
//        int beta = -1;
//        if (i == alpha)
//            beta = 0;
//        else if (j == alpha)
//            beta = 1;
//        
//        std::vector<V> edge(2);
//        edge[0] = vv[i];
//        edge[1] = vv[j];
//        
//        if (beta > -1)
//        {
//            g2 = grad_length(edge,beta,lij);
//        }
//        
//        V g3(0.0);
//        beta = -1;
//        if (i == alpha)
//            beta = 0;
//        else if (j == alpha)
//            beta = 1;
//        else if (k == alpha)
//            beta = 2;
//        
//        std::vector<V> face_l(3);
//        face_l[0] = vv[i];
//        face_l[1] = vv[j];
//        face_l[2] = vv[k];
//        
//        if (beta > -1)
//        {
//            g3 = grad_area(face_l,beta,al);
//        }
//        
//        V g4(0.0);
//        beta = -1;
//        if (i == alpha)
//            beta = 0;
//        else if (j == alpha)
//            beta = 1;
//        else if (l == alpha)
//            beta = 2;
//        
//        std::vector<V> face_k(3);
//        face_k[0] = vv[i];
//        face_k[1] = vv[j];
//        face_k[2] = vv[l];
//        
//        if (beta > -1)
//        {
//            g4 = grad_area(face_k,beta,ak);
//        }
//        
//        return s * (g1/vol + g2/lij - g3/al - g4/ak);
//    }
    
    
    //	template <class VT>
    //	double min_sin_dihedral_angle(std::vector<VT> const & vv, int & i, int & j)
    //	{
    //		double min_s = 2.0;
    //		double v = volume(vv);
    //
    //		if (v < -0.0001) return -1.0;
    //
    //		double a0 = area(vv, 1, 2, 3),
    //			a1 = area(vv, 0, 2, 3),
    //			a2 = area(vv, 0, 1, 3),
    //			a3 = area(vv, 0, 1, 2);
    //		double e01 = length(vv, 0, 1),
    //			e02 = length(vv, 0, 2),
    //			e03 = length(vv, 0, 3),
    //			e12 = length(vv, 1, 2),
    //			e13 = length(vv, 1, 3),
    //			e23 = length(vv, 2, 3);
    //
    //		double s = 1.5 * v * e01 / (a2 * a3);
    //		if (s < min_s)
    //		{
    //			min_s = s;
    //			i = 0;
    //			j = 1;
    //		}
    //
    //		s = 1.5 * v * e02 / (a1 * a3);
    //		if (s < min_s)
    //		{
    //			min_s = s;
    //			i = 0;
    //			j = 2;
    //		}
    //
    //		s = 1.5 * v * e03 / (a1 * a2);
    //		if (s < min_s)
    //		{
    //			min_s = s;
    //			i = 0;
    //			j = 3;
    //		}
    //
    //		s = 1.5 * v * e12 / (a0 * a3);
    //		if (s < min_s)
    //		{
    //			min_s = s;
    //			i = 1;
    //			j = 2;
    //		}
    //
    //		s = 1.5 * v * e13 / (a0 * a2);
    //		if (s < min_s)
    //		{
    //			min_s = s;
    //			i = 1;
    //			j = 3;
    //		}
    //
    //		s = 1.5 * v * e23 / (a0 * a1);
    //		if (s < min_s)
    //		{
    //			min_s = s;
    //			i = 2;
    //			j = 3;
    //		}
    //
    //		return min_s;
    //	}
    
//    template <typename MT>
//    inline typename MT::real_type min_biased_sin_dihedral_angle(std::vector<typename MT::vector3_type> const & vv,
//                                                                int & i, int & j, typename MT::real_type bias)
//    {
//        typedef typename MT::real_type    T;
//        typedef typename MT::vector3_type V;
//        
//        T min_s = 2.0;
//        T v = volume(vv);
//        T a0 = area(vv, 1, 2, 3),
//        a1 = area(vv, 0, 2, 3),
//        a2 = area(vv, 0, 1, 3),
//        a3 = area(vv, 0, 1, 2);
//        T e01 = length(vv, 0, 1),
//        e02 = length(vv, 0, 2),
//        e03 = length(vv, 0, 3),
//        e12 = length(vv, 1, 2),
//        e13 = length(vv, 1, 3),
//        e23 = length(vv, 2, 3);
//        V n0 = normal_direction(vv[1], vv[2], vv[3], vv[0]),
//        n1 = normal_direction(vv[0], vv[2], vv[3], vv[1]),
//        n2 = normal_direction(vv[0], vv[1], vv[3], vv[2]),
//        n3 = normal_direction(vv[0], vv[1], vv[2], vv[3]);
//        
//        T s = 1.5 * v * e01 / (a2 * a3);
//        assert(!MT::is_nan(s));
//        if (dot(n2, n3) > 0) s *= bias;
//        if (s < min_s)
//        {
//            min_s = s;
//            i = 0;
//            j = 1;
//        }
//        
//        s = 1.5 * v * e02 / (a1 * a3);
//        assert(!MT::is_nan(s));
//        if (dot(n1, n3) > 0) s *= bias;
//        if (s < min_s)
//        {
//            min_s = s;
//            i = 0;
//            j = 2;
//        }
//        
//        s = 1.5 * v * e03 / (a1 * a2);
//        assert(!MT::is_nan(s));
//        if (dot(n1, n2) > 0) s *= bias;
//        if (s < min_s)
//        {
//            min_s = s;
//            i = 0;
//            j = 3;
//        }
//        
//        s = 1.5 * v * e12 / (a0 * a3);
//        assert(!MT::is_nan(s));
//        if (dot(n0, n3) > 0) s *= bias;
//        if (s < min_s)
//        {
//            min_s = s;
//            i = 1;
//            j = 2;
//        }
//        
//        s = 1.5 * v * e13 / (a0 * a2);
//        assert(!MT::is_nan(s));
//        if (dot(n0, n2) > 0) s *= bias;
//        if (s < min_s)
//        {
//            min_s = s;
//            i = 1;
//            j = 3;
//        }
//        
//        s = 1.5 * v * e23 / (a0 * a1);
//        assert(!MT::is_nan(s));
//        if (dot(n0, n1) > 0) s *= bias;
//        if (s < min_s)
//        {
//            min_s = s;
//            i = 2;
//            j = 3;
//        }
//        
//        return min_s;
//    }
//    
//    
//    template <typename MT>
//    inline typename MT::vector3_type grad_biased_sin_dihedral_angle(std::vector<typename MT::vector3_type> const & vv,
//                                                                    int i, int j, int alpha,
//                                                                    typename MT::real_type vol,
//                                                                    typename MT::real_type ak,
//                                                                    typename MT::real_type al,
//                                                                    typename MT::real_type lij,
//                                                                    typename MT::real_type bias)
//    {
//        typedef typename MT::real_type    T;
//        typedef typename MT::vector3_type V;
//        
//        int k = -1, l = -1, n = -1;
//        for (n=0; n<4; ++n)
//            if ((n!=i) && (n!=j))
//            {
//                k = n;
//                ++n;
//                break;
//            }
//        for ( ; n<4; ++n)
//            if ((n!=i) && (n!=j))
//            {
//                l = n;
//                break;
//            }
//        
//        V ni = normal_direction(vv[j], vv[k], vv[l], vv[i]),
//        nj = normal_direction(vv[i], vv[k], vv[l], vv[j]);
//        
//        T s = 1.5 * vol * lij / (ak * al);
//        assert(!MT::is_nan(s));
//        
//        V g1 = grad_volume(vv,alpha);
//        
//        V g2(0.0);
//        int beta = -1;
//        if (i == alpha)
//            beta = 0;
//        else if (j == alpha)
//            beta = 1;
//        
//        std::vector<V> edge(2);
//        edge[0] = vv[i];
//        edge[1] = vv[j];
//        
//        if (beta > -1)
//        {
//            g2 = grad_length(edge,beta,lij);
//        }
//        
//        V g3(0.0);
//        beta = -1;
//        if (i == alpha)
//            beta = 0;
//        else if (j == alpha)
//            beta = 1;
//        else if (k == alpha)
//            beta = 2;
//        
//        std::vector<V> face_l(3);
//        face_l[0] = vv[i];
//        face_l[1] = vv[j];
//        face_l[2] = vv[k];
//        
//        if (beta > -1)
//        {
//            g3 = grad_area(face_l,beta,al);
//        }
//        
//        V g4(0.0);
//        beta = -1;
//        if (i == alpha)
//            beta = 0;
//        else if (j == alpha)
//            beta = 1;
//        else if (l == alpha)
//            beta = 2;
//        
//        std::vector<V> face_k(3);
//        face_k[0] = vv[i];
//        face_k[1] = vv[j];
//        face_k[2] = vv[l];
//        
//        if (beta > -1)
//        {
//            g4 = grad_area(face_k,beta,ak);
//        }
//        
//        T b = 1.0;
//        if (dot(ni, nj) > 0)
//            b = bias;
//        
//        return b * s * (g1/vol + g2/lij - g3/al - g4/ak);
//    }
    
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
    inline std::vector<typename MT::real_type> barycentric_coords(typename MT::vector3_type const & v, typename MT::vector3_type const & v0, typename MT::vector3_type const & v1, typename MT::vector3_type const & v2)
    {
        typedef typename MT::real_type      T;
        typedef typename MT::vector3_type   V;
        
        std::vector<T> coords(3);
        double scale = (v1[1] - v2[1])*(v0[0] - v2[0]) + (v2[0] - v1[0])*(v0[1] - v2[1]);
#ifdef DEBUG
        assert(scale != 0.);
#endif
        scale = 1./scale;
        coords[0] = ((v1[1] - v2[1])*(v[0] - v2[0]) + (v2[0] - v1[0])*(v[1] - v2[1])) * scale;
        coords[1] = ((v2[1] - v0[1])*(v[0] - v2[0]) + (v0[0] - v2[0])*(v[1] - v2[1])) * scale;
        coords[2] = 1. - coords[0] - coords[1];        
        return coords;
    }
    
    /**
     * Calculates the barycentric coordinates of a point v in a tetrahedron spanned by the four vertices in verts.
     */
    template <typename MT>
    inline void get_barycentric_coords(const typename MT::vector3_type& v, const std::vector<typename MT::vector3_type>& verts, std::vector<typename MT::real_type> & coords)
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
    
    
    
    //template <class VT>
    //inline void get_barycentric_coords(VT const & v, std::vector<VT> const & verts, std::vector<double> & coords)
    //{
    //    std::vector<CGLA::Vec4d> v4(4);
    //    CGLA::Vec3d c = 0.25 * (verts[0] + verts[1] + verts[2] + verts[3]);
    //#ifndef max
    //    double d = std::max(std::max(length(verts[0]-c),length(verts[1]-c)),
    //                        std::max(length(verts[2]-c),length(verts[3]-c)));
    //#else
    //	double d = max(max(length(verts[0]-c),length(verts[1]-c)),
    //                   max(length(verts[2]-c),length(verts[3]-c)));
    //#endif
    //    for (int i = 0; i < 4; ++i)
    //        v4[i] = CGLA::Vec4d((verts[i][0]-c[0])/d, (verts[i][1]-c[1])/d, (verts[i][2]-c[2])/d, 1.0);
    //    CGLA::Mat4x4d m(v4[0], v4[1], v4[2], v4[3]);
    //    //******* problems here
    //    if (abs(CGLA::determinant(m)) < 1e-6)
    //    {
    //        for (int i = 0; i < 4; ++i)
    //            coords[i] = -0.25;
    //    }
    //    else
    //    {
    //        CGLA::Mat4x4d m_inv = CGLA::invert(CGLA::transpose(m));
    //        CGLA::Vec4d b((v[0]-c[0])/d, (v[1]-c[1])/d, (v[2]-c[2])/d, 1.0), x;
    //        x = m_inv * b;
    //
    //        for (int i = 0; i < 4; ++i)
    //            coords[i] = x[i];
    //    }
    //}
    
    template <typename MT>
    inline typename MT::vector3_type normal_direction(typename MT::vector3_type const & a,
                                                      typename MT::vector3_type const & b,
                                                      typename MT::vector3_type const & c)
    {
        typedef typename MT::vector3_type V;
        V ab = b - a;
        V ac = c - a;
        V n = MT::cross(ab, ac);
        assert(!MT::is_nan(n[0]) && !MT::is_nan(n[1]) && !MT::is_nan(n[2]));
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
        assert(!MT::is_nan(result.length()) || !"????");
        return result;
    }
    
    
//    template<typename MT>
//    inline typename MT::real_type quality(std::vector<typename MT::vector3_type> & vv)
//    {
//        typedef typename MT::real_type      T;
//        
//        T v = Util::volume<MT>(vv);
//        T lrms = rms_length<MT>(vv);
//        
//        T q = 8.48528 * v / (lrms * lrms * lrms);
//        assert(!MT::is_nan(q));
//        return q;
//    }
    
    template<typename MT>
    inline typename MT::real_type quality(typename MT::vector3_type const a, typename MT::vector3_type const b, typename MT::vector3_type const c, typename MT::vector3_type const d)
    {
        typedef typename MT::real_type      T;
        
        T v = Util::signed_volume<MT>(a, b, c, d);
        T lrms = rms_length<MT>(a, b, c, d);
        
        T q = 8.48528 * v / (lrms * lrms * lrms);
        assert(!MT::is_nan(q));
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
     * Calculates the intersection between the line segment |p0 p1| and the plane spanned by the vertices in verts. The intersection point is defined by p0 + t*(p1 - p0) and the function returns t. Returns infinity if it does not intersect.
     */
    template<typename MT>
    typename MT::real_type intersection_ray_plane(const typename MT::vector3_type& p0, const typename MT::vector3_type& p1, const typename MT::vector3_type& v0, const typename MT::vector3_type& v1, const typename MT::vector3_type& v2)
    {
        typedef typename MT::real_type      T;
        typedef typename MT::vector3_type   V;
        
        V normal = normal_direction<MT>(v0, v1, v2);
        
        V ray = p1 - p0;
        T n = MT::dot(normal, v0 - p0);
        T d = MT::dot(normal, ray);
        
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
     * Calculates the intersection between the line segment |p0 p1| and the triangle made up by the vertices in verts. The intersection point is defined by p0 + t*(p1 - p0) and the function returns t. Returns infinity if it does not intersect.
     */
    template<typename MT>
    typename MT::real_type intersection(const typename MT::vector3_type& p0, const typename MT::vector3_type& p1, const std::vector<typename MT::vector3_type>& verts)
    {
        typedef typename MT::real_type      T;
        typedef typename MT::vector3_type   V;
        
        T t = intersection_ray_plane<MT>(p0, p1, verts[0], verts[1], verts[2]);
        if(t < 0.) // The ray goes away from the triangle
        {
            return t;
        }
        V p = p0 + t*(p1 - p0);
        
        std::vector<T> coords = barycentric_coords<MT>(p, verts[0], verts[1], verts[2]);
        if(coords[0] >= 0. && coords[1] >= 0. && coords[2] >= 0.) // The intersection happens inside the triangle.
        {
            return t;
        }
        return INFINITY; // The intersection happens outside the triangle.
    }
    
    
    // Copyright 2001 softSurfer, 2012 Dan Sunday
    // This code may be freely used and modified for any purpose
    // providing that this copyright notice is included with it.
    // SoftSurfer makes no warranty for this code, and cannot be held
    // liable for any real or imagined damage resulting from its use.
    // Users of this code must verify correctness for their application.
    
    // intersect3D_RayTriangle(): find the 3D intersection of a ray with a triangle
    //    Input:  a ray R, and a triangle T
    //    Output: *I = intersection point (when it exists)
    //    Return: -1 = triangle is degenerate (a segment or point)
    //             0 =  disjoint (no intersect)
    //             1 =  intersect in unique point I1
    //             2 =  are in the same plane
    template<typename MT>
    typename MT::real_type intersection(const typename MT::vector3_type& p0, const typename MT::vector3_type& p1, const typename MT::vector3_type& v0, const typename MT::vector3_type& v1, const typename MT::vector3_type& v2)
    {
        typedef typename MT::real_type      T;
        typedef typename MT::vector3_type   V;
        
        // get triangle edge vectors and plane normal
        V u = v1 - v0;
        V v = v2 - v0;
        V n = MT::cross(u, v);              // cross product
        if (n == V(0))             // triangle is degenerate
            return INFINITY;                  // do not deal with this case
        
        V dir = p1 - p0;              // ray direction vector
        V w0 = p0 - v0;
        T a = -MT::dot(n, w0);
        T b = MT::dot(n, dir);
        if (std::abs(b) < EPSILON) {     // ray is  parallel to triangle plane
            if (a == 0) {                // ray lies in triangle plane
                return 0.;
            }
            else {
                return INFINITY;              // ray disjoint from plane
            }
        }
        
        // get intersect point of ray with triangle plane
        T r = a / b;
        if (r < 0.0)                    // ray goes away from triangle
            return r;                   // => no intersect
        // for a segment, also test if (r > 1.0) => no intersect
        
        V I = p0 + r * dir;            // intersect point of ray and plane
        
        // is I inside T?
        T uu = dot(u,u);
        T uv = dot(u,v);
        T vv = dot(v,v);
        V w = I - v0;
        T wu = dot(w,u);
        T wv = dot(w,v);
        T D = uv * uv - uu * vv;
        
        // get and test parametric coords
        float s, t;
        s = (uv * wv - vv * wu) / D;
        if (s < 0.0 || s > 1.0)         // I is outside T
            return INFINITY;
        t = (uv * wu - uu * wv) / D;
        if (t < 0.0 || (s + t) > 1.0)  // I is outside T
            return INFINITY;
        
        return r;                       // I is in T
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

#endif
