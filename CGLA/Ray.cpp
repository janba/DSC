//
// Created by Morten Nobel-Jørgensen on 10/11/14.
// Copyright (c) 2014 ___FULLUSERNAME___. All rights reserved.
//

#include "Ray.h"
#include "Vec4f.h"
#include "Vec4d.h"
#include <cassert>


namespace CGLA {
    Ray::Ray()
            : origin(0), direction(0)
    {
    }

    int max_dim(Vec3d v){
        int maxDim = 0;
        for (int i=1;i<3;i++){
            if (std::fabs(v[i]) > std::fabs(v[maxDim])){
                maxDim = i;
            }
        }
        return maxDim;
    }

    Ray::Ray(Vec3d const &origin, Vec3d const &direction)
            : origin(origin), direction(direction) {
        assert(std::abs(length(direction) - 1.0) < FLT_EPSILON*10);
        // calculate dimension where the ray direction is maximal ∗/
        kz = max_dim(direction);
        kx = kz+1; if (kx == 3) kx = 0;
        ky = kx+1; if (ky == 3) ky = 0;
        // swap kx and ky dimension to preserve winding direction of triangles ∗/
        if (direction[kz] < 0.0f) std::swap(kx,ky);
        // calculate shear constants ∗/
        Sx = direction[kx]/direction[kz];
        Sy = direction[ky]/direction[kz];
        Sz = 1.0f/direction[kz];
    }

    double Ray::distance(Vec3d v) const {
        return length(origin + dot((v-origin), direction) / dot(direction, direction) * direction);
    }

    bool Ray::is_parallel_to_triangle(Vec3d v0, Vec3d v1, Vec3d v2) const {
        Vec3d faceNormal = cross(normalize(v0-v1), normalize(v2-v1));
        double val = fabs(dot(faceNormal, direction));
        return val < FLT_EPSILON;
    }

    bool Ray::intersect_triangle(Vec3d v0, Vec3d v1, Vec3d v2, double &dist) const {
        // calculate vertices relative to ray origin
        const Vec3d A = (v0-origin);
        const Vec3d B = (v1-origin);
        const Vec3d C = (v2-origin);
        // perform shear and scale of vertices ∗/
        const double Ax = A[kx] - Sx*A[kz];
        const double Ay = A[ky] - Sy*A[kz];
        const double Bx = B[kx] - Sx*B[kz];
        const double By = B[ky] - Sy*B[kz];
        const double Cx = C[kx] - Sx*C[kz];
        const double Cy = C[ky] - Sy*C[kz];
        // calculate scaled barycentric coordinates ∗/
        double U = Cx*By - Cy*Bx;
        double V = Ax*Cy - Ay*Cx;
        double W = Bx*Ay - By*Ax;
        // fallback to test against edges using double precision ∗/
        if (U == 0.0f || V == 0.0f || W == 0.0f) {
            double CxBy = Cx * By;
            double CyBx = Cy * Bx;
            U = (float)(CxBy - CyBx);
            double AxCy = Ax * Cy;
            double AyCx = Ay * Cx;
            V = (float)(AxCy - AyCx);
            double BxAy = Bx * Ay;
            double ByAx = By * Ax;
            W = (float)(BxAy - ByAx);
        }
        // Perform edge tests. Moving this test before and at the end of the previous conditional gives higher performance . ∗/
#ifdef BACKFACE_CULLING
  if (U<0.0f || V<0.0f || W<0.0f) return;
#else
        if ((U<0.0f || V<0.0f || W<0.0f) &&
                (U>0.0f || V>0.0f || W>0.0f)) return false;
#endif
        // calculate determinant ∗/
        double det = U+V+W;
        if (det == 0.0f) return false;
        // Calculate scaled z−coordinates of vertices
        //        and use them to calculate the hit distance. ∗/
        const double Az = Sz*A[kz];
        const double Bz = Sz*B[kz];
        const double Cz = Sz*C[kz];
        const double T = U*Az + V*Bz + W*Cz;
#ifdef BACKFACE_CULLING
  if (T < 0.0f || T > hit.t * det)
    return;
#else
        /*static const u_int32_t SIGNMASK = 0x80000000;
        u_int32_t det_sign = ((u_int32_t)det) & SIGNMASK;

        #define ftoi(X) ((u_int32_t)X)
        #define xorf(X,Y) (float)(ftoi(X)^ftoi(Y))

        if ((xorf(T,det_sign) < 0.0f) || (xorf(T,det_sign) > dist * xorf(det, det_sign)))
            return false;   */
#endif
        // normalize U, V, W, and T ∗/
        const double rcpDet = 1.0f/det;
        Vec4d hit;
        hit[0] = U*rcpDet;
        hit[1] = V*rcpDet;
        hit[2] = W*rcpDet;
        hit[3] = T*rcpDet;
        dist = hit[3];
        return true;
    }
}
