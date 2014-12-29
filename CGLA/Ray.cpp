//
// Created by Morten Nobel-Jørgensen on 10/11/14.
// Copyright (c) 2014 ___FULLUSERNAME___. All rights reserved.
//

#include "Ray.h"
#include <cassert>


namespace CGLA {
    Ray::Ray()
            : origin(0), direction(0)
    {
    }

    Ray::Ray(Vec3d const &origin, Vec3d const &direction)
            : origin(origin), direction(direction) {
        assert(std::abs(length(direction) - 1.0) < FLT_EPSILON*10);
    }

    bool Ray::intersect_triangle(Vec3d v0, Vec3d v1, Vec3d v2, double &dist) const {
        // based on RTR 3ed page 750
        Vec3d e1 = v1 - v0;
        Vec3d e2 = v2 - v0;
        Vec3d q = cross(direction, e2);
        double a = dot(e1, q);
        if (a > -FLT_EPSILON && a < FLT_EPSILON) return false;
        double f = 1 / a;
        Vec3d s = origin - v0;
        double u = f * dot(s, q);
        if (u < 0.0) return false;
        Vec3d r = cross(s, e1);
        double v = f * dot(direction, r);
        if (v < 0.0 || u + v > 1.0) return false;
        dist = f * dot(e2, r);
        return true;
    }
}
