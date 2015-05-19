//
// Created by Morten Nobel-Jørgensen on 10/11/14.
// Copyright (c) 2014 ___FULLUSERNAME___. All rights reserved.
//


#pragma once

#include "Vec3d.h"

namespace CGLA {
    class Ray {
        Vec3d origin = Vec3d{0};
        Vec3d direction = Vec3d{0};
        // precalculated values http://jcgt.org/published/0002/01/05/paper.pdf
        int kz;
        int kx;
        int ky;
        double Sx;
        double Sy;
        double Sz;

    public:
        Ray();
        Ray(Vec3d const &origin, Vec3d const &direction);

        Vec3d get_origin() const { return origin; }
        Vec3d get_direction() const { return direction; }

        Vec3d get_point(double f) const { return origin + direction*f; }

        double distance(Vec3d v) const ;

        // Watertight http://jcgt.org/published/0002/01/05/paper.pdf
        bool intersect_triangle(Vec3d v1, Vec3d v2, Vec3d v3, double &dist) const;

        bool is_parallel_to_triangle(Vec3d v0, Vec3d v1, Vec3d v2) const ;
    };
}


