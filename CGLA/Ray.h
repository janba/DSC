//
// Created by Morten Nobel-Jørgensen on 10/11/14.
// Copyright (c) 2014 ___FULLUSERNAME___. All rights reserved.
//


#pragma once

#include "Vec3d.h"

namespace CGLA {
    class Ray {
    public:
        Ray();
        Ray(Vec3d const &origin, Vec3d const &direction);

        Vec3d origin = Vec3d{0};
        Vec3d direction = Vec3d{0};

        // returns negative value if no intersection
        bool intersect_triangle(Vec3d v1, Vec3d v2, Vec3d v3, double &dist) const;
    };
}


