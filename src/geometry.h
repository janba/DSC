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

#include "util.h"

namespace DSC {
    
    class Geometry {
        
    public:
        Geometry()
        {
            
        }
        
        virtual bool is_inside(vec3 p)
        {
            return true;
        }
        
        virtual bool is_all_inside(std::vector<vec3> pos)
        {
            for(vec3& p : pos)
            {
                if(!is_inside(p))
                {
                    return false;
                }
            }
            return true;
        }
    };
    
    class Point : public Geometry {
    protected:
        vec3 point;
        
    public:
        Point(vec3 c) : Geometry(), point(c)
        {
            
        }
        
        virtual bool is_inside(vec3 p) override
        {
            return sqr_length(p - point) < EPSILON;
        }
    };
    
    class Cube : public Point {
    protected:
        vec3 size;
        std::vector<vec3> directions;
        
    public:
        Cube(vec3 c, vec3 s, vec3 x = vec3(1., 0., 0.), vec3 y = vec3(0., 1., 0.)) : Point(c), size(0.5*s)
        {
            directions.push_back(normalize(x));
            directions.push_back(normalize(y));
            directions.push_back(normalize(cross(directions[0], directions[1])));
        }
        
        virtual bool is_inside(vec3 p) override
        {
            for(int i = 0; i < 3; i++)
            {
                real d = dot(p - point, directions[i]);
                if(std::abs(d) > size[i])
                {
                    return false;
                }
            }
            return true;
        }
    };
    
    class Cylinder : public Point {
    protected:
        real sqr_radius, height;
        vec3 up_direction;
    public:
        Cylinder(vec3 c, real r, real h, vec3 up = vec3(0., 1., 0.)) : Point(c), sqr_radius(r*r), height(0.5*h), up_direction(normalize(up))
        {
            
        }
        
        virtual bool is_inside(vec3 p) override
        {
            real d = dot(p - point, up_direction);
            if(std::abs(d) > height)
            {
                return false;
            }
            
            vec3 p_proj = p - up_direction * d;
            return sqr_length(p_proj - point) < sqr_radius;
        }
    };
    
    class Plane : public Point {
    protected:
        vec3 normal;
        
    public:
        Plane(vec3 p, vec3 n) : Point(p), normal(normalize(n))
        {
            
        }
        
        virtual bool is_inside(vec3 p) override
        {
            return std::abs(dot(p - point, normal)) < EPSILON;
        }
    };
    
    class Circle : public Cylinder {
        
    public:
        Circle(vec3 center, real radius, vec3 normal) : Cylinder(center, radius, 2.*EPSILON, normal)
        {
            
        }
    };
    
    class Square : public Cube {
        
    public:
        Square(vec3 center, real width, real height, vec3 width_dir, vec3 height_dir) : Cube(center, vec3(width, height, 2.*EPSILON), width_dir, height_dir)
        {
            
        }
    };
    
}