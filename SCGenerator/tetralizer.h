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
#include "geometry.h"

class Tetralizer
{
    static void tetralize_cube1(int i, int j, int k, int Ni, int Nj, int Nk, std::vector<int>& tets);
    
    static void tetralize_cube2(int i, int j, int k, int Ni, int Nj, int Nk, std::vector<int>& tets);
    
    static void create_tets(int Ni, int Nj, int Nk, std::vector<int>& tets);
    
    static void create_points(const vec3& origin, const vec3& voxel_size, int Ni, int Nj, int Nk, std::vector<vec3>& points);
    
    static void build_boundary_mesh(std::vector<real>& points_boundary, real avg_edge_length, std::vector<int>& faces_boundary, const vec3& size);
    
    static void tetrahedralize_inside(const std::vector<real>& points_interface, const std::vector<int>& faces_interface, std::vector<real>& points_inside, std::vector<int>& tets_inside);
    
    static void tetrahedralize_outside(const std::vector<real>& points_interface, const std::vector<int>&  faces_interface, std::vector<real>& points_boundary, std::vector<int>&  faces_boundary, std::vector<real>& points_outside, std::vector<int>& tets_outside, const vec3& inside_pts);
    
    static void merge_inside_outside(const std::vector<real>& points_interface, const std::vector<int>&  faces_interface, std::vector<real>& points_inside, std::vector<int>&  tets_inside, std::vector<real>& points_outside, std::vector<int>&  tets_outside, std::vector<real>& output_points, std::vector<int>&  output_tets, std::vector<int>&  output_tet_flags);
    
public:
    
    static void tetralize(const vec3& size, real avg_edge_length, const std::vector<vec3>& points_interface, const std::vector<int>& faces_interface, std::vector<vec3>& points, std::vector<int>& tets, std::vector<int>& tet_labels)
    {
        std::vector<real> points_interface_real;
        for (vec3 p : points_interface) {
            points_interface_real.push_back(p[0]);
            points_interface_real.push_back(p[1]);
            points_interface_real.push_back(p[2]);
        }
        
        std::vector<real>    points_boundary;
        std::vector<int>  faces_boundary;
        build_boundary_mesh(points_boundary, avg_edge_length, faces_boundary, size);
        
        std::vector<real> points_inside;
        std::vector<int> tets_inside;
        tetrahedralize_inside(points_interface_real, faces_interface, points_inside, tets_inside);
        
        std::vector<real> points_outside;
        std::vector<int> tets_outside;
        tetrahedralize_outside(points_interface_real, faces_interface, points_boundary, faces_boundary, points_outside, tets_outside, vec3(points_inside[0], points_inside[1], points_inside[2]));

        std::vector<real> points_real;
        merge_inside_outside(points_interface_real, faces_interface, points_inside, tets_inside, points_outside, tets_outside, points_real, tets, tet_labels);
        
        points.resize(points_real.size()/3);
        for (unsigned int i = 0; i < points_real.size()/3; i++) {
            points[i] = vec3(points_real[i*3], points_real[i*3+1], points_real[i*3+2]);
        }
    }
    
    static void tetralize(const vec3& origin, const vec3& voxel_size, int Ni, int Nj, int Nk, const std::vector<int>& voxel_labels, std::vector<vec3>& points, std::vector<int>& tets, std::vector<int>& tet_labels)
    {
        create_points(origin, voxel_size, Ni, Nj, Nk, points);
        create_tets(Ni, Nj, Nk, tets);
        
        for (int l : voxel_labels) {
            for (int i = 0; i < 5; i++) {
                tet_labels.push_back(l);
            }
        }
    }
    
    static void tetralize(const vec3& origin, const vec3& size, real avg_edge_length, std::vector<unsigned int>& labels, std::vector<is_mesh::Geometry*>& geometries, std::vector<vec3>& points, std::vector<int>& tets, std::vector<int>& tet_labels)
    {
        int Ni = std::ceil(size[0]/avg_edge_length);
        int Nj = std::ceil(size[1]/avg_edge_length);
        int Nk = std::ceil(size[2]/avg_edge_length);
        
        create_points(origin, vec3(avg_edge_length), Ni, Nj, Nk, points);
        create_tets(Ni, Nj, Nk, tets);
        
        for (unsigned int i = 0; i < tets.size(); i += 4)
        {
            int label = 0;
            for (auto g = 0; g < labels.size(); g++) {
                vec3 bc = Util::barycenter(points[tets[i]], points[tets[i+1]], points[tets[i+2]], points[tets[i+3]]);
                if(geometries[g]->is_inside(bc))
                {
                    label = labels[g];
                }
            }
            tet_labels.push_back(label);
        }
    }
};
