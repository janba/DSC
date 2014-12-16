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
    
    static void tetrahedralize_inside(const std::vector<double>& points_interface, const std::vector<int>& faces_interface, std::vector<double>& points_inside, std::vector<int>& tets_inside);
    
    static void tetrahedralize_outside(const std::vector<double>& points_interface, const std::vector<int>&  faces_interface, std::vector<double>& points_boundary, std::vector<int>&  faces_boundary, std::vector<double>& points_outside, std::vector<int>& tets_outside, const vec3& inside_pts);
    
    static void merge_inside_outside(const std::vector<double>& points_interface, const std::vector<int>&  faces_interface, std::vector<double>& points_inside, std::vector<int>&  tets_inside, std::vector<double>& points_outside, std::vector<int>&  tets_outside, std::vector<double>& output_points, std::vector<int>&  output_tets, std::vector<int>&  output_tet_flags);
    
public:

    static void build_boundary_mesh(std::vector<double>& points_boundary, double avg_edge_length, std::vector<int>& faces_boundary, const vec3& min, const vec3& max);

    static void tetralize(double padding, double avg_edge_length, const std::vector<vec3>& points_interface, const std::vector<int>& faces_interface, std::vector<vec3>& points, std::vector<int>& tets, std::vector<int>& tet_labels);

    static void tetralize(const vec3& origin, const vec3& voxel_size, int Ni, int Nj, int Nk, const std::vector<int>& voxel_labels, std::vector<vec3>& points, std::vector<int>& tets, std::vector<int>& tet_labels);

    static void tetralize(const vec3& origin, const vec3& size, double avg_edge_length, std::vector<unsigned int>& labels, std::vector<is_mesh::Geometry*>& geometries, std::vector<vec3>& points, std::vector<int>& tets, std::vector<int>& tet_labels);
};
