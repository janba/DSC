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

#include "mesh_io.h"
#include "tetralizer.h"

using namespace std;

const string file_path = string("data/");
const string extension = string(".dsc");

void scale(vector<vec3>& points)
{
    real p_min = INFINITY, p_max = -INFINITY;
    for (vec3 p : points) {
        for (int i = 0; i < 3; i++) {
            p_min = Util::min(p[i], p_min);
            p_max = Util::max(p[i], p_max);
        }
    }
    std::cout << p_min << " to " << p_max << std::endl;
    
    real scale = 0.5*(p_max - p_min);
    for (vec3& p : points) {
        p = (p - vec3(p_min))/scale - vec3(1.);
    }
    
    p_min = INFINITY, p_max = -INFINITY;
    for (vec3 p : points) {
        for (int i = 0; i < 3; i++) {
            p_min = Util::min(p[i], p_min);
            p_max = Util::max(p[i], p_max);
        }
    }
    std::cout << p_min << " to " << p_max << std::endl;
}

void generate_from_obj(const string& input_file_name, const string& output_file_name)
{
    vector<vec3> points;
    vector<int> tets;
    vector<int> tet_labels;
    
    std::vector<vec3> points_interface;
    std::vector<int> faces_interface;
    is_mesh::import_surface_mesh(file_path + input_file_name, points_interface, faces_interface);
    scale(points_interface);
    
    Tetralizer::tetralize(vec3(3.), points_interface, faces_interface, points, tets, tet_labels);
    
    is_mesh::export_tet_mesh(file_path + output_file_name + extension, points, tets, tet_labels);
}

void generate_empty(const string& output_file_name)
{
    vector<vec3> points;
    vector<int>  tets;
    vector<int>  tet_labels;
    
    Tetralizer::tetralize(vec3(3.), 0.1, points, tets, tet_labels);
    is_mesh::export_tet_mesh(file_path + output_file_name + extension, points, tets, tet_labels);
}

int main(int argc, const char * argv[])
{
    if(argc > 1)
    {
        string output_file_name = string(argv[1]);
        
        if(argc == 2) {
            generate_empty(output_file_name);
        }
        else if(argc == 3)
        {
            string input_file_name = string(argv[2]);
            generate_from_obj(input_file_name, output_file_name);
        }
        
        std::cout << "Generated " << output_file_name + extension << std::endl;
    }
    
    return 0;
}

