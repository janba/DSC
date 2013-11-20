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
#include "attributes.h"
#include "object_generator.h"
#include "tetralizer.h"
#include "is_mesh.h"

using namespace std;

const string file_path = string("data/");
const string extension = string(".dsc");

void generate_from_obj(const string& input_file_name, const string& output_file_name)
{
    vector<double> points;
    vector<int> tets;
    vector<int> tet_labels;
    
    std::vector<double> points_interface;
    std::vector<int> faces_interface;
    is_mesh::import_surface_mesh(file_path + input_file_name, points_interface, faces_interface);
    
    Tetralizer::tetralize(points_interface, faces_interface, points, tets, tet_labels);
    
    is_mesh::ISMesh<is_mesh::NodeAttributes, is_mesh::EdgeAttributes, is_mesh::FaceAttributes, is_mesh::TetAttributes> mesh(points, tets, tet_labels);
    
    {
        std::vector<vec3> points;
        std::vector< std::vector<int>> tets;
        mesh.extract_tet_mesh(points, tets);
        is_mesh::export_tet_mesh(file_path + output_file_name + extension, points, tets);
    }
}

void generate_cube()
{
    vector<double> points;
    vector<int>  tets;
    vector<int>  tet_labels;
    
    Tetralizer::tetralize(vec3(1.), 0.05, points, tets);
    
    is_mesh::ISMesh<is_mesh::NodeAttributes, is_mesh::EdgeAttributes, is_mesh::FaceAttributes, is_mesh::TetAttributes> mesh(points, tets, tet_labels);
    
    double size = 0.75;
    ObjectGenerator::create_cube(mesh, vec3(-size/2.), vec3(size), 1);
    
    {
        std::vector<vec3> points;
        std::vector< std::vector<int>> tets;
        mesh.extract_tet_mesh(points, tets);
        is_mesh::export_tet_mesh(file_path + string("cube") + extension, points, tets);
    }
}

void generate_one_cell()
{
    vector<double> points;
    vector<int>  tets;
    vector<int>  tet_labels;
    
    Tetralizer::tetralize(vec3(1.), 1./3., points, tets);
    
    is_mesh::ISMesh<is_mesh::NodeAttributes, is_mesh::EdgeAttributes, is_mesh::FaceAttributes, is_mesh::TetAttributes> mesh(points, tets, tet_labels);
    
    ObjectGenerator::create_cube(mesh, vec3(-1./6.), vec3(1./3.), 1);
    
    {
        std::vector<vec3> points;
        std::vector< std::vector<int>> tets;
        mesh.extract_tet_mesh(points, tets);
        is_mesh::export_tet_mesh(file_path + string("one_cell") + extension, points, tets);
    }
}

int main(int argc, const char * argv[])
{
//    string input_file_name = string("data/") + string(argv[0]);
//    string output_file_name = string("data/") + string(argv[1]);
    
    string input_file_name = string("eight.obj");
    string output_file_name = string("test");
    
    generate_from_obj(input_file_name, output_file_name);
    
    return 0;
}

