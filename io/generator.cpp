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


#include "is_mesh/is_mesh.h"
#include "mesh_io.h"
#include "attributes.h"
#include "object_generator.h"
#include "tetralize.h"
#include "tetralizer.h"

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
    obj_load(file_path + input_file_name, points_interface, faces_interface);
    
    is_mesh::build_tetrahedralization(points_interface, faces_interface, points, tets, tet_labels);
    
    is_mesh::ISMesh<DSC::NodeAttributes, DSC::EdgeAttributes, DSC::FaceAttributes, DSC::TetAttributes> mesh(points, tets, tet_labels);
    
    DSC::export_tet_mesh(mesh, file_path + output_file_name + extension);
}

void generate_cube()
{
    vector<double> points;
    vector<int>  tets;
    vector<int>  tet_labels;
    
    DSC::Tetralizer::tetralize(DSC::vec3(1.), 0.05, points, tets);
    
    is_mesh::ISMesh<DSC::NodeAttributes, DSC::EdgeAttributes, DSC::FaceAttributes, DSC::TetAttributes> mesh(points, tets);
    
    double size = 0.75;
    DSC::ObjectGenerator::create_cube(mesh, DSC::vec3(-size/2.), DSC::vec3(size), 1);
    
    DSC::export_tet_mesh(mesh, file_path + string("cube") + extension);
}

void generate_one_cell()
{
    vector<double> points;
    vector<int>  tets;
    vector<int>  tet_labels;
    
    DSC::Tetralizer::tetralize(DSC::vec3(1.), 1./3., points, tets);
    
    is_mesh::ISMesh<DSC::NodeAttributes, DSC::EdgeAttributes, DSC::FaceAttributes, DSC::TetAttributes> mesh(points, tets);
    
    DSC::ObjectGenerator::create_cube(mesh, DSC::vec3(-1./6.), DSC::vec3(1./3.), 1);
    
    DSC::export_tet_mesh(mesh, file_path + string("one_cell") + extension);
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

