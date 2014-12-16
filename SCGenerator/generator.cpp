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

#ifdef _WIN32
const std::string file_path = "@PROJECT_SOURCE_DIR@/data/";
#else
const std::string file_path = "./data/";
#endif
const string extension = ".dsc";

void generate_from_obj(string input_file_name, string output_file_name, double padding, double avg_edge_length)
{
    vector<vec3> points;
    vector<int> tets;
    vector<int> tet_labels;
    
    std::vector<vec3> points_interface;
    std::vector<int> faces_interface;
    is_mesh::import_surface_mesh(file_path + input_file_name, points_interface, faces_interface);
    // debug
    std::cout<<"Points"<<std::endl;
    for (int i=0;i<points_interface.size();i++){
        std::cout << i<<"\t"<<points_interface[i]<<std::endl;
    }
    std::cout<<"faces_interface"<<std::endl;
    for (int i=0;i<faces_interface.size();i+=3){
        std::cout << i/3<<"\t"<<faces_interface[i]<<" "<<faces_interface[i+1]<<" "<<faces_interface[i+2]<<" "<<std::endl;
    }
    Tetralizer::tetralize(padding, avg_edge_length, points_interface, faces_interface, points, tets, tet_labels);
    
    is_mesh::export_tet_mesh(file_path + output_file_name + extension, points, tets, tet_labels);
}

void generate_from_vg(const string& input_file_name, const string& output_file_name)
{
    vector<vec3> points;
    vector<int> tets;
    vector<int> tet_labels;
    
    vec3 origin, voxel_size;
    int Ni, Nj, Nk;
    std::vector<int> voxel_labels;
    is_mesh::import_voxel_grid(file_path + input_file_name, origin, voxel_size, Ni, Nj, Nk, voxel_labels);
    
    Tetralizer::tetralize(origin, voxel_size, Ni, Nj, Nk, voxel_labels, points, tets, tet_labels);
    
    is_mesh::export_tet_mesh(file_path + output_file_name + extension, points, tets, tet_labels);
}

void generate_from_geo(const string& input_file_name, const string& output_file_name)
{
    vector<vec3> points;
    vector<int> tets;
    vector<int> tet_labels;
    
    vec3 origin;
    vec3 size;
    double discretization;
    std::vector<unsigned int> labels;
    std::vector<is_mesh::Geometry*> geometries;
    
    is_mesh::import_geometry(file_path + input_file_name, origin, size, discretization, labels, geometries);
    
    Tetralizer::tetralize(origin - vec3(3.*discretization), size + vec3(6.*discretization), discretization, labels, geometries, points, tets, tet_labels);
    
    is_mesh::export_tet_mesh(file_path + output_file_name + extension, points, tets, tet_labels);
    
}

int main(int argc, const char * argv[])
{
    if(argc > 2)
    {
        double padding = 0.5;
        double avg_edge_length = 0.5;
        if(argc > 3){
            padding = stod(argv[3]);
        }
        if(argc > 4){
            avg_edge_length = stod(argv[4]);
        }
        string output_file_name = string(argv[1]);
        string input_file_name = string(argv[2]);
        if(input_file_name.compare(input_file_name.size() - 4, 4, ".txt") == 0 || input_file_name.compare(input_file_name.size() - 3, 3, ".vg") == 0)
        {
            generate_from_vg(input_file_name, output_file_name);
            std::cout << "Generated " << output_file_name + extension << " from " << input_file_name << std::endl;
        }
        else if(input_file_name.compare(input_file_name.size() - 4, 4, ".obj") == 0)
        {
            generate_from_obj(input_file_name, output_file_name, padding, avg_edge_length);
            std::cout << "Generated " << output_file_name + extension << " from " << input_file_name << " padding "<<padding<<" avg edge length "<<avg_edge_length<<std::endl;
        }
        else if(input_file_name.compare(input_file_name.size() - 4, 4, ".geo") == 0)
        {
            generate_from_geo(input_file_name, output_file_name);
            std::cout << "Generated " << output_file_name + extension << " from " << input_file_name << std::endl;
        } else {
            std::cout << "Supported input-file formats: txt, vg, obj, geo"<<std::endl;
        }
        
    }
    if (argc <=2) {
        std::cout << R"(
Usage:
SCGenerator [output-file] [input-file] [padding=0.5] [avg_edge_length=0.5]

Supported input-file formats: txt, vg (voxel-grid), obj, geo (is-mesh)

Note that padding and avg edge length are only used for obj-files.
)";
    }
    
    return 0;
}

