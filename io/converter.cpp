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


#include "tetralize.h"
#include "src/mesh_io.h"
#include "is_mesh/is_mesh.h"
#include "src/attributes.h"
#include "src/object_generator.h"

using namespace std;

int main(int argc, const char * argv[])
{
//    string input_file_name = string("data/") + string(argv[0]);
//    string output_file_name = string("data/") + string(argv[1]);
    
    string input_file_name = string("data/eight.obj");
    string output_file_name = string("data/test.dsc");
    DSC::vec3 size(4.);
    
    vector<double> points;
    vector<int>  tets;
    vector<int>  tet_labels;
    vector<DSC::vec3> pts_inside = {DSC::vec3(0.)};
    
    is_mesh::build_tetrahedralization(input_file_name, size, points, tets, tet_labels, pts_inside);

    is_mesh::ISMesh<DSC::NodeAttributes, DSC::EdgeAttributes, DSC::FaceAttributes, DSC::TetAttributes> mesh(points, tets);
    
    DSC::ObjectGenerator::create(mesh, tet_labels);
    
    DSC::export_tet_mesh(mesh, output_file_name);
    
    return 0;
}

