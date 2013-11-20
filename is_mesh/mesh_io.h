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
#include <fstream>

namespace is_mesh {
    
    /**
     * Imports a mesh from a .dsc file.
     */
    void import_tet_mesh(const std::string & filename, std::vector<vec3>& points, std::vector<int>&  tets, std::vector<int>& tet_labels);
    
    /**
     * Imports a surface mesh from an .obj file.
     */
    void import_surface_mesh(const std::string& filename, std::vector<vec3>& points, std::vector<int>& faces);
    
    /**
     * Exports the mesh as a .dsc file.
     */
    void export_tet_mesh(const std::string& filename, const std::vector<vec3>& points, const std::vector<int>& tets, const std::vector<int>& tet_labels);
    
    /**
     * Exports the surface mesh to an .obj file.
     */
    void export_surface_mesh(const std::string& filename, const std::vector<vec3>& points, const std::vector<int>& faces);
    
}
