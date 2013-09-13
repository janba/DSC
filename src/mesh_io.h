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

#include "DSC.h"

#include <fstream>
#include <string>
#include <vector>


/**
 * Exports the mesh as a .dsc file.
 */
inline void export_tet_mesh(DeformableSimplicialComplex<>& dsc, const std::string & filename)
{
	std::vector<vec3> points;
    std::vector< std::vector<int> > tets;
	dsc.extract_tet_mesh(points, tets);
    
    std::ofstream file(filename.data());
    
    for (auto &p : points)
    {
		file << "v " << p[0] << " " << p[1] << " " << p[2] << std::endl;
    }
    
    for (std::vector<int> tet : tets)
    {
		file << "t ";
        for (int i = 0; i < tet.size() - 1; i++)
        {
            file << tet[i];
            file << " ";
        }
        file << tet[tet.size()-1] << std::endl;
    }
    
}

/**
 * Imports a mesh from a .dsc file.
 */
template <typename real>
inline void import_tet_mesh(const std::string & filename, std::vector<real>& points, std::vector<int>&  tets, std::vector<int>&  labels)
{
    
    std::ifstream file(filename.data());
    
    while (!file.eof())
    {
		char c;
		file >> c;
		if (c == 'v')
		{
            real x,y,z; // The (x,y,z) coordinates of a vertex.
            file >> x;
            file >> y;
            file >> z;
            points.push_back(x);
            points.push_back(y);
            points.push_back(z);
		}
		else if (c == 't')
		{
            int v1, v2, v3, v4; // The indeces of the four vertices of a tetrahedron.
            int label; // The label of a tetrahedron.
            file >> v1;
            file >> v2;
            file >> v3;
            file >> v4;
            file >> label;
            
            tets.push_back(v1);
            tets.push_back(v2);
            tets.push_back(v3);
            tets.push_back(v4);
            
            labels.push_back(label);
		}
        c = '\n';
    }
    
    file.close();
}

inline void save_interface(DeformableSimplicialComplex<>& dsc, std::string & filename)
{
	std::vector<vec3> verts;
	std::vector<int> indices;
	dsc.extract_interface(verts, indices);
    
	std::ofstream obj_file;
	obj_file.open(filename.data());
    
	for (unsigned int i = 0; i < verts.size(); ++i)
	{
		obj_file << "v "  <<   verts[i][0] << " " <<   verts[i][1] << " " <<   verts[i][2] << std::endl;
	}
    
	for (unsigned int i = 0; i < indices.size(); ++i)
	{
		if (i%3 == 0) obj_file << "f ";
		obj_file << indices[i];
		if (i%3 == 2) obj_file << std::endl;
		else obj_file << " ";
	}
    
	obj_file.close();
}
