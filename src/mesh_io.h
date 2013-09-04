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

#ifndef MESH_IO_H
#define MESH_IO_H

#include "DSC.h"

#include <fstream>
#include <string>
#include <vector>


/**
 * Exports the mesh as a .dsc file.
 */
template <typename MT>
inline void export_tet_mesh(DeformableSimplicialComplex<MT> & dsc, const std::string & filename)
{
	std::vector<typename MT::vector3_type> points;
	vector< vector<int> > tets;
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
template <typename MT>
inline void import_tet_mesh(const std::string & filename, vector<typename MT::real_type>& points, vector<int>&  tets, vector<int>&  labels)
{
    
    std::ifstream file(filename.data());
    
    while (!file.eof())
    {
		char c;
		file >> c;
		if (c == 'v')
		{
            double x,y,z; // The (x,y,z) coordinates of a vertex.
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

template <typename MT>
inline void save_interface(DeformableSimplicialComplex<MT> & dsc, string & filename)
{
	vector<typename MT::vector3_type> verts;
	vector<int> indices;
	dsc.extract_interface(verts, indices);
    
	ofstream obj_file;
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

/**
 * Document the mesh by screenshot and/or by saving the mesh.
 */
/*template <typename MT>
 inline void document(SimplicialComplex<MT> & dsc, int iteration_nr)
 {
 if (iteration_nr%DOCUMENTATION_GAP == 0)
 {
 if(TAKE_SCREENSHOT)
 {
 cout << "Taking screenshot" << endl;
 string str(details::concat4digits("SS/scr_", iteration_nr/DOCUMENTATION_GAP) + ".bmp");
 SOIL_save_screenshot(str.c_str(), SOIL_SAVE_TYPE_BMP, 0, 0, 640, 640 );
 string str1(details::concat4digits("OBJ/mesh_", iteration_nr/DOCUMENTATION_GAP) + ".obj");
 //save_obj(dsc, str1);
 
 }
 if(SAVE_OBJECT)
 {
 cout << "Saving object" << endl;
 string str(details::concat4digits("OBJ/obj_", iteration_nr/DOCUMENTATION_GAP) + ".obj");
 save_interface(dsc, str);
 }
 }
 }*/

template <typename MT>
inline void save_obj(DeformableSimplicialComplex<MT> & dsc, std::vector<std::string> const & materials, std::string const & path)
{
    typedef typename MT::real_type      T;
    typedef typename MT::vector3_type   V;
    
    std::vector<V> verts;
    std::vector<V> normals;
    std::vector<int> labels;
    std::vector<int> indices;
    std::map<typename DeformableSimplicialComplex<MT>::node_key_type, int> vert_index;
    std::map<typename DeformableSimplicialComplex<MT>::edge_key_type, int> edge_index;
    
    typename DeformableSimplicialComplex<MT>::node_iterator nit = dsc.nodes_begin();
    typename DeformableSimplicialComplex<MT>::node_iterator nn_it = dsc.nodes_end();
    
    while (nit != nn_it)
    {
        if (nit->is_interface())
        {
            verts.push_back(nit->v);
            vert_index[nit.key()] = verts.size();
        }
        ++nit;
    }
    
    typename DeformableSimplicialComplex<MT>::face_iterator fit = dsc.faces_begin();
    typename DeformableSimplicialComplex<MT>::face_iterator ff_it = dsc.faces_end();
    
    while (fit != ff_it)
    {
        if (fit->is_interface())
        {
            int l = enforce_orientation(dsc, fit.key());
            
            std::vector<typename DeformableSimplicialComplex<MT>::node_key_type> nodes(3);
            dsc.vertices(fit.key(), nodes);
            
            indices.push_back(vert_index[nodes[0]]);
            indices.push_back(vert_index[nodes[1]]);
            indices.push_back(vert_index[nodes[2]]);
            
            labels.push_back(l);
        }
        ++fit;
    }
    
    std::ofstream obj_file;
    obj_file.open(path.data());
    
    obj_file << "mtllib dsc_test.mtl" << std::endl;
    
    for (int i = 0; i < verts.size(); ++i)
        obj_file << "v "  <<   verts[i][0] << " " <<   verts[i][1] << " " <<   verts[i][2] << std::endl;
    
    int last_label = -1;
    for (int i = 0; i < indices.size() / 3; ++i)
    {
        int label = labels[i];
        if (label != last_label)
        {
            last_label = label;
            obj_file << "usemtl " << materials[label-1] << std::endl;
        }
        
        obj_file << "f " << indices[3*i  ] << " " <<
        indices[3*i+1] << " " <<
        indices[3*i+2] << std::endl;
    }
    
    obj_file.close();
    
}

// MESH_IO_H
#endif
