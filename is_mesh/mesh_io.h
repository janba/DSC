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

#include "is_mesh.h"
#include <fstream>

namespace is_mesh {
    
    namespace util
    {
        struct edge_key {
            int k1, k2;
            edge_key(int i, int j) : k1(i), k2(j) {}
            bool operator<(const edge_key& k) const
            {
                return k1 < k.k1 || (k1 == k.k1 && k2 < k.k2);
            }
        };
        struct face_key
        {
            int k1, k2, k3;
            face_key(int i, int j, int k) : k1(i), k2(j), k3(k){}
            bool operator<(const face_key& k) const
            {
                //        return k1 < k.k1 || (k1 == k.k1 && k2 < k.k2 || (k2 == k.k2 && k3 < k.k3)) ;
                if (k1 < k.k1) return true;
                if (k1 == k.k1 && k2 < k.k2) return true;
                if (k1 == k.k1 && k2 == k.k2 && k3 < k.k3) return true;
                return false;
            }
        };
        struct tetrahedron_key { int k1, k2, k3, k4; };
        
        template<typename map_type, typename mesh_type>
        inline int create_edge(int i, int j, map_type& edge_map, mesh_type& mesh)
        {
            int a,b;
            if (i <= j) { a = i; b = j; }
            else { a = j; b = i; }
            edge_key key(a, b);
            typename map_type::iterator it = edge_map.find(key);
            if (it == edge_map.end())
            {
                int n = mesh.insert_edge(i, j); //non-sorted
                it = edge_map.insert(std::pair<edge_key,int>(key, n)).first;
            }
            return it->second;
        }
        
        template<typename map_type, typename mesh_type>
        inline int create_face(int i, int j, int k, map_type& face_map, mesh_type& mesh)
        {
            int a[3] = {i, j, k};
            std::sort(a, a+3);
            face_key key(a[0], a[1], a[2]);
            typename map_type::iterator it = face_map.find(key); //lookup in sorted order
            if (it == face_map.end())
            {
                int a = mesh.insert_face(i, j, k); //create in supplied order
                it = face_map.insert(std::pair<face_key,int>(key, a)).first;
            }
            return it->second;
        }
        
    }
    
    template<typename ISMesh>
    bool create(const std::vector<real>& points, const std::vector<int>& tets, ISMesh & mesh)
    {
        std::map<util::edge_key, int> edge_map;
        std::map<util::face_key, int> face_map;
        
        int cnt_nodes = 0;
        for (unsigned int i = 0; i < points.size()/3; ++i)
        {
            real x, y, z;
            x = points[3*i];
            y = points[3*i+1];
            z = points[3*i+2];
            mesh.insert_node(vec3(x,y,z));
            mesh.get(is_mesh::NodeKey(cnt_nodes)).set_destination(vec3(x,y,z));
            ++cnt_nodes;
        }
        
        for (unsigned int j = 0; j < tets.size()/4; ++j)
        {
            int idx[4];
            idx[0] = tets[4*j];
            idx[1] = tets[4*j+1];
            idx[2] = tets[4*j+2];
            idx[3] = tets[4*j+3];
            
            int edges[6];
            edges[0] = util::create_edge(idx[0], idx[1], edge_map, mesh); //edge 01
            edges[1] = util::create_edge(idx[0], idx[2], edge_map, mesh); //edge 02
            edges[2] = util::create_edge(idx[0], idx[3], edge_map, mesh); //edge 03
            edges[3] = util::create_edge(idx[1], idx[2], edge_map, mesh); //edge 12
            edges[4] = util::create_edge(idx[1], idx[3], edge_map, mesh); //edge 13
            edges[5] = util::create_edge(idx[2], idx[3], edge_map, mesh); //edge 23
            
            int faces[4];
            faces[0] = util::create_face(edges[3], edges[5], edges[4], face_map, mesh); //12-23-31
            faces[1] = util::create_face(edges[1], edges[5], edges[2], face_map, mesh); //02-23-30
            faces[2] = util::create_face(edges[0], edges[4], edges[2], face_map, mesh); //01-13-30
            faces[3] = util::create_face(edges[0], edges[3], edges[1], face_map, mesh); //01-12-20
            
            mesh.insert_tetrahedron( faces[0], faces[1], faces[2], faces[3] );
        }
        
        return true;
    }
    
    /**
     * Exports the mesh as a .dsc file.
     */
    template<typename ISMesh>
    inline void export_tet_mesh(const std::string & filename, ISMesh& mesh)
    {
        std::vector<vec3> points;
        std::vector< std::vector<int> > tets;
        mesh.extract_tet_mesh(points, tets);
        
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
    void import_tet_mesh(const std::string & filename, std::vector<real>& points, std::vector<int>&  tets, std::vector<int>& labels);
    
    inline bool import_surface_mesh(const std::string& filename, std::vector<double>& vertices, std::vector<int>& faces)
    {
        std::ifstream ifs(filename.data());
        
        if(ifs)
        {
            while(ifs.good() && !ifs.eof())
            {
                std::string tok;
                ifs >> tok;
                if(tok == "v")
                {
                    float x,y,z;
                    ifs >> x >> y >> z;
                    vertices.push_back(x);
                    vertices.push_back(y);
                    vertices.push_back(z);
                    char line[1000];
                    ifs.getline(line, 998);
                }
                else if(tok == "f")
                {
                    char line[1000];
                    ifs.getline(line, 998);
                    char* pch = strtok(line, " \t");
                    int ctr = 0;
                    while(pch != 0)
                    {
                        int v;
                        sscanf(pch, "%d", &v);
                        faces.push_back(v-1);
                        pch = strtok(0, " \t");
                        ++ctr;
                    }
                }
                else
                {
                    char line[1000];
                    ifs.getline(line, 998);
                }
            }
            
            return true;
        }
        return false;
    }
    
    template <typename ISMesh>
    inline void export_surface_mesh(const std::string& filename, ISMesh& dsc)
    {
        std::vector<vec3> verts;
        std::vector<int> indices;
        dsc.extract_surface_mesh(verts, indices);
        
        std::ofstream obj_file;
        obj_file.open(filename.data());
        
        for (unsigned int i = 0; i < verts.size(); ++i)
        {
            obj_file << "v "  <<   verts[i][0] << " " <<   verts[i][1] << " " <<   verts[i][2] << std::endl;
        }
        
        for (unsigned int i = 0; i < indices.size(); ++i)
        {
            if (i%3 == 0)
            {
                obj_file << "f ";
            }
            obj_file << indices[i];
            if (i%3 == 2)
            {
                obj_file << std::endl;
            }
            else {
                obj_file << " ";
            }
        }
        
        obj_file.close();
    }
    
}
