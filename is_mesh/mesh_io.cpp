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
#include <fstream>

namespace is_mesh {
    
    void scale(std::vector<vec3>& points, real size)
    {
        vec3 p_min(INFINITY), p_max(-INFINITY);
        for (vec3 p : points) {
            for (int i = 0; i < 3; i++) {
                p_min[i] = Util::min(p[i], p_min[i]);
                p_max[i] = Util::max(p[i], p_max[i]);
            }
        }
        
        vec3 center = 0.5*(p_max + p_min);
        real scale = -INFINITY;
        for (int i = 0; i < 3; i++) {
            scale = Util::max(p_max[i] - p_min[i], scale);
        }
        
        for (vec3& p : points) {
            p = size*(p - center)/scale;
        }
    }
    
    void import_tet_mesh(const std::string & filename, std::vector<vec3>& points, std::vector<int>&  tets, std::vector<int>& tet_labels)
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
                points.push_back(vec3(x,y,z));
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
                
                tet_labels.push_back(label);
            }
            c = '\n';
        }
        file.close();
//        scale(points, 3.);
    }
    
    void import_surface_mesh(const std::string& filename, std::vector<vec3>& points, std::vector<int>& faces)
    {
        std::ifstream file(filename.data());
        
        if(file)
        {
            while(!file.eof())
            {
                std::string tok;
                file >> tok;
                if(tok == "v")
                {
                    real x,y,z;
                    file >> x;
                    file >> y;
                    file >> z;
                    points.push_back(vec3(x,y,z));
                    char line[1000];
                    file.getline(line, 998);
                }
                else if(tok == "f")
                {
                    char line[1000];
                    file.getline(line, 998);
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
                    file.getline(line, 998);
                }
            }
            file.close();
        }
//        scale(points, 2.);
    }
    
    void import_voxel_grid(const std::string& filename, vec3& origin, vec3& voxel_size, int& Ni, int& Nj, int& Nk, std::vector<int>& voxels)
    {
        std::ifstream file(filename.data());
        if(file)
        {
            while(!file.eof())
            {
                std::string tok;
                file >> tok;
                if(tok == "n")
                {
                    file >> Ni;
                    file >> Nj;
                    file >> Nk;
                }
                else if(tok == "s")
                {
                    real x, y, z;
                    file >> x;
                    file >> y;
                    file >> z;
                    voxel_size = vec3(x,y,z);
                }
                else if(tok == "o")
                {
                    real x, y, z;
                    file >> x;
                    file >> y;
                    file >> z;
                    origin = vec3(x,y,z);
                }
                else {
                    int l = atoi(tok.c_str());
                    int n;
                    file >> n;
                    for (unsigned int i = 0; i < n; i++) {
                        voxels.push_back(l);
                    }
                }
            }
        }
    }
    
    Geometry* load_geometry(std::ifstream& file)
    {
        std::string tok;
        file >> tok;
        if(tok == "cube")
        {
            real x, y, z;
            file >> x;
            file >> y;
            file >> z;
            vec3 center(x,y,z);
            file >> x;
            file >> y;
            file >> z;
            vec3 size(x,y,z);
            return new Cube(center, size);
        }
        else if(tok == "circle")
        {
            real x, y, z;
            file >> x;
            file >> y;
            file >> z;
            vec3 center(x,y,z);
            real radius;
            file >> radius;
            file >> x;
            file >> y;
            file >> z;
            vec3 normal(x,y,z);
            return new Circle(center, radius, normal);
        }
        else if(tok == "plane")
        {
            real x, y, z;
            file >> x;
            file >> y;
            file >> z;
            vec3 point(x,y,z);
            file >> x;
            file >> y;
            file >> z;
            vec3 normal(x,y,z);
            return new Plane(point, normal);
        }
        return new Geometry();
    }
    
    void import_geometry(const std::string& filename, vec3& origin, vec3& size, real& discretization, std::vector<unsigned int>& labels, std::vector<Geometry*>& geometries)
    {
        std::ifstream file(filename.data());
        if(file)
        {
            while(!file.eof())
            {
                std::string tok;
                file >> tok;
                if(tok == "d")
                {
                    file >> discretization;
                }
                else if(tok == "o")
                {
                    real x, y, z;
                    file >> x;
                    file >> y;
                    file >> z;
                    origin = vec3(x,y,z);
                }
                else if(tok == "s")
                {
                    real x, y, z;
                    file >> x;
                    file >> y;
                    file >> z;
                    size = vec3(x,y,z);
                }
                else if(tok == "c")
                {
                    real x, y, z;
                    file >> x;
                    file >> y;
                    file >> z;
                    vec3 c(x,y,z);
                    file >> x;
                    file >> y;
                    file >> z;
                    vec3 s(x,y,z);
                    geometries.push_back(new Cube(c, s));
                    file >> x;
                    labels.push_back(x);
                }
            }
        }
    }
    
    void export_tet_mesh(const std::string& filename, std::vector<vec3>& points, std::vector<int>& tets, std::vector<int>& tet_labels)
    {
//        scale(points, 3.);
        std::ofstream file(filename.data());
        
        for (auto &p : points)
        {
            file << "v " << p[0] << " " << p[1] << " " << p[2] << std::endl;
        }
        
        for (unsigned int i = 0; i < tet_labels.size(); i++)
        {
            file << "t ";
            file << tets[4*i] << " " << tets[4*i+1] << " " << tets[4*i+2] << " " << tets[4*i+3] << " ";
            file << tet_labels[i] << std::endl;
        }
    }
    
    void export_surface_mesh(const std::string& filename, std::vector<vec3>& points, std::vector<int>& faces)
    {
//        scale(points, 2.);
        std::ofstream obj_file;
        obj_file.open(filename.data());
        
        for (unsigned int i = 0; i < points.size(); ++i)
        {
            obj_file << "v "  <<   points[i][0] << " " <<   points[i][1] << " " <<   points[i][2] << std::endl;
        }
        
        for (unsigned int i = 0; i < faces.size(); ++i)
        {
            if (i%3 == 0)
            {
                obj_file << "f ";
            }
            obj_file << faces[i];
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
