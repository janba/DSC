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

namespace is_mesh {
    
    void import_tet_mesh(const std::string & filename, std::vector<real>& points, std::vector<int>&  tets, std::vector<int>& tet_labels)
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
                
                tet_labels.push_back(label);
            }
            c = '\n';
        }
        
        file.close();
    }
    
    void import_surface_mesh(const std::string& filename, std::vector<real>& points, std::vector<int>& faces)
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
                    points.push_back(x);
                    points.push_back(y);
                    points.push_back(z);
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
        }
    }
    
    void export_tet_mesh(const std::string& filename, const std::vector<vec3>& points, const std::vector<int>& tets, const std::vector<int>& tet_labels)
    {
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
    
    void export_surface_mesh(const std::string& filename, const std::vector<vec3>& points, const std::vector<int>& faces)
    {
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
