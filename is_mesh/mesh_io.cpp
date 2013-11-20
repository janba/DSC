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
    
    void import_tet_mesh(const std::string & filename, std::vector<real>& points, std::vector<int>&  tets, std::vector<int>& labels)
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
    
}
