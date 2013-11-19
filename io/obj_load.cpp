//
//  obj_load.c
//  Converter
//
//  Created by Asger Nyman Christiansen on 19/11/13.
//  Copyright (c) 2013 Asger Nyman Christiansen. All rights reserved.
//

#include "obj_load.h"

#include <fstream>

bool obj_load(const std::string& filename, std::vector<double>& vertices, std::vector<int>& faces)
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
//                if(ctr)
//                    faces.push_back(ctr);
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