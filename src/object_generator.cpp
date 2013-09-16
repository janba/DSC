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

#include "object_generator.h"

namespace DSC {
    
    void ObjectGenerator::fit_mesh_to_object(DeformableSimplicialComplex<>& dsc, const std::vector<vec3>& corners)
    {
        
    }
    
    void ObjectGenerator::label_faces(DeformableSimplicialComplex<>& dsc, const std::vector<vec3>& corners, int label)
    {
//        for (auto fi = dsc.faces_begin(); fi != dsc.faces_end(); fi++)
//        {
//            auto verts = dsc.get_pos(fi.key());
//            if(Util::is_inside(verts[0], corners) &&
//               Util::is_inside(verts[1], corners) &&
//               Util::is_inside(verts[2], corners)) // Face is inside the object.
//            {
//                dsc.update_attributes(*fi, label);
//            }
//        }
    }
    
    void ObjectGenerator::create_object(DeformableSimplicialComplex<>& dsc, const std::vector<vec3>& corners, int label)
    {
//        dsc.init_attributes(); // Maybe not necessary??
        
        fit_mesh_to_object(dsc, corners);
        label_faces(dsc, corners, label);
//        dsc.update_attributes(); // Maybe not necessary??
        
        dsc.fix_complex();
    }
    
}