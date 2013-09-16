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
#include "DSC.h"

namespace DSC {
    
    class ObjectGenerator {
        
        static void fit_mesh_to_object(DeformableSimplicialComplex<>& dsc, const std::vector<vec3>& corners);
        
        static void label_faces(DeformableSimplicialComplex<>& dsc, const std::vector<vec3>& corners, int label);
        
        static void create_object(DeformableSimplicialComplex<>& dsc, const std::vector<vec3>& corners, int label);
        
        
//        static void label_tets(std::vector<int>& tet_labels)
//        {
//            for (int k = 0; k < Nk-1; k++) {
//                for (int j = 0; j < Nj-1; j++) {
//                    for (int i = 0; i < Ni-1; i++)
//                    {
//                        if(i > Ni*3./10. && i < Ni*7./10. && j > Nj*3./10. && j < Nj*7./10. && k > Nk*3./10. && k < Nk*7./10.)
//                        {
//                            for(int t = 0; t < 5; t++)
//                            {
//                                tet_labels.push_back(1);
//                            }
//                        }
//                        else {
//                            for(int t = 0; t < 5; t++)
//                            {
//                                tet_labels.push_back(0);
//                            }
//                        }
//                    }
//                }
//            }
//        }
        
        
        /**
         * Label all tetrahedra according to tet_labels and perform an initial update
         * of flags and attributes of all simplices
         */
        void label_tets(DeformableSimplicialComplex<>& dsc, const std::vector<int>& tet_labels)
        {
            // Label all tetrahedra
            for (auto tit = dsc.tetrahedra_begin(); tit != dsc.tetrahedra_end(); tit++)
            {
                dsc.set_label(tit.key(), tet_labels[tit.key()]);
            }
        }
        
    public:
        
        static void create(DeformableSimplicialComplex<>& dsc, const std::vector<int>& tet_labels)
        {
            // Label all tetrahedra
            for (auto tit = dsc.tetrahedra_begin(); tit != dsc.tetrahedra_end(); tit++)
            {
                dsc.set_label(tit.key(), tet_labels[tit.key()]);
            }
            dsc.fix_complex();
        }
        
        static void create_sphere(DeformableSimplicialComplex<>& dsc, const vec3& center, const real& radius, int label)
        {
//            std::vector<vec3> corners;
//            for (real a = 0; a < 2*M_PI; a += (1./32.)*2*M_PI)
//            {
//                corners.push_back(radius*vec3(std::cos(a), -std::sin(a)) + center);
//            }
//            create_object(dsc, corners, label);
        }
        
        
        
        static void create_box(DeformableSimplicialComplex<>& dsc, const vec3& origin, const vec3& size, int label)
        {
            std::vector<vec3> corners;
            corners.push_back(origin);
            corners.push_back(origin + vec3(size[0], 0., 0.));
            corners.push_back(origin + vec3(0., size[1], 0.));
            corners.push_back(origin + size);
            
            create_object(dsc, corners, label);
        }
    };
    
    
}