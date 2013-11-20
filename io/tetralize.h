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

#include "tetgen.h"
#include "obj_load.h"
#include <vector>

#pragma once

namespace is_mesh {
    
    template<typename real, typename vec3>
    inline void build_boundary_mesh(std::vector<real>& points_boundary, std::vector<int>& faces_boundary, const vec3& size, int n = 3)
    {
        std::vector<std::vector<int> > face_xp_points(n+1),
        face_xm_points(n+1),
        face_yp_points(n+1),
        face_ym_points(n+1),
        face_zp_points(n+1),
        face_zm_points(n+1);
        
        for (int i = 0; i < n+1; ++i)
        {
            face_xp_points[i].resize(n+1);
            face_xm_points[i].resize(n+1);
            face_yp_points[i].resize(n+1);
            face_ym_points[i].resize(n+1);
            face_zp_points[i].resize(n+1);
            face_zm_points[i].resize(n+1);
        }
        
        real x,y,z,
        d = 4.0/static_cast<real>(n);
        int counter = 0;
        
        x = -0.5*size[0];
        for (int iy = 0; iy < n+1; ++iy)
        {
            for (int iz = 0; iz < n+1; ++iz)
            {
                y = -0.5*size[1]+iy*d;
                z = -0.5*size[2]+iz*d;
                points_boundary.push_back(x);
                points_boundary.push_back(y);
                points_boundary.push_back(z);
                face_xm_points[iy][iz] = counter;
                if (iy == 0) face_ym_points[0][iz] = counter;
                if (iy == n) face_yp_points[0][iz] = counter;
                if (iz == 0) face_zm_points[0][iy] = counter;
                if (iz == n) face_zp_points[0][iy] = counter;
                counter++;
            }
        }
        
        x = 0.5*size[0];
        for (int iy = 0; iy < n+1; ++iy)
        {
            for (int iz = 0; iz < n+1; ++iz)
            {
                y = -0.5*size[1]+iy*d;
                z = -0.5*size[2]+iz*d;
                points_boundary.push_back(x);
                points_boundary.push_back(y);
                points_boundary.push_back(z);
                face_xp_points[iy][iz] = counter;
                if (iy == 0) face_ym_points[n][iz] = counter;
                if (iy == n) face_yp_points[n][iz] = counter;
                if (iz == 0) face_zm_points[n][iy] = counter;
                if (iz == n) face_zp_points[n][iy] = counter;
                counter++;
            }
        }
        
        y = -0.5*size[1];
        for (int ix = 1; ix < n; ++ix)
        {
            for (int iz = 0; iz < n+1; ++iz)
            {
                x = -0.5*size[0]+ix*d;
                z = -0.5*size[2]+iz*d;
                points_boundary.push_back(x);
                points_boundary.push_back(y);
                points_boundary.push_back(z);
                face_ym_points[ix][iz] = counter;
                if (iz == 0) face_zm_points[ix][0] = counter;
                if (iz == n) face_zp_points[ix][0] = counter;
                counter++;
            }
        }
        
        y = 0.5*size[1];
        for (int ix = 1; ix < n; ++ix)
        {
            for (int iz = 0; iz < n+1; ++iz)
            {
                x = -0.5*size[0]+ix*d;
                z = -0.5*size[2]+iz*d;
                points_boundary.push_back(x);
                points_boundary.push_back(y);
                points_boundary.push_back(z);
                face_yp_points[ix][iz] = counter;
                if (iz == 0) face_zm_points[ix][n] = counter;
                if (iz == n) face_zp_points[ix][n] = counter;
                counter++;
            }
        }
        
        z = -0.5*size[2];
        for (int ix = 1; ix < n; ++ix)
        {
            for (int iy = 1; iy < n; ++iy)
            {
                x = -0.5*size[0]+ix*d;
                y = -0.5*size[1]+iy*d;
                points_boundary.push_back(x);
                points_boundary.push_back(y);
                points_boundary.push_back(z);
                face_zm_points[ix][iy] = counter;
                counter++;
            }
        }
        
        z = 0.5*size[2];
        for (int ix = 1; ix < n; ++ix)
        {
            for (int iy = 1; iy < n; ++iy)
            {
                x = -0.5*size[0]+ix*d;
                y = -0.5*size[1]+iy*d;
                points_boundary.push_back(x);
                points_boundary.push_back(y);
                points_boundary.push_back(z);
                face_zp_points[ix][iy] = counter;
                counter++;
            }
        }
        
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
            {
                faces_boundary.push_back(face_xm_points[i  ][j  ]);
                faces_boundary.push_back(face_xm_points[i  ][j+1]);
                faces_boundary.push_back(face_xm_points[i+1][j+1]);
                
                faces_boundary.push_back(face_xm_points[i  ][j  ]);
                faces_boundary.push_back(face_xm_points[i+1][j+1]);
                faces_boundary.push_back(face_xm_points[i+1][j  ]);
                
                faces_boundary.push_back(face_xp_points[i  ][j  ]);
                faces_boundary.push_back(face_xp_points[i  ][j+1]);
                faces_boundary.push_back(face_xp_points[i+1][j+1]);
                
                faces_boundary.push_back(face_xp_points[i  ][j  ]);
                faces_boundary.push_back(face_xp_points[i+1][j+1]);
                faces_boundary.push_back(face_xp_points[i+1][j  ]);
                
                faces_boundary.push_back(face_ym_points[i  ][j  ]);
                faces_boundary.push_back(face_ym_points[i  ][j+1]);
                faces_boundary.push_back(face_ym_points[i+1][j+1]);
                
                faces_boundary.push_back(face_ym_points[i  ][j  ]);
                faces_boundary.push_back(face_ym_points[i+1][j+1]);
                faces_boundary.push_back(face_ym_points[i+1][j  ]);
                
                faces_boundary.push_back(face_yp_points[i  ][j  ]);
                faces_boundary.push_back(face_yp_points[i  ][j+1]);
                faces_boundary.push_back(face_yp_points[i+1][j+1]);
                
                faces_boundary.push_back(face_yp_points[i  ][j  ]);
                faces_boundary.push_back(face_yp_points[i+1][j+1]);
                faces_boundary.push_back(face_yp_points[i+1][j  ]);
                
                faces_boundary.push_back(face_zm_points[i  ][j  ]);
                faces_boundary.push_back(face_zm_points[i  ][j+1]);
                faces_boundary.push_back(face_zm_points[i+1][j+1]);
                
                faces_boundary.push_back(face_zm_points[i  ][j  ]);
                faces_boundary.push_back(face_zm_points[i+1][j+1]);
                faces_boundary.push_back(face_zm_points[i+1][j  ]);
                
                faces_boundary.push_back(face_zp_points[i  ][j  ]);
                faces_boundary.push_back(face_zp_points[i  ][j+1]);
                faces_boundary.push_back(face_zp_points[i+1][j+1]);
                
                faces_boundary.push_back(face_zp_points[i  ][j  ]);
                faces_boundary.push_back(face_zp_points[i+1][j+1]);
                faces_boundary.push_back(face_zp_points[i+1][j  ]);
            }
    }
    
    template<typename real>
    inline void tetrahedralize_inside(const std::vector<real>& points_interface, const std::vector<int>& faces_interface, std::vector<real>& points_inside, std::vector<int>& tets_inside)
    {
        tetgenio in, out;
        
        in.firstnumber = 0;
        in.mesh_dim = 3;
        
        
        in.numberofpoints = (int)(points_interface.size()/3);
        //std::vector<T> points;
        //points.resize(points_interface.size() );
        //in.pointlist = &points[0];
        in.pointlist = new real[points_interface.size()];
        
        for (unsigned int i = 0; i < points_interface.size(); ++i)
        {
            in.pointlist[i] = points_interface[i];
        }
        
        in.numberoffacets = (int)(faces_interface.size()/3);
        in.facetlist = new tetgenio::facet[in.numberoffacets];
        in.facetmarkerlist = new int[in.numberoffacets];
        
        for (int i = 0; i < in.numberoffacets; ++i)
        {
            tetgenio::facet* f = &in.facetlist[i];
            f->numberofpolygons = 1;
            f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
            f->numberofholes = 0;
            f->holelist = NULL;
            tetgenio::polygon* p = &(f->polygonlist[0]);
            p->numberofvertices = 3;
            p->vertexlist = new int[p->numberofvertices];
            for (int j = 0; j < p->numberofvertices; ++j)
                p->vertexlist[j] = faces_interface[3*i+j];
        }
        
        char * tetgen_flags = "pq1.5a0.00005YY";
        //
        
        //tetgenbehavior tetbeh = tetgenbehavior();
        //tetrahedralize(&tetbeh, &in, &out);
        tetrahedralize(tetgen_flags, &in, &out);
        
        points_inside.resize(3 * out.numberofpoints);
        for (unsigned int i = 0; i < points_inside.size(); ++i)
            points_inside[i] = out.pointlist[i];
        
        tets_inside.resize(4 * out.numberoftetrahedra);
        for (unsigned int i = 0; i < tets_inside.size(); ++i)
            tets_inside[i] = out.tetrahedronlist[i];
    }
    
    template <typename real, typename vec3>
    inline void tetrahedralize_outside(const std::vector<real>& points_interface, const std::vector<int>&  faces_interface, std::vector<real>& points_boundary, std::vector<int>&  faces_boundary, std::vector<real>& points_outside, std::vector<int>& tets_outside, const vec3& inside_pts)
    {
        tetgenio in, out;
        
        in.firstnumber = 0;
        in.mesh_dim = 3;
        
        in.numberofpoints = (int)(points_interface.size()/3+points_boundary.size()/3);
        in.pointlist = new real[points_interface.size()+points_boundary.size()];
        for (unsigned int i = 0; i < points_interface.size(); ++i)
            in.pointlist[i] = points_interface[i];
        for (unsigned int i = points_interface.size(); i < points_interface.size()+points_boundary.size(); ++i)
            in.pointlist[i] = points_boundary[i-points_interface.size()];
        
        in.numberoffacets = (int)(faces_interface.size()/3+faces_boundary.size()/3);
        in.facetlist = new tetgenio::facet[in.numberoffacets];
        in.facetmarkerlist = new int[in.numberoffacets];
        
        for (int i = 0; i < in.numberoffacets; ++i)
        {
            tetgenio::facet* f = &in.facetlist[i];
            f->numberofpolygons = 1;
            f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
            f->numberofholes = 0;
            f->holelist = NULL;
            tetgenio::polygon* p = &(f->polygonlist[0]);
            p->numberofvertices = 3;
            p->vertexlist = new int[p->numberofvertices];
            for (int j = 0; j < p->numberofvertices; ++j)
            {
                if ((unsigned int)i < faces_interface.size()/3)
                    p->vertexlist[j] = faces_interface[3*i+j];
                else
                    p->vertexlist[j] = faces_boundary[3*i+j-faces_interface.size()]+points_interface.size()/3;
            }
        }
        
        in.numberofholes = 1;
        in.holelist = new real[3*in.numberofholes];
        in.holelist[0] = inside_pts[0];
        in.holelist[1] = inside_pts[1];
        in.holelist[2] = inside_pts[2];
        
        //tetgenbehavior tetbeh = tetgenbehavior();
        //tetrahedralize(&tetbeh, &in, &out);
        
        char * tetgen_flags = "pq1.8a0.005YY";
        tetrahedralize(tetgen_flags, &in, &out);
        
        points_outside.resize(3*out.numberofpoints);
        for (unsigned int i = 0; i < points_outside.size(); ++i)
            points_outside[i] = out.pointlist[i];
        
        tets_outside.resize(4*out.numberoftetrahedra);
        for (unsigned int i = 0; i < tets_outside.size(); ++i)
            tets_outside[i] = out.tetrahedronlist[i];
    }
    
    template<typename real>
    inline void merge_inside_outside(const std::vector<real>& points_interface, const std::vector<int>&  faces_interface, std::vector<real>& points_inside, std::vector<int>&  tets_inside, std::vector<real>& points_outside, std::vector<int>&  tets_outside, std::vector<real>& output_points, std::vector<int>&  output_tets, std::vector<int>&  output_tet_flags)
    {
        int no_interface_points = points_interface.size()/3;
        int no_outside_points = points_outside.size()/3;
        
        output_points.resize(points_outside.size() + points_inside.size() - points_interface.size());
        output_tets.resize(tets_inside.size() + tets_outside.size());
        output_tet_flags.resize(output_tets.size()/3);
        
        unsigned int ip, it;
        for (ip = 0; ip < points_outside.size(); ++ip)
            output_points[ip] = points_outside[ip];
        int i = 0;
        for (; ip < output_points.size(); ++ip, ++i)
            output_points[ip] = points_inside[points_interface.size() + i];
        
        for (it = 0; it < tets_outside.size(); ++it)
        {
            output_tets[it] = tets_outside[it];
            if (it%4 == 0) output_tet_flags[it/4] = 0;
        }
        for (; it < output_tets.size(); ++it)
        {
            if (tets_inside[it-tets_outside.size()] < no_interface_points)
                output_tets[it] = tets_inside[it-tets_outside.size()];
            else
                output_tets[it] = tets_inside[it-tets_outside.size()]-no_interface_points+no_outside_points;
            if (it%4 == 0) output_tet_flags[it/4] = 1;
        }
    }
    
    template <typename real>
    inline void build_tetrahedralization(const std::vector<real>& points_interface, const std::vector<int>& faces_interface, std::vector<real>& points, std::vector<int>& tets, std::vector<int>& tet_labels)
    {
        DSC::vec3 size(4.);
        DSC::vec3 inside_point(0.);
        
        std::vector<double>    points_boundary;
        std::vector<int>  faces_boundary;
        build_boundary_mesh(points_boundary, faces_boundary, size);
        
        std::vector<double> points_inside;
        std::vector<int> tets_inside;
        tetrahedralize_inside(points_interface, faces_interface, points_inside, tets_inside);
        
        std::vector<double> points_outside;
        std::vector<int> tets_outside;
        tetrahedralize_outside(points_interface, faces_interface, points_boundary, faces_boundary, points_outside, tets_outside, inside_point);
        
        merge_inside_outside(points_interface, faces_interface, points_inside, tets_inside, points_outside, tets_outside, points, tets, tet_labels);
    }
    
}
