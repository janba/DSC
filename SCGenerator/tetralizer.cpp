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

#include "tetralizer.h"
#include "tetgen.h"

static char* tetgen_flags = "pq1.5Y";

int get_index(int i, int j, int k, int Ni, int Nj, int Nk)
{
    return i + j*Ni + k*Ni*Nj;
}

void Tetralizer::tetralize_cube1(int i, int j, int k, int Ni, int Nj, int Nk, std::vector<int>& tets)
{
    // First tetrahedron:
    tets.push_back(get_index(i, j, k, Ni, Nj, Nk));
    tets.push_back(get_index(i, j, k+1, Ni, Nj, Nk));
    tets.push_back(get_index(i+1, j, k+1, Ni, Nj, Nk));
    tets.push_back(get_index(i, j+1, k+1, Ni, Nj, Nk));
    
    // Second tetrahedron:
    tets.push_back(get_index(i, j, k, Ni, Nj, Nk));
    tets.push_back(get_index(i, j+1, k+1, Ni, Nj, Nk));
    tets.push_back(get_index(i+1, j+1, k, Ni, Nj, Nk));
    tets.push_back(get_index(i, j+1, k, Ni, Nj, Nk));
    
    // Third tetrahedron:
    tets.push_back(get_index(i, j, k, Ni, Nj, Nk));
    tets.push_back(get_index(i+1, j, k, Ni, Nj, Nk));
    tets.push_back(get_index(i+1, j+1, k, Ni, Nj, Nk));
    tets.push_back(get_index(i+1, j, k+1, Ni, Nj, Nk));
    
    // Fourth tetrahedron:
    tets.push_back(get_index(i+1, j+1, k+1, Ni, Nj, Nk));
    tets.push_back(get_index(i+1, j, k+1, Ni, Nj, Nk));
    tets.push_back(get_index(i, j+1, k+1, Ni, Nj, Nk));
    tets.push_back(get_index(i+1, j+1, k, Ni, Nj, Nk));
    
    // Fifth tetrahedron:
    tets.push_back(get_index(i, j, k, Ni, Nj, Nk));
    tets.push_back(get_index(i+1, j+1, k, Ni, Nj, Nk));
    tets.push_back(get_index(i, j+1, k+1, Ni, Nj, Nk));
    tets.push_back(get_index(i+1, j, k+1, Ni, Nj, Nk));
}

void Tetralizer::tetralize_cube2(int i, int j, int k, int Ni, int Nj, int Nk, std::vector<int>& tets)
{
    // First tetrahedron:
    tets.push_back(get_index(i, j, k, Ni, Nj, Nk));
    tets.push_back(get_index(i, j, k+1, Ni, Nj, Nk));
    tets.push_back(get_index(i+1, j, k, Ni, Nj, Nk));
    tets.push_back(get_index(i, j+1, k, Ni, Nj, Nk));
    
    // Second tetrahedron:
    tets.push_back(get_index(i, j+1, k, Ni, Nj, Nk));
    tets.push_back(get_index(i+1, j+1, k, Ni, Nj, Nk));
    tets.push_back(get_index(i+1, j+1, k+1, Ni, Nj, Nk));
    tets.push_back(get_index(i+1, j, k, Ni, Nj, Nk));
    
    // Third tetrahedron:
    tets.push_back(get_index(i, j+1, k+1, Ni, Nj, Nk));
    tets.push_back(get_index(i, j+1, k, Ni, Nj, Nk));
    tets.push_back(get_index(i+1, j+1, k+1, Ni, Nj, Nk));
    tets.push_back(get_index(i, j, k+1, Ni, Nj, Nk));
    
    // Fourth tetrahedron:
    tets.push_back(get_index(i+1, j, k+1, Ni, Nj, Nk));
    tets.push_back(get_index(i+1, j, k, Ni, Nj, Nk));
    tets.push_back(get_index(i, j, k+1, Ni, Nj, Nk));
    tets.push_back(get_index(i+1, j+1, k+1, Ni, Nj, Nk));
    
    // Fifth tetrahedron:
    tets.push_back(get_index(i, j, k+1, Ni, Nj, Nk));
    tets.push_back(get_index(i+1, j, k, Ni, Nj, Nk));
    tets.push_back(get_index(i, j+1, k, Ni, Nj, Nk));
    tets.push_back(get_index(i+1, j+1, k+1, Ni, Nj, Nk));
}

void Tetralizer::create_tets(int Ni, int Nj, int Nk, std::vector<int>& tets)
{
    for (int k = 0; k < Nk; k++) {
        for (int j = 0; j < Nj; j++) {
            for (int i = 0; i < Ni; i++)
            {
                if((i + j + k)%2 == 0)
                {
                    tetralize_cube1(i, j, k, Ni+1, Nj+1, Nk+1, tets);
                }
                else {
                    tetralize_cube2(i, j, k, Ni+1, Nj+1, Nk+1, tets);
                }
            }
        }
    }
}

void Tetralizer::create_points(const vec3& origin, const vec3& voxel_size, int Ni, int Nj, int Nk, std::vector<vec3>& points)
{
    for (int k = 0; k < Nk+1; k++) {
        for (int j = 0; j < Nj+1; j++) {
            for (int i = 0; i < Ni+1; i++)
            {
                points.push_back(origin + vec3(i*voxel_size[0], j*voxel_size[1], k*voxel_size[2]));
            }
        }
    }
}

void Tetralizer::build_boundary_mesh(std::vector<real>& points_boundary, real d, std::vector<int>& faces_boundary, const vec3& size)
{
    int n = size[0]/d;
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
    
    real x,y,z;
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

void Tetralizer::tetrahedralize_inside(const std::vector<real>& points_interface, const std::vector<int>& faces_interface, std::vector<real>& points_inside, std::vector<int>& tets_inside)
{
    tetgenio in, out;
    
    in.firstnumber = 0;
    in.mesh_dim = 3;
    
    in.numberofpoints = (int)(points_interface.size()/3);
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
    
    tetrahedralize(tetgen_flags, &in, &out);
    
    points_inside.resize(3 * out.numberofpoints);
    for (unsigned int i = 0; i < points_inside.size(); ++i)
    {
        points_inside[i] = out.pointlist[i];
    }
    
    tets_inside.resize(4 * out.numberoftetrahedra);
    for (unsigned int i = 0; i < tets_inside.size(); ++i)
    {
        tets_inside[i] = out.tetrahedronlist[i];
    }
}

void Tetralizer::tetrahedralize_outside(const std::vector<real>& points_interface, const std::vector<int>&  faces_interface, std::vector<real>& points_boundary, std::vector<int>&  faces_boundary, std::vector<real>& points_outside, std::vector<int>& tets_outside, const vec3& inside_pts)
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
    
    in.numberoffacets = static_cast<int>(faces_interface.size()/3+faces_boundary.size()/3);
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
    
    tetrahedralize(tetgen_flags, &in, &out);
    
    points_outside.resize(3*out.numberofpoints);
    for (unsigned int i = 0; i < points_outside.size(); ++i)
        points_outside[i] = out.pointlist[i];
    
    tets_outside.resize(4*out.numberoftetrahedra);
    for (unsigned int i = 0; i < tets_outside.size(); ++i)
    {
        tets_outside[i] = out.tetrahedronlist[i];
    }
}

void Tetralizer::merge_inside_outside(const std::vector<real>& points_interface, const std::vector<int>&  faces_interface, std::vector<real>& points_inside, std::vector<int>&  tets_inside, std::vector<real>& points_outside, std::vector<int>&  tets_outside, std::vector<real>& output_points, std::vector<int>&  output_tets, std::vector<int>&  output_tet_flags)
{
    int no_interface_points = static_cast<int>(points_interface.size()/3);
    int no_outside_points = static_cast<int>(points_outside.size()/3);
    
    output_points.resize(points_outside.size() + points_inside.size() - points_interface.size());
    output_tets.resize(tets_inside.size() + tets_outside.size());
    output_tet_flags.resize(output_tets.size()/4);
    
    unsigned int ip, it;
    for (ip = 0; ip < points_outside.size(); ++ip)
    {
        output_points[ip] = points_outside[ip];
    }
    int i = 0;
    for (; ip < output_points.size(); ++ip, ++i)
    {
        output_points[ip] = points_inside[points_interface.size() + i];
    }
    
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
