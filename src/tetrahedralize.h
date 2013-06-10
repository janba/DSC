#ifndef TETRAHEDRALIZE_H
#define TETRAHEDRALIZE_H

/************************************************************************************************************************
 Marek: I believe we should kill this. We don't want Tetgen to be a part of the DSC release, the input for DSC should be
 a file with a tetrahedral mesh, possibly generate with Brian's scene composer. Same goes for all the TinyXML stuff.
 *************************************************************************************************************************/

#include <vector>

#include "tetgen.h"

#include <CGLA/Vec3d.h>   // only to make HMesh stop complaining?
#include <CGLA/Vec4d.h>
#include <CGLA/Mat4x4d.h>
#include <CGLA/Mat3x3d.h>
//#include <HMesh/FaceCirculator.h>
#include <HMesh/Manifold.h>
#include <HMesh/obj_load.h>
#include <HMesh/x3d_load.h>

#include <cassert>

using namespace std;
using namespace HMesh;

template <typename MT>
inline void init_mesh(Manifold& m,
                      vector<typename MT::real_type>& points_interface,
                      vector<int>&  faces_interface,
                      vector<typename MT::vector3_type>& inside_pts,
                      typename MT::real_type scaling_factor = 1.25)
{
    typedef typename MT::vector3_type V;
    
	V c(0.0f,0.0f,0.0f);
    V c_new(0.0f);
    
	float r = 0.0f;
    
    CGLA::Vec3d cf( c[0], c[1], c[2] );
	bsphere(m, cf, r);    // Kenny: Argh! bsphere is hardwired to Vec3f we must convert types!
	
    std::cout << cf << " " << r << std::endl;
	
    r *= scaling_factor;
    
    VertexAttributeVector<int> touched(0);
	unsigned int i = 0;
	for (VertexIDIterator vi = m.vertices_begin(); vi != m.vertices_end(); ++vi)
	{
        m.pos(*vi) -= cf;   // Kenny: Argh! Manifold is hardwired to Vec3f we must convert types!
        m.pos(*vi) /= r;
        assert(!MT::is_nan(m.pos(*vi)[0]) && !MT::is_nan(m.pos(*vi)[1]) && !MT::is_nan(m.pos(*vi)[2]));
		{
			c_new += V( m.pos(*vi) );
            touched[*vi] = i;
			++i;
			for (int j = 0; j < 3; ++j)
				points_interface.push_back(m.pos(*vi)[j]);
		}
	}
    
	c_new /= m.no_vertices()/2.0;
    assert(!MT::is_nan(c_new[0]) && !MT::is_nan(c_new[1]) && !MT::is_nan(c_new[2]));
	i = 0;
	for (; i < points_interface.size()/3; ++i)
	{
		typename MT::vector3_type p(points_interface[3*i  ],
                                    points_interface[3*i+1],
                                    points_interface[3*i+2]);
		if (points_interface[3*i+1] < 0)
		{
			points_interface[3*i  ] = p[0];
			points_interface[3*i+1] = p[1];
			points_interface[3*i+2] = p[2];
		}
	}
    
	i = 0;
	for (FaceIDIterator fi = m.faces_begin(); fi != m.faces_end(); ++fi, ++i)
	{
		int j = 0;
        Walker w = m.walker(*fi);
		for(; !w.full_circle(); w = w.next())
		{
			if (touched[w.vertex()] != -1)
				faces_interface.push_back(touched[w.vertex()]);
			++j;
		}
	}
}

template<typename T>
inline void build_boundary_mesh(vector<T>& points_boundary,
                                vector<int>&  faces_boundary,
                                int n = 3)
{
	vector<vector<int> > face_xp_points(n+1),
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
    
	T x,y,z,
    d = 4.0/(double)n;
	int counter = 0;
    
	x = -2.0;
	for (int iy = 0; iy < n+1; ++iy)
	{
		for (int iz = 0; iz < n+1; ++iz)
		{
			y = -2.0+iy*d;
			z = -2.0+iz*d;
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
    
	x = 2.0;
	for (int iy = 0; iy < n+1; ++iy)
	{
		for (int iz = 0; iz < n+1; ++iz)
		{
			y = -2.0+iy*d;
			z = -2.0+iz*d;
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
    
	y = -2.0;
	for (int ix = 1; ix < n; ++ix)
	{
		for (int iz = 0; iz < n+1; ++iz)
		{
			x = -2.0+ix*d;
			z = -2.0+iz*d;
			points_boundary.push_back(x);
			points_boundary.push_back(y);
			points_boundary.push_back(z);
			face_ym_points[ix][iz] = counter;
			if (iz == 0) face_zm_points[ix][0] = counter;
			if (iz == n) face_zp_points[ix][0] = counter;
			counter++;
		}
	}
    
	y = 2.0;
	for (int ix = 1; ix < n; ++ix)
	{
		for (int iz = 0; iz < n+1; ++iz)
		{
			x = -2.0+ix*d;
			z = -2.0+iz*d;
			points_boundary.push_back(x);
			points_boundary.push_back(y);
			points_boundary.push_back(z);
			face_yp_points[ix][iz] = counter;
			if (iz == 0) face_zm_points[ix][n] = counter;
			if (iz == n) face_zp_points[ix][n] = counter;
			counter++;
		}
	}
    
	z = -2.0;
	for (int ix = 1; ix < n; ++ix)
	{
		for (int iy = 1; iy < n; ++iy)
		{
			x = -2.0+ix*d;
			y = -2.0+iy*d;
			points_boundary.push_back(x);
			points_boundary.push_back(y);
			points_boundary.push_back(z);
			face_zm_points[ix][iy] = counter;
			counter++;
		}
	}
    
	z = 2.0;
	for (int ix = 1; ix < n; ++ix)
	{
		for (int iy = 1; iy < n; ++iy)
		{
			x = -2.0+ix*d;
			y = -2.0+iy*d;
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


template<typename T>
inline void tetrahedralize_inside(vector<T>& points_interface,
                                  vector<int>&  faces_interface,
                                  vector<T>& points_inside,
                                  vector<int>&  tets_inside)
{
	tetgenio in, out;
    
	in.firstnumber = 0;
	in.mesh_dim = 3;
    
    
	in.numberofpoints = (int)(points_interface.size()/3);
    //std::vector<T> points;
    //points.resize(points_interface.size() );
	//in.pointlist = &points[0];
	in.pointlist = new T[points_interface.size()];
    
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

template <typename MT>
inline void tetrahedralize_outside(vector<typename MT::real_type>& points_interface,
                                   vector<int>&  faces_interface,
                                   vector<typename MT::real_type>& points_boundary,
                                   vector<int>&  faces_boundary,
                                   vector<typename MT::real_type>& points_outside,
                                   vector<int>&  tets_outside,
                                   vector<typename MT::vector3_type>& inside_pts)
{
    typedef typename MT::real_type T;
    
	tetgenio in, out;
    
	in.firstnumber = 0;
	in.mesh_dim = 3;
    
	in.numberofpoints = (int)(points_interface.size()/3+points_boundary.size()/3);
	in.pointlist = new T[points_interface.size()+points_boundary.size()];
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
    
	in.numberofholes = inside_pts.size();
	in.holelist = new T[3*in.numberofholes];
	for (int i = 0; i < in.numberofholes; ++i)
	{
		in.holelist[3*i  ] = inside_pts[i][0];
		in.holelist[3*i+1] = inside_pts[i][1];
		in.holelist[3*i+2] = inside_pts[i][2];
	}
    
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

template<typename T>
inline void merge_inside_outside(vector<T>& points_interface,
                                 vector<int>&  faces_interface,
                                 vector<T>& points_inside,
                                 vector<int>&  tets_inside,
                                 vector<T>& points_outside,
                                 vector<int>&  tets_outside,
                                 vector<T>& output_points,
                                 vector<int>&  output_tets,
                                 vector<int>&  output_tet_flags)
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



template <typename MT>
inline void build_tetrahedralization(string const & filename,
                                     vector<typename MT::real_type> & points,
                                     vector<int>&  tets,
                                     vector<int>&  tet_flags,
                                     vector<typename MT::vector3_type>& inside_pts)
{
    typedef typename MT::real_type T;
    
	Manifold m;
	m.clear();
    
	obj_load(filename, m);
	//x3d_load(filename, m);
    
	vector<T>    points_interface;
	vector<int>  faces_interface;
	vector<T>    points_boundary;
	vector<int>  faces_boundary;
    
	vector<T>     points_inside, points_outside;
	vector<int>  tets_inside, tets_outside;
    
	init_mesh<MT>(m, points_interface, faces_interface, inside_pts);
	build_boundary_mesh<T>(points_boundary, faces_boundary);
    
	tetrahedralize_inside(points_interface, faces_interface, points_inside, tets_inside);
	tetrahedralize_outside<MT>(points_interface, faces_interface,
                               points_boundary, faces_boundary,
                               points_outside, tets_outside,
                               inside_pts);
    
	merge_inside_outside(points_interface, faces_interface,
                         points_inside, tets_inside,
                         points_outside, tets_outside,
                         points, tets, tet_flags);
}


#endif
