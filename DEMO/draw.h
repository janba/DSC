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

#ifndef ___D_DSC__draw__
#define ___D_DSC__draw__


#ifdef WIN32
#include <GL/glut.h>
#include <GLGraphics/SOIL.h>
#else
#include <GLUT/glut.h>
#include <GEL/GLGraphics/SOIL.h>
#endif

const static float ALPHA = 0.2;
const static float BACKGROUND_COLOR[] = {0.7, 0.7, 0.7, 0.};
const static float INVISIBLE[] = {-1., -1., -1.};
const static float DARK_RED[] = {0.66,0.11,0.15, ALPHA};
const static float RED[] = {0.96,0.11,0.15, ALPHA};
const static float YELLOW[] = {0.9,0.9,0., ALPHA};
const static float DARK_BLUE[] = {0.14,0.16,0.88, ALPHA};
const static float BLUE[] = {0.45,0.7,0.9, ALPHA};
const static float GREEN[] = {0.05,1.,0.15, ALPHA};
const static float ORANGE[] = {0.9,0.4,0., ALPHA};
const static float BLACK[] = {0., 0., 0.};
const static float DARK_GRAY[] = {0.5, 0.5, 0.5, ALPHA};
const static float GRAY[] = {0.8, 0.8, 0.8, ALPHA};

const static float POINT_SIZE = 0.5;
const static float LINE_WIDTH = 0.1;

/**
 A painter handles all draw functionality using OpenGL.
 */
template<class MT>
class Painter {
    
    typedef typename MT::real_type      T;
    typedef typename MT::vector3_type   V;
    
public:
    
    Painter()
    {
        
    }
    
    /**
     Saves the current painting to the selected folder.
     */
    static void save_painting(int width, int height, std::string folder = std::string(""), int time_step = -1)
    {
        std::ostringstream s;
        if (folder.length() == 0) {
            s << "scr";
        }
        else {
            s << folder << "/scr";
        }
        
        if (time_step >= 0)
        {
            s << std::string(Util::concat4digits("_", time_step));
        }
        s << ".png";
        int success = SOIL_save_screenshot(s.str().c_str(), SOIL_SAVE_TYPE_PNG, 0, 0, width, height);
        if(!success)
        {
            std::cout << "ERROR: Failed to take screen shot: " << s.str().c_str() << std::endl;
            return;
        }
    }
    
    /**
     Begins drawing.
     */
    static void begin()
    {
        glClearColor(BACKGROUND_COLOR[0], BACKGROUND_COLOR[1], BACKGROUND_COLOR[2],BACKGROUND_COLOR[3]);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    }
    
    /**
     Ends drawing.
     */
    static void end()
    {
        glFlush();
        glutSwapBuffers();
    }
    
    /**
     Draws the simplicial complex.
     */
    static void draw_complex(DeformableSimplicialComplex<MT> *complex)
    {
        draw_faces(complex);
//        draw_nodes(complex);
        draw_edges(complex);
        glDisable(GL_LIGHTING);
        draw_bad_tetrahedra(complex);
        glEnable(GL_LIGHTING);
    }
    
    /**
     Draws the bad tetrahedra.
     */
    static void draw_bad_tetrahedra(DeformableSimplicialComplex<MT> *complex)
    {
        bool low_quality, small_angle;
        std::vector<V> verts;
        glLineWidth(LINE_WIDTH);
        glBegin(GL_TRIANGLES);
        for (auto tit = complex->tetrahedra_begin(); tit != complex->tetrahedra_end(); tit++)
        {
            low_quality = complex->quality(tit.key()) < complex->get_min_tet_quality();
            small_angle = complex->min_dihedral_angle(tit.key()) < complex->get_min_angle();
            if (low_quality || small_angle)
            {
                if (low_quality && small_angle) {
                    glColor4fv(&RED[0]);
                }
                else if(low_quality)
                {
                    glColor4fv(&YELLOW[0]);
                }
                else {
                    glColor4fv(&GREEN[0]);
                }
                
                typename DeformableSimplicialComplex<MT>::simplex_set cl_t;
                complex->closure(tit.key(), cl_t);
                
                for (auto fit = cl_t.faces_begin(); fit != cl_t.faces_end(); fit++)
                {
                    V n = complex->get_normal(*fit);
                    complex->get_pos(*fit, verts);
                    glVertex3dv(&verts[0][0]);
                    glVertex3dv(&verts[1][0]);
                    glVertex3dv(&verts[2][0]);
                    glNormal3dv(&n[0]);
                }
            }
        }
        glEnd();
        
    }
    
    /**
     Draws the faces with the colors defined by the get_face_colors function in the simplicial complex.
     */
    static void draw_faces(DeformableSimplicialComplex<MT> *complex, const float color[] = GRAY)
    {
        glColor3fv(&color[0]);
        glBegin(GL_TRIANGLES);
        for (auto fit = complex->faces_begin(); fit != complex->faces_end(); fit++)
        {
            if (fit->is_interface())
            {
                std::vector<V> verts;
                complex->get_pos(fit.key(), verts);
                V n = complex->get_normal(fit.key());
                
                for (int i = 0; i < 3; ++i)
                {
                    glNormal3f(n[0], n[1], n[2]);
                    glVertex3f(verts[i][0], verts[i][1], verts[i][2]);
                }
            }
        }
        
        glEnd();
    }
    
    static void draw_edges(DeformableSimplicialComplex<MT> *complex, const float color[] = BLACK)
    {
        glColor3fv(&color[0]);
        glLineWidth(LINE_WIDTH);
        V p1, p2;
        glBegin(GL_LINES);
        for(auto eit = complex->edges_begin(); eit != complex->edges_end(); eit++)
        {
            std::vector<V> verts;
            complex->get_pos(eit.key(), verts);
            if(eit->is_interface())
            {
                glVertex3dv(&verts[0][0]);
                glVertex3dv(&verts[1][0]);
            }
        }
        glEnd();
    }

    static void draw_nodes(DeformableSimplicialComplex<MT> *complex, const float color[] = BLACK)
    {
        glColor3fv(&color[0]);
        glPointSize(POINT_SIZE);
        glBegin(GL_POINTS);
        V p;
        for(auto nit = complex->nodes_begin(); nit != complex->nodes_end(); nit++)
        {
            if (nit->is_interface()) {
                p = complex->get_pos(nit.key());
                glVertex3dv(&p[0]);
            }
        }
        glEnd();
    }

    
    /**
     Draws the domain.
     */
//    static void draw_domain(const Domain *domain, CGLA::Vec3d color = GRAY);
//    
//    /**
//     Draws the vertices with the colors defined by the get_vertex_colors function in the simplicial complex.
//     */
//    static void draw_vertices(const SimplicialComplex<MT> *complex);
//
//    /**
//     Draws the edges with the colors defined by the get_edge_colors function in the simplicial complex.
//     */
//    static void draw_edges(const SimplicialComplex<MT> *complex);
//    
//
//    /**
//     Draws the faces with the colors given as input.
//     */
//    static void draw_faces(const SimplicialComplex<MT> *complex, const HMesh::FaceAttributeVector<CGLA::Vec3d> &colors);
//    
//    /**
//     Draws the interface with the color given as input.
//     */
//    static void draw_interface(const SimplicialComplex<MT> *complex, CGLA::Vec3d color = ORANGE);
//    
//    /**
//     Draws the arrows given as input with the color given as input.
//     */
//    static void draw_arrows(const SimplicialComplex<MT> *complex, const HMesh::VertexAttributeVector<CGLA::Vec2d> &arrows, CGLA::Vec3d color = ORANGE);
//
//    /**
//     Draws the lines given as input with the color given as input.
//     */
//    static void draw_lines(const SimplicialComplex<MT> *complex, const HMesh::VertexAttributeVector<CGLA::Vec2d> &lines, CGLA::Vec3d color = GREEN);
};


#endif /* defined(___D_DSC__draw__) */
