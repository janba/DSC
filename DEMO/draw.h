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

#include <GL/glew.h>
#ifdef WIN32
#include <GL/glut.h>
#else
#include <GLUT/glut.h>
#endif

#include <SOIL/SOIL.h>
#include <CGLA/Mat4x4f.h>

#include "DSC.h"

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


inline void _check_gl_error(const char *file, int line)
{
    GLenum err (glGetError());
    
    while(err!=GL_NO_ERROR) {
        std::string error;
        
        switch(err) {
            case GL_INVALID_OPERATION:      error="INVALID_OPERATION";      break;
            case GL_INVALID_ENUM:           error="INVALID_ENUM";           break;
            case GL_INVALID_VALUE:          error="INVALID_VALUE";          break;
            case GL_OUT_OF_MEMORY:          error="OUT_OF_MEMORY";          break;
            case GL_INVALID_FRAMEBUFFER_OPERATION:  error="INVALID_FRAMEBUFFER_OPERATION";  break;
        }
        
        std::cerr << "GL_" << error.c_str() <<" - "<<file<<":"<<line<<std::endl;
        err=glGetError();
    }
}

#define check_gl_error() _check_gl_error(__FILE__,__LINE__)

/**
 A painter handles all draw functionality using OpenGL.
 */
class Painter {
    
    const static unsigned int NULL_LOCATION = -1;
    
    GLuint interface_shader;
    GLuint interface_array, interface_buffer;
    std::vector<DSC::vec3> interface_data;
    
    GLuint interface_position_att, interface_normal_att;
    
    GLuint MVMatrixUniform, MVPMatrixUniform, NormalMatrixUniform, lightPosUniform;
    
public:
    
    Painter(int WIN_SIZE_X, int WIN_SIZE_Y, double r)
    {
        load_shader();
        
        // Generate arrays and buffers
        glGenVertexArrays(1, &interface_array);
        glBindVertexArray(interface_array);
        
        glGenBuffers(1, &interface_buffer);
        glBindBuffer(GL_ARRAY_BUFFER, interface_buffer);
        
        glEnableVertexAttribArray(interface_position_att);
        glEnableVertexAttribArray(interface_normal_att);
        
        // Set up model view projection matrix
        CGLA::Mat4x4f projection = CGLA::perspective_Mat4x4f(53.f, WIN_SIZE_X/float(WIN_SIZE_Y), 0.01*r, 3.*r); // Projection matrix
        CGLA::Mat4x4f view = CGLA::lookAt_Mat4x4f(CGLA::Vec3f(0.3*r, 0.3*r, r), CGLA::Vec3f(0.), CGLA::Vec3f(0., 1., 0.)); // View matrix
        CGLA::Mat4x4f model = CGLA::rotation_Mat4x4f(CGLA::YAXIS, M_PI);
        
        // Send model view projection matrix
        CGLA::Mat4x4f modelViewProjection = projection * view * model;
        glUniformMatrix4fv(MVPMatrixUniform, 1, GL_TRUE, &modelViewProjection[0][0]);
        
        // Send model view matrix
        CGLA::Mat4x4f modelView = view * model;
        glUniformMatrix4fv(MVMatrixUniform, 1, GL_TRUE, &modelView[0][0]);
        
        // Send normal matrix
        CGLA::Mat4x4f normalMatrix = CGLA::invert_ortho(view * model);
        glUniformMatrix4fv(NormalMatrixUniform, 1, GL_FALSE, &normalMatrix[0][0]);
        
        // Send light position
        CGLA::Vec3f light_pos(0., 0.5*r, r);
        glUniform3fv(lightPosUniform, 1, &light_pos[0]);
        
        // Enable states
        glEnable(GL_DEPTH_TEST);
        glDepthMask(GL_TRUE);

//        glEnable(GL_CULL_FACE);
//        glCullFace(GL_BACK);

        glEnable (GL_BLEND);
        glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        
        check_gl_error();
    }
    
private:
    /**
     Loads the shaders.
     */
    void load_shader();
        
    /**
     Draws the bad tetrahedra.
     */
    void draw_bad_tetrahedra(DSC::DeformableSimplicialComplex<> *complex)
    {
        bool low_quality, small_angle;
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
                
                typename DSC::DeformableSimplicialComplex<>::simplex_set cl_t;
                complex->closure(tit.key(), cl_t);
                
                for (auto fit = cl_t.faces_begin(); fit != cl_t.faces_end(); fit++)
                {
                    DSC::vec3 n = complex->get_normal(*fit);
                    auto verts = complex->get_pos(*fit);
                    glVertex3d(static_cast<double>(verts[0][0]), static_cast<double>(verts[0][1]), static_cast<double>(verts[0][2]));
                    glVertex3d(static_cast<double>(verts[1][0]), static_cast<double>(verts[1][1]), static_cast<double>(verts[1][2]));
                    glVertex3d(static_cast<double>(verts[2][0]), static_cast<double>(verts[2][1]), static_cast<double>(verts[2][2]));
                    
                    glNormal3d(static_cast<double>(n[0]), static_cast<double>(n[1]), static_cast<double>(n[2]));
                }
            }
        }
        glEnd();
        
    }
    
    /**
     Draws the faces with the colors defined by the get_face_colors function in the simplicial complex.
     */
    void draw_faces(DSC::DeformableSimplicialComplex<> *complex, const float color[] = GRAY)
    {
        glColor3fv(&color[0]);
        glBegin(GL_TRIANGLES);
        for (auto fit = complex->faces_begin(); fit != complex->faces_end(); fit++)
        {
            if (fit->is_interface())
            {
                auto verts = complex->get_pos(fit.key());
                DSC::vec3 n = complex->get_normal(fit.key());
                
                for (int i = 0; i < 3; ++i)
                {
                    glNormal3f(n[0], n[1], n[2]);
                    glVertex3f(verts[i][0], verts[i][1], verts[i][2]);
                }
            }
        }
        
        glEnd();
    }
    
    void draw_edges(DSC::DeformableSimplicialComplex<> *complex, const float color[] = BLACK)
    {
        glColor3fv(&color[0]);
        glLineWidth(LINE_WIDTH);
        DSC::vec3 p1, p2;
        glBegin(GL_LINES);
        for(auto eit = complex->edges_begin(); eit != complex->edges_end(); eit++)
        {
            auto verts = complex->get_pos(eit.key());
            if(eit->is_interface())
            {
                glVertex3d(static_cast<double>(verts[0][0]), static_cast<double>(verts[0][1]), static_cast<double>(verts[0][2]));
                glVertex3d(static_cast<double>(verts[1][0]), static_cast<double>(verts[1][1]), static_cast<double>(verts[1][2]));
            }
        }
        glEnd();
    }
    
    void draw_nodes(DSC::DeformableSimplicialComplex<> *complex, const float color[] = BLACK)
    {
        glColor3fv(&color[0]);
        glPointSize(POINT_SIZE);
        glBegin(GL_POINTS);
        DSC::vec3 p;
        for(auto nit = complex->nodes_begin(); nit != complex->nodes_end(); nit++)
        {
            if (nit->is_interface()) {
                p = complex->get_pos(nit.key());
                glVertex3dv(&p[0]);
            }
        }
        glEnd();
    }
    
public:
    /**
     Draws the simplicial complex.
     */
    void draw();
    
    /**
     Updates the drawn interface.
     */
    void update_interface(DSC::DeformableSimplicialComplex<>& complex);
    
    
    /**
     Saves the current painting to the selected folder.
     */
    void save_painting(int width, int height, std::string folder = std::string(""), int time_step = -1);
};
