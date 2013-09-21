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
    
    constexpr static float dist = 120.;
    const static unsigned int NULL_LOCATION = -1;
    
    GLuint gouraud_shader;
    
    // Interface variables:
    GLuint interface_array, interface_buffer;
    std::vector<DSC::vec3> interface_data;
    GLuint interface_position_att, interface_normal_att;
    
    // Boundary variables:
    GLuint boundary_array, boundary_buffer;
    std::vector<DSC::vec3> boundary_data;
    GLuint boundary_position_att, boundary_normal_att;
    
    // Design domain variables:
    GLuint domain_shader;
    GLuint domain_array, domain_buffer;
    std::vector<DSC::vec3> domain_data;
    GLuint domain_position_att, domain_normal_att;
    
    // Tetrahedra variables:
    GLuint tetrahedra_array, tetrahedra_buffer;
    std::vector<DSC::vec3> tetrahedra_data;
    GLuint tetrahedra_position_att, tetrahedra_normal_att;
    
    // Uniform variables
    CGLA::Mat4x4f modelViewProjectionMatrix, modelViewMatrix, normalMatrix;
    CGLA::Vec3f light_pos = CGLA::Vec3f(0.f, 0.5*dist, dist);
    CGLA::Vec3f eye_pos = CGLA::Vec3f(0.3*dist, 0.3*dist, dist);
    
public:
    
    Painter(int WIN_SIZE_X, int WIN_SIZE_Y)
    {
        // Initialize uniforms
        CGLA::Mat4x4f projection = CGLA::perspective_Mat4x4f(53.f, WIN_SIZE_X/float(WIN_SIZE_Y), 0.01*dist, 3.*dist); // Projection matrix
        CGLA::Mat4x4f view = CGLA::lookAt_Mat4x4f(eye_pos, CGLA::Vec3f(0.), CGLA::Vec3f(0., 1., 0.)); // View matrix
        CGLA::Mat4x4f model = CGLA::rotation_Mat4x4f(CGLA::YAXIS, M_PI);
        
        modelViewProjectionMatrix = projection * view * model;
        modelViewMatrix = view * model;
        normalMatrix = CGLA::invert_ortho(view * model);
        
        init_gouraud_shader();
        
        init_interface();
        init_boundary();
        init_tetrahedra();
        
        // Enable states
        glEnable(GL_DEPTH_TEST);
        glDepthMask(GL_TRUE);

        glEnable(GL_CULL_FACE);
        
        glEnable(GL_BLEND);
        glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        
        check_gl_error();
    }
    
private:
    
    void init_gouraud_shader();
    
    /**
     Initialize drawing of the interface.
     */
    void init_interface();
    
    /**
     Initialize drawing of the boundary.
     */
    void init_boundary();

    /**
     Initialize drawing of the design domain.
     */
    void init_domain();
    
    void init_tetrahedra();
    
    
    void use_solid_material()
    {
        CGLA::Vec4f ambientMat(0.1, 0.3, 0.1, 1.);
        CGLA::Vec4f diffuseMat(0.5, 0.5, 0.5, 1.);
        CGLA::Vec4f specMat(0.2, 0.2, 0.2, 1.);
        
        GLuint uniform = glGetUniformLocation(gouraud_shader, "ambientMat");
        if (uniform == NULL_LOCATION) {
            std::cerr << "Shader did not contain the 'ambientMat' uniform."<<std::endl;
        }
        glUniform4fv(uniform, 1, &ambientMat[0]);
        
        uniform = glGetUniformLocation(gouraud_shader, "diffuseMat");
        if (uniform == NULL_LOCATION) {
            std::cerr << "Shader did not contain the 'diffuseMat' uniform."<<std::endl;
        }
        glUniform4fv(uniform, 1, &diffuseMat[0]);
        
        uniform = glGetUniformLocation(gouraud_shader, "specMat");
        if (uniform == NULL_LOCATION) {
            std::cerr << "Shader did not contain the 'specMat' uniform."<<std::endl;
        }
        glUniform4fv(uniform, 1, &specMat[0]);
    }
    
    void use_transparent_material()
    {
        CGLA::Vec4f ambientMat(0.4, 0.2, 0.2, 0.1);
        CGLA::Vec4f diffuseMat(0.5, 0.4, 0.4, 0.2);
        CGLA::Vec4f specMat(0.0, 0.0, 0.0, 0.);
        
        GLuint uniform = glGetUniformLocation(gouraud_shader, "ambientMat");
        if (uniform == NULL_LOCATION) {
            std::cerr << "Shader did not contain the 'ambientMat' uniform."<<std::endl;
        }
        glUniform4fv(uniform, 1, &ambientMat[0]);
        
        uniform = glGetUniformLocation(gouraud_shader, "diffuseMat");
        if (uniform == NULL_LOCATION) {
            std::cerr << "Shader did not contain the 'diffuseMat' uniform."<<std::endl;
        }
        glUniform4fv(uniform, 1, &diffuseMat[0]);
        
        uniform = glGetUniformLocation(gouraud_shader, "specMat");
        if (uniform == NULL_LOCATION) {
            std::cerr << "Shader did not contain the 'specMat' uniform."<<std::endl;
        }
        glUniform4fv(uniform, 1, &specMat[0]);
    }
    
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
    void update_interface(DSC::DeformableSimplicialComplex<>& dsc);
    
    /**
     Updates the drawn design domain.
     */
    void update_boundary(DSC::DeformableSimplicialComplex<>& dsc);
    
    void update_design_domain(DSC::DeformableSimplicialComplex<>& dsc);
    
    void update_tetrahedra(DSC::DeformableSimplicialComplex<>& dsc);
    
    /**
     Saves the current painting to the selected folder.
     */
    void save_painting(int width, int height, std::string folder = std::string(""), int time_step = -1);
};
