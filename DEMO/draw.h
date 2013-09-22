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
    
    class GLObject {
        
        GLuint shader;
        
        std::vector<DSC::vec3> data;
        
        GLuint array_id, buffer_id;
        GLuint position_att, normal_att;
        
        CGLA::Vec4f ambient_mat, diffuse_mat, specular_mat;
        
    public:
        
        GLObject(GLuint _shader, const CGLA::Vec4f& ambient_mat = CGLA::Vec4f(1.f), const CGLA::Vec4f& diffuse_mat = CGLA::Vec4f(0.f), const CGLA::Vec4f& specular_mat = CGLA::Vec4f(0.f));
        
        void add_data(std::vector<DSC::vec3> _data);
        
        void draw();
        
    };
    
    constexpr static float dist = 120.;
    
    std::unique_ptr<GLObject> interface, boundary, domain, tetrahedra;
    
    // Uniform variables
    CGLA::Mat4x4f modelViewProjectionMatrix, modelViewMatrix, modelMatrix, normalMatrix;
    CGLA::Vec3f light_pos = CGLA::Vec3f(0.f, 0.5*dist, dist);
    CGLA::Vec3f eye_pos = CGLA::Vec3f(0.3*dist, 0.3*dist, dist);
    
public:
    
    Painter(int WIN_SIZE_X, int WIN_SIZE_Y)
    {
        // Initialize uniforms
        CGLA::Mat4x4f projection = CGLA::perspective_Mat4x4f(53.f, WIN_SIZE_X/float(WIN_SIZE_Y), 0.01*dist, 3.*dist); // Projection matrix
        CGLA::Mat4x4f view = CGLA::lookAt_Mat4x4f(eye_pos, CGLA::Vec3f(0.), CGLA::Vec3f(0., 1., 0.)); // View matrix
        modelMatrix = CGLA::rotation_Mat4x4f(CGLA::YAXIS, M_PI);
        
        modelViewProjectionMatrix = projection * view * modelMatrix;
        modelViewMatrix = view * modelMatrix;
        normalMatrix = CGLA::invert_ortho(view * modelMatrix);
        
        // Initialize shader
        GLuint shader = init_gouraud_shader();
        
        interface = std::unique_ptr<GLObject>(new GLObject(shader, {0.1, 0.3, 0.1, 1.}, {0.5, 0.5, 0.5, 1.}, {0.3, 0.3, 0.3, 1.}));
        boundary = std::unique_ptr<GLObject>(new GLObject(shader, {0.3, 0.3, 0.3, 1.}, {0.3, 0.3, 0.3, 1.}, {0.3, 0.3, 0.3, 1.}));
        domain = std::unique_ptr<GLObject>(new GLObject(shader, {0.3, 0.3, 0.3, 0.3}, {0.3, 0.3, 0.3, 0.3}, {0.3, 0.3, 0.3, 0.3}));
        tetrahedra = std::unique_ptr<GLObject>(new GLObject(shader, {0.3, 0.1, 0.1, 0.1}, {0.6, 0.4, 0.4, 0.2}, {0., 0., 0., 0.}));
        
        // Enable states
        glEnable(GL_DEPTH_TEST);
        glDepthMask(GL_TRUE);

        glEnable(GL_CULL_FACE);
        
        glEnable(GL_BLEND);
        glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        
        check_gl_error();
    }
    
private:
    
    GLuint init_gouraud_shader();
    
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
     Updates the drawn boundary.
     */
    void update_boundary(DSC::DeformableSimplicialComplex<>& dsc);
    
    /**
     Updates the drawn design domain.
     */
    void update_design_domain(DSC::DeformableSimplicialComplex<>& dsc);
    
    /**
     Updates the drawn tetrahedra.
     */
    void update_tetrahedra(DSC::DeformableSimplicialComplex<>& dsc);
    
    /**
     Saves the current painting to the selected folder.
     */
    void save_painting(int width, int height, std::string folder = std::string(""), int time_step = -1);
};
