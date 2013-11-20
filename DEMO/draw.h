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
        
        std::vector<vec3> data;
        
        GLuint array_id, buffer_id;
        GLuint position_att, vector_att;
        
        CGLA::Vec4f ambient_mat, diffuse_mat, specular_mat;
        
    public:
        
        GLObject(GLuint _shader, const CGLA::Vec4f& ambient_mat = CGLA::Vec4f(1.f), const CGLA::Vec4f& diffuse_mat = CGLA::Vec4f(0.f), const CGLA::Vec4f& specular_mat = CGLA::Vec4f(0.f));
        
        void add_data(std::vector<vec3> _data);
        
        void clear_data()
        {
            data.clear();
        }
        
        void draw(GLenum mode = GL_TRIANGLES);
        
    };
    
    int WIDTH, HEIGHT;
    GLuint gouraud_shader, line_shader, wire_shader;
    
    std::unique_ptr<GLObject> interface, wire_frame, domain, low_quality, edges, unmoved;
    
    // Uniform variables
    CGLA::Mat4x4f projectionMatrix, viewMatrix, modelMatrix = CGLA::rotation_Mat4x4f(CGLA::YAXIS, M_PI);
    CGLA::Vec3f center = CGLA::Vec3f(0.);
    
public:
    
    Painter(const vec3& light_pos);
    
private:
    // Create a GLSL program object from vertex and fragment shader files
    GLuint init_shader(const char* vShaderFile, const char* fShaderFile, const char* outputAttributeName, const char* gShaderFile = nullptr);
    
public:
    /**
     Reshape the window.
     */
    void reshape(int width, int height);
    
    /**
     Set the position of the camera/eye.
     */
    void set_view_position(vec3 pos);
    
    /**
     Draws the simplicial complex.
     */
    void draw();
    
    /**
     Updates what to draw.
     */
    void update(DSC::DeformableSimplicialComplex<>& dsc, bool do_wire_frame);
    
    /**
     Saves the current painting to the selected folder.
     */
    void save_painting(std::string folder = std::string(""), int time_step = -1);
    
private:
    /**
     Updates the drawn interface.
     */
    void update_interface(DSC::DeformableSimplicialComplex<>& dsc);
    
    void update_wire_frame(DSC::DeformableSimplicialComplex<>& dsc);
    
    /**
     Updates the drawn edges.
     */
    void update_edges(DSC::DeformableSimplicialComplex<>& dsc);
    
    /**
     Updates the drawn boundary.
     */
    void update_boundary(DSC::DeformableSimplicialComplex<>& dsc);
    
    /**
     Updates the drawn domain.
     */
    void update_domain(DSC::DeformableSimplicialComplex<>& dsc);
    
    /**
     Updates the drawn tetrahedra.
     */
    void update_low_quality(DSC::DeformableSimplicialComplex<>& dsc);
    
    void update_unmoved(DSC::DeformableSimplicialComplex<>& dsc);
};
