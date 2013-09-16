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

#include <SOIL/SOIL.h>

#ifdef WIN32
#include <GL/glut.h>
#else
#include <GLUT/glut.h>
#endif

#include <CGLA/Mat4x4f.h>

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
    
    GLuint shaderProgram;
    
    GLuint VertexArrayID;
    GLuint vertexbuffer;
    
    std::vector<DSC::vec3> vertexdata;
    
    GLuint MVMatrixUniform, MVPMatrixUniform, NormalMatrixUniform, lightPosUniform;
    
    GLuint positionAttribute, normalAttribute;
    
public:
    
    Painter()
    {
        glEnable(GL_DEPTH_TEST);
        glDepthMask(GL_TRUE);
        
        glEnable(GL_CULL_FACE);
        glCullFace(GL_BACK);
        
        glEnable(GL_LIGHTING);
        glEnable(GL_LIGHT0);
        glShadeModel(GL_FLAT);
        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);
        
        glEnable (GL_BLEND);
        glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        
        check_gl_error();
    }
    
    Painter(int WIN_SIZE_X, int WIN_SIZE_Y, double r)
    {
        load_shader();
        
        // Test data
        vertexdata.clear();
        vertexdata.push_back(DSC::vec3(-0.5, -0.5, 0.));
        vertexdata.push_back(DSC::vec3(0., 0., 1.));
        
        vertexdata.push_back(DSC::vec3(0.5, -0.5, 0.));
        vertexdata.push_back(DSC::vec3(0., 0., 1.));
        
        vertexdata.push_back(DSC::vec3(0., 0.5, 0.));
        vertexdata.push_back(DSC::vec3(0., 0., 1.));
        
        vertexdata.push_back(DSC::vec3(-0.5, -0.5, 0.));
        vertexdata.push_back(DSC::vec3(0.3, 0., 1.));
        
        vertexdata.push_back(DSC::vec3(-0.5, 0.5, 0.3));
        vertexdata.push_back(DSC::vec3(0.3, 0., 1.));
        
        vertexdata.push_back(DSC::vec3(0., 0.5, 0.));
        vertexdata.push_back(DSC::vec3(0.3, 0., 1.));
        
        vertexdata.push_back(DSC::vec3(-0.5, -0.5, 0.));
        vertexdata.push_back(DSC::vec3(0., 0.3, 1.));
        
        vertexdata.push_back(DSC::vec3(-0.5, 0.5, 0.3));
        vertexdata.push_back(DSC::vec3(0.3, 0.3, 1.));
        
        vertexdata.push_back(DSC::vec3(-1., 0.5, 0.));
        vertexdata.push_back(DSC::vec3(0., 0.3, 1.));
        
        // Generate arrays and buffers
        glGenVertexArrays(1, &VertexArrayID);
        glBindVertexArray(VertexArrayID);
        
        glGenBuffers(1, &vertexbuffer);
        glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
        glBufferData(GL_ARRAY_BUFFER, sizeof(DSC::vec3)*vertexdata.size(), &vertexdata[0], GL_STATIC_DRAW); // Visualize test data
        
        glEnableVertexAttribArray(positionAttribute);
        glEnableVertexAttribArray(normalAttribute);
        
        glVertexAttribPointer(positionAttribute, 3, GL_DOUBLE, GL_FALSE, 2.*sizeof(DSC::vec3), (const GLvoid *)0); // Visualize test data
        glVertexAttribPointer(normalAttribute, 3, GL_DOUBLE, GL_FALSE, 2.*sizeof(DSC::vec3), (const GLvoid *)sizeof(DSC::vec3)); // Visualize test data
        
        // Set up model view projection matrix
        CGLA::Mat4x4f projection = CGLA::perspective_Mat4x4f(53.f, WIN_SIZE_X/float(WIN_SIZE_Y), 0.01*r, 3.*r); // Projection matrix
        CGLA::Mat4x4f view = CGLA::lookAt_Mat4x4f(CGLA::Vec3f(0., 0., r), CGLA::Vec3f(0.), CGLA::Vec3f(0., 1., 0.)); // View matrix
        CGLA::Mat4x4f model = CGLA::identity_Mat4x4f();
        
        // Send model view projection matrix
        CGLA::Mat4x4f modelViewProjection = projection * view * model;
        glUniformMatrix4fv(MVPMatrixUniform, 1, GL_FALSE, &modelViewProjection[0][0]);
        
        // Send model view matrix
        CGLA::Mat4x4f modelView = view * model;
        glUniformMatrix4fv(MVMatrixUniform, 1, GL_FALSE, &modelView[0][0]);
        
        // Send normal matrix
        CGLA::Mat4x4f normalMatrix = CGLA::transpose(CGLA::invert_ortho(view * model));
        glUniformMatrix4fv(NormalMatrixUniform, 1, GL_FALSE, &normalMatrix[0][0]);
        
        // Send light position
        CGLA::Vec3f light_pos(0., 0.5*r, r);
        glUniform3fv(lightPosUniform, 1, &light_pos[0]);
        
        // Enable states
        glEnable(GL_DEPTH_TEST);
        glDepthMask(GL_TRUE);

//        glEnable(GL_CULL_FACE);
//        glCullFace(GL_BACK);
//        
//        glEnable (GL_BLEND);
//        glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        
        check_gl_error();
    }
    
    /**
     Saves the current painting to the selected folder.
     */
    void save_painting(int width, int height, std::string folder = std::string(""), int time_step = -1)
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
            s << std::string(DSC::Util::concat4digits("_", time_step));
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
    void begin()
    {
        glClearColor(BACKGROUND_COLOR[0], BACKGROUND_COLOR[1], BACKGROUND_COLOR[2],BACKGROUND_COLOR[3]);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    }
    
    /**
     Ends drawing.
     */
    void end()
    {
        glFlush();
        glutSwapBuffers();
    }
    
    /**
     Draws the simplicial complex.
     */
    void draw_complex(DSC::DeformableSimplicialComplex<> *complex)
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
    
    // NEW OPENGL STUFF
    void draw_new()
    {
        glClearColor(BACKGROUND_COLOR[0], BACKGROUND_COLOR[1], BACKGROUND_COLOR[2],BACKGROUND_COLOR[3]);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        glUseProgram(shaderProgram);
        
        if(vertexdata.size() != 0)
        {
            glDrawArrays(GL_TRIANGLES, 0, static_cast<int>(vertexdata.size())/2);
        }
        
        glutSwapBuffers();
        check_gl_error();
    }
    
    void update_interface(DSC::DeformableSimplicialComplex<>& complex, const float color[] = GRAY)
    {
        // Extract interface data
        vertexdata.clear();
        for (auto fit = complex.faces_begin(); fit != complex.faces_end(); fit++)
        {
            if (fit->is_interface())
            {
                auto nodes = complex.get_nodes(fit.key());
                
//                auto verts = complex.get_pos(fit.key());
                DSC::vec3 normal = complex.get_normal(fit.key());
                
                for(auto &n : nodes)
                {
                    vertexdata.push_back(complex.get_pos(n));
                    vertexdata.push_back(normal);
                }
            }
        }
        
        // Send interface data to shader
        glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
        glBufferData(GL_ARRAY_BUFFER, sizeof(DSC::vec3)*vertexdata.size(), &vertexdata[0], GL_STATIC_DRAW);
        
        glEnableVertexAttribArray(positionAttribute);
        glEnableVertexAttribArray(normalAttribute);
        glVertexAttribPointer(positionAttribute, 3, GL_DOUBLE, GL_FALSE, 2.*sizeof(DSC::vec3), (const GLvoid *)0);
        glVertexAttribPointer(normalAttribute, 3, GL_DOUBLE, GL_FALSE, 2.*sizeof(DSC::vec3), (const GLvoid *)sizeof(DSC::vec3));
    }
    
    void load_shader()
    {
        shaderProgram = InitShader("shaders/interface.vert",  "shaders/interface.frag", "fragColour");
        MVMatrixUniform = glGetUniformLocation(shaderProgram, "MVMatrix");
        if (MVMatrixUniform > 10000) {
            std::cerr << "Shader did not contain the 'MVMatrix' uniform."<<std::endl;
        }
        MVPMatrixUniform = glGetUniformLocation(shaderProgram, "MVPMatrix");
        if (MVPMatrixUniform > 10000) {
            std::cerr << "Shader did not contain the 'MVPMatrix' uniform."<<std::endl;
        }
        NormalMatrixUniform = glGetUniformLocation(shaderProgram, "NormalMatrix");
        if (NormalMatrixUniform > 10000) {
            std::cerr << "Shader did not contain the 'NormalMatrix' uniform."<<std::endl;
        }
        lightPosUniform = glGetUniformLocation(shaderProgram, "lightPos");
        if (lightPosUniform > 10000) {
            std::cerr << "Shader did not contain the 'lightPos' uniform."<<std::endl;
        }
        positionAttribute = glGetAttribLocation(shaderProgram, "position");
        if (positionAttribute > 10000) {
            std::cerr << "Shader did not contain the 'position' attribute." << std::endl;
        }
        normalAttribute = glGetAttribLocation(shaderProgram, "normal");
        if (normalAttribute > 10000) {
            std::cerr << "Shader did not contain the 'normal' attribute." << std::endl;
        }
    }
    
    // Create a NULL-terminated string by reading the provided file
	char* readShaderSource(const char* shaderFile)
	{
		FILE *filePointer;
		char *content = NULL;
        
		int count=0;
        
		if (shaderFile != NULL) {
			filePointer = fopen(shaderFile,"rt");
            
			if (filePointer != NULL) {
                
				fseek(filePointer, 0, SEEK_END);
				count = static_cast<int>(ftell(filePointer));
				rewind(filePointer);
                
				if (count > 0) {
					content = (char *)malloc(sizeof(char) * (count+1));
					count = static_cast<int>(fread(content,sizeof(char),count,filePointer));
					content[count] = '\0';
				}
				fclose(filePointer);
			}
		}
		return content;
	}
    
	// Create a GLSL program object from vertex and fragment shader files
	GLuint InitShader(const char* vShaderFile, const char* fShaderFile, const char* outputAttributeName)
    {
		struct Shader {
			const char*  filename;
			GLenum       type;
			GLchar*      source;
		}  shaders[2] = {
			{ vShaderFile, GL_VERTEX_SHADER, NULL },
			{ fShaderFile, GL_FRAGMENT_SHADER, NULL }
		};
        
		GLuint program = glCreateProgram();
        
		for ( int i = 0; i < 2; ++i ) {
			Shader& s = shaders[i];
			s.source = readShaderSource( s.filename );
			if ( shaders[i].source == NULL ) {
				std::cerr << "Failed to read " << s.filename << std::endl;
				exit( EXIT_FAILURE );
			}
			GLuint shader = glCreateShader( s.type );
			glShaderSource( shader, 1, (const GLchar**) &s.source, NULL );
			glCompileShader( shader );
            
			GLint  compiled;
			glGetShaderiv( shader, GL_COMPILE_STATUS, &compiled );
			if ( !compiled ) {
				std::cerr << s.filename << " failed to compile:" << std::endl;
				GLint  logSize;
				glGetShaderiv( shader, GL_INFO_LOG_LENGTH, &logSize );
				char* logMsg = new char[logSize];
				glGetShaderInfoLog( shader, logSize, NULL, logMsg );
				std::cerr << logMsg << std::endl;
				delete [] logMsg;
                
				exit( EXIT_FAILURE );
			}
            
			delete [] s.source;
            
			glAttachShader( program, shader );
		}
        
		/* Link output */
		glBindFragDataLocation(program, 0, outputAttributeName);
        
		/* link  and error check */
		glLinkProgram(program);
        
		GLint  linked;
		glGetProgramiv( program, GL_LINK_STATUS, &linked );
		if ( !linked ) {
			std::cerr << "Shader program failed to link" << std::endl;
			GLint  logSize;
			glGetProgramiv( program, GL_INFO_LOG_LENGTH, &logSize);
			char* logMsg = new char[logSize];
			glGetProgramInfoLog( program, logSize, NULL, logMsg );
			std::cerr << logMsg << std::endl;
			delete [] logMsg;
            
			exit( EXIT_FAILURE );
		}
		
		/* use program object */
		glUseProgram(program);
        
		return program;
	}
    
};
