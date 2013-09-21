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
//
#include "draw.h"


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

void Painter::init_gouraud_shader()
{
    // Load the interface shader
    gouraud_shader = InitShader("shaders/gouraud.vert",  "shaders/gouraud.frag", "fragColour");
    
    // Send uniforms to the shader
    GLuint MVMatrixUniform = glGetUniformLocation(gouraud_shader, "MVMatrix");
    if (MVMatrixUniform == NULL_LOCATION) {
        std::cerr << "Shader did not contain the 'MVMatrix' uniform."<<std::endl;
    }
    glUniformMatrix4fv(MVMatrixUniform, 1, GL_TRUE, &modelViewMatrix[0][0]);
    
    GLuint MVPMatrixUniform = glGetUniformLocation(gouraud_shader, "MVPMatrix");
    if (MVPMatrixUniform == NULL_LOCATION) {
        std::cerr << "Shader did not contain the 'MVPMatrix' uniform."<<std::endl;
    }
    glUniformMatrix4fv(MVPMatrixUniform, 1, GL_TRUE, &modelViewProjectionMatrix[0][0]);
    
    GLuint NormalMatrixUniform = glGetUniformLocation(gouraud_shader, "NormalMatrix");
    if (NormalMatrixUniform == NULL_LOCATION) {
        std::cerr << "Shader did not contain the 'NormalMatrix' uniform."<<std::endl;
    }
    glUniformMatrix4fv(NormalMatrixUniform, 1, GL_FALSE, &normalMatrix[0][0]);
    
    GLuint lightPosUniform = glGetUniformLocation(gouraud_shader, "lightPos");
    if (lightPosUniform == NULL_LOCATION) {
        std::cerr << "Shader did not contain the 'lightPos' uniform."<<std::endl;
    }
    glUniform3fv(lightPosUniform, 1, &light_pos[0]);
    check_gl_error();
    
    GLuint eyePosUniform = glGetUniformLocation(gouraud_shader, "eyePos");
    if (eyePosUniform == NULL_LOCATION) {
        std::cerr << "Shader did not contain the 'eyePos' uniform."<<std::endl;
    }
    glUniform3fv(eyePosUniform, 1, &eye_pos[0]);
    
}

void Painter::init_interface()
{
    // Generate arrays and buffers for visualising the interface
    glGenVertexArrays(1, &interface_array);
    glBindVertexArray(interface_array);
    
    glGenBuffers(1, &interface_buffer);
    glBindBuffer(GL_ARRAY_BUFFER, interface_buffer);
    check_gl_error();
    
    // Initialize shader attributes
    interface_position_att = glGetAttribLocation(gouraud_shader, "position");
    if (interface_position_att == NULL_LOCATION) {
        std::cerr << "Shader did not contain the 'position' attribute." << std::endl;
    }
    interface_normal_att = glGetAttribLocation(gouraud_shader, "normal");
    if (interface_normal_att == NULL_LOCATION) {
        std::cerr << "Shader did not contain the 'normal' attribute." << std::endl;
    }
    
    glEnableVertexAttribArray(interface_position_att);
    glEnableVertexAttribArray(interface_normal_att);
    check_gl_error();
}

void Painter::init_boundary()
{
    // Generate arrays and buffers for visualising the boundary
    glGenVertexArrays(1, &boundary_array);
    glBindVertexArray(boundary_array);
    
    glGenBuffers(1, &boundary_buffer);
    glBindBuffer(GL_ARRAY_BUFFER, boundary_buffer);
    
    // Initialize shader attributes
    boundary_position_att = glGetAttribLocation(gouraud_shader, "position");
    if (boundary_position_att == NULL_LOCATION) {
        std::cerr << "Shader did not contain the 'position' attribute." << std::endl;
    }
    boundary_normal_att = glGetAttribLocation(gouraud_shader, "normal");
    if (boundary_normal_att == NULL_LOCATION) {
        std::cerr << "Shader did not contain the 'normal' attribute." << std::endl;
    }
    
    glEnableVertexAttribArray(boundary_position_att);
    glEnableVertexAttribArray(boundary_normal_att);
    check_gl_error();
}

void Painter::init_tetrahedra()
{
    // Generate arrays and buffers for visualising the boundary
    glGenVertexArrays(1, &tetrahedra_array);
    glBindVertexArray(tetrahedra_array);
    
    glGenBuffers(1, &tetrahedra_buffer);
    glBindBuffer(GL_ARRAY_BUFFER, tetrahedra_buffer);
    
    // Initialize shader attributes
    tetrahedra_position_att = glGetAttribLocation(gouraud_shader, "position");
    if (tetrahedra_position_att == NULL_LOCATION) {
        std::cerr << "Shader did not contain the 'position' attribute." << std::endl;
    }
    tetrahedra_normal_att = glGetAttribLocation(gouraud_shader, "normal");
    if (tetrahedra_normal_att == NULL_LOCATION) {
        std::cerr << "Shader did not contain the 'normal' attribute." << std::endl;
    }
    
    glEnableVertexAttribArray(tetrahedra_position_att);
    glEnableVertexAttribArray(tetrahedra_normal_att);
    check_gl_error();
}

void Painter::init_domain()
{
    // Generate arrays and buffers for visualising the boundary
    glGenVertexArrays(1, &domain_array);
    glBindVertexArray(domain_array);
    
    glGenBuffers(1, &domain_buffer);
    glBindBuffer(GL_ARRAY_BUFFER, domain_buffer);
    
    // Initialize shader attributes
    domain_position_att = glGetAttribLocation(gouraud_shader, "position");
    if (boundary_position_att == NULL_LOCATION) {
        std::cerr << "Shader did not contain the 'position' attribute." << std::endl;
    }
    domain_normal_att = glGetAttribLocation(gouraud_shader, "normal");
    if (domain_normal_att == NULL_LOCATION) {
        std::cerr << "Shader did not contain the 'normal' attribute." << std::endl;
    }
    
    glEnableVertexAttribArray(domain_position_att);
    glEnableVertexAttribArray(domain_normal_att);
    check_gl_error();
}

void Painter::save_painting(int width, int height, std::string folder, int time_step)
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

void Painter::draw()
{
    glClearColor(BACKGROUND_COLOR[0], BACKGROUND_COLOR[1], BACKGROUND_COLOR[2],BACKGROUND_COLOR[3]);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    glUseProgram(gouraud_shader);
    use_solid_material();
    
    glCullFace(GL_BACK);
    if(interface_data.size() != 0)
    {
        glBindBuffer(GL_ARRAY_BUFFER, interface_buffer);
        
        glVertexAttribPointer(interface_position_att, 3, GL_DOUBLE, GL_FALSE, 2.*sizeof(DSC::vec3), (const GLvoid *)0);
        glVertexAttribPointer(interface_normal_att, 3, GL_DOUBLE, GL_FALSE, 2.*sizeof(DSC::vec3), (const GLvoid *)sizeof(DSC::vec3));
        
        glDrawArrays(GL_TRIANGLES, 0, static_cast<int>(interface_data.size())/2);
    }
    
    if(tetrahedra_data.size() != 0)
    {
        use_transparent_material();
        glDisable(GL_CULL_FACE);
        glBindBuffer(GL_ARRAY_BUFFER, tetrahedra_buffer);
        
        glVertexAttribPointer(tetrahedra_position_att, 3, GL_DOUBLE, GL_FALSE, 2.*sizeof(DSC::vec3), (const GLvoid *)0);
        glVertexAttribPointer(tetrahedra_normal_att, 3, GL_DOUBLE, GL_FALSE, 2.*sizeof(DSC::vec3), (const GLvoid *)sizeof(DSC::vec3));
        
        glDrawArrays(GL_TRIANGLES, 0, static_cast<int>(tetrahedra_data.size())/2);
        glEnable(GL_CULL_FACE);
        use_solid_material();
    }
    
    glCullFace(GL_FRONT);
    if(boundary_data.size() != 0)
    {
        glBindBuffer(GL_ARRAY_BUFFER, boundary_buffer);
        
        glVertexAttribPointer(boundary_position_att, 3, GL_DOUBLE, GL_FALSE, 2.*sizeof(DSC::vec3), (const GLvoid *)0);
        glVertexAttribPointer(boundary_normal_att, 3, GL_DOUBLE, GL_FALSE, 2.*sizeof(DSC::vec3), (const GLvoid *)sizeof(DSC::vec3));
        
        glDrawArrays(GL_TRIANGLES, 0, static_cast<int>(boundary_data.size())/2);
    }
    
    if(domain_data.size() != 0)
    {
        glBindBuffer(GL_ARRAY_BUFFER, domain_buffer);
        
        glVertexAttribPointer(domain_position_att, 3, GL_DOUBLE, GL_FALSE, 2.*sizeof(DSC::vec3), (const GLvoid *)0);
        glVertexAttribPointer(domain_normal_att, 3, GL_DOUBLE, GL_FALSE, 2.*sizeof(DSC::vec3), (const GLvoid *)sizeof(DSC::vec3));
        
        glDrawArrays(GL_TRIANGLES, 0, static_cast<int>(domain_data.size())/2);
    }
    
    glutSwapBuffers();
    check_gl_error();
}

void Painter::update_interface(DSC::DeformableSimplicialComplex<>& dsc)
{
    interface_data.clear();
    for (auto fit = dsc.faces_begin(); fit != dsc.faces_end(); fit++)
    {
        if (fit->is_interface())
        {
            auto nodes = dsc.get_nodes(fit.key());
            DSC::vec3 normal = dsc.get_normal(fit.key());
            
            for(auto &n : nodes)
            {
                interface_data.push_back(dsc.get_pos(n));
                interface_data.push_back(normal);
            }
        }
    }
    glBindBuffer(GL_ARRAY_BUFFER, interface_buffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(DSC::vec3)*interface_data.size(), &interface_data[0], GL_STATIC_DRAW);
}

void Painter::update_boundary(DSC::DeformableSimplicialComplex<>& dsc)
{
    boundary_data.clear();
    for (auto fit = dsc.faces_begin(); fit != dsc.faces_end(); fit++)
    {
        if (fit->is_boundary())
        {
            auto nodes = dsc.get_nodes(fit.key());
            DSC::vec3 normal = -dsc.get_normal(fit.key());
            
            for(auto &n : nodes)
            {
                boundary_data.push_back(dsc.get_pos(n));
                boundary_data.push_back(normal);
            }
        }
    }
    glBindBuffer(GL_ARRAY_BUFFER, boundary_buffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(DSC::vec3)*boundary_data.size(), &boundary_data[0], GL_STATIC_DRAW);
}


void Painter::update_design_domain(DSC::DeformableSimplicialComplex<>& dsc)
{
    domain_data.clear();
    const DSC::DesignDomain *domain = dsc.get_design_domain();
    
    glBindBuffer(GL_ARRAY_BUFFER, domain_buffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(DSC::vec3)*domain_data.size(), &domain_data[0], GL_STATIC_DRAW);
}

void Painter::update_tetrahedra(DSC::DeformableSimplicialComplex<>& dsc)
{
    tetrahedra_data.clear();
    for (auto tit = dsc.tetrahedra_begin(); tit != dsc.tetrahedra_end(); tit++)
    {
        bool low_quality = dsc.quality(tit.key()) < dsc.get_min_tet_quality();
        bool small_angle = dsc.min_dihedral_angle(tit.key()) < dsc.get_min_angle();
        if(low_quality || small_angle)
        {
            
            typename DSC::DeformableSimplicialComplex<>::simplex_set cl_t;
            dsc.closure(tit.key(), cl_t);
            
            for (auto fit = cl_t.faces_begin(); fit != cl_t.faces_end(); fit++)
            {
                auto nodes = dsc.get_nodes(*fit);
                DSC::vec3 normal = dsc.get_normal(*fit);
                
                for(auto &n : nodes)
                {
                    tetrahedra_data.push_back(dsc.get_pos(n));
                    tetrahedra_data.push_back(normal);
                }
            }
        }
    }
    glBindBuffer(GL_ARRAY_BUFFER, tetrahedra_buffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(DSC::vec3)*tetrahedra_data.size(), &tetrahedra_data[0], GL_STATIC_DRAW);
}

