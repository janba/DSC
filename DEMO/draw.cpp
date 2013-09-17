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

void Painter::load_shader()
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
    
    glUseProgram(shaderProgram);
    
    if(vertexdata.size() != 0)
    {
        glDrawArrays(GL_TRIANGLES, 0, static_cast<int>(vertexdata.size())/2);
    }
    
    glutSwapBuffers();
    check_gl_error();
}

void Painter::update_interface(DSC::DeformableSimplicialComplex<>& complex)
{
    // Extract interface data
    vertexdata.clear();
    for (auto fit = complex.faces_begin(); fit != complex.faces_end(); fit++)
    {
        if (fit->is_interface())
        {
            auto nodes = complex.get_nodes(fit.key());
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

