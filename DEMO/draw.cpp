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

Painter::GLObject::GLObject(GLuint _shader, const CGLA::Vec4f& ambient_mat_, const CGLA::Vec4f& diffuse_mat_, const CGLA::Vec4f& specular_mat_) : shader(_shader), ambient_mat(ambient_mat_), diffuse_mat(diffuse_mat_), specular_mat(specular_mat_)
{
    // Generate arrays and buffers for visualising the interface
    glGenVertexArrays(1, &array_id);
    glBindVertexArray(array_id);
    
    glGenBuffers(1, &buffer_id);
    glBindBuffer(GL_ARRAY_BUFFER, buffer_id);
    
    // Initialize shader attributes
    position_att = glGetAttribLocation(shader, "position");
    if (position_att == NULL_LOCATION) {
        std::cerr << "Shader did not contain the 'position' attribute." << std::endl;
    }
    vector_att = glGetAttribLocation(shader, "vector");
    if (vector_att == NULL_LOCATION) {
        std::cerr << "Shader did not contain the 'vector' attribute." << std::endl;
    }
    
    glEnableVertexAttribArray(position_att);
    glEnableVertexAttribArray(vector_att);
    check_gl_error();
}

void Painter::GLObject::add_data(std::vector<vec3> _data)
{
    data = std::vector<vec3>(_data);
    glBindBuffer(GL_ARRAY_BUFFER, buffer_id);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vec3)*data.size(), &data[0], GL_STATIC_DRAW);
}

void Painter::GLObject::draw(GLenum mode)
{
    if(data.size() != 0)
    {
        glUseProgram(shader);
        GLuint uniform = (GLuint) glGetUniformLocation(shader, "ambientMat");
        if (uniform == NULL_LOCATION) {
            std::cerr << "Shader did not contain the 'ambientMat' uniform."<<std::endl;
        }
        glUniform4fv(uniform, 1, &ambient_mat[0]);
        
        uniform = (GLuint) glGetUniformLocation(shader, "diffuseMat");
        if (uniform == NULL_LOCATION) {
            std::cerr << "Shader did not contain the 'diffuseMat' uniform."<<std::endl;
        }
        glUniform4fv(uniform, 1, &diffuse_mat[0]);
        
        uniform = (GLuint) glGetUniformLocation(shader, "specMat");
        if (uniform == NULL_LOCATION) {
            std::cerr << "Shader did not contain the 'specMat' uniform."<<std::endl;
        }
        glUniform4fv(uniform, 1, &specular_mat[0]);
        
        glUseProgram(shader);
        glBindBuffer(GL_ARRAY_BUFFER, buffer_id);
        
        glVertexAttribPointer(position_att, 3, GL_DOUBLE, GL_FALSE, 2*sizeof(vec3), (const GLvoid *)0);
        glVertexAttribPointer(vector_att, 3, GL_DOUBLE, GL_FALSE, 2*sizeof(vec3), (const GLvoid *)sizeof(vec3));
        
        glDrawArrays(mode, 0, static_cast<int>(data.size())/2);
        check_gl_error();
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
GLuint Painter::init_shader(const char* vShaderFile, const char* fShaderFile, const char* outputAttributeName, const char* gShaderFile)
{
    struct Shader {
        const char*  filename;
        GLenum       type;
        GLchar*      source;
    };
    std::vector<Shader> shaders = { { vShaderFile, GL_VERTEX_SHADER, NULL }, { fShaderFile, GL_FRAGMENT_SHADER, NULL } };
    if(gShaderFile)
    {
        shaders.push_back({ gShaderFile, GL_GEOMETRY_SHADER, NULL });
    }
    
    GLuint program = glCreateProgram();
    
    for ( int i = 0; i < shaders.size(); ++i ) {
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

Painter::Painter(const vec3& light_pos)
{    
    // Initialize shader
    wire_shader = init_shader("shaders/gouraud.vert", "shaders/wire.frag", "fragColour", "shaders/wire.geom");
    line_shader = init_shader("shaders/line.vert", "shaders/line.frag", "fragColour", "shaders/line.geom");
    gouraud_shader = init_shader("shaders/gouraud.vert",  "shaders/gouraud.frag", "fragColour");
    
    // Send light position uniform to the shader
    glUseProgram(gouraud_shader);
    CGLA::Vec3f lp(light_pos);
    GLuint lightPosUniform = glGetUniformLocation(gouraud_shader, "lightPos");
    if (lightPosUniform == NULL_LOCATION) {
        std::cerr << "Shader did not contain the 'lightPos' uniform."<<std::endl;
    }
    glUniform3fv(lightPosUniform, 1, &lp[0]);
    
    glUseProgram(wire_shader);
    lightPosUniform = glGetUniformLocation(wire_shader, "lightPos");
    if (lightPosUniform == NULL_LOCATION) {
        std::cerr << "Shader did not contain the 'lightPos' uniform."<<std::endl;
    }
    glUniform3fv(lightPosUniform, 1, &lp[0]);
    
    CGLA::Vec4f wire_col(0.25,0.5,0.6, 1.);
    GLuint wireColUniform = glGetUniformLocation(wire_shader, "wireCol");
    if (wireColUniform == NULL_LOCATION) {
        std::cerr << "Shader did not contain the 'wireCol' uniform."<<std::endl;
    }
    glUniform4fv(wireColUniform, 1, &wire_col[0]);
    
    check_gl_error();
    
    interface = std::unique_ptr<GLObject>(new GLObject(gouraud_shader, {0.15f,0.4f,0.5f, 1.f}, {0.2f, 0.3f, 0.4f, 1.f}, {0.2f, 0.3f, 0.4f, 1.f}));
    wire_frame = std::unique_ptr<GLObject>(new GLObject(wire_shader, {0.15f,0.4f,0.5f, 1.f}, {0.2f, 0.3f, 0.4f, 1.f}, {0.2f, 0.3f, 0.4f, 1.f}));
    domain = std::unique_ptr<GLObject>(new GLObject(gouraud_shader, {0.1f, 0.1f, 0.3f, 1.f}, {0.2f, 0.2f, 0.3f, 1.f}, {0.f, 0.f, 0.f, 1.f}));
    low_quality = std::unique_ptr<GLObject>(new GLObject(gouraud_shader, {0.3f, 0.1f, 0.1f, 0.1f}, {0.6f, 0.4f, 0.4f, 0.2f}, {0.f, 0.f, 0.f, 0.f}));
    edges = std::unique_ptr<GLObject>(new GLObject(line_shader, {0.f,0.f,0.f, 1.f}, {0.f, 0.f, 0.f, 0.f}, {0.f, 0.f, 0.f, 0.f}));
    unmoved = std::unique_ptr<GLObject>(new GLObject(line_shader, {0.2f, 0.2f, 0.7f, 1.f}, {0.f, 0.f, 0.f, 0.f}, {0.f, 0.f, 0.f, 0.f}));
    
    // Enable states
    glEnable(GL_DEPTH_TEST);
    glDepthMask(GL_TRUE);
    
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    
    glEnable(GL_BLEND);
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    check_gl_error();
}

void Painter::save_painting(std::string folder, int time_step)
{
    draw();
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
    int success = SOIL_save_screenshot(s.str().c_str(), SOIL_SAVE_TYPE_PNG, 0, 0, WIDTH, HEIGHT);
    if(!success)
    {
        std::cout << "ERROR: Failed to take screen shot: " << s.str().c_str() << std::endl;
        return;
    }
}

void Painter::reshape(int width, int height)
{
    WIDTH = width;
    HEIGHT = height;
    projectionMatrix = CGLA::perspective_Mat4x4f(53.f, width/float(height), 1., 1000000.);
    glViewport(0, 0, width, height);
    CGLA::Mat4x4f modelViewProjectionMatrix = projectionMatrix * viewMatrix * modelMatrix;
    
    glUseProgram(gouraud_shader);
    GLuint MVPMatrixUniform = glGetUniformLocation(gouraud_shader, "MVPMatrix");
    if (MVPMatrixUniform == NULL_LOCATION) {
        std::cerr << "Shader did not contain the 'MVPMatrix' uniform."<<std::endl;
    }
    glUniformMatrix4fv(MVPMatrixUniform, 1, GL_TRUE, &modelViewProjectionMatrix[0][0]);
    
    glUseProgram(line_shader);
    GLuint PMatrixUniform = glGetUniformLocation(line_shader, "PMatrix");
    if (PMatrixUniform == NULL_LOCATION) {
        std::cerr << "Shader did not contain the 'PMatrix' uniform."<<std::endl;
    }
    glUniformMatrix4fv(PMatrixUniform, 1, GL_TRUE, &projectionMatrix[0][0]);
    
    glUseProgram(wire_shader);
    MVPMatrixUniform = glGetUniformLocation(wire_shader, "MVPMatrix");
    if (MVPMatrixUniform == NULL_LOCATION) {
        std::cerr << "Shader did not contain the 'MVPMatrix' uniform."<<std::endl;
    }
    glUniformMatrix4fv(MVPMatrixUniform, 1, GL_TRUE, &modelViewProjectionMatrix[0][0]);
    
    GLuint scaleUniform = glGetUniformLocation(wire_shader, "WScale");
    if (scaleUniform == NULL_LOCATION) {
        std::cerr << "Shader did not contain the 'WScale' uniform."<<std::endl;
    }
    glUniform2f(scaleUniform, width, height);
    check_gl_error();
}

void Painter::set_view_position(vec3 pos)
{
    viewMatrix = CGLA::lookAt_Mat4x4f(CGLA::Vec3f(pos), center, CGLA::Vec3f(0., 1., 0.));
    CGLA::Mat4x4f modelViewMatrix = viewMatrix * modelMatrix;
    CGLA::Mat4x4f normalMatrix = CGLA::invert_ortho(modelViewMatrix);
    CGLA::Mat4x4f modelViewProjectionMatrix = projectionMatrix * modelViewMatrix;
    
    glUseProgram(gouraud_shader);
    GLuint MVMatrixUniform = glGetUniformLocation(gouraud_shader, "MVMatrix");
    if (MVMatrixUniform == NULL_LOCATION) {
        std::cerr << "Shader did not contain the 'MVMatrix' uniform."<<std::endl;
    }
    glUniformMatrix4fv(MVMatrixUniform, 1, GL_TRUE, &modelViewMatrix[0][0]);
    
    GLuint NormalMatrixUniform = glGetUniformLocation(gouraud_shader, "NormalMatrix");
    if (NormalMatrixUniform == NULL_LOCATION) {
        std::cerr << "Shader did not contain the 'NormalMatrix' uniform."<<std::endl;
    }
    glUniformMatrix4fv(NormalMatrixUniform, 1, GL_FALSE, &normalMatrix[0][0]);
    
    GLuint MVPMatrixUniform = glGetUniformLocation(gouraud_shader, "MVPMatrix");
    if (MVPMatrixUniform == NULL_LOCATION) {
        std::cerr << "Shader did not contain the 'MVPMatrix' uniform."<<std::endl;
    }
    glUniformMatrix4fv(MVPMatrixUniform, 1, GL_TRUE, &modelViewProjectionMatrix[0][0]);
    
    glUseProgram(line_shader);
    MVMatrixUniform = glGetUniformLocation(line_shader, "MVMatrix");
    if (MVMatrixUniform == NULL_LOCATION) {
        std::cerr << "Shader did not contain the 'MVMatrix' uniform."<<std::endl;
    }
    glUniformMatrix4fv(MVMatrixUniform, 1, GL_TRUE, &modelViewMatrix[0][0]);
    
    NormalMatrixUniform = glGetUniformLocation(line_shader, "NormalMatrix");
    if (NormalMatrixUniform == NULL_LOCATION) {
        std::cerr << "Shader did not contain the 'NormalMatrix' uniform."<<std::endl;
    }
    glUniformMatrix4fv(NormalMatrixUniform, 1, GL_FALSE, &normalMatrix[0][0]);
    
    glUseProgram(wire_shader);
    MVMatrixUniform = glGetUniformLocation(wire_shader, "MVMatrix");
    if (MVMatrixUniform == NULL_LOCATION) {
        std::cerr << "Shader did not contain the 'MVMatrix' uniform."<<std::endl;
    }
    glUniformMatrix4fv(MVMatrixUniform, 1, GL_TRUE, &modelViewMatrix[0][0]);
    
    MVPMatrixUniform = glGetUniformLocation(wire_shader, "MVPMatrix");
    if (MVPMatrixUniform == NULL_LOCATION) {
        std::cerr << "Shader did not contain the 'MVPMatrix' uniform."<<std::endl;
    }
    glUniformMatrix4fv(MVPMatrixUniform, 1, GL_TRUE, &modelViewProjectionMatrix[0][0]);
    
    NormalMatrixUniform = glGetUniformLocation(wire_shader, "NormalMatrix");
    if (NormalMatrixUniform == NULL_LOCATION) {
        std::cerr << "Shader did not contain the 'NormalMatrix' uniform."<<std::endl;
    }
    glUniformMatrix4fv(NormalMatrixUniform, 1, GL_FALSE, &normalMatrix[0][0]);
    
    check_gl_error();
}

void Painter::draw()
{
    glClearColor(1., 1., 1., 0.);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    interface->draw();
    wire_frame->draw();
    
    domain->draw();
    
    glDisable(GL_CULL_FACE);
    edges->draw(GL_POINTS);
    unmoved->draw(GL_POINTS);
    low_quality->draw();
    glEnable(GL_CULL_FACE);
    
    check_gl_error();
}

void Painter::switch_display_type()
{
    display_type = static_cast<DISPLAY_TYPE>((display_type+1)%DISPLAY_TYPE_END);
}

void Painter::update(DSC::DeformableSimplicialComplex& dsc)
{
    interface->clear_data();
    wire_frame->clear_data();
    edges->clear_data();
    domain->clear_data();
    low_quality->clear_data();
    unmoved->clear_data();
    switch (display_type) {
        case INTERFACE:
            update_interface(dsc);
            break;
        case WIRE_FRAME:
            update_wire_frame(dsc);
            break;
        case BOUNDARY:
            update_interface(dsc);
            update_domain(dsc);
            break;
        case EDGES:
            update_wire_frame(dsc);
            update_edges(dsc);
            break;
        case LOW_QUALITY:
            update_interface(dsc);
            update_low_quality(dsc);
            break;
        case UNMOVED:
            update_interface(dsc);
            update_unmoved(dsc);
            break;
            
        default:
            break;
    }
    
}

void Painter::update_interface(DSC::DeformableSimplicialComplex& dsc)
{
    std::vector<vec3> data;
    for (auto fit : dsc.get_is_mesh().faces())
    {
        if (fit->is_interface())
        {
            auto verts = dsc.get_is_mesh().get_pos(dsc.get_is_mesh().get_sorted_nodes(fit.key()));
            vec3 normal = Util::normal_direction(verts[0], verts[1], verts[2]);
            
            for(auto &p : verts)
            {
                data.push_back(p);
                data.push_back(normal);
            }
        }
    }
    interface->add_data(data);
}

void Painter::update_wire_frame(DSC::DeformableSimplicialComplex& dsc)
{
    std::vector<vec3> data;
    for (auto fit : dsc.get_is_mesh().faces())
    {
        if (fit->is_interface())
        {
            auto verts = dsc.get_is_mesh().get_pos(dsc.get_is_mesh().get_sorted_nodes(fit.key()));
            vec3 normal = Util::normal_direction(verts[0], verts[1], verts[2]);
            
            for(auto &p : verts)
            {
                data.push_back(p);
                data.push_back(normal);
            }
        }
    }
    wire_frame->add_data(data);
}

void Painter::update_edges(DSC::DeformableSimplicialComplex& dsc)
{
    auto & mesh = dsc.get_is_mesh();
    std::vector<vec3> data;
    for (auto eit : dsc.get_is_mesh().edges())
    {
        if(!eit->is_interface())
        {
            auto nids = eit->node_keys();
            vec3 vector = mesh.get(nids[1]).get_pos() - mesh.get(nids[0]).get_pos();
            
            data.push_back(mesh.get(nids[0]).get_pos());
            data.push_back(vector);
        }
    }
    edges->add_data(data);
}

void Painter::update_domain(DSC::DeformableSimplicialComplex& dsc)
{
    auto & mesh = dsc.get_is_mesh();
    std::vector<vec3> data;
    for (auto fit : mesh.faces())
    {
        if(fit->is_boundary())
        {
            auto verts = mesh.get_pos(mesh.get_sorted_nodes(fit.key()));
            std::swap(verts[0], verts[2]);
            vec3 normal = Util::normal_direction(verts[0], verts[1], verts[2]);
            for (auto &p : verts) {
                data.push_back(p);
                data.push_back(normal);
            }
        }
    }
    domain->add_data(data);
}

void Painter::update_low_quality(DSC::DeformableSimplicialComplex& dsc)
{
    auto & mesh = dsc.get_is_mesh();
    std::vector<vec3> data;
    for (auto tit : mesh.tetrahedra())
    {
        if(tit->quality() < dsc.get_min_tet_quality())
        {
            for (auto f : tit->face_keys())
            {
                auto nodes = mesh.get_sorted_nodes(f, tit.key());
                vec3 normal = -dsc.get_normal(f);
                
                for(auto &n : nodes)
                {
                    data.push_back(dsc.get_is_mesh().get(n).get_pos());
                    data.push_back(normal);
                }
            }
        }
    }
    for (auto fit : mesh.faces())
    {
        if(fit->quality() < dsc.get_min_face_quality())
        {
            auto nodes = mesh.get_nodes(fit.key());
            vec3 normal = dsc.get_normal(fit.key());
            
            for(auto &n : nodes)
            {
                data.push_back(mesh.get(n).get_pos());
                data.push_back(normal);
            }
        }
    }
    low_quality->add_data(data);
}

void Painter::update_unmoved(DSC::DeformableSimplicialComplex& dsc) {
    std::vector<vec3> data;
    data.push_back(vec3(0.));
    data.push_back(vec3(20., 0., 0.));
    data.push_back(vec3(0.));
    data.push_back(vec3(0., 20., 0.));
    data.push_back(vec3(0.));
    data.push_back(vec3(0., 0., 20.));

    for (auto nit : dsc.get_is_mesh().nodes())
    {
        vec3 vector = nit->get_destination() - nit->get_pos();
        if(vector.length() > EPSILON)
        {
            data.push_back(nit->get_pos());
            data.push_back(vector);
        }
    }
    unmoved->add_data(data);
}

