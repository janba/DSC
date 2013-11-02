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

void Painter::GLObject::add_data(std::vector<DSC::vec3> _data)
{
    data = std::vector<DSC::vec3>(_data);
    glBindBuffer(GL_ARRAY_BUFFER, buffer_id);
    glBufferData(GL_ARRAY_BUFFER, sizeof(DSC::vec3)*data.size(), &data[0], GL_STATIC_DRAW);
}

void Painter::GLObject::draw(GLenum mode)
{
    if(data.size() != 0)
    {
        glUseProgram(shader);
        GLuint uniform = glGetUniformLocation(shader, "ambientMat");
        if (uniform == NULL_LOCATION) {
            std::cerr << "Shader did not contain the 'ambientMat' uniform."<<std::endl;
        }
        glUniform4fv(uniform, 1, &ambient_mat[0]);
        
        uniform = glGetUniformLocation(shader, "diffuseMat");
        if (uniform == NULL_LOCATION) {
            std::cerr << "Shader did not contain the 'diffuseMat' uniform."<<std::endl;
        }
        glUniform4fv(uniform, 1, &diffuse_mat[0]);
        
        uniform = glGetUniformLocation(shader, "specMat");
        if (uniform == NULL_LOCATION) {
            std::cerr << "Shader did not contain the 'specMat' uniform."<<std::endl;
        }
        glUniform4fv(uniform, 1, &specular_mat[0]);
        
        glUseProgram(shader);
        glBindBuffer(GL_ARRAY_BUFFER, buffer_id);
        
        glVertexAttribPointer(position_att, 3, GL_DOUBLE, GL_FALSE, 2.*sizeof(DSC::vec3), (const GLvoid *)0);
        glVertexAttribPointer(vector_att, 3, GL_DOUBLE, GL_FALSE, 2.*sizeof(DSC::vec3), (const GLvoid *)sizeof(DSC::vec3));
        
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

Painter::Painter(const DSC::vec3& light_pos)
{    
    // Initialize shader
    line_shader = init_shader("shaders/line.vert", "shaders/line.frag", "fragColour", "shaders/line.geom");
    gouraud_shader = init_shader("shaders/gouraud.vert",  "shaders/gouraud.frag", "fragColour");
    
    // Send light position uniform to the shader
    CGLA::Vec3f lp(light_pos);
    GLuint lightPosUniform = glGetUniformLocation(gouraud_shader, "lightPos");
    if (lightPosUniform == NULL_LOCATION) {
        std::cerr << "Shader did not contain the 'lightPos' uniform."<<std::endl;
    }
    glUniform3fv(lightPosUniform, 1, &lp[0]);
    check_gl_error();
    
    interface = std::unique_ptr<GLObject>(new GLObject(gouraud_shader, {0.1, 0.3, 0.1, 1.}, {0.5, 0.5, 0.5, 1.}, {0.3, 0.3, 0.3, 1.}));
    domain = std::unique_ptr<GLObject>(new GLObject(gouraud_shader, {0.1, 0.1, 0.3, 1.}, {0.2, 0.2, 0.3, 1.}, {0., 0., 0., 1.}));
    low_quality = std::unique_ptr<GLObject>(new GLObject(gouraud_shader, {0.3, 0.1, 0.1, 0.1}, {0.6, 0.4, 0.4, 0.2}, {0., 0., 0., 0.}));
    edges = std::unique_ptr<GLObject>(new GLObject(line_shader, {0.7, 0.2, 0.2, 1.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}));
    unmoved = std::unique_ptr<GLObject>(new GLObject(line_shader, {0.2, 0.2, 0.7, 1.}, {0., 0., 0., 0.}, {0., 0., 0., 0.}));
    
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
        s << std::string(DSC::Util::concat4digits("_", time_step));
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
    projectionMatrix = CGLA::perspective_Mat4x4f(53.f, width/float(height), 1., 300.);
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
    check_gl_error();
}

void Painter::set_view_position(DSC::vec3 pos)
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
    check_gl_error();
}

void Painter::draw()
{
    glClearColor(0.7, 0.7, 0.7, 0.);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    interface->draw();
    
    domain->draw();
    
    glDisable(GL_CULL_FACE);
    edges->draw(GL_POINTS);
    unmoved->draw(GL_POINTS);
    low_quality->draw();
    glEnable(GL_CULL_FACE);
    
    check_gl_error();
}

void Painter::update(DSC::DeformableSimplicialComplex<>& dsc)
{
    update_interface(dsc);
//    update_unmoved(dsc);
    update_edges(dsc);
//    update_domain(dsc);
//    update_low_quality(dsc);
}

void Painter::update_interface(DSC::DeformableSimplicialComplex<>& dsc)
{
    std::vector<DSC::vec3> data;
    for (auto fit = dsc.faces_begin(); fit != dsc.faces_end(); fit++)
    {
        if (fit->is_interface())
        {
            auto verts = dsc.get_pos(dsc.get_sorted_nodes(fit.key()));
            DSC::vec3 normal = DSC::Util::normal_direction(verts[0], verts[1], verts[2]);
            
            for(auto &p : verts)
            {
                data.push_back(p);
                data.push_back(normal);
            }
        }
    }
    interface->add_data(data);
}

void Painter::update_edges(DSC::DeformableSimplicialComplex<>& dsc)
{
    std::vector<DSC::vec3> data;
    for (auto eit = dsc.edges_begin(); eit != dsc.edges_end(); eit++)
    {
        auto nids = dsc.get_nodes(eit.key());
        if(eit->is_interface())
        {
            DSC::vec3 vector = dsc.get_pos(nids[1]) - dsc.get_pos(nids[0]);
            
            data.push_back(dsc.get_pos(nids[0]));
            data.push_back(vector);
        }
    }
    edges->add_data(data);
}

bool is_boundary(DSC::DeformableSimplicialComplex<>& dsc, const DSC::DeformableSimplicialComplex<>::face_key& fid, std::vector<DSC::vec3>& verts)
{
    if(dsc.get(fid).is_boundary())
    {
        std::swap(verts[0], verts[2]);
        return true;
    }
    auto *design_domain = dsc.get_design_domain();
    if(design_domain)
    {
        for (auto &plane : design_domain->get_planes()) {
            bool boundary = true;
            for (auto &p : verts) {
                if (DSC::Util::is_inside(p, plane.p, plane.n)) {
                    boundary = false;
                    break;
                }
            }
            if(boundary)
            {
                auto apices = dsc.get_nodes(dsc.get_tets(fid)) - dsc.get_nodes(fid);
                if(DSC::Util::is_inside(dsc.get_pos(apices[0]), plane.p, plane.n) != DSC::Util::is_inside(dsc.get_pos(apices[1]), plane.p, plane.n))
                {
                    if(dot(dsc.get_normal(fid), plane.n) > 0.)
                    {
                        std::swap(verts[0], verts[2]);
                    }
                    return true;
                }
            }
        }
    }
    return false;
}

void Painter::update_domain(DSC::DeformableSimplicialComplex<>& dsc)
{
    std::vector<DSC::vec3> data;
    for (auto fit = dsc.faces_begin(); fit != dsc.faces_end(); fit++)
    {
        auto verts = dsc.get_pos(dsc.get_sorted_nodes(fit.key()));
        if(is_boundary(dsc, fit.key(), verts))
        {
            DSC::vec3 normal = DSC::Util::normal_direction(verts[0], verts[1], verts[2]);
            for (auto &p : verts) {
                data.push_back(p);
                data.push_back(normal);
            }
        }
    }
    domain->add_data(data);
}

void Painter::update_low_quality(DSC::DeformableSimplicialComplex<>& dsc)
{
    std::vector<DSC::vec3> data;
    for (auto tit = dsc.tetrahedra_begin(); tit != dsc.tetrahedra_end(); tit++)
    {
        if(dsc.quality(tit.key()) < dsc.get_min_tet_quality())
        {            
            for (auto f : dsc.get_faces(tit.key()))
            {
                auto nodes = dsc.get_sorted_nodes(f, tit.key());
                DSC::vec3 normal = -dsc.get_normal(f);
                
                for(auto &n : nodes)
                {
                    data.push_back(dsc.get_pos(n));
                    data.push_back(normal);
                }
            }
        }
    }
    for (auto fit = dsc.faces_begin(); fit != dsc.faces_end(); fit++)
    {
        if(dsc.quality(fit.key()) < dsc.get_min_face_quality())
        {
            auto nodes = dsc.get_nodes(fit.key());
            DSC::vec3 normal = dsc.get_normal(fit.key());
            
            for(auto &n : nodes)
            {
                data.push_back(dsc.get_pos(n));
                data.push_back(normal);
            }
        }
    }
    low_quality->add_data(data);
}

void Painter::update_unmoved(DSC::DeformableSimplicialComplex<>& dsc)
{
    std::vector<DSC::vec3> data;
    for (auto nit = dsc.nodes_begin(); nit != dsc.nodes_end(); nit++)
    {
        DSC::vec3 vector = nit->get_destination() - nit->get_pos();
        if(vector.length() > DSC::EPSILON)
        {
            data.push_back(nit->get_pos());
            data.push_back(vector);
        }
    }
    unmoved->add_data(data);
}

