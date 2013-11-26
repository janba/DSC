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

#include "user_interface.h"
#include "rotate_function.h"
#include "average_function.h"
#include "normal_function.h"

#include "mesh_io.h"

using namespace DSC;

void display_(){
    UI::get_instance()->display();
}

void keyboard_(unsigned char key, int x, int y){
    UI::get_instance()->keyboard(key, x, y);
}

void reshape_(int width, int height){
    UI::get_instance()->reshape(width, height);
}

void visible_(int v){
    UI::get_instance()->visible(v);
}

void animate_(){
    UI::get_instance()->animate();
}

UI* UI::instance = NULL;

UI::UI(int &argc, char** argv)
{
    instance = this;

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_3_2_CORE_PROFILE | GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH | GLUT_MULTISAMPLE);
    glutCreateWindow("");
    
    glutDisplayFunc(display_);
    glutKeyboardFunc(keyboard_);
	glutIgnoreKeyRepeat(true);
    glutVisibilityFunc(visible_);
    glutReshapeFunc(reshape_);
	glutIdleFunc(animate_);
    
	glewExperimental = GL_TRUE;  // See http://www.opengl.org/wiki/OpenGL_Loading_Library
	GLint GlewInitResult = glewInit();
	if (GlewInitResult != GLEW_OK) {
		printf("ERROR: %s\n", glewGetErrorString(GlewInitResult));
	}
    check_gl_error(); // Catches a GL_INVALID_ENUM error. See http://www.opengl.org/wiki/OpenGL_Loading_Library
    
    // Read input
    std::string motion = "";
    
    if(argc > 1)
    {
        QUIT_ON_COMPLETION = true;
        CONTINUOUS = true;
        RECORD = true;
        
        for(int i = 0; i < argc; ++i)
        {
            std::string str(argv[i]);
            if (str == "nu") {
                VELOCITY = std::atof(argv[i+1]);
            }
            else if (str == "delta") {
                DISCRETIZATION = std::atof(argv[i+1]);
            }
            else if (str == "alpha") {
                ACCURACY = std::atof(argv[i+1]);
            }
            else if (str == "model") {
                model_file_name = argv[i+1];
            }
            else if (str == "motion") {
                motion = argv[i+1];
            }
        }
    }
    
    update_title();
	glutReshapeWindow(WIN_SIZE_X, WIN_SIZE_Y);
    check_gl_error();
    
    painter = std::unique_ptr<Painter>(new Painter(light_pos));
    load_model();
    keyboard(*motion.data(), 0, 0);
}

void UI::load_model()
{
    dsc = nullptr;
    std::vector<vec3> points;
    std::vector<int>  tets;
    std::vector<int>  tet_labels;
    is_mesh::import_tet_mesh(get_data_file_path(model_file_name + extension), points, tets, tet_labels);
    
    dsc = std::unique_ptr<DeformableSimplicialComplex<>>(new DeformableSimplicialComplex<>(DISCRETIZATION, points, tets, tet_labels));
    dsc->set_design_domain(new DesignDomain(DesignDomain::CUBE, vec3(50.)));
    dsc->scale(vec3(20.));
    painter->update(*dsc);
}

void UI::update_title()
{
    std::ostringstream oss;
    oss << "3D DSC\t";
    if(vel_fun)
    {
        oss << vel_fun->get_name();
        oss << ", Time step " << vel_fun->get_time_step();
    }
    oss << " (Nu = " << VELOCITY << ", Delta = " << DISCRETIZATION << ", Alpha = " << ACCURACY << ")";
    std::string str(oss.str());
    glutSetWindowTitle(str.c_str());
}

void UI::display()
{
    if (glutGet(GLUT_WINDOW_WIDTH) != WIN_SIZE_X || glutGet(GLUT_WINDOW_HEIGHT) != WIN_SIZE_Y) {
        return;
    }
    GLfloat timeValue = glutGet(GLUT_ELAPSED_TIME)*0.0002;
    vec3 ep = vec3( eye_pos[0] * sinf(timeValue), eye_pos[1] * cosf(timeValue) , eye_pos[2] * cosf(timeValue));
    painter->set_view_position(ep);
    painter->draw();
    glutSwapBuffers();
    update_title();
    check_gl_error();
}

void UI::reshape(int width, int height)
{
    WIN_SIZE_X = width;
    WIN_SIZE_Y = height;
    painter->reshape(width, height);
}

void UI::animate()
{
    if(vel_fun && CONTINUOUS)
    {
        vel_fun->take_time_step(*dsc);
        painter->update(*dsc);
        if(RECORD)
        {
            painter->set_view_position(camera_pos);
            painter->save_painting(basic_log->get_path(), vel_fun->get_time_step());
        }
        
        basic_log->write_timestep(*vel_fun, *dsc);
        if (vel_fun->is_motion_finished(*dsc))
        {
            stop();
            if (QUIT_ON_COMPLETION) {
                exit(0);
            }
        }
    }
    glutPostRedisplay();
}

void UI::keyboard(unsigned char key, int x, int y) {
    switch(key) {
        case '\033':
            stop();
            exit(0);
            break;
        case '0':
            stop();
            break;
        case '1':
            stop();
            vel_fun = std::unique_ptr<VelocityFunc<>>(new RotateFunc(VELOCITY, ACCURACY));
            start();
            break;
        case '2':
            stop();
            vel_fun = std::unique_ptr<VelocityFunc<>>(new AverageFunc(VELOCITY, ACCURACY));
            start();
            break;
        case '3':
            stop();
            vel_fun = std::unique_ptr<VelocityFunc<>>(new NormalFunc(VELOCITY, ACCURACY));
            start();
            break;
        case ' ':
            if(!CONTINUOUS)
            {
                std::cout << "MOTION STARTED" << std::endl;
            }
            else {
                std::cout << "MOTION PAUSED" << std::endl;
            }
            CONTINUOUS = !CONTINUOUS;
            break;
        case 'm':
            std::cout << "MOVE" << std::endl;
            if(vel_fun)
            {
                vel_fun->take_time_step(*dsc);
                painter->update(*dsc);
            }
            else {
                dsc->deform();
                painter->update(*dsc);
            }
            break;
        case 'r':
            std::cout << "RELOAD MODEL" << std::endl;
            load_model();
            break;
        case 't':
            if(dsc)
            {
                std::cout << "TEST DSC" << std::endl;
                dsc->test_flip23_flip32();
                painter->update(*dsc);
                dsc->test_split_collapse();
                painter->update(*dsc);
                dsc->test_flip44();
                painter->update(*dsc);
                dsc->test_flip22();
                painter->update(*dsc);
                if(vel_fun)
                {
                    std::cout << "TEST VELOCITY FUNCTION" << std::endl;
                    vel_fun->test(*dsc);
                    painter->update(*dsc);
                }
            }
            break;
        case '\t':
            if(dsc)
            {
                painter->switch_display_type();
                painter->update(*dsc);
            }
            break;
        case 's':
            if(dsc)
            {
                std::cout << "TAKING SCREEN SHOT" << std::endl;
                painter->set_view_position(camera_pos);
                painter->save_painting("LOG");
            }
            break;
        case 'e':
            if(dsc)
            {
                std::cout << "EXPORTING MESH" << std::endl;
                std::string filename("data/mesh.dsc");
                std::vector<vec3> points;
                std::vector<int> tets;
                std::vector<int> tet_labels;
                dsc->extract_tet_mesh(points, tets, tet_labels);
                is_mesh::export_tet_mesh(filename, points, tets, tet_labels);
            }
            break;
        case 'i':
            if(dsc)
            {
                std::cout << "EXPORTING SURFACE MESH" << std::endl;
                std::string filename("data/mesh.obj");
                std::vector<vec3> points;
                std::vector<int> faces;
                dsc->extract_surface_mesh(points, faces);
                is_mesh::export_surface_mesh(filename, points, faces);
            }
            break;
        case '+':
            if(dsc)
            {
                VELOCITY = std::min(VELOCITY + 1., 100.);
                if(vel_fun)
                {
                    vel_fun->set_velocity(VELOCITY);
                }
                update_title();
            }
            break;
        case '-':
            if(dsc)
            {
                VELOCITY = std::max(VELOCITY - 1., 0.);
                if(vel_fun)
                {
                    vel_fun->set_velocity(VELOCITY);
                }
                update_title();
            }
            break;
        case '.':
            if(dsc)
            {
                DISCRETIZATION = std::min(DISCRETIZATION + 0.5, 100.);
                dsc->set_discretization(DISCRETIZATION);
                update_title();
            }
            break;
        case ',':
            if(dsc)
            {
                DISCRETIZATION = std::max(DISCRETIZATION - 0.5, 1.);
                dsc->set_discretization(DISCRETIZATION);
                update_title();
            }
            break;
        case '<':
            if(dsc)
            {
                ACCURACY = std::min(ACCURACY + 1., 100.);
                if(vel_fun)
                {
                    vel_fun->set_accuracy(ACCURACY);
                }
                update_title();
            }
            break;
        case '>':
            if(dsc)
            {
                ACCURACY = std::max(ACCURACY - 1., 1.);
                if(vel_fun)
                {
                    vel_fun->set_accuracy(ACCURACY);
                }
                update_title();
            }
            break;
    }
}

void UI::visible(int v)
{
    if(v==GLUT_VISIBLE)
        glutIdleFunc(animate_);
    else
        glutIdleFunc(0);
}

void UI::stop()
{
    if(basic_log)
    {
        basic_log->write_message("MOTION STOPPED");
        basic_log->write_log(*dsc);
        if(vel_fun)
        {
            basic_log->write_log(*vel_fun);
            basic_log->write_timings(*vel_fun);
            std::string filename = basic_log->get_path() + std::string("/mesh.obj");
            std::vector<vec3> points;
            std::vector<int> tets;
            std::vector<int> tet_labels;
            dsc->extract_tet_mesh(points, tets, tet_labels);
            is_mesh::export_tet_mesh(filename, points, tets, tet_labels);
        }
    }
    basic_log = nullptr;
    vel_fun = nullptr;
}

void UI::start()
{
    basic_log = std::unique_ptr<Log>(new Log(create_log_path()));
    if(RECORD && vel_fun)
    {
        painter->set_view_position(camera_pos);
        painter->save_painting(basic_log->get_path(), vel_fun->get_time_step());
    }
    
    if(vel_fun)
    {
        basic_log->write_message(vel_fun->get_name().c_str());
        basic_log->write_log(*vel_fun);
    }
    basic_log->write_log(*dsc);
    
    update_title();
	glutReshapeWindow(WIN_SIZE_X, WIN_SIZE_Y);
    glutPostRedisplay();
}
