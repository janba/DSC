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

#include <iostream>
#include <iomanip>
#include <ctime>
#include <chrono>

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
#ifdef _WIN32
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH | GLUT_MULTISAMPLE);
#else
    glutInitDisplayMode(GLUT_3_2_CORE_PROFILE | GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH | GLUT_MULTISAMPLE);
#endif
    glutCreateWindow("");
    
    glutDisplayFunc(display_);
    glutKeyboardFunc(keyboard_);
	glutIgnoreKeyRepeat(true);
    glutVisibilityFunc(visible_);
    glutReshapeFunc(reshape_);
	glutIdleFunc(animate_);

#ifndef __APPLE__
	glewExperimental = GL_TRUE;  // See http://www.opengl.org/wiki/OpenGL_Loading_Library
	GLint GlewInitResult = glewInit();
	if (GlewInitResult != GLEW_OK) {
		printf("ERROR: %s\n", glewGetErrorString(GlewInitResult));
	}
    check_gl_error(); // Catches a GL_INVALID_ENUM error. See http://www.opengl.org/wiki/OpenGL_Loading_Library
#endif
    
    // Read input
    std::string motion = "";
    double discretization = 2.5;
    double velocity = 5.;
    double accuracy = 0.25;
    
    if(argc == 2)
    {
        model_file_name = std::string(argv[1]);
    }
    else if(argc > 2)
    {
        for(int i = 0; i < argc; ++i)
        {
            std::string str(argv[i]);
            if (str == "nu") {
                velocity = std::atof(argv[i+1]);
            }
            else if (str == "delta") {
                discretization = std::atof(argv[i+1]);
            }
            else if (str == "alpha") {
                accuracy = std::atof(argv[i+1]);
            }
            else if (str == "model") {
                model_file_name = argv[i+1];
            }
            else if (str == "motion") {
                motion = argv[i+1];
            }
        }
    }
    painter = std::unique_ptr<Painter>(new Painter(light_pos));
    load_model(model_file_name, discretization);
    
    if(motion.empty())
    {
        vel_fun = std::unique_ptr<VelocityFunc>(new VelocityFunc(velocity, accuracy, 500));
        start("");
    }
    else {
        keyboard(*motion.data(), 0, 0);
    }
    
	glutReshapeWindow(WIN_SIZE_X, WIN_SIZE_Y);
    check_gl_error();
}

void UI::load_model(const std::string& file_name, double discretization)
{
    std::cout << "\nLoading " << obj_path + file_name + ".dsc" << std::endl;
    dsc = nullptr;
    std::vector<vec3> points;
    std::vector<int>  tets;
    std::vector<int>  tet_labels;
    is_mesh::import_tet_mesh(obj_path + file_name + ".dsc", points, tets, tet_labels);
    
    dsc = std::unique_ptr<DeformableSimplicialComplex>(new DeformableSimplicialComplex(points, tets, tet_labels));
    
    vec3 p_min(INFINITY), p_max(-INFINITY);
    for (auto & nit : dsc->get_is_mesh().nodes()) {
        auto pos = nit.get_pos();
        for (int i = 0; i < 3; i++) {
            p_min[i] = Util::min(pos[i], p_min[i]);
            p_max[i] = Util::max(pos[i], p_max[i]);
        }
    }
    
    vec3 size = p_max - p_min;
    double var = Util::max(Util::max(size[0], size[1]), size[2]);
    double dist = 1.2*var;
    eye_pos = {dist, var, dist};
    camera_pos = {var, var, -dist};
    light_pos = {0., 0., dist};
    
    painter->update(*dsc);
    std::cout << "Loading done" << std::endl << std::endl;
}

void UI::update_title()
{
    std::ostringstream oss;
    oss << "3D DSC\t" << vel_fun->get_name() << ", Time step " << vel_fun->get_time_step();
    oss << " (Nu = " << vel_fun->get_velocity() << ", Delta = " << dsc->get_avg_edge_length() << ", Alpha = " << vel_fun->get_accuracy() << ")";
    std::string str(oss.str());
    glutSetWindowTitle(str.c_str());
}

void UI::display()
{
    if (glutGet(GLUT_WINDOW_WIDTH) != WIN_SIZE_X || glutGet(GLUT_WINDOW_HEIGHT) != WIN_SIZE_Y) {
        return;
    }
    GLfloat timeValue = glutGet(GLUT_ELAPSED_TIME)*0.0002;
    vec3 ep( eye_pos[0] * sinf(timeValue), eye_pos[1] * cosf(timeValue) , eye_pos[2] * cosf(timeValue));
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
    if(CONTINUOUS)
    {
        std::cout << "\n***************TIME STEP " << vel_fun->get_time_step() + 1 <<  " START*************\n" << std::endl;
        vel_fun->take_time_step(*dsc);
        painter->update(*dsc);
        if(RECORD && basic_log)
        {
            painter->set_view_position(camera_pos);
            painter->save_painting(basic_log->get_path(), vel_fun->get_time_step());
            basic_log->write_timestep(*vel_fun, *dsc);
        }
        if (vel_fun->is_motion_finished(*dsc))
        {
            stop();
            if (QUIT_ON_COMPLETION) {
                exit(0);
            }
        }
        std::cout << "\n***************TIME STEP " << vel_fun->get_time_step() <<  " STOP*************\n" << std::endl;
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
            QUIT_ON_COMPLETION = false;
            RECORD = false;
            vel_fun = std::unique_ptr<VelocityFunc>(new VelocityFunc(vel_fun->get_velocity(), vel_fun->get_accuracy(), 500));
            start("");
            break;
        case '1':
            stop();
            QUIT_ON_COMPLETION = true;
            RECORD = true;
            vel_fun = std::unique_ptr<VelocityFunc>(new RotateFunc(vel_fun->get_velocity(), vel_fun->get_accuracy()));
            start("rotate");
            break;
        case '2':
            stop();
            QUIT_ON_COMPLETION = true;
            RECORD = true;
            vel_fun = std::unique_ptr<VelocityFunc>(new AverageFunc(vel_fun->get_velocity(), vel_fun->get_accuracy()));
            start("smooth");
            break;
        case '3':
            stop();
            QUIT_ON_COMPLETION = true;
            RECORD = true;
            vel_fun = std::unique_ptr<VelocityFunc>(new NormalFunc(vel_fun->get_velocity(), vel_fun->get_accuracy()));
            start("expand");
            break;
        case ' ':
            if(!CONTINUOUS)
            {
                std::cout << "MOTION STARTED" << std::endl;
                if(RECORD && basic_log)
                {
                    painter->set_view_position(camera_pos);
                    painter->save_painting(basic_log->get_path(), vel_fun->get_time_step());
                }
            }
            else {
                std::cout << "MOTION PAUSED" << std::endl;
            }
            CONTINUOUS = !CONTINUOUS;
            break;
        case 'm':
            std::cout << "MOVE" << std::endl;
            vel_fun->take_time_step(*dsc);
            painter->update(*dsc);
            break;
        case 'r':
            std::cout << "RELOAD MODEL" << std::endl;
            load_model(model_file_name, dsc->get_avg_edge_length());
            break;
        case 't':
            std::cout << "TEST VELOCITY FUNCTION" << std::endl;
            vel_fun->test(*dsc);
            painter->update(*dsc);
            break;
        case '\t':
            painter->switch_display_type();
            painter->update(*dsc);
            break;
        case 's':
            std::cout << "TAKING SCREEN SHOT" << std::endl;
            painter->set_view_position(camera_pos);
            painter->save_painting("LOG");
            break;
        case 'e':
        {
            std::cout << "EXPORTING MESH" << std::endl;
            std::string filename("data/mesh.dsc");
            std::vector<vec3> points;
            std::vector<int> tets;
            std::vector<int> tet_labels;
            dsc->get_is_mesh().extract_tet_mesh(points, tets, tet_labels);
            is_mesh::export_tet_mesh(filename, points, tets, tet_labels);
        }
            break;
        case 'i':
        {
            std::cout << "EXPORTING SURFACE MESH" << std::endl;
            std::string filename("data/mesh.obj");
            std::vector<vec3> points;
            std::vector<int> faces;
            dsc->get_is_mesh().extract_surface_mesh(points, faces);
            is_mesh::export_surface_mesh(filename, points, faces);
        }
            break;
        case '+':
        {
            double velocity = std::min(vel_fun->get_velocity() + 1., 100.);
            vel_fun->set_velocity(velocity);
            update_title();
        }
            break;
        case '-':
        {
            double velocity = std::max(vel_fun->get_velocity() - 1., 0.);
            vel_fun->set_velocity(velocity);
            update_title();
        }
            break;
        case '.':
        {
            double discretization = std::min(dsc->get_avg_edge_length() + 0.5, 100.);
            dsc->set_avg_edge_length(discretization);
            update_title();
        }
            break;
        case ',':
        {
            double discretization = std::max(dsc->get_avg_edge_length() - 0.5, 1.);
            dsc->set_avg_edge_length(discretization);
            update_title();
        }
            break;
        case '<':
        {
            double accuracy = std::min(vel_fun->get_accuracy() + 1., 100.);
            vel_fun->set_accuracy(accuracy);
            update_title();
        }
            break;
        case '>':
        {
            double accuracy = std::max(vel_fun->get_accuracy() - 1., 1.);
            vel_fun->set_accuracy(accuracy);
            update_title();
        }
        break;
        case 'x':
        {
            if (dsc->get_is_mesh().get_subdomain()){
                dsc->clear_subdomain();
            } else {
                dsc->set_subdomain(std::make_shared<is_mesh::Sphere>(vec3{0.0,0.0,0.0},0.06));
            }
            painter->update(*dsc);
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
    if(RECORD && basic_log)
    {
        basic_log->write_message("MOTION STOPPED");
        basic_log->write_log(*dsc);
        basic_log->write_log(*vel_fun);
        basic_log->write_timings(*vel_fun);
        
        std::vector<vec3> points;
        std::vector<int> faces;
        std::vector<int> tets;
        std::vector<int> tet_labels;
        dsc->get_is_mesh().extract_tet_mesh(points, tets, tet_labels);
        is_mesh::export_tet_mesh(basic_log->get_path() + std::string("/mesh.dsc"), points, tets, tet_labels);
        points.clear();
        dsc->get_is_mesh().extract_surface_mesh(points, faces);
        is_mesh::export_surface_mesh(basic_log->get_path() + std::string("/mesh.obj"), points, faces);
        basic_log = nullptr;
    }
    
    CONTINUOUS = false;
    painter->update(*dsc);
    update_title();
    glutPostRedisplay();
}

void UI::start(const std::string& log_folder_name)
{
    if(RECORD)
    {
        basic_log = std::unique_ptr<Log>(new Log(log_path + log_folder_name));
        painter->set_view_position(camera_pos);
        painter->save_painting(log_path, vel_fun->get_time_step());
        basic_log->write_message(vel_fun->get_name().c_str());
        basic_log->write_log(*vel_fun);
        basic_log->write_log(*dsc);
    }
    
    painter->update(*dsc);
    update_title();
    glutPostRedisplay();
}
