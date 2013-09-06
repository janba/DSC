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

//#include "tetrahedralize.h"
#include "mesh_io.h"
#include "tetralizer.h"


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
	WIN_SIZE_X = 1000;
    WIN_SIZE_Y = 1000;

    glutInit(&argc, argv);
    glutInitWindowSize(WIN_SIZE_X,WIN_SIZE_Y);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH | GLUT_MULTISAMPLE);
    glutCreateWindow("");
    
    glutDisplayFunc(display_);
    glutKeyboardFunc(keyboard_);
	glutIgnoreKeyRepeat(true);
    glutVisibilityFunc(visible_);
    glutReshapeFunc(reshape_);
	glutIdleFunc(animate_);
    
	glewExperimental = GL_TRUE;  // Added because of http://openglbook.com/glgenvertexarrays-access-violationsegfault-with-glew/
    
	GLint GlewInitResult = glewInit();
	if (GlewInitResult != GLEW_OK) {
		printf("ERROR: %s\n", glewGetErrorString(GlewInitResult));
	}
        
    painter = new Painter<GELTypes>();
    vel_fun = nullptr;
    dsc = nullptr;
    
    if(argc > 1)
    {
        QUIT_ON_COMPLETION = true;
        CONTINUOUS = true;
        RECORD = true;
        
        int motion;
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
            else if (str == "motion") {
                motion = std::atoi(argv[i+1]);
            }
        }
        
        switch (motion) {
            case 1:
                motion1();
                break;
            case 2:
                motion2();
                break;
            case 3:
                motion3();
                break;
            default:
                break;
        }
    }
    else {
        VELOCITY = 5.;
        DISCRETIZATION = 2.5;
        ACCURACY = 1.;
        
        CONTINUOUS = false;
        RECORD = true;
        QUIT_ON_COMPLETION = true;
    }
    update_title();
    check_gl_error();
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
    draw();
    update_title();
    check_gl_error();
}

void UI::reshape(int width, int height)
{
    WIN_SIZE_X = width;
    WIN_SIZE_Y = height;
	glViewport(0, 0, WIN_SIZE_X, WIN_SIZE_Y);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(53, WIN_SIZE_X/float(WIN_SIZE_Y), 0.01*r, 3.*r);
    glMatrixMode(GL_MODELVIEW);
}

void UI::animate()
{
    if(vel_fun && CONTINUOUS)
    {
        vel_fun->take_time_step(*dsc);
        basic_log->write_timestep(vel_fun, dsc);
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
            motion1();
            break;
        case '2':
            motion2();
            break;
        case '3':
            motion3();
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
            if(vel_fun)
            {
                std::cout << "MOVE" << std::endl;
                vel_fun->take_time_step(*dsc);
            }
            break;
        case 't':
            if(vel_fun)
            {
                std::cout << "TEST" << std::endl;
                vel_fun->test(*dsc);
            }
            break;
        case '\t':
            if(dsc)
            {
                switch_display_type();
            }
            break;
        case 's':
            if(dsc)
            {
                std::cout << "TAKING SCREEN SHOT" << std::endl;
                Painter<GELTypes>::save_painting(WIN_SIZE_X, WIN_SIZE_Y, "LOG");
            }
            break;
        case 'e':
            if(dsc)
            {
                std::cout << "EXPORTING MESH" << std::endl;
                std::string filepath("data/mesh.dsc");
                export_tet_mesh<GELTypes>(*dsc, filepath);
            }
            break;
        case '+':
            if(!vel_fun)
            {
                VELOCITY = std::min(VELOCITY + 1., 100.);
                update_title();
            }
            break;
        case '-':
            if(!vel_fun)
            {
                VELOCITY = std::max(VELOCITY - 1., 0.);
                update_title();
            }
            break;
        case '.':
            if(!vel_fun)
            {
                DISCRETIZATION = std::min(DISCRETIZATION + 0.5, 100.);
                update_title();
            }
            break;
        case ',':
            if(!vel_fun)
            {
                DISCRETIZATION = std::max(DISCRETIZATION - 0.5, 1.);
                update_title();
            }
            break;
        case '<':
            if(!vel_fun)
            {
                ACCURACY = std::min(ACCURACY + 1., 100.);
                update_title();
            }
            break;
        case '>':
            if(!vel_fun)
            {
                ACCURACY = std::max(ACCURACY - 1., 1.);
                update_title();
            }
            break;
    }
    draw();
}

void UI::visible(int v)
{
    if(v==GLUT_VISIBLE)
        glutIdleFunc(animate_);
    else
        glutIdleFunc(0);
}

void UI::draw()
{
    Painter<GELTypes>::begin();
    if (dsc)
    {
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        auto center = dsc->get_center();
        gluLookAt(0., 0., -r, center[0], center[1], center[2], 0., 1., 0.);
        
        Painter<GELTypes>::draw_complex(dsc);
        if(vel_fun && RECORD && CONTINUOUS)
        {
            Painter<GELTypes>::save_painting(WIN_SIZE_X, WIN_SIZE_Y, basic_log->get_path(), vel_fun->get_time_step());
        }
    }
    Painter<GELTypes>::end();
}

void UI::stop()
{
    if(vel_fun)
    {
        draw();
        basic_log->write_message("MOTION STOPPED");
        basic_log->write_log(dsc);
        basic_log->write_log(vel_fun);
        basic_log->write_timings(vel_fun);

        delete vel_fun;
        delete dsc;
        delete basic_log;
        vel_fun = nullptr;
        dsc = nullptr;
    }
}

void UI::motion1()
{
    stop();
    // Build the Simplicial Complex
    std::vector<double> points;
    std::vector<int>  tets;
    std::vector<int>  tet_labels;
    
//    import_tet_mesh<GELTypes>(get_data_file_path("armadillo.dsc").data(), points, tets, tet_labels);
    Tetralizer<GELTypes> tetralizer(50., 50., 50., DISCRETIZATION);
    tetralizer.tetralize(points, tets, tet_labels);
    
    dsc = new DeformableSimplicialComplex<GELTypes>(DISCRETIZATION, points, tets, tet_labels);
    vel_fun = new RotateFunc<GELTypes>(VELOCITY, ACCURACY);
    
    basic_log = new Log(create_log_path());
    
    basic_log->write_message(vel_fun->get_name().c_str());
    basic_log->write_log(dsc);
    basic_log->write_log(vel_fun);
    
    update_title();
    reshape(WIN_SIZE_X, WIN_SIZE_Y);
}

void UI::motion2()
{
    stop();
    // Build the Simplicial Complex
    std::vector<double> points;
    std::vector<int>  tets;
    std::vector<int>  tet_labels;
    import_tet_mesh<GELTypes>(get_data_file_path("armadillo.dsc").data(), points, tets, tet_labels);
    
    dsc = new DeformableSimplicialComplex<GELTypes>(DISCRETIZATION, points, tets, tet_labels);
    vel_fun = new AverageFunc<GELTypes>(VELOCITY, ACCURACY);
    
    basic_log = new Log(create_log_path());
    
    basic_log->write_message(vel_fun->get_name().c_str());
    basic_log->write_log(dsc);
    basic_log->write_log(vel_fun);
    
    update_title();
    reshape(WIN_SIZE_X, WIN_SIZE_Y);
}

void UI::motion3()
{
    stop();
    // Build the Simplicial Complex
    std::vector<double> points;
    std::vector<int>  tets;
    std::vector<int>  tet_labels;
    import_tet_mesh<GELTypes>(get_data_file_path("armadillo.dsc").data(), points, tets, tet_labels);
    
//    std::vector<GELTypes::vector3_type> pts_inside = {GELTypes::vector3_type( 0.0f, 0.0f, 0.0f)};
//    build_tetrahedralization<GELTypes>(get_data_file_path("armadillo-very-simple.obj"), points, tets, tet_labels, pts_inside);
    
    dsc = new DeformableSimplicialComplex<GELTypes>(DISCRETIZATION, points, tets, tet_labels);
    vel_fun = new NormalFunc<GELTypes>(VELOCITY, ACCURACY);
    
    basic_log = new Log(create_log_path());
    
    basic_log->write_message(vel_fun->get_name().c_str());
    basic_log->write_log(dsc);
    basic_log->write_log(vel_fun);
    
    update_title();
    reshape(WIN_SIZE_X, WIN_SIZE_Y);
}

