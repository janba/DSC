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
#include "tetrahedralize.h"
#include "mesh_io.h"
#include "tetralizer.h"

void _check_gl_error(const char *file, int line)
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

void motion_(int x, int y){
    UI::get_instance()->motion(x,y);
}

void mouse_(int button, int state, int x, int y)
{
    UI::get_instance()->mouse(button, state, x, y);
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
	glutMouseFunc(mouse_);
	glutMotionFunc(motion_);
	glutIdleFunc(animate_);
    
    glEnable(GL_MULTISAMPLE);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glEnable(GL_POINT_SMOOTH);
    glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
    
	glEnable(GL_DEPTH_TEST);
    glDepthMask(GL_TRUE);
    
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    
	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
    
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glShadeModel(GL_FLAT);
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);
    
	glEnable (GL_BLEND);
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
        
    vel_fun = nullptr;
    dsc = nullptr;
    
    if(argc > 1)
    {
        QUIT_ON_COMPLETION = true;
        CONTINUOUS = true;
        RECORD = true;
        
        Util::ArgExtracter ext(argc, argv);
        ext.extract("nu", VELOCITY);
        ext.extract("delta", DISCRETIZATION);
        ext.extract("alpha", ACCURACY);
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
    
    if (vel_fun) {
        view_ctrl->reshape(WIN_SIZE_X,WIN_SIZE_Y);
    }
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

void UI::mouse(int button, int state, int x, int y)
{
    if (dsc) {
        CGLA::Vec2i pos(x,y);
        if (state == GLUT_DOWN)
        {
            if (button == GLUT_LEFT_BUTTON)
                view_ctrl->grab_ball(GLGraphics::ROTATE_ACTION,pos);
            else if (button == GLUT_MIDDLE_BUTTON)
                view_ctrl->grab_ball(GLGraphics::ZOOM_ACTION,pos);
            else if (button == GLUT_RIGHT_BUTTON)
                view_ctrl->grab_ball(GLGraphics::PAN_ACTION,pos);
        }
        else if (state == GLUT_UP)
            view_ctrl->release_ball();
    }
}

void UI::motion(int x, int y)
{
    if (dsc) {
        CGLA::Vec2i pos(x,y);
        view_ctrl->roll_ball(pos);
    }
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
        view_ctrl->set_gl_modelview();
        
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
        delete view_ctrl;
        vel_fun = nullptr;
        dsc = nullptr;
    }
}

void UI::motion1()
{
    typedef typename GELTypes::vector3_type V;
    stop();
    // Build the Simplicial Complex
    std::vector<double> points;
    std::vector<int>  tets;
    std::vector<int>  tet_labels;
    std::vector<V> pts_inside(1);
    pts_inside[0] = V( 0.0f, 0.0f, 0.0f);
    
//    import_tet_mesh<GELTypes>(get_data_file_path("armadillo.dsc").data(), points, tets, tet_labels);
//    build_tetrahedralization<GELTypes>(get_data_file_path("armadillo-very-simple.obj"), points, tets, tet_labels, pts_inside);
    Tetralizer<GELTypes> tetralizer(50., 50., 50., DISCRETIZATION);
    tetralizer.tetralize(points, tets, tet_labels);
    
    dsc = new DeformableSimplicialComplex<GELTypes>(DISCRETIZATION, points, tets, tet_labels);
    vel_fun = new RotateFunc<GELTypes>(VELOCITY, ACCURACY);
    
    basic_log = new Log(create_log_path());
    
    basic_log->write_message(vel_fun->get_name().c_str());
    basic_log->write_log(dsc);
    basic_log->write_log(vel_fun);
    
    view_ctrl = new GLGraphics::GLViewController(WIN_SIZE_X, WIN_SIZE_Y, (CGLA::Vec3f)dsc->get_center(), r);
    view_ctrl->set_view_param(CGLA::Vec3f(0.,0., -r), (CGLA::Vec3f)dsc->get_center(), CGLA::Vec3f(0.,1.,0.));
    
    update_title();
    reshape(WIN_SIZE_X, WIN_SIZE_Y);
}

void UI::motion2()
{
    typedef typename GELTypes::vector3_type V;
    stop();
    // Build the Simplicial Complex
    vector<double> points;
    vector<int>  tets;
    vector<int>  tet_labels;
    vector<V> pts_inside(1);
    pts_inside[0] = V( 0.0f, 0.0f, 0.0f);
    
    import_tet_mesh<GELTypes>(get_data_file_path("armadillo.dsc").data(), points, tets, tet_labels);
//    build_tetrahedralization<GELTypes>(get_data_file_path("armadillo-very-simple.obj"), points, tets, tet_labels, pts_inside);
    
    dsc = new DeformableSimplicialComplex<GELTypes>(DISCRETIZATION, points, tets, tet_labels);
    vel_fun = new AverageFunc<GELTypes>(VELOCITY, ACCURACY);
    
    basic_log = new Log(create_log_path());
    
    basic_log->write_message(vel_fun->get_name().c_str());
    basic_log->write_log(dsc);
    basic_log->write_log(vel_fun);
    
    view_ctrl = new GLGraphics::GLViewController(WIN_SIZE_X, WIN_SIZE_Y, (CGLA::Vec3f)dsc->get_center(), r);
    view_ctrl->set_view_param(CGLA::Vec3f(0.,0., -r), (CGLA::Vec3f)dsc->get_center(), CGLA::Vec3f(0.,1.,0.));
    
    update_title();
    reshape(WIN_SIZE_X, WIN_SIZE_Y);
}

void UI::motion3()
{
    typedef typename GELTypes::vector3_type V;
    stop();
    // Build the Simplicial Complex
    vector<double> points;
    vector<int>  tets;
    vector<int>  tet_labels;
    vector<V> pts_inside(1);
    pts_inside[0] = V( 0.0f, 0.0f, 0.0f);
    
    import_tet_mesh<GELTypes>(get_data_file_path("armadillo.dsc").data(), points, tets, tet_labels);
//    build_tetrahedralization<GELTypes>(get_data_file_path("armadillo-very-simple.obj"), points, tets, tet_labels, pts_inside);
    
    dsc = new DeformableSimplicialComplex<GELTypes>(DISCRETIZATION, points, tets, tet_labels);
    vel_fun = new NormalFunc<GELTypes>(VELOCITY, ACCURACY);
    
    basic_log = new Log(create_log_path());
    
    basic_log->write_message(vel_fun->get_name().c_str());
    basic_log->write_log(dsc);
    basic_log->write_log(vel_fun);
    
    view_ctrl = new GLGraphics::GLViewController(WIN_SIZE_X, WIN_SIZE_Y, (CGLA::Vec3f)dsc->get_center(), r);
    view_ctrl->set_view_param(CGLA::Vec3f(0.,0., -r), (CGLA::Vec3f)dsc->get_center(), CGLA::Vec3f(0.,1.,0.));
    
    update_title();
    reshape(WIN_SIZE_X, WIN_SIZE_Y);
}

