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

#include "DSC.h"
#include "velocity_function.h"
#include "log.h"
#include "draw.h"

/**
 A default application which utilizes OpenGL, GLEW and GLUT for visualization. Three sample velocity functions (rotation, smoothing and expansion) can be applies to a model specified by the model_file_name variable or as input variable. See https://github.com/asny/DSC/wiki/DEMO-instructions for details on how to use this DEMO application. See https://github.com/asny/DSC/wiki/Instructions for instructions on how to build your own application which uses the implementation of the DSC method.
 */
class UI
{
    std::unique_ptr<DSC::VelocityFunc<>> vel_fun;
    std::unique_ptr<DSC::DeformableSimplicialComplex<>> dsc;
    std::unique_ptr<Log> basic_log;
    std::unique_ptr<Painter> painter;
    
    std::string model_file_name = "armadillo";
    
    vec3 eye_pos = {70., 30., 70.};
    vec3 camera_pos = {30., 30., 70.};
    vec3 light_pos = {0., 0., 70.};
    
    vec3 center = vec3(0.);
    vec3 size = vec3(40.);
    
    int WIN_SIZE_X = 700.;
    int WIN_SIZE_Y = 700;
    
    bool CONTINUOUS = false;
    bool RECORD = false;
    bool QUIT_ON_COMPLETION = false;
    
    static UI* instance;
    
#ifdef _WIN32
    const std::string obj_path = "data\\";
    const std::string log_path = "LOG\\";
#else
    const std::string obj_path = "./data/";
    const std::string log_path = "./LOG/";
#endif
    const std::string extension = ".dsc";
    
public:
    
    UI(int &argc, char** argv);
    
    static UI* get_instance()
    {
        return instance;
    }
    
    void display();
    
    void animate();
    
    void reshape(int width, int height);
    
    void visible(int v);
    
    void mouse(int button, int state, int x, int y);
    
    void motion(int x, int y);
    
    /**
     The keyboard is used for all inputs. See https://github.com/asny/DSC/wiki/DEMO-instructions for up-to-date instructions on how to use the DEMO application.
     */
    void keyboard(unsigned char key, int x, int y);
    
private:
    
    void set_size(vec3 s)
    {
        size = s;
        real var = Util::max(Util::max(s[0], s[1]), s[2]);
        real dist = 1.2*var;
        eye_pos = {dist, var, dist};
        camera_pos = {var, var, -dist};
        light_pos = {0., 0., dist};
    }
    
    /**
     Loads the .dsc file specified by the model_file_name variable.
     */
    void load_model(const std::string& file_name, real discretization);
    
    /**
     Updates the window title.
     */
    void update_title();
    
    /**
     Starts the motion.
     */
    void start(const std::string& log_folder_name);
    
    /**
     Stops the motion and deletes the DSC object.
     */
    void stop();
};
