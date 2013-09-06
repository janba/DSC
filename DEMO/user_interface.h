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

#ifndef ___D_DSC__DSC_UI__
#define ___D_DSC__DSC_UI__

#include <GL/glew.h>

#ifdef WIN32
#include <GL/glut.h>
#else
#include <GLUT/glut.h>
#endif

#include "GEL_types.h"
//#include "OTB_types.h"
#include "GLM_types.h"

#include "DSC.h"
#include "velocity_function.h"
#include "log.h"
#include "draw.h"

/**
 A default user interface which utilize OpenGL, GLEW and GLUT. At least some of the motion functions should be overridden.
 */
class UI
{
protected:
    float r = 150.;
    
    VelocityFunc<GELTypes> *vel_fun;
    DeformableSimplicialComplex<GELTypes> *dsc;
    
    Log *basic_log;
    
    int WIN_SIZE_X;
    int WIN_SIZE_Y;
    
    bool CONTINUOUS;
    bool RECORD;
    bool QUIT_ON_COMPLETION;
    
    double VELOCITY;
    double DISCRETIZATION;
    double ACCURACY;
    static UI* instance;
    
#ifdef WIN32
    const std::string obj_path = "@PROJECT_SOURCE_DIR@/data/";
#else
    const std::string obj_path = "./data/";
#endif
    
public:
    
    UI(int &argc, char** argv);
    
    static UI* get_instance()
    {
        return instance;
    }
    
    virtual std::string create_log_path()
    {
        std::ostringstream s;
        s << "LOG/delta" << DISCRETIZATION << "_nu" << VELOCITY << "_alpha" << ACCURACY;
        return s.str();
    }
    
    std::string get_data_file_path(std::string const & file)
    {
        std::cout << obj_path + file << std::endl;
        return obj_path + file;
    }
    
    virtual void display();
    
    virtual void animate();
    
    virtual void reshape(int width, int height);
    
    virtual void visible(int v);
    
    void mouse(int button, int state, int x, int y);
    
    void motion(int x, int y);
    
    /**
     The keyboard is used for all inputs.
     The workflow is to select parameters (discretization, velocity and accuracy), then the type of motion (the type of velocity function) and finally start the motion.
     A complete list of options are:
     
     *** SELECT PARAMETERS ***
     ,:         Decreases discretization by 0.5 to a minimum of 1.
     .:         Increases discretization by 0.5 to a maximum of 100.
     -:         Decreases velocity by 1 to a minimum of 1.
     +:         Increases velocity by 1 to a maximum of 100.
     >:         Decreases accuracy by 1 to a minimum of 1.
     <:         Increases accuracy by 1 to a maximum of 100.
     
     *** START/STOP MOTION ***
     SPACE:     Starts/pauses the current motion.
     0:         Stops the current motion.
     ESCAPE:    Stops the current motion and exits the application
     m:         Moves the interface vertices one time step according to the current velocity function.
     
     *** MISCELLANEOUS ***
     t:         Performs a test on the current velocity function.
     s:         Takes a screen shot.
     e:         Export the simplicial complex to a .dsc file.
     TAB:       Switches the display type.
     
     *** SELECT MOTION ***
     1:         Selects motion type 1.
     2:         Selects motion type 2.
     3:         Selects motion type 3.
     4:         Selects motion type 4.
     5:         Selects motion type 5.
     6:         Selects motion type 6.
     */
    virtual void keyboard(unsigned char key, int x, int y);
    
    void motion1();
    
    void motion2();
    
    void motion3();
    
    /**
     Updates the window title.
     */
    void update_title();
    
    /**
     Switches between different types of display if implemented.
     */
    void switch_display_type()
    {
        
    }
    
    /**
     Draws the simplicial complex.
     */
    virtual void draw();
    
    /**
     Stops the motion and deletes the DSC object.
     */
    virtual void stop();
};
#endif /* defined(___D_DSC__DSC_UI__) */
