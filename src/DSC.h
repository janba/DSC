//
//  DSC.h
//  3D_DSC
//
//  Created by Asger Nyman Christiansen on 2/14/13.
//  Copyright (c) 2013 Asger Nyman Christiansen. All rights reserved.
//

#ifndef _D_DSC_DSC_h
#define _D_DSC_DSC_h

#include "velocity_function.h"
#include "log.h"

/**
 The fundamental DSC class which contains the simplicial complex and the velocity function. It also handles logging and keeps track of the time step number.
 */
template<class MT>
class DSC {
    
    VelocityFunc<MT> *vel_fun;
    SimplicialComplex<MT> *complex;
    Log *basic_log;
    
    int time_step;
    int MAX_TIME_STEPS;
    
public:
    
    DSC(VelocityFunc<MT> *vel_fun_, SimplicialComplex<MT> *complex_, Log *log_, int max_time_steps = 500):
        vel_fun(vel_fun_), complex(complex_), basic_log(log_), time_step(0), MAX_TIME_STEPS(max_time_steps)
    {
        basic_log->write_message(vel_fun->get_name().c_str());
        basic_log->write_log(complex);
        basic_log->write_log(vel_fun);
    }
    
    ~DSC()
    {
        basic_log->write_message("MOTION STOPPED");
        basic_log->write_log(complex);
        basic_log->write_log(vel_fun);
        basic_log->write_timings(vel_fun);
        delete basic_log;
        delete complex;
        delete vel_fun;
    }
    
    /**
     Returns the simplicial complex. Primarily used to be able to draw the simplicial complex.
     */
    SimplicialComplex<MT>* get_complex()
    {
        return complex;
    }
    
    /**
     Returns the velocity function. Primarily used to be able to draw some aspects of the velocity function.
     */
    const VelocityFunc<MT>* get_velocity_func()
    {
        return vel_fun;
    }
    
    /**
     Returns the current time step.
     */
    int get_time_step()
    {
        return time_step;
    }
    
    /**
     Returns the path where the log is saved.
     */
    std::string get_log_path()
    {
        return basic_log->get_path();
    }
    
    /**
     Returns the title of the deformation.
     */
    std::string get_title()
    {
        std::ostringstream oss;
        oss << vel_fun->get_name();
        oss << ", Time step " << time_step;
        return oss.str();
    }
    
    /**
     Takes one time step thereby deforming the simplicial complex according to the velocity function.
     */
    bool take_time_step()
    {
        vel_fun->reset_times();
        vel_fun->deform(complex);
        basic_log->write_timestep(time_step, vel_fun, complex);
        
        time_step++;
        
        return (time_step == MAX_TIME_STEPS) | vel_fun->is_motion_finished(complex);
    }
    
    /**
     Calls the optional test function in the velocity function class.
     */
    void test()
    {
        vel_fun->test(complex);
    }

};

#endif
