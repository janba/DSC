//
//  log.h
//  3D_DSC
//
//  Created by Asger Nyman Christiansen on 9/10/12.
//  Copyright (c) 2012 DTU Informatics. All rights reserved.
//

#ifndef LOG_h
#define LOG_h

#include "util.h"
#include "simplicial_complex.h"
#include "velocity_function.h"

#include <fstream>
#include <string>

/**
 A class for logging information.
 */

class Log {
    
    std::ofstream log;
    std::string path;
    
public:
    /**
     Creates a log at the path given as input.
     */
    Log(std::string path);
    
    ~Log()
    {        
        log.close();
    }
    
    /**
     Returns the path where the log is saved.
     */
    std::string get_path()
    {
        return path;
    }
    
protected:
    /**
     Write a variable with name and value to the log.
     */
    void write_variable(const char* name, double value);
    
    /**
     Write a variable with name, value and change in value to the log.
     */
    void write_variable(const char* name, double value, double change);
    
    /**
     Write a variable with name, value and unit of the value to the log.
     */
    void write_variable(const char* name, double value, const char* unit);
    
    /**
     Write a variable with name and values to the log.
     */
    void write_variable(const char* name, const std::vector<double>& values);
    
    void write_variable(const char* name, const std::vector<int>& values);
    
public:
    
    /**
     Write a message to the terminal and the log.
     */
    virtual void write_message(const char* message);
    
    /**
     Write the time step number, timings and additional time step information to the log.
     */
    template<class MT>
    void write_timestep(int time_step, const VelocityFunc<MT> *vel_fun, SimplicialComplex<MT> *complex)
    {
        //    std::cout << "\n\n*** Time step #" << vel_fun->get_time_step() << " ***" << std::endl;
        log << std::endl << "*** Time step #" << time_step << " ***" << std::endl;
        log << std::endl;
        write_variable("Compute time", vel_fun->get_compute_time(), "s");
        write_variable("Deform time", vel_fun->get_deform_time(), "s");
        write_variable("Total time", vel_fun->get_compute_time() + vel_fun->get_deform_time(), "s");
        
        
        write_variable("Min quality", complex->min_quality());
        
        std::vector<int> hist;
        double min_a, max_a;
        complex->calc_dihedral_angles(hist, min_a, max_a);
        write_variable("Min dih. angle", min_a, "degrees");
        write_variable("Max dih. angle", max_a, "degrees");
    }
    
    /**
     Writes simplicial complex information to the log.
     */
    template<class MT>
    void write_log(SimplicialComplex<MT> *complex)
    {
        write_message("SIMPLICIAL COMPLEX INFO");
//        write_variable("Size X\t", complex->get_size_x());
//        write_variable("Size Y\t", complex->get_size_y());
//        write_variable("Avg edge length", complex->get_avg_edge_length());
//        write_variable("Min deformation", complex->get_min_deformation());
//        write_variable("Total volume", complex->get_volume());
        
        int total, object;
        complex->count_nodes(total, object);
        write_variable("#nodes\t", total);
        write_variable("#obj nodes", object);
        
        complex->count_edges(total, object);
        write_variable("#edges\t", total);
        write_variable("#obj edges", object);
        
        complex->count_faces(total, object);
        write_variable("#faces\t", total);
        write_variable("#obj faces", object);
        
        complex->count_tetrahedra(total, object);
        write_variable("#tetrahedra", total);
        write_variable("#obj tets", object);
        
        write_variable("Min quality", complex->min_quality());
        
        std::vector<int> hist;
        double min_a, max_a;
        complex->calc_dihedral_angles(hist, min_a, max_a);
        write_variable("Min dih. angle", min_a, "degrees");
        write_variable("Max dih. angle", max_a, "degrees");
        write_variable("hist", hist);
    }
    
    /**
     Writes velocity function information to the log.
     */
    template<class MT>
    void write_log(const VelocityFunc<MT> *vel_fun)
    {
        write_message("VELOCITY FUNCTION INFO");
        write_variable("Velocity", vel_fun->get_velocity());
        write_variable("Accuracy", vel_fun->get_accuracy());
        
    }
    
    /**
     Writes timings to the log.
     */
    template<class MT>
    void write_timings(const VelocityFunc<MT> *vel_fun)
    {
        double deform_time = vel_fun->get_total_deform_time();
        double compute_time = vel_fun->get_total_compute_time();
        write_message("TIMINGS");
        write_variable("Total time", deform_time + compute_time, "s");
        write_variable("Compute time", compute_time, "s");
        write_variable("Deform time", deform_time , "s");
    }
};



#endif
