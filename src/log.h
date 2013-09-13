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

#include "util.h"
#include "DSC.h"
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
    void write_variable(const char* name, real value);
    
    /**
     Write a variable with name, value and change in value to the log.
     */
    void write_variable(const char* name, real value, real change);
    
    /**
     Write a variable with name, value and unit of the value to the log.
     */
    void write_variable(const char* name, real value, const char* unit);
    
    /**
     Write a variable with name and values to the log.
     */
    void write_variable(const char* name, const std::vector<real>& values);
    
    void write_variable(const char* name, const std::vector<int>& values);
    
public:
    
    /**
     Write a message to the terminal and the log.
     */
    virtual void write_message(const char* message);
    
    /**
     Write the time step number, timings and additional time step information to the log.
     */
    void write_timestep(const VelocityFunc *vel_fun, DeformableSimplicialComplex<> *complex)
    {
        //    std::cout << "\n\n*** Time step #" << vel_fun->get_time_step() << " ***" << std::endl;
        log << std::endl << "*** Time step #" << vel_fun->get_time_step() << " ***" << std::endl;
        log << std::endl;
        write_variable("Compute time", vel_fun->get_compute_time(), "s");
        write_variable("Deform time", vel_fun->get_deform_time(), "s");
        write_variable("Total time", vel_fun->get_compute_time() + vel_fun->get_deform_time(), "s");
        
        write_variable("Min quality", complex->min_quality());
        
        std::vector<int> hist;
        real min_a, max_a;
        complex->get_dihedral_angles(hist, min_a, max_a);
        write_variable("Min dih. angle", min_a, "degrees");
        write_variable("Max dih. angle", max_a, "degrees");
    }
    
    /**
     Writes simplicial complex information to the log.
     */
    void write_log(DeformableSimplicialComplex<> *complex)
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
        
        std::vector<int> hist;
        real min_a, max_a;
        complex->get_dihedral_angles(hist, min_a, max_a);
        write_variable("Min dih. angle", min_a, "degrees");
        write_variable("Max dih. angle", max_a, "degrees");
        write_variable("DAhist", hist);
        
        complex->get_qualities(hist, min_a);
        write_variable("Min quality", min_a);
        write_variable("Qhist", hist);
        
    }
    
    /**
     Writes velocity function information to the log.
     */
    void write_log(const VelocityFunc *vel_fun)
    {
        write_message("VELOCITY FUNCTION INFO");
        write_variable("Velocity", vel_fun->get_velocity());
        write_variable("Accuracy", vel_fun->get_accuracy());
    }
    
    /**
     Writes timings to the log.
     */
    void write_timings(const VelocityFunc *vel_fun)
    {
        real deform_time = vel_fun->get_total_deform_time();
        real compute_time = vel_fun->get_total_compute_time();
        write_message("TIMINGS");
        write_variable("Total time", deform_time + compute_time, "s");
        write_variable("Compute time", compute_time, "s");
        write_variable("Deform time", deform_time , "s");
    }
};
