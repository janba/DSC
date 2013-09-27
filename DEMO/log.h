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

#ifdef WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#include <sys/types.h>
#endif

/**
 A class for logging information.
 */
class Log {
    
    std::string path;
    std::ofstream log;
    
public:
    /**
     Creates a log at the path given as input.
     */
    Log(std::string path_)
    {
        std::string temp;
        int error;
        int i = 0;
        do {
            temp = DSC::Util::concat4digits(path_ + "_test",i);
#ifdef WIN32
            error = _mkdir(temp.c_str());
#else
            error = mkdir(temp.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
#endif
            i++;
        } while (error != 0 && i < 1000);
        path = temp;
        log.open(path + "/log.txt");
    }
    
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
    
private:
    /**
     Write a variable with name and value to the log.
     */
    void write_variable(const char* name, DSC::real value);
    
    /**
     Write a variable with name, value and change in value to the log.
     */
    void write_variable(const char* name, DSC::real value, DSC::real change);
    
    /**
     Write a variable with name, value and unit of the value to the log.
     */
    void write_variable(const char* name, DSC::real value, const char* unit);
    
    /**
     Write a variable with name and values to the log.
     */
    void write_variable(const char* name, const std::vector<DSC::real>& values);
    
    void write_variable(const char* name, const std::vector<int>& values);
    
public:
    
    /**
     Write a message to the terminal and the log.
     */
    void write_message(const char* message);
    
    /**
     Write the time step number, timings and additional time step information to the log.
     */
    void write_timestep(const DSC::VelocityFunc<> *vel_fun, DSC::DeformableSimplicialComplex<> *complex);
    
    /**
     Writes simplicial complex information to the log.
     */
    void write_log(DSC::DeformableSimplicialComplex<> *complex);
    
    /**
     Writes velocity function information to the log.
     */
    void write_log(const DSC::VelocityFunc<> *vel_fun);
    
    /**
     Writes timings to the log.
     */
    void write_timings(const DSC::VelocityFunc<> *vel_fun);
};
