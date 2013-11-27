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

#include "log.h"

#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#include <sys/types.h>
#endif

using namespace DSC;

Log::Log(const std::string& path_)
{
    std::string temp;
    int error;
    int i = 0;
    do {
        temp = Util::concat4digits(path_ + "_test",i);
#ifdef _WIN32
        error = _mkdir(temp.c_str());
#else
        error = mkdir(temp.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
#endif
        i++;
    } while (error != 0 && i < 1000);
    path = temp;
    log.open(path + "/log.txt");
}

void Log::write_variable(const std::string& name, real value)
{
    log << "\t" << name << "\t:\t" << value << std::endl;
}

void Log::write_variable(const std::string& name, real value, real change)
{
    log << "\t" << name << "\t:\t" << value << "\t\tChange\t:\t" << change << std::endl;
}

void Log::write_variable(const std::string& name, real value, const std::string& unit)
{
    log << "\t" << name << "\t:\t" << value << " " << unit << std::endl;
}

void Log::write_variable(const std::string& name, const std::vector<real>& values)
{
    if(values.size() > 0)
    {
        log << "\t" << name << " = [";
        for (auto val = values.begin(); val != values.end(); val++)
        {
            log << *val;
            if(val + 1 != values.end())
            {
                log << ",\t";
            }
        }
        log << "];" << std::endl << std::endl;
    }
}

void Log::write_variable(const std::string& name, const std::vector<int>& values)
{
    if(values.size() > 0)
    {
        log << "\t" << name << " = [";
        for (auto val = values.begin(); val != values.end(); val++)
        {
            log << *val;
            if(val + 1 != values.end())
            {
                log << ",\t";
            }
        }
        log << "];" << std::endl << std::endl;
    }
}

void Log::write_message(const std::string& message)
{
    log << std::endl  << "*** " << message << " ***" << std::endl;
    std::cout << "*** " << message << " ***" << std::endl;
}

void Log::write_timestep(const VelocityFunc<>& vel_fun, DeformableSimplicialComplex<>& dsc)
{
    log << std::endl << "*** Time step #" << vel_fun.get_time_step() << " ***" << std::endl;
    log << std::endl;
    write_variable("Compute time", vel_fun.get_compute_time(), "s");
    write_variable("Deform time", vel_fun.get_deform_time(), "s");
    write_variable("Total time", vel_fun.get_compute_time() + vel_fun.get_deform_time(), "s");
    
    write_variable("Min quality", dsc.min_quality());
    
    std::vector<int> hist;
    real min_a, max_a;
    dsc.get_dihedral_angles(hist, min_a, max_a);
    write_variable("Min dih. angle", min_a, "degrees");
    write_variable("Max dih. angle", max_a, "degrees");
}

void Log::write_log(DeformableSimplicialComplex<>& dsc)
{
    write_message("SIMPLICIAL COMPLEX INFO");
    write_variable("Discretization", dsc.get_avg_edge_length());
    
    int total, object;
    dsc.count_nodes(total, object);
    write_variable("#nodes\t", total);
    write_variable("#obj nodes", object);
    
    dsc.count_edges(total, object);
    write_variable("#edges\t", total);
    write_variable("#obj edges", object);
    
    dsc.count_faces(total, object);
    write_variable("#faces\t", total);
    write_variable("#obj faces", object);
    
    dsc.count_tetrahedra(total, object);
    write_variable("#tetrahedra", total);
    write_variable("#obj tets", object);
    
    std::vector<int> hist;
    real min_a, max_a;
    dsc.get_dihedral_angles(hist, min_a, max_a);
    write_variable("Min dih. angle", min_a, "degrees");
    write_variable("Max dih. angle", max_a, "degrees");
    write_variable("DAhist", hist);
    
    dsc.get_qualities(hist, min_a);
    write_variable("Min quality", min_a);
    write_variable("Qhist", hist);
    
}

void Log::write_log(const VelocityFunc<>& vel_fun)
{
    write_message("VELOCITY FUNCTION INFO");
    write_variable("Velocity", vel_fun.get_velocity());
    write_variable("Accuracy", vel_fun.get_accuracy());
}

void Log::write_timings(const VelocityFunc<>& vel_fun)
{
    real deform_time = vel_fun.get_total_deform_time();
    real compute_time = vel_fun.get_total_compute_time();
    write_message("TIMINGS");
    write_variable("Total time", deform_time + compute_time, "s");
    write_variable("Compute time", compute_time, "s");
    write_variable("Deform time", deform_time , "s");
}
