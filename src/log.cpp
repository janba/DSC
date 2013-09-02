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


#ifdef WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#include <sys/types.h>
#endif

Log::Log(std::string path_)
{
    std::string temp;
    int error;
    int i = 0;
    do {
        temp = Util::concat4digits(path_ + "_test",i);
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

void Log::write_message(const char* message)
{
    log << std::endl  << "*** " << message << " ***" << std::endl;
    std::cout << "*** " << message << " ***" << std::endl;
}

//void Log::write_timestep(int time_step, const VelocityFunc *vel_fun)
//{
////    std::cout << "\n\n*** Time step #" << vel_fun->get_time_step() << " ***" << std::endl;
//    log << std::endl << "*** Time step #" << time_step << " ***" << std::endl;
//    log << std::endl;
//    write_variable("Compute time", vel_fun->get_compute_time(), "s");
//    write_variable("Deform time", vel_fun->get_deform_time(), "s");
//    write_variable("Total time", vel_fun->get_compute_time() + vel_fun->get_deform_time(), "s");
//}

void Log::write_variable(const char* name, double value)
{
    log << "\t" << name << "\t:\t" << value << std::endl;
}

void Log::write_variable(const char* name, double value, double change)
{
    log << "\t" << name << "\t:\t" << value << "\t\tChange\t:\t" << change << std::endl;
}

void Log::write_variable(const char* name, double value, const char* unit)
{
    log << "\t" << name << "\t:\t" << value << " " << unit << std::endl;
}

void Log::write_variable(const char* name, const std::vector<double>& values)
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

void Log::write_variable(const char* name, const std::vector<int>& values)
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

//void Log::write_log(const VelocityFunc *vel_fun)
//{
//    write_message("VELOCITY FUNCTION INFO");
//    write_variable("Velocity", vel_fun->get_velocity());
//    write_variable("Accuracy", vel_fun->get_accuracy());
//    
//}
//
//void Log::write_timings(const VelocityFunc *vel_fun)
//{
//    double deform_time = vel_fun->get_total_deform_time();
//    double compute_time = vel_fun->get_total_compute_time();
//    write_message("TIMINGS");
//    write_variable("Total time", deform_time + compute_time, "s");
//    write_variable("Compute time", compute_time, "s");
//    write_variable("Deform time", deform_time , "s");
//}

//void Log::write_log(const SimplicialComplex *complex)
//{
//    write_message("SIMPLICIAL COMPLEX INFO");
//    write_variable("Size X\t", complex->get_size_x());
//    write_variable("Size Y\t", complex->get_size_y());
//    write_variable("Avg edge length", complex->get_avg_edge_length());
//    write_variable("Min deformation", complex->get_min_deformation());
//    write_variable("Total volume", complex->get_volume());
//    write_variable("#vertices", complex->get_no_vertices());
//    write_variable("#edges\t", complex->get_no_halfedges());
//    write_variable("#faces\t", complex->get_no_faces());
//}
