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

/**
 An abstract class which a specific velocity function should enherit from.
 */
template<class MT>
class VelocityFunc
{
    int time_step;
    int MAX_TIME_STEPS;
    
    real compute_time;
    real deform_time;
    
    real total_compute_time;
    real total_deform_time;
    
protected:
    real VELOCITY; // Determines the distance each interface vertex moves at each iteration.
    real ACCURACY; // Determines the accuracy of the final result.
    
    std::vector<vec3> pos_old;
    
    /**
     Creates a velocity function which is applied to the simplicial complex defined by the first input parameter. The velocity parameter determines the velocity of the function.
     */
    VelocityFunc(real velocity, real accuracy, int max_time_steps):
        time_step(0), MAX_TIME_STEPS(max_time_steps), VELOCITY(velocity), ACCURACY(accuracy),
        compute_time(0.), deform_time(0.), total_compute_time(0.), total_deform_time(0.)
    {
        
    }
    
public:
    virtual ~VelocityFunc()
    {
        pos_old.clear();
    }
    
    /**
     Returns the name of the velocity function.
     */
    virtual std::string get_name() const = 0;
    
    /**
     Returns the current time step.
     */
    int get_time_step() const
    {
        return time_step;
    }
    
    /**
     Returns the velocity.
     */
    real get_velocity() const
    {
        return VELOCITY;
    }
    
    /**
     Returns the accuracy.
     */
    real get_accuracy() const
    {
        return ACCURACY;
    }
    
    /**
     Returns the time it took to deform the interface in this time step.
     */
    real get_deform_time() const
    {
        return deform_time;
    }
    
    /**
     Returns the time it took to compute the new positions of the interface in this time step.
     */
    real get_compute_time() const
    {
        return compute_time;
    }
    
    /**
     Returns the total time it took to deform the interface.
     */
    real get_total_deform_time() const
    {
        return total_deform_time;
    }
    
    /**
     Returns the total time it took to compute the new positions of the interface.
     */
    real get_total_compute_time() const
    {
        return total_compute_time;
    }
    
protected:
    /**
     Updates the time it took to compute new positions for the interface vertices.
     */
    void update_compute_time(const std::chrono::time_point<std::chrono::system_clock>& start_time)
    {
        std::chrono::duration<real> t = std::chrono::system_clock::now() - start_time;
        compute_time += t.count();
        total_compute_time += t.count();
    }
    /**
     Updates the time it took to deform the interface.
     */
    void update_deform_time(const std::chrono::time_point<std::chrono::system_clock>& start_time)
    {
        std::chrono::duration<real> t = std::chrono::system_clock::now() - start_time;
        deform_time += t.count();
        total_deform_time += t.count();
    }
    
    /**
     Computes the motion of each interface vertex and stores the new position in new_pos in the simplicial complex class.
     */
    virtual void deform(DeformableSimplicialComplex<MT>& dsc) = 0;
    
public:
    /**
     Returns wether the motion has finished.
     */
    virtual bool is_motion_finished(const DeformableSimplicialComplex<MT>& dsc)
    {
//        std::vector<typename MT::vector3_type> pos = complex->get_design_variable_positions();
//        for (auto p = pos.begin(); p != pos.end(); p++)
//        {
//            bool match = false;
//            for (int i = 0; i+1 < pos_old.size(); i += 2)
//            {
//                if (min_dist(pos_old[i], pos_old[i+1], *p) < ACCURACY)
//                {
//                    match = true;
//                    break;
//                }
//            }
//            if (!match) {
//                std::cout << "Stopping criteria: Position " << *p << " has moved." << std::endl;
//                pos_old = complex->get_interface_edge_positions();
//                return false;
//            }
//        }
//        pos_old = complex->get_interface_edge_positions();
//        return true;
        return time_step == MAX_TIME_STEPS;
    }
    
    /**
     Takes one time step thereby deforming the simplicial complex according to the velocity function.
     */
    void take_time_step(DeformableSimplicialComplex<MT>& dsc)
    {
        compute_time = 0.;
        deform_time = 0.;
        
        deform(dsc);
        
        time_step++;
    }
    
    /**
     An optional test function which can be used to test some aspect of the velocity function.
     */
    virtual void test(DeformableSimplicialComplex<MT>& dsc)
    {
        
    }
    
};
