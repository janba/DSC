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

#include <chrono>
#include <ctime>
#include "DSC.h"

namespace DSC {
    
    /**
     An abstract class which a specific velocity function should inherit from.
     */
    class VelocityFunc
    {
        double compute_time = 0.;
        double deform_time = 0.;
        
        double total_compute_time = 0.;
        double total_deform_time = 0.;
        
    protected:
        int time_step = 0;
        int MAX_TIME_STEPS;
        
        double VELOCITY; // Determines the distance each interface vertex moves at each iteration.
        double ACCURACY; // Determines the accuracy of the final result.
        
        std::vector<vec3> pos_old;
        
    public:
        /**
        Creates a velocity function which is applied to the simplicial complex defined by the first input parameter. The velocity parameter determines the velocity of the function.
        */
        VelocityFunc(double velocity, double accuracy, int max_time_steps);

        virtual ~VelocityFunc();

        /**
        Returns the name of the velocity function.
        */
        virtual std::string get_name() const;

        /**
        Returns the current time step.
        */
        int get_time_step() const;

        virtual void set_max_time_steps(int max_time_steps);

        /**
        Returns the velocity.
        */
        double get_velocity() const;

        virtual void set_velocity(double vel);

        /**
        Returns the accuracy.
        */
        double get_accuracy() const;

        virtual void set_accuracy(double acc);

        /**
        Returns the time it took to deform the interface in this time step.
        */
        double get_deform_time() const;

        /**
        Returns the time it took to compute the new positions of the interface in this time step.
        */
        double get_compute_time() const;

        /**
        Returns the total time it took to deform the interface.
        */
        double get_total_deform_time() const;

        /**
        Returns the total time it took to compute the new positions of the interface.
        */
        double get_total_compute_time() const;
        
    protected:
        /**
        Updates the time it took to compute new positions for the interface vertices.
        */
        void update_compute_time(const std::chrono::time_point<std::chrono::system_clock>& start_time);

        /**
        Updates the time it took to deform the interface.
        */
        void update_deform_time(const std::chrono::time_point<std::chrono::system_clock>& start_time);
        
        /**
         Computes the motion of each interface vertex and stores the new position in new_pos in the simplicial complex class.
         */
        virtual void deform(DeformableSimplicialComplex& dsc);
        
    public:
        /**
         Returns whether the motion has finished.
         */
        virtual bool is_motion_finished(DeformableSimplicialComplex& dsc);
        
        /**
         Takes one time step thereby deforming the simplicial complex according to the velocity function.
         */
        void take_time_step(DeformableSimplicialComplex& dsc);
        
        /**
         An optional test function which can be used to test some aspect of the velocity function.
         */
        virtual void test(DeformableSimplicialComplex& dsc);
        
    };
    
}
