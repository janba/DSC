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
     An abstract class which a specific velocity function should enherit from.
     */
    template <typename DeformableSimplicialComplex = DeformableSimplicialComplex<>>
    class VelocityFunc
    {
        real compute_time = 0.;
        real deform_time = 0.;
        
        real total_compute_time = 0.;
        real total_deform_time = 0.;
        
    protected:
        int time_step = 0;
        int MAX_TIME_STEPS;
        
        real VELOCITY; // Determines the distance each interface vertex moves at each iteration.
        real ACCURACY; // Determines the accuracy of the final result.
        
        std::vector<vec3> pos_old;
        
    public:
        /**
         Creates a velocity function which is applied to the simplicial complex defined by the first input parameter. The velocity parameter determines the velocity of the function.
         */
        VelocityFunc(real velocity, real accuracy, int max_time_steps): MAX_TIME_STEPS(max_time_steps)
        {
            set_velocity(velocity);
            set_accuracy(accuracy);
        }
        
        virtual ~VelocityFunc()
        {
            pos_old.clear();
        }
        
        /**
         Returns the name of the velocity function.
         */
        virtual std::string get_name() const
        {
            return std::string("NO MOTION");
        }
        
        /**
         Returns the current time step.
         */
        int get_time_step() const
        {
            return time_step;
        }
        
        virtual void set_max_time_steps(int max_time_steps)
        {
            MAX_TIME_STEPS = max_time_steps;
        }
        
        /**
         Returns the velocity.
         */
        real get_velocity() const
        {
            return VELOCITY;
        }
        
        virtual void set_velocity(real vel)
        {
            VELOCITY = vel;
        }
        
        /**
         Returns the accuracy.
         */
        real get_accuracy() const
        {
            return ACCURACY;
        }
        
        virtual void set_accuracy(real acc)
        {
            ACCURACY = acc;
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
        virtual void deform(DeformableSimplicialComplex& dsc)
        {
            auto init_time = std::chrono::system_clock::now();
            
            dsc.deform();
            
            update_deform_time(init_time);
        }
        
    public:
        /**
         Returns wether the motion has finished.
         */
        virtual bool is_motion_finished(DeformableSimplicialComplex& dsc)
        {
            if(time_step < MAX_TIME_STEPS)
            {
                for (auto nit = dsc.nodes_begin(); nit != dsc.nodes_end(); nit++)
                {
                    if(dsc.is_movable(nit.key()))
                    {
                        bool match = false;
                        for (int i = 0; i+2 < pos_old.size(); i += 3)
                        {
                            if (Util::distance_point_triangle<real>(nit->get_pos(), pos_old[i], pos_old[i+1], pos_old[i+2]) < ACCURACY)
                            {
                                match = true;
                                break;
                            }
                        }
                        if (!match) {
                            std::cout << "Stopping criteria: Position " << nit->get_pos() << " has moved." << std::endl;
                            pos_old = dsc.get_interface_face_positions();
                            return false;
                        }
                    }
                }
                pos_old = dsc.get_interface_face_positions();
            }
            return true;
        }
        
        /**
         Takes one time step thereby deforming the simplicial complex according to the velocity function.
         */
        void take_time_step(DeformableSimplicialComplex& dsc)
        {
            compute_time = 0.;
            deform_time = 0.;
            
            deform(dsc);
            
            time_step++;
        }
        
        /**
         An optional test function which can be used to test some aspect of the velocity function.
         */
        virtual void test(DeformableSimplicialComplex& dsc)
        {
            dsc.validity_check();
            
            dsc.test_flip23_flip32();
            dsc.test_split_collapse();
            dsc.test_flip44();
            dsc.test_flip22();
        }
        
    };
    
}
