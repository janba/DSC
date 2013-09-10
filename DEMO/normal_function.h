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

#include "velocity_function.h"

/**
 A velocity function which moves the interface vertices in the normal direction.
 */
template<class MT>
class NormalFunc: public VelocityFunc<MT> {
    
    
public:
    /**
     Creates a velocity function which moves the interface vertices in the normal direction.
     */
    NormalFunc(double velocity, double accuracy, int max_time_steps = 500):
        VelocityFunc<MT>(velocity/10., accuracy, max_time_steps)
    {
        
    }
    
    /**
     Returns the name of the velocity function.
     */
    virtual std::string get_name() const
    {
        return std::string("NORMAL MOTION");
    }
    
    /**
     Computes the motion of each interface vertex and stores the new position in new_pos in the simplicial complex class.
     */
    virtual void deform(DeformableSimplicialComplex<MT>& dsc)
    {
        auto init_time = std::chrono::system_clock::now();
        vec3 new_pos;
        for(auto nit = dsc.nodes_begin(); nit != dsc.nodes_end(); nit++)
        {
            if(nit->is_interface() && !nit->is_crossing())
            {
                new_pos = dsc.get_pos(nit.key()) + VelocityFunc<MT>::VELOCITY * dsc.get_normal(nit.key());
                dsc.set_destination(nit.key(), new_pos);
            }
        }
        VelocityFunc<MT>::update_compute_time(init_time);
        init_time = std::chrono::system_clock::now();
        
        dsc.deform();
        
        VelocityFunc<MT>::update_deform_time(init_time);
    }
};
