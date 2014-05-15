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
 A velocity function which moves the interface vertices towards the average of their neighbouring interface vertices, i.e. a constant smoothing of the interface.
 */
class AverageFunc: public DSC::VelocityFunc<> {
    
public:
    /**
     Creates a velocity function which smooths the interface.
     */
    AverageFunc(real velocity, real accuracy, int max_time_steps = 500):
        VelocityFunc<>(velocity, accuracy, max_time_steps)
    {
        
    }
    
    /**
     Returns the name of the velocity function.
     */
    virtual std::string get_name() const
    {
        return std::string("AVERAGE MOTION");
    }
    
    /**
     Computes the motion of each interface vertex and stores the new position in new_pos in the simplicial complex class.
     */
    virtual void deform(DSC::DeformableSimplicialComplex<>& dsc)
    {
        auto init_time = std::chrono::system_clock::now();
        vec3 new_pos, p;
        for(auto nit = dsc.nodes_begin(); nit != dsc.nodes_end(); nit++)
        {
            if(dsc.is_movable(nit.key()))
            {
                p = nit->get_pos();
                new_pos = p + 0.1*VELOCITY * (dsc.get_barycenter(nit.key(), true) - p);
                dsc.set_destination(nit.key(), new_pos);
            }
        }
        update_compute_time(init_time);
        VelocityFunc::deform(dsc);
    }
};
