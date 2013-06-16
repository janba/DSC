//
//  normal_function.h
//  2D_DSC
//
//  Created by Asger Nyman Christiansen on 2/14/13.
//  Copyright (c) 2013 Asger Nyman Christiansen. All rights reserved.
//

#ifndef ___D_DSC__normal_func__
#define ___D_DSC__normal_func__

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
        VelocityFunc<MT>(velocity/100., accuracy/100., max_time_steps)
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
        typedef typename MT::vector3_type V;
        clock_t init_time = clock();
        for(auto nit = dsc.nodes_begin(); nit != dsc.nodes_end(); nit++)
        {
            if(nit->is_interface())
            {
                V p = dsc.get_pos(nit.key());
                V p_new = p + VelocityFunc<MT>::VELOCITY * dsc.get_normal(nit.key());
                dsc.set_destination(nit.key(), p_new);
            }
        }
        VelocityFunc<MT>::update_compute_time(clock() - init_time);
        init_time = clock();
        
        dsc.deform();
        
        VelocityFunc<MT>::update_deform_time(clock() - init_time);
    }
};

#endif /* defined(___D_DSC__normal_func__) */
