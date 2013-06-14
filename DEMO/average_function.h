//
//  average_function.h
//  2D_DSC
//
//  Created by Asger Nyman Christiansen on 2/14/13.
//  Copyright (c) 2013 Asger Nyman Christiansen. All rights reserved.
//

#ifndef ___D_DSC__average_func__
#define ___D_DSC__average_func__

#include "velocity_function.h"

/**
 A velocity function which moves the interface vertices towards the average of their neighbouring interface vertices, i.e. a constant smoothing of the interface.
 */
template<class MT>
class AverageFunc: public VelocityFunc<MT> {
    
public:
    /**
     Creates a velocity function which smooths the interface.
     */
    AverageFunc(double velocity, double accuracy): VelocityFunc<MT>(velocity/100., accuracy/100.)
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
    virtual void deform(DeformableSimplicialComplex<MT> *complex)
    {
        typedef typename MT::vector3_type V;
        clock_t init_time = clock();
        for(auto nit = complex->nodes_begin(); nit != complex->nodes_end(); nit++)
        {
            if(nit->is_interface())
            {
                V p = complex->get_pos(nit.key());
                V p_new = p + VelocityFunc<MT>::VELOCITY * (complex->get_barycenter(nit.key(), true) - p);
                complex->set_destination(nit.key(), p_new);
            }
        }
        VelocityFunc<MT>::update_compute_time(clock() - init_time);
        init_time = clock();
        
        complex->deform();
        
        VelocityFunc<MT>::update_deform_time(clock() - init_time);
    }
};

#endif /* defined(___D_DSC__average_func__) */
