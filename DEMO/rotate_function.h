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

#ifndef ___D_DSC__rotate_func__
#define ___D_DSC__rotate_func__

#include "velocity_function.h"


/**
 A rotating velocity function.
 */
template <class MT>
class RotateFunc: public VelocityFunc<MT>
{
    
    
public:
    /**
     Creates a rotating velocity function.
     */
    RotateFunc(double velocity, double accuracy, int max_time_steps = 500):
        VelocityFunc<MT>(velocity/5., accuracy, max_time_steps)
    {
        
    }
    
    /**
     Returns the name of the velocity function.
     */
    virtual std::string get_name() const
    {
        return std::string("ROTATE MOTION");
    }
    
    /**
     Computes the motion of each interface vertex and stores the new position in new_pos in the simplicial complex class.
     */
    virtual void deform(DeformableSimplicialComplex<MT>& dsc)
    {
        typedef typename MT::vector3_type V;
        typedef typename MT::vector4_type V4;
        clock_t init_time = clock();
        
        V center = dsc.get_center();
        typename MT::matrix4x4_type mrot = MT::get_rotation_matrix(MT::get_z_axis(), VelocityFunc<MT>::VELOCITY);
        V p;
        V new_pos;
        for(auto nit = dsc.nodes_begin(); nit != dsc.nodes_end(); nit++)
        {
            if(nit->is_interface() && !nit->is_crossing())
            {
                V a = dsc.get_pos(nit.key()) - center;
                V4 a4(a[0], a[1], a[2], 1.);
                V4 b = mrot * a4;
                new_pos = center + V(b[0], b[1], b[2]);
                dsc.set_destination(nit.key(), new_pos);
            }
        }
        VelocityFunc<MT>::update_compute_time(clock() - init_time);
        init_time = clock();
        
        dsc.deform();
        
        VelocityFunc<MT>::update_deform_time(clock() - init_time);
    }
    
    /**
     Returns wether the motion has finished.
     */
    virtual bool is_motion_finished()
    {
        return false;
    }
    
    virtual void test(DeformableSimplicialComplex<MT>& dsc)
    {
        std::vector<typename DeformableSimplicialComplex<MT>::edge_key> edges;
        for (auto eit = dsc.edges_begin(); eit != dsc.edges_end(); eit++)
        {
            if (eit->is_interface())
            {
                edges.push_back(eit.key());
            }
        }
                
//        for (auto e : edges) {
//            dsc.split_edge(e);
//        }
//        
//        for (auto e : edges) {
//            if(dsc.exists(e))
//            {
//                dsc.collapse(e);
//                break;
//            }
//        }
        
//        for (auto e : edges) {
//            if(dsc.exists(e) && dsc.is_flippable(e))
//            {
//                dsc.flip(e);
//            }
//        }
//        dsc.validity_check();
    }
    
};

#endif /* defined(___D_DSC__rotate_func__) */
