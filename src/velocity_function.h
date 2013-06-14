//
//  velocity_function.h
//  3D_DSC
//
//  Created by Asger Nyman Christiansen on 9/27/12.
//  Copyright (c) 2012 DTU Informatics. All rights reserved.
//

#ifndef _DDSC_velocity_func_h
#define _DDSC_velocity_func_h

#include "DSC.h"
#include "ctime"

/**
 An abstract class which a specific velocity function should enherit from.
 */
template<class MT>
class VelocityFunc
{
    int time_step;
    int MAX_TIME_STEPS;
    
    double compute_time;
    double deform_time;
    
    double total_compute_time;
    double total_deform_time;
    
protected:
    DeformableSimplicialComplex<MT> *dsc;
    
    double VELOCITY; // Determines the distance each interface vertex moves at each iteration.
    double ACCURACY; // Determines the accuracy of the final result.
    
    std::vector<typename MT::vector3_type> pos_old;
    
    /**
     Creates a velocity function which is applied to the simplicial complex defined by the first input parameter. The velocity parameter determines the velocity of the function.
     */
    VelocityFunc(DeformableSimplicialComplex<MT> *dsc_, double velocity, double accuracy, int max_time_steps):
        dsc(dsc_), time_step(0), MAX_TIME_STEPS(max_time_steps), VELOCITY(velocity), ACCURACY(accuracy),
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
    double get_velocity() const
    {
        return VELOCITY;
    }
    
    /**
     Returns the accuracy.
     */
    double get_accuracy() const
    {
        return ACCURACY;
    }
    
    /**
     Returns the time it took to deform the interface in this time step.
     */
    double get_deform_time() const
    {
        return deform_time;
    }
    
    /**
     Returns the time it took to compute the new positions of the interface in this time step.
     */
    double get_compute_time() const
    {
        return compute_time;
    }
    
    /**
     Returns the total time it took to deform the interface.
     */
    double get_total_deform_time() const
    {
        return total_deform_time;
    }
    
    /**
     Returns the total time it took to compute the new positions of the interface.
     */
    double get_total_compute_time() const
    {
        return total_compute_time;
    }
    
protected:
    /**
     Updates the time it took to compute new positions for the interface vertices.
     */
    void update_compute_time(const clock_t& compute_time_)
    {
        double t = (double)(compute_time_) / ((double)CLOCKS_PER_SEC);
        compute_time += t;
        total_compute_time += t;
    }
    /**
     Updates the time it took to deform the interface.
     */
    void update_deform_time(const clock_t& deform_time_)
    {
        double t = (double)(deform_time_) / ((double)CLOCKS_PER_SEC);
        deform_time += t;
        total_deform_time += t;
    }
    
public:
    /**
     Computes the motion of each interface vertex and stores the new position in new_pos in the simplicial complex class.
     */
    virtual void deform() = 0;
    
    /**
     Returns wether the motion has finished.
     */
    virtual bool is_motion_finished()
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
        return false;
    }
    
    /**
     Takes one time step thereby deforming the simplicial complex according to the velocity function.
     */
    bool take_time_step()
    {
        compute_time = 0.;
        deform_time = 0.;
        
        deform();
        
        time_step++;
        
        return (time_step == MAX_TIME_STEPS) | is_motion_finished();
    }
    
    /**
     An optional test function which can be used to test some aspect of the velocity function.
     */
    virtual void test()
    {
        
    }
    
};

#endif
