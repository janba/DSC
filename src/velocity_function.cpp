#include "velocity_function.h"

DSC::VelocityFunc::VelocityFunc(double velocity, double accuracy, int max_time_steps) : MAX_TIME_STEPS(max_time_steps) {
    set_velocity(velocity);
    set_accuracy(accuracy);
}

DSC::VelocityFunc::~VelocityFunc() {
    pos_old.clear();
}

std::string DSC::VelocityFunc::get_name() const {
    return std::string("NO MOTION");
}

int DSC::VelocityFunc::get_time_step() const {
    return time_step;
}

void DSC::VelocityFunc::set_max_time_steps(int max_time_steps) {
    MAX_TIME_STEPS = max_time_steps;
}

double DSC::VelocityFunc::get_velocity() const {
    return VELOCITY;
}

void DSC::VelocityFunc::set_velocity(double vel) {
    VELOCITY = vel;
}

double DSC::VelocityFunc::get_accuracy() const {
    return ACCURACY;
}

void DSC::VelocityFunc::set_accuracy(double acc) {
    ACCURACY = acc;
}

double DSC::VelocityFunc::get_deform_time() const {
    return deform_time;
}

double DSC::VelocityFunc::get_compute_time() const {
    return compute_time;
}

double DSC::VelocityFunc::get_total_deform_time() const {
    return total_deform_time;
}

double DSC::VelocityFunc::get_total_compute_time() const {
    return total_compute_time;
}

void DSC::VelocityFunc::update_compute_time(const std::chrono::time_point<std::chrono::system_clock>& start_time) {
    std::chrono::duration<double> t = std::chrono::system_clock::now() - start_time;
    compute_time += t.count();
    total_compute_time += t.count();
}

void DSC::VelocityFunc::update_deform_time(const std::chrono::time_point<std::chrono::system_clock>& start_time) {
    std::chrono::duration<double> t = std::chrono::system_clock::now() - start_time;
    deform_time += t.count();
    total_deform_time += t.count();
}

void DSC::VelocityFunc::deform(DeformableSimplicialComplex& dsc) {
    auto init_time = std::chrono::system_clock::now();

    dsc.deform(deform_time_steps);

    update_deform_time(init_time);
}

bool DSC::VelocityFunc::is_motion_finished(DeformableSimplicialComplex& dsc) {
    if(time_step < MAX_TIME_STEPS)
    {
        for (auto & nit : dsc.nodes())
        {
            if(dsc.is_movable(nit.key()))
            {
                bool match = false;
                for (int i = 0; i+2 < pos_old.size(); i += 3)
                {
                    if (Util::distance_point_triangle(nit.get_pos(), pos_old[i], pos_old[i+1], pos_old[i+2]) < ACCURACY)
                    {
                        match = true;
                        break;
                    }
                }
                if (!match) {
#ifdef DEBUG
                    std::cout << "Stopping criteria: Position " << nit.get_pos() << " has moved." << std::endl;
#endif
                    pos_old = dsc.get_interface_face_positions();
                    return false;
                }
            }
        }
        pos_old = dsc.get_interface_face_positions();
    }
    return true;
}

void DSC::VelocityFunc::take_time_step(DeformableSimplicialComplex& dsc) {
    compute_time = 0.;
    deform_time = 0.;

    deform(dsc);

    time_step++;
}

void DSC::VelocityFunc::test(DeformableSimplicialComplex& dsc) {
    dsc.get_is_mesh().validity_check();

    dsc.test_flip23_flip32();
    dsc.test_split_collapse();
    dsc.test_flip44();
    dsc.test_flip22();
}