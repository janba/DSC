//
//  sample.cpp
//  3D_DSC
//
//  Created by Asger Nyman Christiansen on 2/14/13.
//  Copyright (c) 2013 Asger Nyman Christiansen. All rights reserved.
//

#include "rotate_function.h"
#include "average_function.h"
#include "normal_function.h"
#include "user_interface.h"
#include "tetrahedralize.h"

class DemoUI : public UI
{    
    
public:
    
    DemoUI(int &argc, char** argv): UI(argc, argv)
    {

	}
    
protected:
    
    virtual void motion1()
    {
        typedef typename GELTypes::vector3_type V;
        stop();
        // Build the Simplicial Complex
        vector<double> points;
        vector<int>  tets;
        vector<int>  tet_labels;
        vector<V> pts_inside(1);
        pts_inside[0] = V( 0.0f, 0.0f, 0.0f);
        
        //import_tet_mesh(get_data_file_path("bunny.dsc").data(), points, tets, tet_labels);
        build_tetrahedralization<GELTypes>(get_data_file_path("blob.obj"), points, tets, tet_labels, pts_inside);
        
        DeformableSimplicialComplex<GELTypes> *complex = new DeformableSimplicialComplex<GELTypes>(points, tets, tet_labels);
        
        VelocityFunc<GELTypes> *vel_fun = new RotateFunc<GELTypes>(VELOCITY, ACCURACY);
        Log *log = new Log(create_log_path());
        dsc = new DSC<GELTypes>(vel_fun, complex, log);
        
        view_ctrl = new GLGraphics::GLViewController(WIN_SIZE_X, WIN_SIZE_Y, (CGLA::Vec3f)complex->get_center(), r);
        
        update_title();
        reshape(WIN_SIZE_X, WIN_SIZE_Y);
    }
    
    virtual void motion2()
    {        
        typedef typename GELTypes::vector3_type V;
        stop();
        // Build the Simplicial Complex
        vector<double> points;
        vector<int>  tets;
        vector<int>  tet_labels;
        vector<V> pts_inside(1);
        pts_inside[0] = V( 0.0f, 0.0f, 0.0f);
        
        //import_tet_mesh(get_data_file_path("bunny.dsc").data(), points, tets, tet_labels);
        build_tetrahedralization<GELTypes>(get_data_file_path("blob.obj"), points, tets, tet_labels, pts_inside);
        
        DeformableSimplicialComplex<GELTypes> *complex = new DeformableSimplicialComplex<GELTypes>(points, tets, tet_labels);
        
        VelocityFunc<GELTypes> *vel_fun = new AverageFunc<GELTypes>(VELOCITY, ACCURACY);
        Log *log = new Log(create_log_path());
        dsc = new DSC<GELTypes>(vel_fun, complex, log);
        
        view_ctrl = new GLGraphics::GLViewController(WIN_SIZE_X, WIN_SIZE_Y, (CGLA::Vec3f)complex->get_center(), r);
        
        update_title();
        reshape(WIN_SIZE_X, WIN_SIZE_Y);
    }
    
    virtual void motion3()
    {
        typedef typename GELTypes::vector3_type V;
        stop();
        // Build the Simplicial Complex
        vector<double> points;
        vector<int>  tets;
        vector<int>  tet_labels;
        vector<V> pts_inside(1);
        pts_inside[0] = V( 0.0f, 0.0f, 0.0f);
        
        //import_tet_mesh(get_data_file_path("bunny.dsc").data(), points, tets, tet_labels);
        build_tetrahedralization<GELTypes>(get_data_file_path("armadillo-very-simple.obj"), points, tets, tet_labels, pts_inside);
        
        DeformableSimplicialComplex<GELTypes> *complex = new DeformableSimplicialComplex<GELTypes>(points, tets, tet_labels);
        
        VelocityFunc<GELTypes> *vel_fun = new NormalFunc<GELTypes>(VELOCITY, ACCURACY);
        Log *log = new Log(create_log_path());
        dsc = new DSC<GELTypes>(vel_fun, complex, log);
        
        view_ctrl = new GLGraphics::GLViewController(WIN_SIZE_X, WIN_SIZE_Y, (CGLA::Vec3f)complex->get_center(), r);
        view_ctrl->set_view_param(CGLA::Vec3f(0.,0., -r), (CGLA::Vec3f)complex->get_center(), CGLA::Vec3f(0.,1.,0.));
        
        update_title();
        reshape(WIN_SIZE_X, WIN_SIZE_Y);
    }
};


int main(int argc, char** argv)
{
    DemoUI ui(argc, argv);
    glutMainLoop();
    return 0;
}