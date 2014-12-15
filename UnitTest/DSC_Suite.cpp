//
// Created by Morten Nobel-Jørgensen on 12/12/14.
// Copyright (c) 2014 Asger Nyman Christiansen. All rights reserved.
//

#include "DSC_Suite.h"
#include "mesh_io.h"
#include "tetralizer.h"
#include <vector>

bool equal(vec3 v1, vec3 v2){

    for (int i=0;i<3;i++){
        if (abs(v1[i]-v2[i])>0.01f){
            return false;
        }
    }
    return true;
}

int tetGenTest(void) {
//    int* p = NULL;
//    TINYTEST_ASSERT(!p);
//    TINYTEST_ASSERT(!printf(""));
    vector<vec3> points;
    vector<int> tets;
    vector<int> tet_labels;

    std::vector<vec3> points_interface;
    std::vector<int> faces_interface;
    is_mesh::import_surface_mesh("data/blob.obj", points_interface, faces_interface);

    Tetralizer::tetralize(vec3(3.), 0.5, points_interface, faces_interface, points, tets, tet_labels);

    is_mesh::export_tet_mesh( "data/blob-test.dsc", points, tets, tet_labels);

    vector<vec3> points2;
    vector<int> tets2;
    vector<int> tet_labels2;
    is_mesh::import_tet_mesh( "data/blob-test.dsc", points2, tets2, tet_labels2);

    TINYTEST_ASSERT(points.size() == points2.size());
    TINYTEST_ASSERT(tets.size() == tets2.size());
    TINYTEST_ASSERT(tet_labels.size() == tet_labels2.size());

    for (int i=0;i<points.size();i++){
        TINYTEST_ASSERT(equal(points[i], points2[i]));
    }
    for (int i=0;i<tets.size();i++){
        TINYTEST_ASSERT(tets[i] == tets2[i]);
        TINYTEST_ASSERT(tets[i] < points.size());

    }
    for (int i=0;i<tet_labels.size();i++){
        TINYTEST_ASSERT(tet_labels[i] == tet_labels2[i]);
    }

    return 1; // Always return a value different than 0 at test end.
}