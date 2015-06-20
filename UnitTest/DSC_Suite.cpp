//
// Created by Morten Nobel-Jørgensen on 12/12/14.
// Copyright (c) 2014 Asger Nyman Christiansen. All rights reserved.
//

#include "DSC_Suite.h"
#include "mesh_io.h"
#include "tetralizer.h"
#include "is_mesh.h"
#include "query.h"
#include <vector>
#include <chrono>

bool equal(vec3 v1, vec3 v2){

    for (int i=0;i<3;i++){
        if (abs(v1[i]-v2[i])>0.01f){
            return false;
        }
    }
    return true;
}

int build_boundary_mesh_test(void){
    std::vector<double> points_boundary;
    double avg_edge_length = 0.5;
    std::vector<int> faces_boundary;
    const vec3 min{-2,-3,-4};
    const vec3 max{5};
    Tetralizer::build_boundary_mesh(points_boundary, avg_edge_length,faces_boundary,min, max);

    std::vector<vec3> pos;
    for (int i=0;i<points_boundary.size();i+=3){
        pos.push_back(vec3{points_boundary[i],points_boundary[i+1],points_boundary[i+2]});
    }
    for (auto & i : faces_boundary){
        i++;
    }

    is_mesh::export_surface_mesh("data/output.obj", pos, faces_boundary);
    return 1;
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

    Tetralizer::tetralize(0.5f, 0.5f, points_interface, faces_interface, points, tets, tet_labels);

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

int connectedTest() {
    using namespace is_mesh;
    vector<vec3> points2;
    vector<int> tets2;
    vector<int> tet_labels2;
    import_tet_mesh( "data/blob-test.dsc", points2, tets2, tet_labels2);

    ISMesh mesh(points2, tets2, tet_labels2);

    int count = 0;
    int total = 0;
    TetrahedronKey someKey;
    for (auto & t : mesh.tetrahedra()) {
        if (t.label()){
            someKey = t.key();
            count++;
        }
        total++;
    }
    Query q{&mesh};
    auto res = q.connected<TetrahedronKey>(someKey, [&](TetrahedronKey k){return mesh.get(k).label()==1;});
    cout << "Finished "<< count<<"/"<<total<<endl;
    TINYTEST_ASSERT(res.size() == count);



    return 1;
}

int qualityTest() {
    using namespace is_mesh;
    typedef std::chrono::high_resolution_clock Clock;
    vector<vec3> points2;
    vector<int> tets2;
    vector<int> tet_labels2;
    import_tet_mesh( "data/blob-test.dsc", points2, tets2, tet_labels2);

    DSC::DeformableSimplicialComplex dsc(points2, tets2, tet_labels2);
    auto t1 = Clock::now();
    double q = dsc.min_quality();
    auto t2 = Clock::now();
    double avgEdgeLen = dsc.compute_avg_edge_length();
    auto t3 = Clock::now();

    cout << "min_quality "<<std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << "ms"<< endl;
    cout << "compute_avg_edge_length "<<std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count() << "ms"<< endl;
    TINYTEST_ASSERT(abs(avgEdgeLen-6.45694)<0.00001);
    TINYTEST_ASSERT(abs(q-0.00212245)<0.00001);

    std::cout << "Min quality: " <<q<<endl;
    std::cout.flush();

    return 1;
}





int forEachTest(){
    using namespace is_mesh;
    vector<vec3> points2;
    vector<int> tets2;
    vector<int> tet_labels2;
    import_tet_mesh( "data/blob-test.dsc", points2, tets2, tet_labels2);

    ISMesh mesh(points2, tets2, tet_labels2);

    auto fnSmooth = [](Node& node, int threadid){
        node.set_pos(node.smart_laplacian());
    };
    typedef std::chrono::high_resolution_clock Clock;
    auto t1 = Clock::now();

    mesh.for_each_par<Node>(fnSmooth);
    auto t2 = Clock::now();
    mesh.for_each_par_sp<Node>(0.1, 0, fnSmooth);
    auto t3 = Clock::now();
    // naive
    for (auto & node : mesh.nodes()){
        node.set_pos(node.smart_laplacian());
    }
    auto t4 = Clock::now();

    cout << "For each "<<std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << "ms"<< endl;
    cout << "For each sp " << std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count() << "ms"<< endl;
    cout << "Sequential " << std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3).count() << "ms"<< endl;
    return 1;
}

int rayTest() {
    using namespace is_mesh;
    using namespace CGLA;
    vector<vec3> points2;
    vector<int> tets2;
    vector<int> tet_labels2;
    import_tet_mesh( "data/cantilever-test.dsc", points2, tets2, tet_labels2);

    ISMesh mesh(points2, tets2, tet_labels2);
    Ray ray(Vec3d(1000,0.00,0.00),Vec3d(-1,0,0));
    int faceCount = 0;
    int faceInterfaces = 0;
    int faceIntersections = 0;
    for (Face &f:mesh.faces()){
        if (f.is_interface()){
            faceInterfaces++;
            auto pos = mesh.get_pos(f.node_keys());
            double dist;
            if (ray.intersect_triangle(pos[0], pos[1], pos[2], dist)){
                faceIntersections++;
            }
        }
        faceCount++;
    }
    cout <<"faces "<<faceCount<<" face interfaces "<<faceInterfaces<<" ray interface intersections "<< faceIntersections<<endl;

    Query q(&mesh);

    auto iter = q.raycast_faces(ray, QueryType::Interface);

    int count = 0;
    for (auto p : iter){
        cout <<"Collision point "<<endl;
        count++;
    }
    TINYTEST_ASSERT(2 == count);

    for (auto i = iter.begin();i != iter.end();++i){
        cout <<"Collision point "<<i.collision_point()<<endl;
    }

    return 1;
}

template<typename ptr_type>
bool allDifferent(ptr_type** ptr, int size){
    for (int i=0;i<size;i++){
        for (int j=i+1;j<size;j++) {
            if (ptr[i]==ptr[j]){
                return false;
            }
        }
    }
    return true;
}

int testMeshNavigation(){
    using namespace is_mesh;
    vec3 origin{0,0,0};
    vec3 voxel_size{1,1,1};
    std::vector<int> voxel_labels;
    std::vector<vec3> points;
    std::vector<int> tets;
    std::vector<int> tet_labels;

    Tetralizer::tetralize(origin, voxel_size, 2,2,2, voxel_labels, points, tets, tet_labels);

    is_mesh::ISMesh mesh(points, tets, tet_labels);
    NodeKey centerNodeKey;
    for (auto & n : mesh.nodes()){
        if (n.get_pos() == vec3(1,1,1)){
            centerNodeKey = n.key();
            // label single upper tet (with positive center values)
            for (auto neighborgTets : mesh.get_tets(centerNodeKey)){
                is_mesh::Tetrahedron & tet = mesh.get(neighborgTets);
                auto center = tet.get_center();
                if (center[0]>1.1 && center[1]>1.1 && center[2]>1.1){
                    mesh.set_label(neighborgTets, 1);

                    break;
                }
            }
            // label single upper tet (with negative center values)
            for (auto neighborgTets : mesh.get_tets(centerNodeKey)){
                is_mesh::Tetrahedron & tet = mesh.get(neighborgTets);
                auto center = tet.get_center();
                if (center[0] < 0.9 && center[1] < 0.9 && center[2] < 0.9){
                    mesh.set_label(neighborgTets, 1);
                    break;
                }
            }
            break;
        }
    }
    auto & node = mesh.get(centerNodeKey);
    TINYTEST_ASSERT(node.tets().size() == 8);
    TINYTEST_ASSERT(node.faces().size() == 12);
    TINYTEST_ASSERT(node.edges().size() == 6);
    TINYTEST_ASSERT(allDifferent(node.tets().data(), 8));
    TINYTEST_ASSERT(allDifferent(node.faces().data(), 12));
    TINYTEST_ASSERT(allDifferent(node.edges().data(), 6));

    auto & edge = *mesh.get(centerNodeKey).edges()[0];
    TINYTEST_ASSERT(edge.tets().size() == 4);
    TINYTEST_ASSERT(edge.faces().size() == 4);
    TINYTEST_ASSERT(edge.nodes().size() == 2);
    TINYTEST_ASSERT(allDifferent(edge.tets().data(), 4));
    TINYTEST_ASSERT(allDifferent(edge.faces().data(), 4));
    TINYTEST_ASSERT(allDifferent(edge.nodes().data(), 2));

    auto & face = *mesh.get(centerNodeKey).faces()[0];
    TINYTEST_ASSERT(face.tets().size() == 2);
    TINYTEST_ASSERT(face.edges().size() == 3);
    TINYTEST_ASSERT(face.nodes().size() == 3);
    TINYTEST_ASSERT(allDifferent(face.tets().data(), 2));
    TINYTEST_ASSERT(allDifferent(face.edges().data(), 3));
    TINYTEST_ASSERT(allDifferent(face.nodes().data(), 3));

    auto & tet = *mesh.get(centerNodeKey).tets()[0];
    TINYTEST_ASSERT(tet.faces().size() == 4);
    TINYTEST_ASSERT(tet.edges().size() == 6);
    TINYTEST_ASSERT(tet.nodes().size() == 4);
    TINYTEST_ASSERT(allDifferent(tet.faces().data(), 4));
    TINYTEST_ASSERT(allDifferent(tet.edges().data(), 6));
    TINYTEST_ASSERT(allDifferent(tet.nodes().data(), 4));

    return 1;
}

int testNumberOfClusters(){
    using namespace is_mesh;
    vec3 origin{0,0,0};
    vec3 voxel_size{1,1,1};
    std::vector<int> voxel_labels;
    std::vector<vec3> points;
    std::vector<int> tets;
    std::vector<int> tet_labels;

    Tetralizer::tetralize(origin, voxel_size, 2,2,2, voxel_labels, points, tets, tet_labels);

    is_mesh::ISMesh mesh(points, tets, tet_labels);
    NodeKey centerNodeKey;
    for (auto & n : mesh.nodes()){
        if (n.get_pos() == vec3(1,1,1)){
            centerNodeKey = n.key();
            // label single upper tet (with positive center values)
            for (auto neighborgTets : mesh.get_tets(centerNodeKey)){
                is_mesh::Tetrahedron & tet = mesh.get(neighborgTets);
                auto center = tet.get_center();
                if (center[0]>1.1 && center[1]>1.1 && center[2]>1.1){
                    mesh.set_label(neighborgTets, 1);

                    break;
                }
            }
            // label single upper tet (with negative center values)
            for (auto neighborgTets : mesh.get_tets(centerNodeKey)){
                is_mesh::Tetrahedron & tet = mesh.get(neighborgTets);
                auto center = tet.get_center();
                if (center[0] < 0.9 && center[1] < 0.9 && center[2] < 0.9){
                    mesh.set_label(neighborgTets, 1);
                    break;
                }
            }
            break;
        }
    }

    int labels = 0;
    for (auto & tet : mesh.tetrahedra()){
        if (tet.label() == 1){
            labels++;
        }
    }
    TINYTEST_ASSERT(labels==2);

    Node &node = mesh.get(centerNodeKey);
    TINYTEST_ASSERT(node.is_crossing());
    TINYTEST_ASSERT(!node.is_boundary());
    TINYTEST_ASSERT(node.is_interface());

    TINYTEST_ASSERT(node.get_number_of_neighbour_tet_clusters(1)==2);


    for (auto & tet : mesh.tetrahedra()){
        int label = tet.label();
        label = (label+1)%2;
        mesh.set_label(tet.key(), label);
    }

    TINYTEST_ASSERT(node.get_number_of_neighbour_tet_clusters(1)==1);

    return 1;
}