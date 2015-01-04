//
// Created by Morten Nobel-Jørgensen on 19/11/14.
// Copyright (c) 2014 ___FULLUSERNAME___. All rights reserved.
//


#pragma once

#include "DSC.h"
#include "CGLA/Ray.h"
#include <vector>

/**
* Raytracing in ISMesh datastructure.
* Assumes that the boundary surface of the ISMesh is convex
* Assumes that the set of boundary faces are static (if not call rebuild_boundary_cache() after changes).
*/
namespace is_mesh {
    enum class QueryType {
        All,
        Interface,
        Boundary
    };

    class QueryResultIterator {
        FaceKey first_intersection;
        double dist;
        CGLA::Ray ray;
        QueryType query_type;
        SimplexSet<TetrahedronKey> tetrahedron_key;
        FaceKey face_key;
        ISMesh *mesh;
        void next();
    public:

        QueryResultIterator();
        QueryResultIterator(FaceKey const &first_intersection, double dist, CGLA::Ray const &ray, QueryType const &query_type, ISMesh *mesh);
        QueryResultIterator(const QueryResultIterator& iter);

        FaceKey operator*();
        QueryResultIterator &operator++();
        bool operator!=(const QueryResultIterator & other);
    };

    class QueryResult {
        FaceKey first_intersection;
        double dist;
        CGLA::Ray ray;
        QueryType query_type;
        ISMesh *mesh;
    public:
        QueryResult():mesh(nullptr){}
        QueryResult(FaceKey const &first_intersection, double dist, CGLA::Ray const &ray, QueryType const &query_type, ISMesh *mesh);

        QueryResultIterator begin();
        QueryResultIterator end();
    };

    class Query {
        ISMesh *mesh;
        std::vector<FaceKey> *boundary_faces = nullptr;

    public:
        Query(ISMesh *mesh);
        ~Query();
        QueryResult raycast_faces(CGLA::Ray ray, QueryType queryType = QueryType::All);

        // Must be called explicit prior to raycasting if the boundary triangles of the ISMesh has changed
        void rebuild_boundary_cache();

        // returns the nodes within max distance to a given from_node
        // includes the from_node in the result
        std::vector<NodeKey> neighborhood(NodeKey from_node, double max_distance);

        // return the list of edges where both nodes are contained in the nodeKeys
        std::vector<EdgeKey> edges(std::vector<NodeKey> nodeKeys);

        std::vector<FaceKey> faces(std::vector<EdgeKey> edgeKeys);

        std::vector<TetrahedronKey> tetrahedra(std::vector<FaceKey> faceKeys);
    };
}



