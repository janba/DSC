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
        is_mesh::FaceKey first_intersection;
        double dist;
        CGLA::Ray ray;
        QueryType query_type;
        is_mesh::SimplexSet<is_mesh::TetrahedronKey> tetrahedron_key;
        is_mesh::FaceKey face_key;
        is_mesh::ISMesh *mesh;
        void next();
    public:

        QueryResultIterator();
        QueryResultIterator(is_mesh::FaceKey const &first_intersection, double dist, CGLA::Ray const &ray, QueryType const &query_type, is_mesh::ISMesh *mesh);
        QueryResultIterator(const QueryResultIterator& iter);

        is_mesh::FaceKey operator*();
        QueryResultIterator &operator++();
        bool operator!=(const QueryResultIterator & other);
    };

    class QueryResult {
        is_mesh::FaceKey first_intersection;
        double dist;
        CGLA::Ray ray;
        QueryType query_type;
        is_mesh::ISMesh *mesh;
    public:
        QueryResult(){}
        QueryResult(is_mesh::FaceKey const &first_intersection, double dist, CGLA::Ray const &ray, QueryType const &query_type, is_mesh::ISMesh *mesh);

        QueryResultIterator begin();
        QueryResultIterator end();
    };

    class Query {
        is_mesh::ISMesh *mesh;
        std::vector<is_mesh::FaceKey> boundary_faces;

    public:
        Query(is_mesh::ISMesh *mesh);
        QueryResult raycast_faces(CGLA::Ray ray, QueryType queryType = QueryType::All);

        // Must be called explicit prior to raycasting if the boundary triangles of the ISMesh has changed
        void rebuild_boundary_cache();
    };
}



