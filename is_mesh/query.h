//
// Created by Morten Nobel-Jørgensen on 19/11/14.
// Copyright (c) 2014 ___FULLUSERNAME___. All rights reserved.
//


#pragma once

#include "DSC.h"
#include "CGLA/Ray.h"
#include "attribute_vector.h"
#include <vector>
#include <set>

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
        double dist;
        CGLA::Ray ray;
        QueryType query_type;
        SimplexSet<TetrahedronKey> tetrahedron_key;
        FaceKey face_key;
        ISMesh *mesh;
        void next();
    public:

        QueryResultIterator();
        QueryResultIterator(FaceKey const &first_boundary_intersection, double dist, CGLA::Ray const &ray, QueryType const &query_type, ISMesh *mesh);
        QueryResultIterator(const QueryResultIterator& iter);
        CGLA::Vec3d collision_point();
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
        std::set<NodeKey> neighborhood(NodeKey from_node, double max_distance);

        // returns the nodes within max distance to a given point
        std::set<NodeKey> neighborhood(vec3 from, double max_distance);

        // return the list of edges where both nodes are contained in the nodeKeys
        std::set<EdgeKey> edges(const std::set<NodeKey> nodeKeys);

        // return the list of faces where both nodes are contained in the edgeKeys
        std::set<FaceKey> faces(const std::set<EdgeKey> edgeKeys);

        // return the list of tets where both nodes are contained in the faceKeys
        std::set<TetrahedronKey> tetrahedra(const std::set<FaceKey> faceKeys);

        // exclude 'hanging' nodes, edges and faces. (such as nodes with only one edge etc)
        void filter_subset(std::set<NodeKey> &nodes, std::set<EdgeKey> &edges, std::set<FaceKey> &faces, std::set<TetrahedronKey> &tets);

        template <typename K>
        SimplexSet<K> connected(K initialKey, std::function<bool(K k)> includeKey);

        std::set<NodeKey> nodes(is_mesh::Geometry *geometry);
    };

    template <typename K>
    inline SimplexSet<K> Query::connected(K initialKey, std::function<bool(K k)> includeKey) {
        AttributeVector<K, char> explored;
        SimplexSet<K> res;
        std::vector<K> q{initialKey};
        while (!q.empty()){
            K qKey = q.back();
            q.pop_back();

            if (explored[qKey]){
                continue;
            }
            explored[qKey] = true;

            if (includeKey(qKey)){
                res.push_back(qKey);

                const auto & simplex = mesh->get(qKey);
                auto boundary = simplex.get_boundary();
                for (auto b : boundary){
                    for (auto bb : mesh->get(b).get_co_boundary()){
                        if (explored[bb]){
                            continue;
                        }
                        q.push_back(bb);
                    }
                }
            }
        }

        return res;
    }

}



