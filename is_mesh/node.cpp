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


#include "node.h"
#include "is_mesh.h"

namespace is_mesh
{
    Node::Node(ISMesh *owner) noexcept
            : Simplex<Key, EdgeKey>(owner)
    {
    }

    Node::Node(ISMesh *owner,vec3 _p) noexcept
            : p(_p), p_new(_p), Simplex<Key, EdgeKey>(owner)
    {
    }


    Node::Node(Node&& other) noexcept
            : p (other.p), p_new(other.p_new),
              crossing(other.crossing),
              boundary(other.boundary),
              interface(other.interface),
              Simplex<Key, EdgeKey>(std::move(other))
    {
    }

    Node &Node::operator=(Node&& other) noexcept {
        if (this != &other){
            p = other.p;
            p_new = other.p_new;
            crossing = other.crossing;
            boundary = other.boundary;
            interface = other.interface;
            ((Simplex<Key, EdgeKey>*)this)->operator=(std::move(other));
        }
        return *this;
    }

    const SimplexSet<EdgeKey> &Node::edge_keys() const noexcept{
        return get_co_boundary();
    }

    const vec3 &Node::get_pos() const {
        return p;
    }

    const vec3 &Node::get_destination() const {
        return p_new;
    }

    void Node::set_pos(const vec3& p_) {
#ifdef DEBUG
        if (p_ != p){
            assert(!is_boundary());
        }
#endif
        p = p_;
    }

    void Node::set_destination(const vec3& p_) {
        if (is_boundary()){
            assert(p_ == p);
        }
        p_new = p_;
    }

    bool Node::is_crossing() const noexcept {
        return crossing;
    }

    bool Node::is_boundary() const noexcept {
        return boundary;
    }

    bool Node::is_interface() const noexcept {
        return interface;
    }

    void Node::set_crossing(bool b) {
        crossing = b;
    }

    void Node::set_boundary(bool b) {
        boundary = b;
    }

    void Node::set_interface(bool b) {
        interface = b;
    }

    int Node::get_number_of_neighbour_tet_clusters(int label){
        auto tets = tet_keys();
        int clusters = 0;

        while (tets.size()>0){
            auto tetKey = tets.pop_back();
            if (m_mesh->get(tetKey).label() == label){
                SimplexSet<TetrahedronKey> boundary{{tetKey}};
                do {
                    tetKey = boundary.pop_back();
                    if (m_mesh->get(tetKey).label() == label){
                        auto faces = m_mesh->get(tetKey).face_keys();
                        auto neighbourTets = m_mesh->get_tets(faces) - tetKey;
                        for (auto nTet : neighbourTets){
                            if (tets.contains(nTet)){
                                tets -= nTet;
                                boundary.push_back(nTet);
                            }
                        }
                    }
                } while (boundary.size()>0);
                clusters++;
            }
        }
        return clusters;
    }

    vec3 Node::smart_laplacian(double alpha) const {
        auto k = key();
        SimplexSet<TetrahedronKey> tids = tet_keys();
        SimplexSet<FaceKey> fids = m_mesh->get_faces(tids) - face_keys();

        vec3 avg_pos = m_mesh->get_barycenter(m_mesh->get_nodes(fids));
        return p + alpha * (avg_pos - p);
    }

    NodeKey Node::key() const noexcept {
        long index = ((char*)this - m_mesh->m_node_kernel.data())/sizeof(util::kernel_element<NodeKey, Node> );
        assert(index >= 0);
        assert(index < m_mesh->m_node_kernel.capacity());
        return NodeKey((unsigned int) index);
    }

    const vec3 Node::get_center() const {
        return p;
    }

    std::vector<Edge *> Node::edges() const {
        std::vector<Edge *> res;
        for (auto k : edge_keys()){
            res.push_back(&m_mesh->get(k));
        }
        return res;
    }

    std::vector<Face *> Node::faces() const {
        std::vector<Face *> res;
        for (auto key : face_keys()){
            res.push_back(&m_mesh->get(key));
        }
        return res;
    }

    std::vector<Tetrahedron *> Node::tets() const {
        std::vector<Tetrahedron *> res;
        for (auto key : tet_keys()){
            res.push_back(&m_mesh->get(key));
        }
        return res;
    }

    SimplexSet<FaceKey> Node::face_keys() const {
        SimplexSet<FaceKey> resKey;
        for (auto edge : edges()){
            resKey += edge->face_keys();
        }
        return resKey;
    }

    SimplexSet<TetrahedronKey> Node::tet_keys() const {
        SimplexSet<TetrahedronKey> resKey;
        for (auto edge : faces()){
            resKey += edge->tet_keys();
        }
        return resKey;
    }
}