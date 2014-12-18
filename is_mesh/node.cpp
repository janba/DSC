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

namespace is_mesh
{
    Node::Node(ISMesh *owner) : Simplex<Key, EdgeKey>(owner) {

    }

    Node::Node(ISMesh *owner,vec3 _p) : p(_p), p_new(_p), Simplex<Key, EdgeKey>(owner) {

    }


    Node::Node(Node&& other) : p (other.p), p_new(other.p_new),
                               flags(std::move(other.flags)),Simplex<Key, EdgeKey>(std::move(other)) {
    }

    Node &Node::operator=(Node&& other) {
        if (this != &other){
            p = std::move(other.p);
            p_new = std::move(other.p_new);
            flags = std::move(other.flags);
            ((Simplex<Key, EdgeKey>*)this)->operator=(std::move(other));
        }
        return *this;
    }

    const SimplexSet<EdgeKey> &Node::edge_keys() const {
        return get_co_boundary();
    }

    const vec3 &Node::get_pos() const {
        return p;
    }

    const vec3 &Node::get_destination() const {
        return p_new;
    }

    void Node::set_pos(const vec3& p_) {
        p = p_;
    }

    void Node::set_destination(const vec3& p_) {
        p_new = p_;
    }

    bool Node::is_crossing() const {
        return flags[2];
    }

    bool Node::is_boundary() const {
        return flags[1];
    }

    bool Node::is_interface() const {
        return flags[0];
    }

    void Node::set_crossing(bool b) {
        flags[2] = b;
    }

    void Node::set_boundary(bool b) {
        flags[1] = b;
    }

    void Node::set_interface(bool b) {
        flags[0] = b;
    }
}