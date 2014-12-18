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

#include "edge.h"

namespace is_mesh
{
    Edge::Edge(ISMesh *owner) : Simplex<NodeKey, FaceKey>(owner) {

    }

    Edge::Edge(Edge&& other) : flags(std::move(other.flags)), Simplex<NodeKey, FaceKey>(std::move(other)) {

    }

    Edge &Edge::operator=(Edge&& other) {
        if (this != &other){
            flags = std::move(other.flags);
            ((Simplex<NodeKey, FaceKey>*)this)->operator=(std::move(other));
        }
        return *this;
    }

    const SimplexSet<NodeKey> &Edge::node_keys()  const {
        return get_boundary();
    }

    const SimplexSet<FaceKey> &Edge::face_keys()  const {
        return get_co_boundary();
    }

    bool Edge::is_crossing() {
        return flags[2];
    }

    bool Edge::is_boundary() {
        return flags[1];
    }

    bool Edge::is_interface() {
        return flags[0];
    }

    void Edge::set_crossing(bool b) {
        flags[2] = b;
    }

    void Edge::set_boundary(bool b) {
        flags[1] = b;
    }

    void Edge::set_interface(bool b) {
        flags[0] = b;
    }
}