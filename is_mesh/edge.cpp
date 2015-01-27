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
#include "is_mesh.h"

namespace is_mesh
{
    Edge::Edge(ISMesh *owner) noexcept : Simplex<NodeKey, FaceKey>(owner) {

    }

    Edge::Edge(Edge&& other) noexcept
            : flags(std::move(other.flags)), Simplex<NodeKey, FaceKey>(std::move(other))
    {
    }

    Edge &Edge::operator=(Edge&& other) noexcept {
        if (this != &other){
            std::swap(flags, other.flags);
            ((Simplex<NodeKey, FaceKey>*)this)->operator=(std::move(other));
        }
        return *this;
    }

    const SimplexSet<NodeKey> &Edge::node_keys()  const noexcept {
        return get_boundary();
    }

    const SimplexSet<FaceKey> &Edge::face_keys()  const noexcept {
        return get_co_boundary();
    }

    bool Edge::is_crossing() noexcept {
        return flags[2];
    }

    bool Edge::is_boundary() noexcept {
        return flags[1];
    }

    bool Edge::is_interface() noexcept {
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

    double Edge::length() {
        const SimplexSet<NodeKey> & nids = node_keys();
        return Util::length(m_mesh->get(nids[0]).get_pos() - m_mesh->get(nids[1]).get_pos());
    }

    double Edge::length_destination() {
        const SimplexSet<NodeKey> & nids = node_keys();
        return Util::length(m_mesh->get(nids[0]).get_destination() - m_mesh->get(nids[1]).get_destination());
    }

    EdgeKey Edge::key() const noexcept {
        long index = ((char*)this - m_mesh->m_edge_kernel->data())/sizeof(util::kernel_element<EdgeKey, Edge>);
        assert(index >= 0);
        assert(index < m_mesh->m_edge_kernel->capacity());
        return EdgeKey((unsigned int) index);
    }
}