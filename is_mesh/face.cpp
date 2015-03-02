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


#include "face.h"
#include "is_mesh.h"

#include <algorithm>
#include <utility>
#include <cmath>

namespace is_mesh
{
    Face::Face(ISMesh *owner) noexcept : Simplex<EdgeKey, TetrahedronKey>(owner) {

    }

    Face::Face(Face&& other) noexcept
            : interface(other.interface), Simplex<EdgeKey, TetrahedronKey>(std::move(other)) {}

    Face &Face::operator=(Face&& other) noexcept {
        if (this != &other){
            std::swap(interface, other.interface);
            ((Simplex<EdgeKey, TetrahedronKey>*)this)->operator=(std::move(other));
        }
        return *this;
    }

    const SimplexSet<EdgeKey> &Face::edge_keys()  const noexcept{
        return get_boundary();
    }

    const SimplexSet<TetrahedronKey> &Face::tet_keys()  const noexcept{
        return get_co_boundary();
    }

    bool Face::is_boundary() noexcept {
        return tet_keys().size() < 2;
    }

    bool Face::is_interface() noexcept {
        return interface;
    }

    void Face::set_interface(bool b) {
        interface = b;
    }

    double Face::area() const {
        const SimplexSet<NodeKey>& nids = node_keys();
        return Util::area(m_mesh->get(nids[0]).get_pos(), m_mesh->get(nids[1]).get_pos(), m_mesh->get(nids[2]).get_pos());
    }

    const SimplexSet<NodeKey> Face::node_keys() const noexcept{
        const SimplexSet<EdgeKey>& eids = edge_keys();
        SimplexSet<NodeKey> nids = m_mesh->get_nodes(eids[0]);
        nids += m_mesh->get_nodes(eids[1]);
        return nids;
    }

    double Face::area_destination() const {
        const SimplexSet<NodeKey>& nids = node_keys();
        return Util::area(m_mesh->get(nids[0]).get_destination(), m_mesh->get(nids[1]).get_destination(), m_mesh->get(nids[2]).get_destination());
    }

    double Face::min_angle() const {
        SimplexSet<NodeKey> nids = node_keys();
        return Util::min_angle(m_mesh->get(nids[0]).get_pos(), m_mesh->get(nids[1]).get_pos(), m_mesh->get(nids[2]).get_pos());
    }

    double Face::max_angle() const {
        SimplexSet<NodeKey> nids = node_keys();
        return Util::max_angle(m_mesh->get(nids[0]).get_pos(), m_mesh->get(nids[1]).get_pos(), m_mesh->get(nids[2]).get_pos());
    }

    double Face::quality() const {
        SimplexSet<NodeKey> nids = node_keys();
        auto angles = Util::cos_angles(m_mesh->get(nids[0]).get_pos(), m_mesh->get(nids[1]).get_pos(), m_mesh->get(nids[2]).get_pos());
        double worst_a = -INFINITY;
        for(double a : angles)
        {
            worst_a = std::max(worst_a, fabs(a));
        }
        return 1. - worst_a;
    }


    EdgeKey Face::longest_edge() const {
        double max_l = -INFINITY;
        EdgeKey max_e;
        for(auto e : edge_keys())
        {
            double l = m_mesh->get(e).length();
            if(l > max_l)
            {
                max_l = l;
                max_e = e;
            }
        }
        return max_e;
    }

    vec3 Face::barycenter() const {
        auto nids = node_keys();
        return Util::barycenter(m_mesh->get(nids[0]).get_pos(), m_mesh->get(nids[1]).get_pos(), m_mesh->get(nids[2]).get_pos());
    }

    FaceKey Face::key()  const noexcept {
        long index = ((char*)this - m_mesh->m_face_kernel.data())/sizeof(util::kernel_element<FaceKey, Face>);
        assert(index >= 0);
        assert(index < m_mesh->m_face_kernel.capacity());
        return FaceKey((unsigned int) index);
    }

    vec3 Face::get_center() const {
        return barycenter();
    }
}