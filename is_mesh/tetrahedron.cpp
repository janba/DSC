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


#include "tetrahedron.h"
#include "is_mesh.h"

namespace is_mesh
{
    Tetrahedron::Tetrahedron(ISMesh *owner) noexcept
            : Simplex<FaceKey, Key>(owner)
    {
    }

    Tetrahedron::Tetrahedron(Tetrahedron&& other) noexcept
            : l(other.l), Simplex<FaceKey, Key>(std::move(other)) {}

    Tetrahedron &Tetrahedron::operator=(Tetrahedron&& other) noexcept {
        if (this != &other){
            l = other.l;
            ((Simplex<FaceKey, Key>*)this)->operator=(std::move(other));
        }
        return *this;
    }

    const SimplexSet<FaceKey> &Tetrahedron::face_keys() const noexcept{
        return get_boundary();
    }

    int Tetrahedron::label() {
        return l;
    }

    void Tetrahedron::label(unsigned int _label) {
        l = _label;
    }

    double Tetrahedron::volume() {
        SimplexSet<NodeKey> nids = node_keys();
        return Util::volume(m_mesh->get(nids[0]).get_pos(), m_mesh->get(nids[1]).get_pos(), m_mesh->get(nids[2]).get_pos(), m_mesh->get(nids[3]).get_pos());
    }

    const SimplexSet<NodeKey> Tetrahedron::node_keys() const noexcept{
        const SimplexSet<FaceKey>& fids = face_keys();
        SimplexSet<NodeKey> nids = m_mesh->get_nodes(fids[0]);
        nids += m_mesh->get_nodes(fids[1]);
        return nids;

    }

    double Tetrahedron::volume_destination() {
        SimplexSet<NodeKey> nids = node_keys();
        return Util::volume(m_mesh->get(nids[0]).get_destination(), m_mesh->get(nids[1]).get_destination(), m_mesh->get(nids[2]).get_destination(), m_mesh->get(nids[3]).get_destination());
    }

    vec3 Tetrahedron::barycenter() {
        SimplexSet<NodeKey> nids = node_keys();
        return Util::barycenter(m_mesh->get(nids[0]).get_pos(), m_mesh->get(nids[1]).get_pos(), m_mesh->get(nids[2]).get_pos(), m_mesh->get(nids[3]).get_pos());
    }

    vec3 Tetrahedron::barycenter_destination() {
        SimplexSet<NodeKey> nids = node_keys();
        return Util::barycenter(m_mesh->get(nids[0]).get_destination(), m_mesh->get(nids[1]).get_destination(), m_mesh->get(nids[2]).get_destination(), m_mesh->get(nids[3]).get_destination());
    }

    double Tetrahedron::quality() {
        SimplexSet<NodeKey> nids = node_keys();
        return fabs(Util::quality(m_mesh->get(nids[0]).get_pos(), m_mesh->get(nids[1]).get_pos(), m_mesh->get(nids[2]).get_pos(), m_mesh->get(nids[3]).get_pos()));
    }

    TetrahedronKey Tetrahedron::key() const noexcept {
        long index = ((char*)this - m_mesh->m_tetrahedron_kernel.data())/sizeof(util::kernel_element<TetrahedronKey, Tetrahedron>);
        assert(index >= 0);
        assert(index < m_mesh->m_tetrahedron_kernel.capacity());
        return TetrahedronKey((unsigned int) index);
    }
}