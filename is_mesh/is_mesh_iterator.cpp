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

#include "is_mesh_iterator.h"

namespace is_mesh {

    NodeIterator::NodeIterator(const kernel<Node, NodeKey> *m_node_kernel) : m_node_kernel(m_node_kernel) {
    }

    typename kernel<Node, NodeKey>::iterator NodeIterator::begin() const {
        return m_node_kernel->begin();
    }

    typename kernel<Node, NodeKey>::iterator NodeIterator::end() const {
        return m_node_kernel->end();
    }

    EdgeIterator::EdgeIterator(const kernel<Edge, EdgeKey> *m_edge_kernel) : m_edge_kernel(m_edge_kernel) {
    }

    typename kernel<Edge, EdgeKey>::iterator EdgeIterator::begin() const {
        return m_edge_kernel->begin();
    }

    typename kernel<Edge, EdgeKey>::iterator EdgeIterator::end() const {
        return m_edge_kernel->end();
    }

    FaceIterator::FaceIterator(const kernel<Face, FaceKey> *m_face_kernel) : m_face_kernel(m_face_kernel) {
    }

    typename kernel<Face, FaceKey>::iterator FaceIterator::begin() const {
        return m_face_kernel->begin();
    }

    typename kernel<Face, FaceKey>::iterator FaceIterator::end() const {
        return m_face_kernel->end();
    }

    TetrahedronIterator::TetrahedronIterator(const kernel<Tetrahedron, TetrahedronKey> *m_tetrahedron_kernel)
            : m_tetrahedron_kernel(m_tetrahedron_kernel) {
    }

    typename kernel<Tetrahedron, TetrahedronKey>::iterator TetrahedronIterator::begin() const {
        return m_tetrahedron_kernel->begin();
    }

    typename kernel<Tetrahedron, TetrahedronKey>::iterator TetrahedronIterator::end() const {
        return m_tetrahedron_kernel->end();
    }
}