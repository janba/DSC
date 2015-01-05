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

    NodeIterator::NodeIterator(const kernel<NodeKey,Node> *m_node_kernel, const std::vector<NodeKey> *subset)
            : m_node_kernel(m_node_kernel),subset(subset) {
    }

    typename kernel<NodeKey,Node>::iterator NodeIterator::begin() const {
        if (subset){
            return kernel<NodeKey,Node>::iterator(m_node_kernel, subset);
        }
        return m_node_kernel->begin();
    }

    typename kernel<NodeKey,Node>::iterator NodeIterator::end() const {
        return m_node_kernel->end();
    }

    EdgeIterator::EdgeIterator(const kernel<EdgeKey,Edge> *m_edge_kernel, const std::vector<EdgeKey> *subset)
            : m_edge_kernel(m_edge_kernel),subset(subset) {
    }

    typename kernel<EdgeKey,Edge>::iterator EdgeIterator::begin() const {
        if (subset){
            return kernel<EdgeKey,Edge>::iterator(m_edge_kernel, subset);
        }
        return m_edge_kernel->begin();
    }

    typename kernel<EdgeKey,Edge>::iterator EdgeIterator::end() const {
        return m_edge_kernel->end();
    }

    FaceIterator::FaceIterator(const kernel<FaceKey,Face> *m_face_kernel,const std::vector<FaceKey> *subset)
            : m_face_kernel(m_face_kernel),subset(subset) {
    }

    typename kernel<FaceKey,Face>::iterator FaceIterator::begin() const {
        if (subset){
            return kernel<FaceKey,Face>::iterator(m_face_kernel, subset);
        }
        return m_face_kernel->begin();
    }

    typename kernel<FaceKey,Face>::iterator FaceIterator::end() const {
        return m_face_kernel->end();
    }

    TetrahedronIterator::TetrahedronIterator(const kernel<TetrahedronKey,Tetrahedron> *m_tetrahedron_kernel,const std::vector<TetrahedronKey> *subset)
            : m_tetrahedron_kernel(m_tetrahedron_kernel),subset(subset) {
    }

    typename kernel<TetrahedronKey,Tetrahedron>::iterator TetrahedronIterator::begin() const {
        if (subset){
            return kernel<TetrahedronKey,Tetrahedron>::iterator(m_tetrahedron_kernel, subset);
        }
        return m_tetrahedron_kernel->begin();
    }

    typename kernel<TetrahedronKey,Tetrahedron>::iterator TetrahedronIterator::end() const {
        return m_tetrahedron_kernel->end();
    }
}