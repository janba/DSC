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

#pragma once

#include "kernel.h"
#include "key.h"
#include "simplex.h"
#include "simplex_set.h"

namespace is_mesh {
    class NodeIterator {
        const kernel<Node, NodeKey> *m_node_kernel;
    public:
        NodeIterator(const kernel<Node, NodeKey> *m_node_kernel);

        typename kernel<Node, NodeKey>::iterator begin() const;

        typename kernel<Node, NodeKey>::iterator end() const;
    };

    class EdgeIterator {
        const kernel<Edge, EdgeKey>*                  m_edge_kernel;
    public:
        EdgeIterator(const kernel<Edge, EdgeKey> *m_edge_kernel);

        typename kernel<Edge, EdgeKey>::iterator begin() const;

        typename kernel<Edge, EdgeKey>::iterator end() const;
    };

    class FaceIterator {
        const kernel<Face, FaceKey>*                  m_face_kernel;
    public:
        FaceIterator(const kernel<Face, FaceKey> *m_face_kernel);

        typename kernel<Face, FaceKey>::iterator begin() const;

        typename kernel<Face, FaceKey>::iterator end() const;

    };

    class TetrahedronIterator {
        const kernel<Tetrahedron, TetrahedronKey>*           m_tetrahedron_kernel;
    public:
        TetrahedronIterator(const kernel<Tetrahedron, TetrahedronKey> *m_tetrahedron_kernel);

        typename kernel<Tetrahedron, TetrahedronKey>::iterator begin() const;

        typename kernel<Tetrahedron, TetrahedronKey>::iterator end() const;

    };
}