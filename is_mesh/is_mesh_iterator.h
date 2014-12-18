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
        const kernel<NodeKey, Node> *m_node_kernel;
    public:
        NodeIterator(const kernel<NodeKey, Node> *m_node_kernel);

        typename kernel<NodeKey, Node>::iterator begin() const;

        typename kernel<NodeKey, Node>::iterator end() const;
    };

    class EdgeIterator {
        const kernel<EdgeKey, Edge>*                  m_edge_kernel;
    public:
        EdgeIterator(const kernel<EdgeKey, Edge> *m_edge_kernel);

        typename kernel<EdgeKey, Edge>::iterator begin() const;

        typename kernel<EdgeKey, Edge>::iterator end() const;
    };

    class FaceIterator {
        const kernel<FaceKey, Face>*                  m_face_kernel;
    public:
        FaceIterator(const kernel<FaceKey, Face> *m_face_kernel);

        typename kernel<FaceKey, Face>::iterator begin() const;

        typename kernel<FaceKey, Face>::iterator end() const;

    };

    class TetrahedronIterator {
        const kernel<TetrahedronKey, Tetrahedron>*           m_tetrahedron_kernel;
    public:
        TetrahedronIterator(const kernel<TetrahedronKey, Tetrahedron> *m_tetrahedron_kernel);

        typename kernel<TetrahedronKey, Tetrahedron>::iterator begin() const;

        typename kernel<TetrahedronKey, Tetrahedron>::iterator end() const;

    };
}