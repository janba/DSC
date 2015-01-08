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
#include "node.h"
#include "edge.h"
#include "tetrahedron.h"
#include "face.h"

namespace is_mesh {
    class NodeIterator {
        const kernel<NodeKey, Node> *m_node_kernel;
        const std::vector<NodeKey> *subset;
#ifdef DEBUG
        kernel<NodeKey, Node>::iterator first_end_iter;
#endif
    public:
        NodeIterator(const kernel<NodeKey, Node> *m_node_kernel, const std::vector<NodeKey> *subset);
        ~NodeIterator();

        typename kernel<NodeKey, Node>::iterator begin() const;

        typename kernel<NodeKey, Node>::iterator end() const;
    };

    class EdgeIterator {
        const kernel<EdgeKey, Edge>* m_edge_kernel;
        const std::vector<EdgeKey> *subset;
#ifdef DEBUG
        kernel<EdgeKey, Edge>::iterator first_end_iter;
#endif
    public:
        EdgeIterator(const kernel<EdgeKey, Edge> *m_edge_kernel, const std::vector<EdgeKey> *subset);
        ~EdgeIterator();

        typename kernel<EdgeKey, Edge>::iterator begin() const;

        typename kernel<EdgeKey, Edge>::iterator end() const;
    };

    class FaceIterator {
        const kernel<FaceKey, Face>* m_face_kernel;
        const std::vector<FaceKey> *subset;
#ifdef DEBUG
        kernel<FaceKey, Face>::iterator first_end_iter;
#endif
    public:
        FaceIterator(const kernel<FaceKey, Face> *m_face_kernel, const std::vector<FaceKey> *subset);
        ~FaceIterator();

        typename kernel<FaceKey, Face>::iterator begin() const;

        typename kernel<FaceKey, Face>::iterator end() const;

    };

    class TetrahedronIterator {
        const kernel<TetrahedronKey, Tetrahedron>* m_tetrahedron_kernel;
        const std::vector<TetrahedronKey> *subset;
#ifdef DEBUG
        kernel<TetrahedronKey, Tetrahedron>::iterator first_end_iter;
#endif
    public:
        TetrahedronIterator(const kernel<TetrahedronKey, Tetrahedron> *m_tetrahedron_kernel, const std::vector<TetrahedronKey> *subset);
        ~TetrahedronIterator();

        typename kernel<TetrahedronKey, Tetrahedron>::iterator begin() const;

        typename kernel<TetrahedronKey, Tetrahedron>::iterator end() const;

    };
}
