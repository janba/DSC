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

#include "simplex.h"

namespace is_mesh
{
    class Tetrahedron;
    class Edge;
    class Node;

    class Face : public Simplex<EdgeKey, TetrahedronKey>
    {
        bool interface = false;
    public:
        Face(ISMesh *owner) noexcept;

        Face(Face&& other) noexcept;

        Face& operator=(Face&& other) noexcept;

        const SimplexSet<EdgeKey> & edge_keys() const noexcept;
        const SimplexSet<TetrahedronKey> & tet_keys() const noexcept;
        const SimplexSet<NodeKey> node_keys() const noexcept;

        std::vector<Tetrahedron*> tets() const;
        std::vector<Edge*> edges() const;
        std::vector<Node*> nodes() const;

        vec3 get_center() const;

        bool is_boundary() noexcept;

        bool is_interface() noexcept;

        double area() const;

        double area_destination() const;

        double min_angle() const;

        double max_angle() const;

        double quality() const;

        EdgeKey longest_edge() const;

        vec3 barycenter() const;

        FaceKey key() const noexcept;
    private:
        void set_interface(bool b);

        friend class ISMesh;
    };
}