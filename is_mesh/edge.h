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
    class Edge : public Simplex<NodeKey, FaceKey>
    {
        std::bitset<3> flags;
    public:
        Edge(ISMesh *owner);

        Edge(Edge&& other) noexcept;

        Edge& operator=(Edge&& other) noexcept;

        const SimplexSet<NodeKey> & node_keys()  const noexcept;

        const SimplexSet<FaceKey> & face_keys() const noexcept;

        bool is_crossing() noexcept;

        bool is_boundary() noexcept;

        bool is_interface() noexcept;

        double length();

        double length_destination();

        EdgeKey key();
    private:
        void set_crossing(bool b);

        void set_boundary(bool b);

        void set_interface(bool b);

        friend class ISMesh;
    };
}