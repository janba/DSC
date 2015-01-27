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
    class Node : public Simplex<Key, EdgeKey>
    {
        vec3 p;
        vec3 p_new;
        std::bitset<3> flags;
    public:
        Node(ISMesh *owner) noexcept;
        Node(ISMesh *owner, vec3 _p) noexcept;

        Node(Node&& other) noexcept;

        Node& operator=(Node&& other) noexcept;

        const SimplexSet<EdgeKey> & edge_keys() const noexcept;

        const vec3 & get_pos() const;

        const vec3 & get_destination() const;

        void set_pos(const vec3& p_);

        void set_destination(const vec3& p_);

        bool is_crossing() const noexcept;

        bool is_boundary() const noexcept;

        bool is_interface() const noexcept;

        NodeKey key();
    private:
        void set_crossing(bool b);

        void set_boundary(bool b);

        void set_interface(bool b);

        friend class ISMesh;
    };
}