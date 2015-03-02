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
    class Tetrahedron : public Simplex<FaceKey, Key>
    {
        unsigned int l = 0;

        // modify through ISMesh::set_label
        void label(unsigned int _label);
    public:
        Tetrahedron(ISMesh *owner) noexcept;
        Tetrahedron(ISMesh *owner, int l) noexcept;

        Tetrahedron(Tetrahedron&& other) noexcept;

        Tetrahedron& operator=(Tetrahedron&& other) noexcept;

        const SimplexSet<FaceKey> & face_keys() const noexcept;

        const SimplexSet<NodeKey> node_keys() const noexcept;

        vec3 get_center() const;

        int label() const;

        double volume() const;

        double volume_destination() const;

        vec3 barycenter() const;
        vec3 barycenter_destination() const;

        double quality() const;

        TetrahedronKey key() const noexcept;

        friend class ISMesh;
    };
}