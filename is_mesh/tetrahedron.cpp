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

namespace is_mesh
{
    Tetrahedron::Tetrahedron(ISMesh *owner) : Simplex<FaceKey, Key>(owner) {

    }

    Tetrahedron::Tetrahedron(Tetrahedron&& other)
            : l(other.l), Simplex<FaceKey, Key>(std::move(other)) {}

    Tetrahedron &Tetrahedron::operator=(Tetrahedron&& other) {
        if (this != &other){
            l = other.l;
            ((Simplex<FaceKey, Key>*)this)->operator=(std::move(other));
        }
        return *this;
    }

    const SimplexSet<FaceKey> &Tetrahedron::face_keys()  const {
        return get_boundary();
    }

    int Tetrahedron::label() {
        return l;
    }

    void Tetrahedron::label(unsigned int _label) {
        l = _label;
    }
}