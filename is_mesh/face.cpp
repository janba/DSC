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


#include "face.h"

namespace is_mesh
{
    Face::Face(ISMesh *owner) : Simplex<EdgeKey, TetrahedronKey>(owner) {

    }

    Face::Face(Face&& other) : flags(std::move(other.flags)), Simplex<EdgeKey, TetrahedronKey>(std::move(other)) {}

    Face &Face::operator=(Face&& other) {
        if (this != &other){
            flags = std::move(other.flags);
            ((Simplex<EdgeKey, TetrahedronKey>*)this)->operator=(std::move(other));
        }
        return *this;
    }

    const SimplexSet<EdgeKey> &Face::edge_keys()  const {
        return get_boundary();
    }

    const SimplexSet<TetrahedronKey> &Face::tet_keys()  const {
        return get_co_boundary();
    }

    bool Face::is_boundary() {
        return flags[1];
    }

    bool Face::is_interface() {
        return flags[0];
    }

    void Face::set_boundary(bool b) {
        flags[1] = b;
    }

    void Face::set_interface(bool b) {
        flags[0] = b;
    }
}