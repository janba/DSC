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

#include "simplex_set.h"
#include "util.h"

namespace is_mesh
{
    class ISMesh;

    ///////////////////////////////////////////////////////////////////////////////
    // S I M P L E X   B A S E   C L A S S
    ///////////////////////////////////////////////////////////////////////////////
    /**
     * Base class for all simplex classes
     */
    template<typename boundary_key_type, typename co_boundary_key_type>
    class Simplex
    {
        SimplexSet<boundary_key_type>* m_boundary = nullptr;
        SimplexSet<co_boundary_key_type>* m_co_boundary = nullptr;
    protected:
        ISMesh *m_mesh = nullptr;
    public:

        Simplex(ISMesh * owner) : m_mesh {owner} {
            m_boundary = new SimplexSet<boundary_key_type>();
            m_co_boundary = new SimplexSet<co_boundary_key_type>();
        }
        
        Simplex(const Simplex& s):m_mesh {s.m_mesh }
        {
            m_boundary = new SimplexSet<boundary_key_type>(*s.m_boundary);
            m_co_boundary = new SimplexSet<co_boundary_key_type>(*s.m_co_boundary);
        }
        
        Simplex(Simplex&& s) noexcept
        {
            std::swap(m_boundary, s.m_boundary);
            std::swap(m_co_boundary, s.m_co_boundary);
            std::swap(m_mesh, s.m_mesh );
        }

        Simplex<boundary_key_type, co_boundary_key_type>& operator=(Simplex<boundary_key_type, co_boundary_key_type>&& other) noexcept
        {
            if (this != &other){
                std::swap(m_boundary, other.m_boundary);
                std::swap(m_co_boundary, other.m_co_boundary);
                std::swap(m_mesh , other.m_mesh );
            }
            return *this;
        }
        
        ~Simplex()
        {
            delete m_boundary;
            delete m_co_boundary;
        }
        
    public:
        
        const SimplexSet<co_boundary_key_type>& get_co_boundary() const noexcept
        {
            return *m_co_boundary;
        }
        const SimplexSet<boundary_key_type>& get_boundary() const noexcept
        {
            return *m_boundary;
        }
        
        void add_co_face(const co_boundary_key_type& key)
        {
            *m_co_boundary += key;
        }
        
        void add_face(const boundary_key_type& key)
        {
            *m_boundary += key;
        }
        
        void remove_co_face(const co_boundary_key_type& key)
        {
            *m_co_boundary -= key;
        }
        
        void remove_face(const boundary_key_type& key)
        {
            *m_boundary -= key;
        }
    };
}
