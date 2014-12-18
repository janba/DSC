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
        ISMesh *owner = nullptr;
    public:
        
        Simplex(ISMesh *owner):owner{owner}
        {
            m_boundary = new SimplexSet<boundary_key_type>();
            m_co_boundary = new SimplexSet<co_boundary_key_type>();
        }
        
        Simplex(const Simplex& s):owner{s.owner}
        {
            m_boundary = new SimplexSet<boundary_key_type>(*s.m_boundary);
            m_co_boundary = new SimplexSet<co_boundary_key_type>(*s.m_co_boundary);
        }
        
        Simplex(Simplex&& s)
        {
            std::swap(m_boundary, s.m_boundary);
            std::swap(m_co_boundary, s.m_co_boundary);
            std::swap(owner, s.owner);
        }

        Simplex<boundary_key_type, co_boundary_key_type>& operator=(Simplex<boundary_key_type, co_boundary_key_type>&& other){
            if (this != &other){
                std::swap(m_boundary, other.m_boundary);
                std::swap(m_co_boundary, other.m_co_boundary);
                std::swap(owner, other.owner);
            }
            return *this;
        }
        
        ~Simplex()
        {
            delete m_boundary;
            delete m_co_boundary;
        }
        
    public:
        
        const SimplexSet<co_boundary_key_type>& get_co_boundary() const
        {
            return *m_co_boundary;
        }
        const SimplexSet<boundary_key_type>& get_boundary() const
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
    
    ///////////////////////////////////////////////////////////////////////////////
    ///  N O D E
    ///////////////////////////////////////////////////////////////////////////////
    class Node : public Simplex<Key, EdgeKey>
    {
        vec3 p;
        vec3 p_new;
        std::bitset<3> flags;
    public:
        Node(ISMesh *owner);
        Node(ISMesh *owner,vec3 _p);

        Node(Node&& other);

        Node& operator=(Node&& other);

        const SimplexSet<EdgeKey> & edge_keys() const;

        /**
        * Returns the position of the node.
        */
        const vec3& get_pos() const
        {
            return p;
        }

        /**
        * Returns the destination of the node.
        */
        const vec3& get_destination() const
        {
            return p_new;
        }

        void set_pos(const vec3& p_)
        {
            p = p_;
        }

        void set_destination(const vec3& p_)
        {
            p_new = p_;
        }

        bool is_crossing() const
        {
            return flags[2];
        }

        bool is_boundary() const
        {
            return flags[1];
        }

        bool is_interface() const
        {
            return flags[0];
        }
    private:
        void set_crossing(bool b)
        {
            flags[2] = b;
        }

        void set_boundary(bool b)
        {
            flags[1] = b;
        }

        void set_interface(bool b)
        {
            flags[0] = b;
        }
        friend class ISMesh;
    };
    
    ///////////////////////////////////////////////////////////////////////////////
    ///  E D G E
    ///////////////////////////////////////////////////////////////////////////////
    class Edge : public Simplex<NodeKey, FaceKey>
    {
        std::bitset<3> flags;
    public:
        Edge(ISMesh *owner);

        Edge(Edge&& other);

        Edge& operator=(Edge&& other);

        const SimplexSet<NodeKey> & node_keys()  const;

        const SimplexSet<FaceKey> & face_keys() const;

        bool is_crossing()
        {
            return flags[2];
        }

        bool is_boundary()
        {
            return flags[1];
        }

        bool is_interface()
        {
            return flags[0];
        }
    private:
        void set_crossing(bool b)
        {
            flags[2] = b;
        }

        void set_boundary(bool b)
        {
            flags[1] = b;
        }

        void set_interface(bool b)
        {
            flags[0] = b;
        }
        friend class ISMesh;
    };
    
    ///////////////////////////////////////////////////////////////////////////////
    //  F A C E
    ///////////////////////////////////////////////////////////////////////////////
    class Face :  public Simplex<EdgeKey, TetrahedronKey>
    {
        std::bitset<2> flags;
    public:
        Face(ISMesh *owner);

        Face(Face&& other);

        Face& operator=(Face&& other);

        const SimplexSet<EdgeKey> & edge_keys() const;

        const SimplexSet<TetrahedronKey> & tet_keys() const;

        bool is_boundary()
        {
            return flags[1];
        }

        bool is_interface()
        {
            return flags[0];
        }
    private:
        void set_boundary(bool b)
        {
            flags[1] = b;
        }

        void set_interface(bool b)
        {
            flags[0] = b;
        }
        friend class ISMesh;
    };
    
    ///////////////////////////////////////////////////////////////////////////////
    // T E T R A H E D R O N
    ///////////////////////////////////////////////////////////////////////////////
    class Tetrahedron : public Simplex<FaceKey, Key>
    {
        unsigned int l = 0;
    public:
        Tetrahedron(ISMesh *owner);

        Tetrahedron(Tetrahedron&& other);

        Tetrahedron& operator=(Tetrahedron&& other);

        const SimplexSet<FaceKey> & face_keys() const;
        int label()
        {
            return l;
        }

        void label(unsigned int _label)
        {
            l = _label;
        }
        friend class ISMesh;
    };
}
