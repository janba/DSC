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

namespace is_mesh
{
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
        
    public:
        
        Simplex()
        {
            m_boundary = new SimplexSet<boundary_key_type>();
            m_co_boundary = new SimplexSet<co_boundary_key_type>();
        }
        
        Simplex(const Simplex& s)
        {
            m_boundary = new SimplexSet<boundary_key_type>(*s.m_boundary);
            m_co_boundary = new SimplexSet<co_boundary_key_type>(*s.m_co_boundary);
        }
        
        Simplex(Simplex&& s)
        {
            m_boundary = s.m_boundary;
            m_co_boundary = s.m_co_boundary;
            s.m_boundary = nullptr;
            s.m_co_boundary = nullptr;
        }

        Simplex<boundary_key_type, co_boundary_key_type>& operator=(Simplex<boundary_key_type, co_boundary_key_type>&& other){
            if (this != &other){
                std::swap(m_boundary, other.m_boundary);
                std::swap(m_co_boundary, other.m_co_boundary);
            }
            return *this;
        }
        
        ~Simplex()
        {
            if(m_boundary)
            {
                delete m_boundary;
            }
            if(m_co_boundary)
            {
                delete m_co_boundary;
            }
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
    template<typename NodeTraits>
    class Node : public NodeTraits, public Simplex<Key, EdgeKey>
    {
    public:
        typedef NodeTraits  type_traits;
        
        Node() : Simplex<Key, EdgeKey>()
        {
            
        }
        Node(const type_traits & t) : type_traits(t), Simplex<Key, EdgeKey>()
        {
            
        }

        Node(Node&& other)
        :NodeTraits(std::move(other)), Simplex<Key, EdgeKey>(std::move(other))
        {}

        Node& operator=(Node&& other){
            if (this != &other){
                ((NodeTraits*)this)->operator=(std::move(other));
                ((Simplex<Key, EdgeKey>*)this)->operator=(std::move(other));
            }
            return *this;
        }

        const SimplexSet<EdgeKey>& edge_keys(){
            return get_co_boundary();
        }

    };
    
    ///////////////////////////////////////////////////////////////////////////////
    ///  E D G E
    ///////////////////////////////////////////////////////////////////////////////
    template<typename EdgeTraits>
    class Edge : public EdgeTraits, public Simplex<NodeKey, FaceKey>
    {
    public:
        typedef EdgeTraits type_traits;
        
        Edge() : Simplex<NodeKey, FaceKey>()
        {
            
        }
        Edge(const type_traits & t) : type_traits(t), Simplex<NodeKey, FaceKey>()
        {
            
        }

        Edge(Edge&& other)
        :EdgeTraits(std::move(other)), Simplex<NodeKey, FaceKey>(std::move(other))
        {

        }

        Edge& operator=(Edge&& other){
            if (this != &other){
                ((EdgeTraits*)this)->operator=(std::move(other));
                ((Simplex<NodeKey, FaceKey>*)this)->operator=(std::move(other));
            }
            return *this;
        }

        const SimplexSet<NodeKey>& node_keys(){
            return get_boundary();
        }

        const SimplexSet<FaceKey>& face_keys(){
            return get_co_boundary();
        }
    };
    
    ///////////////////////////////////////////////////////////////////////////////
    //  F A C E
    ///////////////////////////////////////////////////////////////////////////////
    template<typename FaceTraits>
    class Face : public FaceTraits, public Simplex<EdgeKey, TetrahedronKey>
    {
    public:
        typedef FaceTraits type_traits;
        
        Face() : Simplex<EdgeKey, TetrahedronKey>()
        {
            
        }
        Face(const type_traits & t) : type_traits(t), Simplex<EdgeKey, TetrahedronKey>()
        {
            
        }

        Face(Face&& other)
        : FaceTraits(std::move(other)), Simplex<EdgeKey, TetrahedronKey>(std::move(other))
        {}

        Face& operator=(Face&& other){
            if (this != &other){
                ((FaceTraits*)this)->operator=(std::move(other));
                ((Simplex<EdgeKey, TetrahedronKey>*)this)->operator=(std::move(other));
            }
            return *this;
        }

        const SimplexSet<EdgeKey>& edge_keys(){
            return get_boundary();
        }

        const SimplexSet<TetrahedronKey>& tet_keys(){
            return get_co_boundary();
        }
    };
    
    ///////////////////////////////////////////////////////////////////////////////
    // T E T R A H E D R O N
    ///////////////////////////////////////////////////////////////////////////////
    template<typename TetrahedronTraits>
    class Tetrahedron : public TetrahedronTraits, public Simplex<FaceKey, Key>
    {
    public:
        typedef TetrahedronTraits  type_traits;
        
        Tetrahedron() : Simplex<FaceKey, Key>()
        {
            
        }
        Tetrahedron(const type_traits & t) : type_traits(t), Simplex<FaceKey, Key>()
        {
            
        }

        Tetrahedron(Tetrahedron&& other)
        :TetrahedronTraits(std::move(other)), Simplex<FaceKey, Key>(std::move(other))
        {}

        Tetrahedron& operator=(Tetrahedron&& other){
            if (this != &other){
                ((TetrahedronTraits*)this)->operator=(std::move(other));
                ((Simplex<FaceKey, Key>*)this)->operator=(std::move(other));
            }
            return *this;
        }

        const SimplexSet<FaceKey>& face_keys(){
            return get_boundary();
        }
    };
}
