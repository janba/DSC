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
    public:
        typedef SimplexSet<boundary_key_type>       boundary_list;
        typedef SimplexSet<co_boundary_key_type>    co_boundary_list;
        
    protected:
        boundary_list* m_boundary = nullptr;
        co_boundary_list* m_co_boundary = nullptr;
        int              		m_label;       //used in coloring - to identify connected components.
        
    public:
        
        Simplex() : m_label(0)
        {
            m_boundary    = new boundary_list();
            m_co_boundary = new co_boundary_list();
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
        
        void set_label(const int& l)
        {
            m_label = l;
        }
        
        int get_label()
        {
            return m_label;
        }
        
        void reset_label()
        {
            m_label = 0;
        }
        
        //copy constructor - needed because simplices are stored in STL-like containers, ie. the kernel.
        Simplex(const Simplex& s) : Simplex()
        {
            m_label      = s.m_label;
            std::copy(s.m_boundary->begin(), s.m_boundary->end(), m_boundary->begin());
            std::copy(s.m_co_boundary->begin(), s.m_co_boundary->end(), m_co_boundary->begin());
        }
        
    public:
        
        co_boundary_list* get_co_boundary() const
        {
            return m_co_boundary;
        }
        boundary_list* get_boundary() const
        {
            return m_boundary;
        }
        co_boundary_list* get_co_boundary()
        {
            return m_co_boundary;
        }
        boundary_list* get_boundary()
        {
            return m_boundary;
        }
        
        void add_co_face(co_boundary_key_type key)
        {
            if(!m_co_boundary->contains(key))
            {
                m_co_boundary->push_back(key);
            }
        }
        
        void add_face(boundary_key_type key)
        {
            if(!m_boundary->contains(key))
            {
                m_boundary->push_back(key);
            }
        }
        
        void remove_co_face(const co_boundary_key_type& key)
        {
            auto iter = std::find(m_co_boundary->begin(), m_co_boundary->end(), key);
            assert(iter != m_co_boundary->end());
            m_co_boundary->erase(iter);
        }
        
        void remove_face(const boundary_key_type& key)
        {
            auto iter = std::find(m_boundary->begin(), m_boundary->end(), key);
            assert(iter != m_boundary->end());
            m_boundary->erase(iter);
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
        
        void invert_orientation()
        {
            assert(m_boundary->size() == 2);
            std::swap((*m_boundary)[0], (*m_boundary)[1]);
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
        
        void invert_orientation()
        {
            assert(m_boundary->size() == 3);
            std::swap((*m_boundary)[1], (*m_boundary)[2]);
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
        
        void invert_orientation()
        {
            assert(m_boundary->size() == 4);
            std::swap((*m_boundary)[2], (*m_boundary)[3]);
        }
    };
}
