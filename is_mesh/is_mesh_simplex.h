#pragma once

#include <set>
#include <is_mesh/is_mesh_key_type.h>

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
        typedef std::vector<boundary_key_type> boundary_type;
        typedef          std::set<co_boundary_key_type>         set_type;
        
        typedef          std::vector<boundary_key_type>*                             boundary_list;
        typedef          set_type*                              co_boundary_set;
        typedef typename std::vector<boundary_key_type>::iterator                    boundary_iterator;
        typedef typename set_type::iterator                     co_boundary_iterator;
        
    protected:
        std::vector<boundary_key_type>* m_boundary = nullptr;
        set_type* m_co_boundary = nullptr;
        bool             		m_is_compact;
        int              		m_label;       //used in coloring - to identify connected components.
        
    public:
        
        Simplex() : m_is_compact(false), m_label(0)
        {
            m_boundary    = new std::vector<boundary_key_type>();
            m_co_boundary = new std::set<co_boundary_key_type>();
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
        Simplex(const Simplex& s) //: m_co_boundary(0), m_boundary(0)
        {
            m_is_compact = s.m_is_compact;
            m_label      = s.m_label;
            m_boundary   = 0;
            m_co_boundary= 0;
            if (s.m_boundary != 0)
            {
                m_boundary = new std::vector<boundary_key_type>();
                std::copy(s.m_boundary->begin(), s.m_boundary->end(), m_boundary->begin());
            }
            if (s.m_co_boundary != 0)
            {
                m_co_boundary = new set_type();
                typename set_type::iterator sb, se;
                sb = s.m_co_boundary->begin();
                se = s.m_co_boundary->end();
                while (sb != se)
                {
                    m_co_boundary->insert(*sb);
                    ++sb;
                }
                //std::copy(s.m_co_boundary->begin(), s.m_co_boundary->end(), m_co_boundary->begin());
            }
        }
        
    public:
        
        set_type* get_co_boundary() const
        {
            return m_co_boundary;
        }
        boundary_list get_boundary() const
        {
            return m_boundary;
        }
        set_type* get_co_boundary()
        {
            return m_co_boundary;
        }
        boundary_list get_boundary()
        {
            return m_boundary;
        }
        void set_co_boundary_set(set_type* co_boundary)
        {
            m_co_boundary = co_boundary;
        }
        void set_boundary_list(boundary_list boundary)
        {
            m_boundary = boundary;
        }
        void add_co_face(co_boundary_key_type n)
        {
            m_co_boundary->insert(n);
        }
        void add_face(boundary_key_type n)
        {
            m_boundary->push_back(n);
        }
        void remove_co_face(co_boundary_key_type const & n)
        {
            m_co_boundary->erase(n);
        }
        
        bool is_compact(){ return m_is_compact; }
        void set_compact(bool c) { m_is_compact = c; }
    };
    
    ///////////////////////////////////////////////////////////////////////////////
    ///  N O D E
    ///////////////////////////////////////////////////////////////////////////////
    template<typename NodeTraits, typename Mesh>
    class Node : public NodeTraits, public Simplex<Key, EdgeKey>
    {
    public:
        typedef Simplex<Key, EdgeKey>         Simplex;
        typedef          NodeTraits                                                         type_traits;
        typedef typename Mesh::edge_type                                              co_boundary_type;
        
        Node() : Simplex()
        {
            
        }
        Node(const type_traits & t) : type_traits(t), Simplex()
        {
            
        }
    };
    
    ///////////////////////////////////////////////////////////////////////////////
    ///  E D G E
    ///////////////////////////////////////////////////////////////////////////////
    template<typename EdgeTraits, typename Mesh>
    class Edge : public EdgeTraits, public Simplex<NodeKey, FaceKey>
    {
    public:
        typedef Simplex<NodeKey, FaceKey> Simplex;
        typedef          EdgeTraits                                                       type_traits;
        typedef typename Mesh::node_type                                            boundary_type;
        typedef typename Mesh::face_type                                            co_boundary_type;
        
        Edge() : Simplex()
        {
            
        }
        Edge(const type_traits & t) : type_traits(t), Simplex()
        {
            
        }
    };
    
    ///////////////////////////////////////////////////////////////////////////////
    //  F A C E
    ///////////////////////////////////////////////////////////////////////////////
    template<typename FaceTraits, typename Mesh>
    class Face : public FaceTraits, public Simplex<EdgeKey, TetrahedronKey>
    {
    public:
        typedef Simplex<EdgeKey, TetrahedronKey>   Simplex;
        typedef          FaceTraits                                       type_traits;
        typedef typename Mesh::edge_type                            boundary_type;
        typedef typename Mesh::tetrahedron_type                     co_boundary_type;
        
        Face() : Simplex()
        {
            Simplex::set_compact(true);
        }
        Face(const type_traits & t) : type_traits(t), Simplex()
        {
            Simplex::set_compact(true);
        }
    };
    
    ///////////////////////////////////////////////////////////////////////////////
    // T E T R A H E D R O N
    ///////////////////////////////////////////////////////////////////////////////
    template<typename TetrahedronTraits, typename Mesh>
    class Tetrahedron : public TetrahedronTraits, public Simplex<FaceKey, Key>
    {
    public:
        typedef Simplex<FaceKey, Key>        Simplex;
        typedef          TetrahedronTraits                                type_traits;
        typedef typename Mesh::face_type                            boundary_type;
        
        Tetrahedron() : Simplex()
        {
            Simplex::set_compact(true);
        }
        Tetrahedron(const type_traits & t) : type_traits(t), Simplex()
        {
            Simplex::set_compact(true);
        }
    };
}
