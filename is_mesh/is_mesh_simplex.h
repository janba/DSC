#pragma once

#include <set>
#include <is_mesh/is_mesh_key_type.h>

namespace is_mesh
{
    ///////////////////////////////////////////////////////////////////////////////
    // S I M P L E X   B A S E   C L A S S
    ///////////////////////////////////////////////////////////////////////////////
    namespace util
    {
        /**
         * Base class for all simplex classes
         */
        template<typename mesh_type_, typename key_type_>
        class Simplex
        {
        public:
            typedef          mesh_type_                             mesh_type;
            typedef          key_type_                              key_type;
            typedef typename mesh_type::node_type                   node_type;
            typedef typename mesh_type::edge_type                   edge_type;
            typedef typename mesh_type::face_type                   face_type;
            typedef typename mesh_type::tetrahedron_type            tetrahedron_type;
            
            typedef typename key_traits<key_type::dim-1>::key_type  boundary_key_type;
            typedef typename key_traits<key_type::dim+1>::key_type  co_boundary_key_type;
            
            typedef          std::vector<boundary_key_type>         list_type;
            typedef          std::set<co_boundary_key_type>         set_type;
            
            typedef          list_type*                             boundary_list;
            typedef          set_type*                              co_boundary_set;
            typedef typename list_type::iterator                    boundary_iterator;
            typedef typename set_type::iterator                     co_boundary_iterator;
            
            static const unsigned int dim = key_type::dim;
            
        protected:
            list_type*       		m_boundary;
            set_type*     		 	m_co_boundary;
            bool             		m_is_compact;
            int              		m_label;       //used in coloring - to identify connected components.
            
        public:
            
            Simplex() : m_boundary(0), m_co_boundary(0), m_is_compact(false), m_label(0)
            { }
            
            ~Simplex() { }
            
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
                    m_boundary = new list_type();
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
            list_type* get_boundary() const
            {
                return m_boundary;
            }
            set_type* get_co_boundary()
            {
                return m_co_boundary;
            }
            list_type* get_boundary()
            {
                return m_boundary;
            }
            void set_co_boundary_set(set_type* co_boundary)
            {
                m_co_boundary = co_boundary;
            }
            void set_boundary_list(list_type* boundary)
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
    }
    
    
    ///////////////////////////////////////////////////////////////////////////////
    ///  N O D E
    ///////////////////////////////////////////////////////////////////////////////
    template<typename NodeTraits, typename Mesh>
    class Node : public NodeTraits, public util::Simplex<Mesh, typename Mesh::node_key_type>
    {
    public:
        typedef util::Simplex<Mesh, typename Mesh::node_key_type>         Simplex;
        typedef          NodeTraits                                                         type_traits;
        typedef typename Mesh::edge_type                                              co_boundary_type;
        
        Node()
        {
            Simplex::m_co_boundary = new typename Simplex::set_type();
        }
        Node(const type_traits & t) : type_traits(t)
        {
            Simplex::m_co_boundary = new typename Simplex::set_type();
        }
        ~Node()
        {
            delete Simplex::m_co_boundary;
        }
    };
    
    ///////////////////////////////////////////////////////////////////////////////
    ///  E D G E
    ///////////////////////////////////////////////////////////////////////////////
    template<typename EdgeTraits, typename Mesh>
    class Edge : public EdgeTraits, public util::Simplex<Mesh, typename Mesh::edge_key_type>
    {
    public:
        typedef util::Simplex<Mesh, typename Mesh::edge_key_type> Simplex;
        typedef          EdgeTraits                                                       type_traits;
        typedef typename Simplex::node_type                                            boundary_type;
        typedef typename Simplex::face_type                                            co_boundary_type;
        
        Edge()
        {
            Simplex::m_co_boundary = new typename Simplex::set_type();
            Simplex::m_boundary    = new typename Simplex::list_type();
        }
        Edge(const type_traits & t) : type_traits(t)
        {
            Simplex::m_co_boundary = new typename Simplex::set_type();
            Simplex::m_boundary    = new typename Simplex::list_type();
        }
        ~Edge()
        {
            delete Simplex::m_co_boundary;
            delete Simplex::m_boundary;
        }
    };
    
    ///////////////////////////////////////////////////////////////////////////////
    //  F A C E
    ///////////////////////////////////////////////////////////////////////////////
    template<typename FaceTraits, typename Mesh>
    class Face : public FaceTraits, public util::Simplex<Mesh, typename Mesh::face_key_type>
    {
    public:
        typedef util::Simplex<Mesh, typename Mesh::face_key_type>   Simplex;
        typedef          FaceTraits                                       type_traits;
        typedef typename Simplex::edge_type                            boundary_type;
        typedef typename Simplex::tetrahedron_type                     co_boundary_type;

        Face()
        {
            Simplex::m_co_boundary = new typename Simplex::set_type();
            Simplex::m_boundary    = new typename Simplex::list_type();
            Simplex::set_compact(true);
        }
        Face(const type_traits & t) : type_traits(t)
        {
            Simplex::m_co_boundary = new typename Simplex::set_type();
            Simplex::m_boundary    = new typename Simplex::list_type();
            Simplex::set_compact(true);
        }
        ~Face()
        {
            delete Simplex::m_co_boundary;
            delete Simplex::m_boundary;
        }
    };
    
    ///////////////////////////////////////////////////////////////////////////////
    // T E T R A H E D R O N
    ///////////////////////////////////////////////////////////////////////////////
    template<typename TetrahedronTraits, typename Mesh>
    class Tetrahedron : public TetrahedronTraits, public util::Simplex<Mesh, typename Mesh::tetrahedron_key_type>
    {
    public:
        typedef util::Simplex<Mesh, typename Mesh::tetrahedron_key_type>        Simplex;
        typedef          TetrahedronTraits                                type_traits;
        typedef typename Simplex::face_type                            boundary_type;
        
        Tetrahedron() 
        { 
            Simplex::m_boundary    = new typename Simplex::list_type();
            Simplex::set_compact(true);
        }
        Tetrahedron(const type_traits & t) : type_traits(t) 
        { 
            Simplex::m_boundary    = new typename Simplex::list_type();
            Simplex::set_compact(true);
        }
        ~Tetrahedron() 
        { 
            delete Simplex::m_boundary;
        }
    };
}
