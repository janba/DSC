#ifndef PROGRESSIVE_MESH_INCIDENCE_SIMPLICIAL_SIMPLEX_H
#define PROGRESSIVE_MESH_INCIDENCE_SIMPLICIAL_SIMPLEX_H

// Disclaimer / copyrights / stuff

#include <iterator>
#include <vector>
#include <set>

#include <is_mesh/is_mesh.h>
#include <is_mesh/is_mesh_key_type.h>

#include <is_mesh/io/is_mesh_xml_write.h>

namespace OpenTissue
{
  namespace is_mesh
  {
    // Forward declarations
    template<typename NT,typename TT, typename ET, typename FT> class t4mesh;

    ///////////////////////////////////////////////////////////////////////////////
    // S I M P L E X   B A S E   C L A S S
    ///////////////////////////////////////////////////////////////////////////////
    namespace util
    {
      /**
       * Base class for all simplex classes
       */
      template<typename simplex_type_, typename mesh_type_, typename key_type_>
      class Simplex
      {
      public:
        typedef          simplex_type_                          simplex_type;
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

       // friend mesh_type;
        friend void xml_write(std::string const &, mesh_type const &);


        static const unsigned int dim = key_type::dim;

      protected:
        list_type*       		m_boundary;
        set_type*     		 	m_co_boundary;
        bool             		m_is_compact;
        int              		m_label;       //used in coloring - to identify connected components.

      public:
        template<int n>
        void init()
        {
          if (n > 0)
            m_boundary = new list_type();
          if (n < 3)
            m_co_boundary = new set_type();
        }

		template<int n>
        void destroy() 
        {
          if (n > 0)
            delete m_boundary;
          if (n < 3)
            delete m_co_boundary;
        }

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
            init<3>(); //init m_boundary
            std::copy(s.m_boundary->begin(), s.m_boundary->end(), m_boundary->begin());
          }
          if (s.m_co_boundary != 0)
          {
            init<0>(); //init m_co_boundary
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
        Simplex& operator=(const Simplex& s) //asignment
        {
          m_is_compact = s.m_is_compact;
          m_label      = s.m_label;
          m_boundary   = 0;
          m_co_boundary= 0;
          if (s.m_boundary != 0)
          {
            init<3>(); //init m_boundary
            m_boundary->insert(m_boundary->begin(), s.m_boundary->begin(), s.m_boundary->end());
          }
          if (s.m_co_boundary != 0)
          {
            init<0>(); //init m_boundary
            m_co_boundary->insert(s.m_co_boundary->begin(), s.m_co_boundary->end());
          }
          return *this;
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
      }; //Simplex
    } //namespace util


    ///////////////////////////////////////////////////////////////////////////////
    ///  N O D E
    ///////////////////////////////////////////////////////////////////////////////
    template<
        typename NodeTraits
      , typename TetrahedronTraits
      , typename EdgeTraits
      , typename FaceTraits
    >
    class Node :
        public NodeTraits
      , public util::Simplex<
            Node<NodeTraits,TetrahedronTraits,EdgeTraits,FaceTraits>
          , t4mesh<NodeTraits,TetrahedronTraits,EdgeTraits,FaceTraits>
          , typename t4mesh<NodeTraits,TetrahedronTraits,EdgeTraits,FaceTraits>::node_key_type
        >
    {
    public:
      typedef util::Simplex<
            Node<NodeTraits,TetrahedronTraits,EdgeTraits,FaceTraits>
          , t4mesh<NodeTraits,TetrahedronTraits,EdgeTraits,FaceTraits>
          , typename t4mesh<NodeTraits,TetrahedronTraits,EdgeTraits,FaceTraits>::node_key_type
        >                                                                                 base_class;
      typedef          NodeTraits                                                         type_traits;
      typedef typename base_class::edge_type                                              co_boundary_type;

      friend class t4mesh<NodeTraits,TetrahedronTraits,EdgeTraits,FaceTraits>;

    private:
    public:
      Node() 
      {
        base_class::m_co_boundary = new typename base_class::set_type();
	  }
      Node(const type_traits & t) : type_traits(t)
      {
        base_class::m_co_boundary = new typename base_class::set_type();
      }
      ~Node() 
      { 
	    delete base_class::m_co_boundary;
      }
    };

    ///////////////////////////////////////////////////////////////////////////////
    ///  E D G E
    ///////////////////////////////////////////////////////////////////////////////
    template<
        typename NodeTraits
      , typename TetrahedronTraits
      , typename EdgeTraits
      , typename FaceTraits
    >
    class Edge :
        public EdgeTraits
      , public util::Simplex<
            Edge<NodeTraits,TetrahedronTraits,EdgeTraits,FaceTraits>
          , t4mesh<NodeTraits,TetrahedronTraits,EdgeTraits,FaceTraits>
          , typename t4mesh<NodeTraits,TetrahedronTraits,EdgeTraits,FaceTraits>::edge_key_type
        >
    {
    private:
    public:
      typedef util::Simplex<
            Edge<NodeTraits,TetrahedronTraits,EdgeTraits,FaceTraits>
          , t4mesh<NodeTraits,TetrahedronTraits,EdgeTraits,FaceTraits>
          , typename t4mesh<NodeTraits,TetrahedronTraits,EdgeTraits,FaceTraits>::edge_key_type
        >                                                                               base_class;
      typedef          EdgeTraits                                                       type_traits;
      typedef typename base_class::node_type                                            boundary_type;
      typedef typename base_class::face_type                                            co_boundary_type;

      friend class t4mesh<NodeTraits,TetrahedronTraits,EdgeTraits,FaceTraits>;
    private:
    public:
      Edge() 
      { 
      	base_class::m_co_boundary = new typename base_class::set_type();
      	base_class::m_boundary    = new typename base_class::list_type();
      }
      Edge(const type_traits & t) : type_traits(t) 
      { 
      	base_class::m_co_boundary = new typename base_class::set_type();
      	base_class::m_boundary    = new typename base_class::list_type();
      }
      ~Edge() 
      { 
      	delete base_class::m_co_boundary;
      	delete base_class::m_boundary;
      }
    };

    ///////////////////////////////////////////////////////////////////////////////
    //  F A C E
    ///////////////////////////////////////////////////////////////////////////////
    template<
        typename NodeTraits
      , typename TetrahedronTraits
      , typename EdgeTraits
      , typename FaceTraits
    >
    class Face :
        public FaceTraits
      , public util::Simplex<
            Face<NodeTraits,TetrahedronTraits,EdgeTraits,FaceTraits>
          , t4mesh<NodeTraits,TetrahedronTraits,EdgeTraits,FaceTraits>
          , typename t4mesh<NodeTraits,TetrahedronTraits,EdgeTraits,FaceTraits>::face_key_type
        >
    {
    private:
    public:
      typedef util::Simplex<
            Face<NodeTraits,TetrahedronTraits,EdgeTraits,FaceTraits>
          , t4mesh<NodeTraits,TetrahedronTraits,EdgeTraits,FaceTraits>
          , typename t4mesh<NodeTraits,TetrahedronTraits,EdgeTraits,FaceTraits>::face_key_type
        >                                                               base_class;
      typedef          FaceTraits                                       type_traits;
      typedef typename base_class::edge_type                            boundary_type;
      typedef typename base_class::tetrahedron_type                     co_boundary_type;

      friend class t4mesh<NodeTraits,TetrahedronTraits,EdgeTraits,FaceTraits>;
    private:
    public:
      Face() 
      { 
      	base_class::m_co_boundary = new typename base_class::set_type();
      	base_class::m_boundary    = new typename base_class::list_type();
      	base_class::set_compact(true);
      }
      Face(const type_traits & t) : type_traits(t) 
      { 
      	base_class::m_co_boundary = new typename base_class::set_type();
      	base_class::m_boundary    = new typename base_class::list_type();
      	base_class::set_compact(true); 
      }
      ~Face() 
      { 
      	delete base_class::m_co_boundary;
      	delete base_class::m_boundary; 
      }
    };

    ///////////////////////////////////////////////////////////////////////////////
    // T E T R A H E D R O N
    ///////////////////////////////////////////////////////////////////////////////
    template<
        typename NodeTraits
      , typename TetrahedronTraits
      , typename EdgeTraits
      , typename FaceTraits
    >
    class Tetrahedron :
        public TetrahedronTraits
      , public util::Simplex<
            Tetrahedron<NodeTraits,TetrahedronTraits,EdgeTraits,FaceTraits>
          , t4mesh<NodeTraits, TetrahedronTraits, EdgeTraits, FaceTraits>
          , typename t4mesh<NodeTraits,TetrahedronTraits,EdgeTraits,FaceTraits>::tetrahedron_key_type
        >
    {
    private:
    public:
      typedef util::Simplex<
            Tetrahedron<NodeTraits,TetrahedronTraits,EdgeTraits,FaceTraits>
          , t4mesh<NodeTraits,TetrahedronTraits,EdgeTraits,FaceTraits>
          , typename t4mesh<NodeTraits,TetrahedronTraits,EdgeTraits,FaceTraits>::tetrahedron_key_type
        >                                                               base_class;
      typedef          TetrahedronTraits                                type_traits;
      typedef typename base_class::face_type                            boundary_type;

      friend class t4mesh<NodeTraits,TetrahedronTraits,EdgeTraits,FaceTraits>;
    private:
    public:
      Tetrahedron() 
      { 
      	base_class::m_boundary    = new typename base_class::list_type();
      	base_class::set_compact(true); 
      }
      Tetrahedron(const type_traits & t) : type_traits(t) 
      { 
      	base_class::m_boundary    = new typename base_class::list_type();
      	base_class::set_compact(true); 
      }
      ~Tetrahedron() 
      { 
      	delete base_class::m_boundary;
	  }
    };
  }//namespace is_mesh
} //namespace OpenTissue}
#endif //PROGRESSIVE_MESH_INCIDENCE_SIMPLICIAL_SIMPLEX_H
