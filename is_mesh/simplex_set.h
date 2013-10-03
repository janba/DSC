
#pragma once

#include <is_mesh/is_mesh_utils.h>

namespace is_mesh
{
    
    /**
     * This class holds simplices of dimension 0 to 3, and provides basic set
     * operations upon a set of simplices.
     * Only handles are stored in the simplex set, actual data must be looked up elsewhere.
     */
    template<typename node_key_t_, typename edge_key_t_, typename face_key_t_, typename tetrahedron_key_t_>
    class simplex_set
    {
    public:
        typedef          node_key_t_                                  node_key_type;
        typedef          edge_key_t_                                  edge_key_type;
        typedef          face_key_t_                                  face_key_type;
        typedef          tetrahedron_key_t_                           tetrahedron_key_type;
        typedef simplex_set<node_key_type, edge_key_type, face_key_type, tetrahedron_key_type>      simplex_set_type;
        
        typedef          std::set<node_key_type>                      node_set;
        typedef          std::set<edge_key_type>                      edge_set;
        typedef          std::set<face_key_type>                      face_set;
        typedef          std::set<tetrahedron_key_type>               tetrahedron_set;
        typedef typename node_set::iterator                           node_set_iterator;
        typedef typename edge_set::iterator                           edge_set_iterator;
        typedef typename face_set::iterator                           face_set_iterator;
        typedef typename tetrahedron_set::iterator                    tetrahedron_set_iterator;
        
        typedef typename node_set::size_type                          size_type;
        
        template<int d> struct iterator_type{};
        
    private:
        node_set         m_nodes;
        edge_set         m_edges;
        face_set         m_faces;
        tetrahedron_set  m_tetrahedra;
        
        // alters a so that only elements in both a and b are in a.
        template<typename set_type>
        void intersection_helper(set_type& a, set_type const & b)
        {
            if (b.empty())
            {
                a.clear();
                return;
            }
            typename set_type::size_type size = (a.size() > b.size())?a.size():b.size();
            std::vector<typename set_type::value_type> res(size);
            typename std::vector<typename set_type::value_type>::iterator new_last;
            new_last = std::set_intersection(a.begin(), a.end(), b.begin(), b.end(), res.begin());
            a.clear();
            a.insert(res.begin(), new_last);
        }
        
        //removes all elements from a that are also in b
        template<typename set_type>
        void difference_helper(set_type& a, set_type const & b)
        {
            if (b.empty()) return;
            typename set_type::size_type size = (a.size() > b.size())?a.size():b.size();
            std::vector<typename set_type::value_type> res(size);
            typename std::vector<typename set_type::value_type>::iterator new_last;
            new_last = std::set_difference(a.begin(), a.end(), b.begin(), b.end(), res.begin());
            a.clear();
            a.insert(res.begin(), new_last);
        }
        
        template<typename set_type>
        void union_helper(set_type& a, set_type const & b)
        {
            if (b.empty()) return;
            std::vector<typename set_type::value_type> res(a.size() + b.size());
            std::merge(a.begin(), a.end(), b.begin(), b.end(), res.begin());
            typename std::vector<typename set_type::value_type>::iterator new_last = std::unique(res.begin(), res.end());
            a.clear();
            a.insert(res.begin(), new_last);
        }
    public:
        
        //TODO constructor
        
        //TODO copy constructor
        
        //TODO assignment
        
        
        template<typename iterator>
        void insert_nodes(iterator begin, iterator end)
        {
            m_nodes.insert(begin, end);
        }
        
        template<typename iterator>
        void insert_edges(iterator begin, iterator end)
        {
            m_edges.insert(begin, end);
        }
        
        template<typename iterator>
        void insert_faces(iterator begin, iterator end)
        {
            m_faces.insert(begin, end);
        }
        
        template<typename iterator>
        void insert_tetrahedron(iterator begin, iterator end)
        {
            m_tetrahedra.insert(begin, end);
        }
        
        node_set_iterator insert(node_key_type const & k) { return m_nodes.insert(k).first; }
        edge_set_iterator insert(edge_key_type const & k) { return m_edges.insert(k).first; }
        face_set_iterator insert(face_key_type const & k) { return m_faces.insert(k).first; }
        tetrahedron_set_iterator insert(tetrahedron_key_type const & k) { return m_tetrahedra.insert(k).first; }
        
        template<typename _first, typename _last>
        void insert(_first first, _last last)
        {
            while(first != last)
            {
                insert(*first);
                ++first;
            }
        }
        
        void filter(util::node_simplex_tag)
        {
            m_edges.clear();
            m_faces.clear();
            m_tetrahedra.clear();
        }
        void filter(util::edge_simplex_tag)
        {
            m_nodes.clear();
            m_faces.clear();
            m_tetrahedra.clear();
        }
        void filter(util::face_simplex_tag)
        {
            m_nodes.clear();
            m_edges.clear();
            m_tetrahedra.clear();
        }
        void filter(util::tetrahedron_simplex_tag)
        {
            m_nodes.clear();
            m_edges.clear();
            m_faces.clear();
        }
        
        void difference(const simplex_set_type& s)
        {
            difference_helper(m_nodes, s.m_nodes);
            difference_helper(m_edges, s.m_edges);
            difference_helper(m_faces, s.m_faces);
            difference_helper(m_tetrahedra, s.m_tetrahedra);
        }
        
        void intersection(const simplex_set_type& s)
        {
            //find intersection between all four sets, one at a time
            intersection_helper(m_nodes, s.m_nodes);
            intersection_helper(m_edges, s.m_edges);
            intersection_helper(m_faces, s.m_faces);
            intersection_helper(m_tetrahedra, s.m_tetrahedra);
        }
        
        //the union operation - cannot be called union as this is a reserved word.
        void add(const simplex_set_type& s)
        {
            union_helper(m_nodes, s.m_nodes);
            union_helper(m_edges, s.m_edges);
            union_helper(m_faces, s.m_faces);
            union_helper(m_tetrahedra, s.m_tetrahedra);
        }
        
        bool contains(node_key_type const & k) { return m_nodes.count(k) > 0; }
        bool contains(edge_key_type const & k) { return m_edges.count(k) > 0; }
        bool contains(face_key_type const & k) { return m_faces.count(k) > 0; }
        bool contains(tetrahedron_key_type const & k) { return m_tetrahedra.count(k) > 0; }
        
        void erase(node_key_type const & k) { m_nodes.erase(k); }
        void erase(edge_key_type const & k) { m_edges.erase(k); }
        void erase(face_key_type const & k) { m_faces.erase(k); }
        void erase(tetrahedron_key_type const & k) { m_tetrahedra.erase(k); }
        
        void clear()
        {
            m_nodes.clear();
            m_edges.clear();
            m_faces.clear();
            m_tetrahedra.clear();
        }
        
        bool empty() const
        {
            return m_nodes.empty() && m_edges.empty() && m_edges.empty() && m_tetrahedra.empty();
        }
        
        size_type size() const
        {
            return m_nodes.size() + m_edges.size() + m_faces.size() + m_tetrahedra.size();
        }
        
        size_type size_nodes() const
        {
            return m_nodes.size();
        }
        size_type size_edges() const
        {
            return m_edges.size();
        }
        size_type size_faces() const
        {
            return m_faces.size();
        }
        size_type size_tetrahedra() const
        {
            return m_tetrahedra.size();
        }
        
        node_set_iterator nodes_begin() { return m_nodes.begin(); }
        node_set_iterator nodes_end()   { return m_nodes.end(); }
        edge_set_iterator edges_begin() { return m_edges.begin(); }
        edge_set_iterator edges_end()   { return m_edges.end(); }
        face_set_iterator faces_begin() { return m_faces.begin(); }
        face_set_iterator faces_end()   { return m_faces.end(); }
        tetrahedron_set_iterator tetrahedra_begin() { return m_tetrahedra.begin(); }
        tetrahedron_set_iterator tetrahedra_end()   { return m_tetrahedra.end(); }
        
        node_set_iterator begin(util::node_simplex_tag) { return nodes_begin(); }
        edge_set_iterator begin(util::edge_simplex_tag) { return edges_begin(); }
        face_set_iterator begin(util::face_simplex_tag) { return faces_begin(); }
        tetrahedron_set_iterator begin(util::tetrahedron_simplex_tag) { return tetrahedra_begin(); }
        
        node_set_iterator end(util::node_simplex_tag) { return nodes_end(); }
        edge_set_iterator end(util::edge_simplex_tag) { return edges_end(); }
        face_set_iterator end(util::face_simplex_tag) { return faces_end(); }
        tetrahedron_set_iterator end(util::tetrahedron_simplex_tag) { return tetrahedra_end(); }
    };
    
    template<> template<>
    struct simplex_set<NodeKey,EdgeKey,FaceKey,TetrahedronKey>::iterator_type<0>
    {
        typedef node_set_iterator           iterator;
    };
    
    template<> template<>
    struct simplex_set<NodeKey,EdgeKey,FaceKey,TetrahedronKey>::iterator_type<1>
    {
        typedef edge_set_iterator           iterator;
    };
    
    template<> template<>
    struct simplex_set<NodeKey,EdgeKey,FaceKey,TetrahedronKey>::iterator_type<2>
    {
        typedef face_set_iterator           iterator;
    };
    
    template<> template<>
    struct simplex_set<NodeKey,EdgeKey,FaceKey,TetrahedronKey>::iterator_type<3>
    {
        typedef tetrahedron_set_iterator    iterator;
    };
    
}
