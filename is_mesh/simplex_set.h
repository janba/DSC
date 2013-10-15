
#pragma once

namespace is_mesh
{
    
    template<typename key_type>
    class SimplexSet : public std::vector<key_type>
    {
    public:
        using std::vector<key_type>::vector;
        
        bool contains(const key_type& k) const
        {
            return std::find(this->begin(), this->end(), k) != this->end();
        }
        
        SimplexSet<key_type>& operator+=(const SimplexSet<key_type>& set)
        {
            for (auto &k : set) {
                *this += k;
            }
            return *this;
        }
        
        SimplexSet<key_type>& operator+=(const key_type& key)
        {
            if(!contains(key))
            {
                this->push_back(key);
            }
            return *this;
        }
        
        SimplexSet<key_type>& operator-=(const SimplexSet<key_type>& set)
        {
            for (auto &k : set) {
                *this -= k;
            }
            return *this;
        }
        
        SimplexSet<key_type>& operator-=(const key_type& key)
        {
            auto iter = std::find(this->begin(), this->end(), key);
            if(iter != this->end())
            {
                this->erase(iter);
            }
            return *this;
        }
    };
    
    template<typename key_type>
    bool operator==(const SimplexSet<key_type>& A, const SimplexSet<key_type>& B)
    {
        if(A.size() == B.size())
        {
            for (auto k : A)
            {
                if(!B.contains(k))
                {
                    return false;
                }
            }
            return true;
        }
        return false;
    }
    
    /**
     *  Returns the union of the two sets A and B.
     */
    template<typename key_type>
    SimplexSet<key_type> operator+(const SimplexSet<key_type>& A, const SimplexSet<key_type>& B)
    {
        SimplexSet<key_type> C;
        return (C += A) += B;
    }
    
    /**
     *  Returns the union of the two sets A and B.
     */
    template<typename key_type>
    SimplexSet<key_type>&& operator+(SimplexSet<key_type>&& A, const SimplexSet<key_type>& B)
    {
        return std::move(A += B);
    }
    
    /**
     *  Returns the union of the two sets A and B.
     */
    template<typename key_type>
    SimplexSet<key_type>&& operator+(SimplexSet<key_type>&& A, SimplexSet<key_type>&& B)
    {
        return std::move(A) + B;
    }

    /**
     *  Returns the union of the two sets A and B.
     */
    template<typename key_type>
    SimplexSet<key_type>&& operator+(const SimplexSet<key_type>& A, SimplexSet<key_type>&& B)
    {
        return std::move(B) + A;
    }
    
    /**
     *  Returns set A without the elements in set B.
     */
    template<typename key_type>
    SimplexSet<key_type> operator-(const SimplexSet<key_type>& A, const SimplexSet<key_type>& B)
    {
        SimplexSet<key_type> C;
        return (C += A) -= B;
    }
    
    /**
     *  Returns set A without the elements in set B.
     */
    template<typename key_type>
    SimplexSet<key_type>&& operator-(SimplexSet<key_type>&& A, const SimplexSet<key_type>& B)
    {
        return std::move(A -= B);
    }
    
    /**
     *  Returns set A without the element key.
     */
    template<typename key_type>
    SimplexSet<key_type> operator-(const SimplexSet<key_type>& A, const key_type& key)
    {
        SimplexSet<key_type> B = {key};
        return A - B;
    }
    
    /**
     *  Returns set A without the element key.
     */
    template<typename key_type>
    SimplexSet<key_type>&& operator-(SimplexSet<key_type>&& A, const key_type& key)
    {
        return std::move(A -= key);
    }
    
    /**
     *  Returns the intersection of the two sets.
     */
    template<typename key_type>
    SimplexSet<key_type> operator&(const SimplexSet<key_type>& set1, const SimplexSet<key_type>& set2)
    {
        return (set1 + set2) - (set1 | set2);
    }
    
    /**
     *  Returns the difference between the two sets, ie. the elements which are in either set1 or set2 but not both.
     */
    template<typename key_type>
    SimplexSet<key_type> operator|(const SimplexSet<key_type>& set1, const SimplexSet<key_type>& set2)
    {
        return (set1 - set2) + (set2 - set1);
    }
    
    inline void simplex_set_test()
    {
        std::cout << "Testing simplex set class: ";
        SimplexSet<int> A = {1,3,9,4};
        SimplexSet<int> B = {1,7,5,3,10};
        
        SimplexSet<int> U = {1,3,9,4,7,5,10};
        assert((A+B) == U);
        
        SimplexSet<int> C = {9,4};
        assert((A-B) == C);
        
        SimplexSet<int> D = {9,4,7,5,10};
        assert((A|B) == D);
        
        SimplexSet<int> I = {1,3};
        assert((A&B) == I);
        
        A -= 3;
        A += 9;
        A += 11;
        SimplexSet<int> E = {1,9,4,11};
        assert(A == E);
        
        std::cout << "PASSED" << std::endl;
    }
    
    
    /**
     * This class holds simplices of dimension 0 to 3, and provides basic set
     * operations upon a set of simplices.
     * Only handles are stored in the simplex set, actual data must be looked up elsewhere.
     */
    template<typename node_key_type, typename edge_key_type, typename face_key_type, typename tetrahedron_key_type>
    class simplex_set
    {
    public:
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
        
        typename node_set::iterator insert(node_key_type const & k) { return m_nodes.insert(k).first; }
        typename edge_set::iterator insert(edge_key_type const & k) { return m_edges.insert(k).first; }
        typename face_set::iterator insert(face_key_type const & k) { return m_faces.insert(k).first; }
        typename tetrahedron_set::iterator insert(tetrahedron_key_type const & k) { return m_tetrahedra.insert(k).first; }
        
        template<typename _first, typename _last>
        void insert(_first first, _last last)
        {
            while(first != last)
            {
                insert(*first);
                ++first;
            }
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
        
        void clear_nodes()
        {
            m_nodes.clear();
        }
        
        void clear_edges()
        {
            m_edges.clear();
        }
        
        void clear_faces()
        {
            m_faces.clear();
        }
        
        void clear_tetrahedra()
        {
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
        
        node_set get_nodes()
        {
            return m_nodes;
        }
        
        edge_set get_edges()
        {
            return m_edges;
        }
        
        face_set get_faces()
        {
            return m_faces;
        }
        
        tetrahedron_set get_tets()
        {
            return m_tetrahedra;
        }
        
        node_set_iterator nodes_begin() { return m_nodes.begin(); }
        node_set_iterator nodes_end()   { return m_nodes.end(); }
        edge_set_iterator edges_begin() { return m_edges.begin(); }
        edge_set_iterator edges_end()   { return m_edges.end(); }
        face_set_iterator faces_begin() { return m_faces.begin(); }
        face_set_iterator faces_end()   { return m_faces.end(); }
        tetrahedron_set_iterator tetrahedra_begin() { return m_tetrahedra.begin(); }
        tetrahedron_set_iterator tetrahedra_end()   { return m_tetrahedra.end(); }
    };
    
}
