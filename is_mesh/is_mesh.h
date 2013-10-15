
#pragma once

#include <algorithm>
#include <functional>
#include <queue>
#include <map>
#include <list>

#include <is_mesh/kernel.h>
#include <is_mesh/simplex.h>
#include <is_mesh/simplex_set.h>

namespace is_mesh
{
    /**
     * A data structure for managing a Simplicial Complex. Based on the work
     * by de Floriani and Hui, the Incidence Simplicial.
     * The complex is specialixed for 3-dimensional simplices, and can only
     * store 0-, 1-, 2-, nad 3-simplices.
     * Simplices are explicitly stored using a different memory kernel for
     * each type.h
     */
    template<typename NodeTraits, typename TetrahedronTraits, typename EdgeTraits, typename FaceTraits>
    class t4mesh
    {
    public:
        typedef Node<NodeTraits>                                         node_type;
        typedef Edge<EdgeTraits>                                         edge_type;
        typedef Face<FaceTraits>                                         face_type;
        typedef Tetrahedron<TetrahedronTraits>                           tetrahedron_type;
        
        typedef          simplex_set<NodeKey, EdgeKey, FaceKey, TetrahedronKey>            simplex_set_type;
        
    public:
        typedef          kernel<node_type, NodeKey>                            node_kernel_type;
        typedef          kernel<edge_type, EdgeKey>                            edge_kernel_type;
        typedef          kernel<face_type, FaceKey>                            face_kernel_type;
        typedef          kernel<tetrahedron_type, TetrahedronKey>              tetrahedron_kernel_type;
        
    public:
        
        typedef typename node_kernel_type::iterator                                  node_iterator;
        typedef typename edge_kernel_type::iterator                                  edge_iterator;
        typedef typename face_kernel_type::iterator                                  face_iterator;
        typedef typename tetrahedron_kernel_type::iterator                           tetrahedron_iterator;
        
        //expose some types from the kernel
        typedef typename node_kernel_type::size_type                                 size_type;
        
    public:
        node_kernel_type*                  m_node_kernel;
        edge_kernel_type*                  m_edge_kernel;
        face_kernel_type*                  m_face_kernel;
        tetrahedron_kernel_type*           m_tetrahedron_kernel;
        
    public:
        
        /**
         *
         */
        node_type & lookup_simplex(NodeKey const & k)
        {
            return m_node_kernel->find(k);
        }
        
        /**
         *
         */
        edge_type & lookup_simplex(EdgeKey const & k)
        {
            return m_edge_kernel->find(k);
        }
        
        /**
         *
         */
        face_type & lookup_simplex(FaceKey const & k)
        {
            return m_face_kernel->find(k);
        }
        
        /**
         *
         */
        tetrahedron_type & lookup_simplex(TetrahedronKey const & k)
        {
            return m_tetrahedron_kernel->find(k);
        }
        
        /**
         * Marek
         * Helper function for boundary and closure.
         */
        void boundary_helper(TetrahedronKey const & k, simplex_set_type & set)
        {
            auto s_boundary = lookup_simplex(k).get_boundary();
            auto it = s_boundary->begin();
            while (it != s_boundary->end())
            {
                set.insert(*it);
                boundary_helper(*it, set);
                ++it;
            }
        }
        
        void boundary_helper(FaceKey const & k, simplex_set_type & set)
        {
            auto s_boundary = lookup_simplex(k).get_boundary();
            auto it = s_boundary->begin();
            while (it != s_boundary->end())
            {
                set.insert(*it);
                boundary_helper(*it, set);
                ++it;
            }
        }
        
        void boundary_helper(EdgeKey const & k, simplex_set_type & set)
        {
            auto s_boundary = lookup_simplex(k).get_boundary();
            auto it = s_boundary->begin();
            while (it != s_boundary->end())
            {
                set.insert(*it);
                ++it;
            }
        }
        
        void boundary_helper(NodeKey const & k, simplex_set_type & set) {}
        
        
        void closure_helper(simplex_set_type & input_set, simplex_set_type & set)
        {
            for (auto &tid : input_set.get_tets())
            {
                set.insert(tid);
            }
            
            for (auto tid : set.get_tets())
            {
                auto t_boundary = lookup_simplex(tid).get_boundary();
                auto it = t_boundary->begin();
                while (it != t_boundary->end())
                {
                    set.insert(*it);
                    ++it;
                }
            }
            
            typename simplex_set_type::face_set_iterator fit = input_set.faces_begin();
            
            while (fit != input_set.faces_end())
            {
                set.insert(*fit);
                ++fit;
            }
            
            fit = set.faces_begin();
            while (fit != set.faces_end())
            {
                auto f_boundary = lookup_simplex(*fit).get_boundary();
                auto it = f_boundary->begin();
                while (it != f_boundary->end())
                {
                    set.insert(*it);
                    ++it;
                }
                ++fit;
            }
            
            typename simplex_set_type::edge_set_iterator eit = input_set.edges_begin();
            
            while (eit != input_set.edges_end())
            {
                set.insert(*eit);
                ++eit;
            }
            
            eit = set.edges_begin();
            while (eit != set.edges_end())
            {
                auto e_boundary = lookup_simplex(*eit).get_boundary();
                auto it = e_boundary->begin();
                while (it != e_boundary->end())
                {
                    set.insert(*it);
                    ++it;
                }
                ++eit;
            }
            
            typename simplex_set_type::node_set_iterator nit = input_set.nodes_begin();
            
            while (nit != input_set.nodes_end())
            {
                set.insert(*nit);
                ++nit;
            }
        }
        
        /**
         * Marek
         * Helper function performing a single transposition on simplex's boundary list.
         */
        template<typename key_type>
        void invert_orientation(key_type const & k)
        {
            if (k.dim == 0) return;
            
            auto boundary = lookup_simplex(k).get_boundary();
            auto it = boundary->begin();
            
            ++it;
            
            auto temp = *it;
            *it = *(boundary->begin());
            *(boundary->begin()) = temp;
        }
        
        /**
         * Marek
         * Helper function for orient_face[...] methods.
         * fk must be a face of sk, dim(sk) = dimension, dim(fk) = dimension-1.
         */
        void orient_face_helper(const TetrahedronKey& tid, const FaceKey& fid, bool consistently)
        {
            auto face_boundary = lookup_simplex(fid).get_boundary();
            std::vector<EdgeKey> new_face_boundary(face_boundary->size());
            
            unsigned char f_index = 0, i = 0;
            
            for (auto f : *lookup_simplex(tid).get_boundary())
            {
                if (f == fid)
                {
                    f_index = i+1;
                }
                else
                {
                    EdgeKey eid;
                    bool res = get_intersection(fid, f, eid);
                    assert(res || !"Two faces of the same simplex do not intersect?!");
                    new_face_boundary[i] = eid;
                    ++i;
                }
            }
            
            assert(f_index > 0 || !"fid is not a face of tid");
            
            face_boundary->clear();
            for (auto e : new_face_boundary)
            {
                face_boundary->push_back(e);
            }
            
            f_index %= 2;
            if ((f_index == 0 && consistently) || (f_index == 1 && !consistently))
            {
                invert_orientation(fid);
            }
        }
        
        /**
         * Marek
         * Helper function for orient_face[...] methods.
         * fk must be a face of sk, dim(sk) = dimension, dim(fk) = dimension-1.
         */
        void orient_face_helper(const FaceKey& fid, const EdgeKey& eid, bool consistently)
        {
            auto simplex_boundary = lookup_simplex(fid).get_boundary();
            auto face_boundary = lookup_simplex(eid).get_boundary();
            
            auto sb_it = simplex_boundary->begin();
            auto fb_it = face_boundary->begin();
            
            std::vector<NodeKey> new_face_boundary(face_boundary->size());
            
            unsigned char f_index = 0, i = 0;
            
            while (sb_it != simplex_boundary->end())
            {
                if (*sb_it == eid)
                {
                    f_index = i+1;
                }
                else
                {
                    NodeKey ek;
                    bool res = get_intersection(eid, *sb_it, ek);
                    assert(res || !"Two faces of the same simplex do not intersect?!");
                    new_face_boundary[i] = ek;
                    ++fb_it;
                    ++i;
                }
                ++sb_it;
            }
            
            assert(f_index > 0 || !"fk is not a face of sk");
            
            i = 0;
            /**/
            face_boundary->clear();
            for (; i < new_face_boundary.size(); ++i)
                face_boundary->push_back(new_face_boundary[i]);
            
            f_index %= 2;
            if ((f_index == 0 && consistently) ||
                (f_index == 1 && !consistently))
                invert_orientation(eid);
        }
        
        
        void orient_coface_helper(const FaceKey& sk, const TetrahedronKey& cfk, bool consistently)
        {
            assert(sk.dim == (cfk.dim - 1) || !"sk is not a boundary face of cfk.");
            assert(sk.dim < 3 || !"No simplices of dimension more than three.");
            assert(sk.dim > 0 || !"Vertices are not oriented.");
            
            auto coface_boundary = lookup_simplex(cfk).get_boundary();
            
            std::map<EdgeKey, FaceKey> face_to_simplex;
            
            for (auto cfb_it : *coface_boundary)
            {
                if (cfb_it != sk)
                {
                    EdgeKey k;
                    bool res = get_intersection(sk, cfb_it, k);
                    assert(res || !"Two faces of the same simplex do not intersect?!");
                    face_to_simplex[k] = cfb_it;
                }
            }
            
            auto cfb_it = coface_boundary->begin();
            *cfb_it = sk;
            ++cfb_it;
            
            for (auto sb_it : *lookup_simplex(sk).get_boundary())
            {
                *cfb_it = face_to_simplex[sb_it];
                ++cfb_it;
            }
            
            if (!consistently)
            {
                invert_orientation(cfk);
            }
        }
        
        /**
         * Marek
         * Helper function for orient_coface[...] methods.
         * sk must be a face of cfk, dim(sk) = dimension, dim(cfk) = dimension+1.
         */
        void orient_coface_helper(const EdgeKey& sk, const FaceKey& cfk, bool consistently)
        {
            assert(sk.dim == (cfk.dim - 1) || !"sk is not a boundary face of cfk.");
            assert(sk.dim < 3 || !"No simplices of dimension more than three.");
            assert(sk.dim > 0 || !"Vertices are not oriented.");
            
            auto coface_boundary = lookup_simplex(cfk).get_boundary();
            
            std::map<NodeKey, EdgeKey> face_to_simplex;
            
            for (auto cfb_it : *coface_boundary)
            {
                if (cfb_it != sk)
                {
                    NodeKey k;
                    bool res = get_intersection(sk, cfb_it, k);
                    assert(res || !"Two faces of the same simplex do not intersect?!");
                    face_to_simplex[k] = cfb_it;
                }
            }
            
            auto cfb_it = coface_boundary->begin();
            *cfb_it = sk;
            ++cfb_it;
            
            for (auto sb_it : *lookup_simplex(sk).get_boundary())
            {
                *cfb_it = face_to_simplex[sb_it];
                ++cfb_it;
            }
            
            if (!consistently)
            {
                invert_orientation(cfk);
            }
        }
        
    private:
        
        /**
         *
         */
        template<typename simplex_type>
        void reset_label(simplex_type & s)
        {
            s.reset_label();
            for(auto co_bound_it : *s.get_co_boundary())
            {
                auto& simplex = lookup_simplex(co_bound_it);
                if (simplex.get_label() != 0)
                {
                    reset_co_label(simplex);
                }
            }
        }
        
        /**
         *
         */
        template<typename simplex_type>
        void reset_co_label(simplex_type & s)
        {
            s.reset_label();
            for(auto bound_it : *s.get_boundary())
            {
                auto& simplex = lookup_simplex(bound_it);
                if (simplex.get_label() != 0) reset_label(simplex);
            }
        }
        
        /**
         * Returns whether a face (f) is in the boundary of a simplex (s).
         */
        template<typename key_type_face, typename key_type_simplex>
        bool in_boundary(key_type_face const & f, key_type_simplex const & s)
        {
            if (f.dim >= s.dim)
            {
                //cannot be in boundary as dimensions are wrong
                return false;
            }
            auto b_list = lookup_simplex(s).get_boundary();
            if (f.dim+1 < s.dim)
            {
                //too far apart.. need to call recursively
                bool in_bound = false;
                auto bound_it = b_list->begin();
                ++bound_it;
                for( ; bound_it != b_list->end() ; ++bound_it)
                {
                    in_bound = in_bound || in_boundary(f, *(bound_it));
                    if (in_bound) break;
                }
                return in_bound;
            }
            //this is the right dimensions.. search the boundary list
            auto bound_it = b_list->begin();
            for( ; bound_it != b_list->end() ; ++bound_it)
            {
                if (f == *bound_it) return true;
            }
            return false;
        }
        
        template<typename key_type_face>
        bool in_boundary(key_type_face const & f, NodeKey const & s)
        { //no simplex is in the boundary of a node...
            return false;
        }
        
        /**
         * Helper for the star operation. Searches up in the simplex hieracy.
         * Assumes that mesh is compressed. Should work in a non-compressed state too.
         * Template should traverse nodes and edges.
         *
         * @param s The original simplex that we are calulating the star of
         * @param t The current simplex that we are searching through
         * @param set The result is returned here.
         */
        template<typename simplex_key_s, typename simplex_key_t>
        void star_helper(simplex_key_s const & s, simplex_key_t const & t,  simplex_set_type & set)
        {
            assert(s.dim < t.dim || !"Star traversed a wrong dimension simplex");
            
            auto& simplex = lookup_simplex(t);
            set.insert(t); //add ourself to the set
            simplex.set_label(1);
            
            //iterate up - dim lower than 3
            star_helper_recurse_up(s, t, set);
            //iterate down through the structure in a recursive manner
            star_helper_recurse_down(s, t, set);
        }
        
        /**
         *
         */
        template<typename simplex_key_s>
        void star_helper_recurse_up(simplex_key_s const & s, TetrahedronKey const & t, simplex_set_type & set) {}
        
        /**
         *
         */
        template<typename simplex_key_s, typename simplex_key_t>
        void star_helper_recurse_up(simplex_key_s const & s, simplex_key_t const & t, simplex_set_type & set)
        {
            for(auto co_b_it : *lookup_simplex(t).get_co_boundary())
            {
                //only recurse if it has not been visited previously - no need to check boundary relations as a co-face
                //of t will have s on it's boundary if t has it on it's boundary - and it does
                if(lookup_simplex(co_b_it).get_label() == 0)
                {
                    //recurse through children
                    star_helper(s, co_b_it, set);
                }
            }
        }
        
        
        template<typename N, typename M>
        void star_helper_recurse_down(N const &, M const &, simplex_set_type &)
        { }
        
        
        /**
         *
         */
        void star_helper_recurse_down(EdgeKey const & s,
                                      TetrahedronKey const & t,
                                      simplex_set_type & set)
        {
            star_helper_recurse_down_(s, t, set);
        }
        
        /**
         *
         */
        void star_helper_recurse_down(NodeKey const & s,
                                      TetrahedronKey const & t,
                                      simplex_set_type & set)
        {
            star_helper_recurse_down_(s, t, set);
        }
        
        /**
         *
         */
        void star_helper_recurse_down(NodeKey const & s,
                                      FaceKey const & t,
                                      simplex_set_type & set)
        {
            star_helper_recurse_down_(s, t, set);
        }
        
        /**
         *
         */
        template<typename simplex_key_s, typename simplex_key_t>
        void star_helper_recurse_down_(simplex_key_s const & s, simplex_key_t const & t,  simplex_set_type & set)
        {
            for(auto b_it : *lookup_simplex(t).get_boundary())
            {
                //only recurse if the simplex is a co-face* of s and if it has not been visited previously
                if(lookup_simplex(b_it).get_label() == 0 && in_boundary(s, b_it))
                {
                    //recurse through children
                    star_helper(s, b_it, set);
                }
            }
        }
        
    public:
        
        /**
         *
         */
        t4mesh() {
            m_node_kernel = new node_kernel_type();
            m_edge_kernel = new edge_kernel_type();
            m_face_kernel = new face_kernel_type();
            m_tetrahedron_kernel = new tetrahedron_kernel_type();
        }
        
        /**
         *
         */
        ~t4mesh()
        {
            delete m_tetrahedron_kernel;
            delete m_face_kernel;
            delete m_edge_kernel;
            delete m_node_kernel;
        }
        
        /**
         *
         */
        node_iterator nodes_begin() { return m_node_kernel->begin(); }
        
        /**
         *
         */
        node_iterator nodes_end() { return m_node_kernel->end(); }
        
        /**
         *
         */
        typename edge_kernel_type::iterator edges_begin() { return m_edge_kernel->begin(); }
        
        /**
         *
         */
        typename edge_kernel_type::iterator edges_end() { return m_edge_kernel->end(); }
        
        /**
         *
         */
        typename face_kernel_type::iterator faces_begin() { return m_face_kernel->begin(); }
        
        /**
         *
         */
        typename face_kernel_type::iterator faces_end() { return m_face_kernel->end(); }
        
        /**
         *
         */
        typename tetrahedron_kernel_type::iterator tetrahedra_begin() { return m_tetrahedron_kernel->begin(); }
        
        /**
         *
         */
        typename tetrahedron_kernel_type::iterator tetrahedra_end() { return m_tetrahedron_kernel->end(); }
        
        /**
         * Inserts a node into the mesh. Trivial.
         */
        NodeKey insert_node()
        {
            node_iterator node = m_node_kernel->create();
            return node.key();
        }
        
        /**
         * Inserts an edge into the mesh. Updates the co-boundary of the boundary nodes with the newly created edge.
         * Leaves the closure of the edge in an uncompressed state.
         */
        EdgeKey insert_edge(NodeKey node1, NodeKey node2)
        {
            edge_iterator edge = m_edge_kernel->create();
            //add the new simplex to the co-boundary relation of the boundary simplices
            m_node_kernel->find(node1).add_co_face(edge.key());
            m_node_kernel->find(node2).add_co_face(edge.key());
            //set the boundary relation
            edge->add_face(node1);
            edge->add_face(node2);
            return edge.key();
        }
        
        /**
         * Inserts a face into the mesh. Updates the co-boundary of the boundary faces with the newly created face.
         * Leaves the closure of the face in an uncompressed state.
         */
        FaceKey insert_face(EdgeKey edge1, EdgeKey edge2, EdgeKey edge3)
        {
            face_iterator face = m_face_kernel->create();
            //update relations
            m_edge_kernel->find(edge1).add_co_face(face.key());
            m_edge_kernel->find(edge2).add_co_face(face.key());
            m_edge_kernel->find(edge3).add_co_face(face.key());
            face->add_face(edge1);
            face->add_face(edge2);
            face->add_face(edge3);
            return face.key();
        }
        
        /**
         * Inserts a tetrahedron into the mesh. Updates the co-boundary of the boundary edges with the newly created tetrahedron.
         * Leaves the closure of the tetrahedron in an uncompressed state.
         */
        TetrahedronKey insert_tetrahedron(FaceKey face1, FaceKey face2, FaceKey face3, FaceKey face4)
        {
            auto tetrahedron = m_tetrahedron_kernel->create();
            //update relations
            m_face_kernel->find(face1).add_co_face(tetrahedron.key());
            m_face_kernel->find(face2).add_co_face(tetrahedron.key());
            m_face_kernel->find(face3).add_co_face(tetrahedron.key());
            m_face_kernel->find(face4).add_co_face(tetrahedron.key());
            tetrahedron->add_face(face1);
            tetrahedron->add_face(face2);
            tetrahedron->add_face(face3);
            tetrahedron->add_face(face4);
            return tetrahedron.key();
        }
        
        /**
         *
         */
        void remove(const NodeKey& nid)
        {
            auto& node = lookup_simplex(nid);
            for(auto eid : *node.get_co_boundary())
            {
                lookup_simplex(eid).remove_face(nid);
            }
            m_node_kernel->erase(nid);
        }
        
        /**
         *
         */
        void remove(const EdgeKey& eid)
        {
            auto& edge = lookup_simplex(eid);
            for(auto fid : *edge.get_co_boundary())
            {
                lookup_simplex(fid).remove_face(eid);
            }
            for(auto nid : *edge.get_boundary())
            {
                lookup_simplex(nid).remove_co_face(eid);
            }
            m_edge_kernel->erase(eid);
        }
        
        /**
         *
         */
        void remove(const FaceKey& fid)
        {
            auto& face = lookup_simplex(fid);
            for(auto tid : *face.get_co_boundary())
            {
                lookup_simplex(tid).remove_face(fid);
            }
            for(auto eid : *face.get_boundary())
            {
                lookup_simplex(eid).remove_co_face(fid);
            }
            m_face_kernel->erase(fid);
        }
        
        /**
         *
         */
        void remove(const TetrahedronKey& tid)
        {
            auto& tet = lookup_simplex(tid);
            for(auto fid : *tet.get_boundary())
            {
                lookup_simplex(fid).remove_co_face(tid);
            }
            m_tetrahedron_kernel->erase(tid);
        }
        
        void merge(const NodeKey& key1, const NodeKey& key2)
        {
            auto& simplex = lookup_simplex(key2);
            for(auto cob : *simplex.get_co_boundary())
            {
                lookup_simplex(cob).add_face(key1);
            }
            
            lookup_simplex(key1).merge(simplex);
            remove(key2);
        }
        
        template<typename key_type>
        void merge(const key_type& key1, const key_type& key2)
        {
            auto& simplex = lookup_simplex(key2);
            for(auto cob : *simplex.get_co_boundary())
            {
                lookup_simplex(cob).add_face(key1);
            }
            for(auto bou : *simplex.get_boundary())
            {
                lookup_simplex(bou).add_co_face(key1);
            }
            
            lookup_simplex(key1).merge(simplex);
            remove(key2);
        }
        
        size_type size_nodes() { return m_node_kernel->size(); }
        size_type size_edges() { return m_edge_kernel->size(); }
        size_type size_faces() { return m_face_kernel->size(); }
        size_type size_tetrahedra() { return m_tetrahedron_kernel->size(); }
        size_type size() { return size_nodes() + size_edges() + size_faces() + size_tetrahedra(); }
        
        /**
         * Returns the restricted star of a simplex.
         */
        void star(NodeKey const & key, simplex_set_type & s_set)
        {
            for( auto co_bit : *lookup_simplex(key).get_co_boundary())
            {
                star_helper(key, co_bit, s_set);
            }
            //now that we are done, remember to reset all labels
            //reset_label(simplex);
            //reset for each level - exept tets - so we get the entire component reset
            for (simplex_set_type::edge_set_iterator eit = s_set.edges_begin(); eit != s_set.edges_end(); ++eit)
                reset_label(lookup_simplex(*eit));
            for (simplex_set_type::face_set_iterator fit = s_set.faces_begin(); fit != s_set.faces_end(); ++fit)
                reset_label(lookup_simplex(*fit));
            for (simplex_set_type::tetrahedron_set_iterator tit = s_set.tetrahedra_begin(); tit != s_set.tetrahedra_end(); ++tit)
                lookup_simplex(*tit).reset_label();
        }
        
        
        void star(EdgeKey const & key, simplex_set_type & s_set)
        {
            for(auto co_bit : *lookup_simplex(key).get_co_boundary())
            {
                star_helper(key, co_bit, s_set);
            }
            //now that we are done, remember to reset all labels
            //reset_label(simplex);
            //reset for each level - exept tets - so we get the entire component reset
            for (simplex_set_type::face_set_iterator fit = s_set.faces_begin(); fit != s_set.faces_end(); ++fit)
                reset_label(lookup_simplex(*fit));
            for (simplex_set_type::tetrahedron_set_iterator tit = s_set.tetrahedra_begin(); tit != s_set.tetrahedra_end(); ++tit)
                lookup_simplex(*tit).reset_label();
        }
        
        void star(FaceKey const & f, simplex_set_type & s_set)
        {
            for(auto co_bit : *lookup_simplex(f).get_co_boundary())
            {
                s_set.insert(co_bit);
            }
        }
        
        void star(TetrahedronKey const &, simplex_set_type &)
        { /* do nothing */ }
        
        /**
         * Marek
         */
        void star(simplex_set_type & set, simplex_set_type & result_set)
        {
            simplex_set_type::node_set_iterator nit = set.nodes_begin();
            while (nit != set.nodes_end())
            {
                simplex_set_type st_n;
                star(*nit, st_n);
                result_set.add(st_n);
                ++nit;
            }
            
            simplex_set_type::edge_set_iterator eit = set.edges_begin();
            while (eit != set.edges_end())
            {
                if (!result_set.contains(*eit))
                {
                    simplex_set_type st_e;
                    star(*eit, st_e);
                    result_set.add(st_e);
                }
                ++eit;
            }
            
            simplex_set_type::face_set_iterator fit = set.faces_begin();
            while (fit != set.faces_end())
            {
                if (!result_set.contains(*fit))
                {
                    simplex_set_type st_f;
                    star(*fit, st_f);
                    result_set.add(st_f);
                }
                ++fit;
            }
            
            simplex_set_type::tetrahedron_set_iterator tit = set.tetrahedra_begin();
            while (tit != set.tetrahedra_end())
            {
                result_set.insert(*tit);
                ++tit;
            }
        }
        
        /**
         * Marek
         * Boundary of a simplex.
         */
        template<typename key_type>
        void boundary(key_type const & k, simplex_set_type & result_set)
        {
            boundary_helper(k, result_set);
        }
        
        /**
         * Marek
         * Boundary of a set of tetrahedra.
         */
        void boundary(simplex_set_type & tetrahedra, simplex_set_type & result_set)
        {
            std::map<FaceKey, char> face_occurrences;
            for (auto pit = tetrahedra.tetrahedra_begin(); pit != tetrahedra.tetrahedra_end(); ++pit)
            {
                for (auto bit : *lookup_simplex(*pit).get_boundary())
                {
                    face_occurrences[bit]++;
                }
            }
            
            simplex_set_type boundary_faces;
            
            for (auto foit = face_occurrences.begin(); foit != face_occurrences.end(); ++foit)
            {
                if (foit->second == 1)
                {
                    boundary_faces.insert(foit->first);
                }
            }
            
            closure(boundary_faces, result_set);
        }
        
        /**
         * Marek
         * Closure of a simplex.
         */
        template<typename key_type>
        void closure(key_type const & k, simplex_set_type & result_set)
        {
            boundary_helper(k, result_set);
            result_set.insert(k);
        }
        
        /**
         * Marek.
         * Closure of a simplex set.
         */
        void closure(simplex_set_type & input_set, simplex_set_type & result_set)
        {
            closure_helper(input_set, result_set);
        }
        
        /**
         * Marek
         * Induces consistent orientations on all faces of the tetrahedra with ID tid.
         */
        void orient_faces_consistently(const TetrahedronKey& tid)
        {
            for (auto it : *lookup_simplex(tid).get_boundary())
            {
                orient_face_helper(tid, it, true);
            }
        }
        
        /**
         * Marek
         * Induces opposite orientations on all faces of a simplex sk.
         */
        template<typename key_type>
        void orient_faces_oppositely(key_type const & sk)
        {
            for (auto it : *lookup_simplex(sk).get_boundary())
            {
                orient_face_helper(sk, it, false);
            }
        }
        
        /**
         * Marek
         * Sets the orientation of the simplex sk so that is opposite to the orientation of its boundary face fk.
         */
        template<typename key_type_simplex
        , typename key_type_coface>
        void orient_coface_oppositely(key_type_simplex const & fk, key_type_coface const & sk)
        {
            orient_coface_helper(fk, sk, false);
        }
        
        /**
         * Marek
         * Provided that simplices k1 and k2 are (dimension-1)-adjacent, finds the shared (dimension-1)-simplex k.
         */
        template<typename key_type_simplex
        , typename key_type_face>
        bool get_intersection(key_type_simplex const & k1, key_type_simplex const & k2, key_type_face & k)
        {
            assert(k1 != k2 || !"The same key for both input simplices");
            assert(k1.dim > 0 || !"Cannot intersect two vertices");
            
            for (auto bi1 : *lookup_simplex(k1).get_boundary())
            {
                for (auto bi2 : *lookup_simplex(k2).get_boundary())
                {
                    if (bi1 == bi2)
                    {
                        k = bi1;
                        return true;
                    }
                }
            }
            
            return false;
        }
        
        void link(TetrahedronKey const & k, simplex_set_type & result){}
        
        void link(FaceKey const & f, simplex_set_type & result)
        {
            simplex_set_type st_f, cl_f;
            star(f, st_f);
            closure(st_f, result);
            closure(f, cl_f);
            result.difference(cl_f);
            result.difference(st_f);
            result.clear_edges();
            result.clear_faces();
            result.clear_tetrahedra();
        }
                
        void link(EdgeKey const & e, simplex_set_type & result)
        {
            simplex_set_type st_e, cl_e, temp;
            star(e, st_e);
            closure(st_e, temp);
            closure(e, cl_e);
            NodeKey n1, n2;
            simplex_set_type::node_set_iterator nit = cl_e.nodes_begin();
            n1 = *nit;  ++nit;  n2 = *nit;
            temp.difference(st_e);
            temp.difference(cl_e);
            simplex_set_type::edge_set_iterator eit = temp.edges_begin();
            while (eit != temp.edges_end())
            {
                auto ebnd = lookup_simplex(*eit).get_boundary();
                auto ebit = ebnd->begin();
                if (*ebit == n1 || *ebit == n2)
                {
                    ++eit;
                    continue;
                }
                ++ebit;
                if (*ebit == n1 || *ebit == n2)
                {
                    ++eit;
                    continue;
                }
                result.insert(*eit);
                ++eit;
            }
            temp.clear_edges();
            temp.clear_faces();
            temp.clear_tetrahedra();
            result.add(temp);
        }
        
        void link(NodeKey const & n, simplex_set_type & result)
        {
            simplex_set_type st_n;
            star(n, st_n);
            st_n.insert(n);
            closure(st_n, result);
            result.difference(st_n);
        }
        
        /**
         *
         */
        bool exists(TetrahedronKey const & t)
        {
            return m_tetrahedron_kernel->is_valid(t);
        }
        
        /**
         *
         */
        bool exists(FaceKey const & f)
        {
            return m_face_kernel->is_valid(f);
        }
        
        /**
         *
         */
        bool exists(EdgeKey const & e)
        {
            return m_edge_kernel->is_valid(e);
        }
        
        /**
         *
         */
        bool exists(NodeKey const & n)
        {
            return m_node_kernel->is_valid(n);
        }
        
        /**
         *
         */
        void garbage_collect()
        {
            m_node_kernel->garbage_collect();
            m_edge_kernel->garbage_collect();
            m_face_kernel->garbage_collect();
            m_tetrahedron_kernel->garbage_collect();
        }
        
    };
}
