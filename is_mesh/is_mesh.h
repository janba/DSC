
#pragma once

#include <algorithm>
#include <functional>
#include <queue>
#include <map>
#include <list>

#include <is_mesh/kernel.h>
#include <is_mesh/is_mesh_simplex.h>
#include <is_mesh/simplex_set.h>

namespace is_mesh
{
    
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
     * Simplex typebinding traits
     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
    
    template<typename mesh_type, int dim>
    struct simplex_traits{};
    
    template<typename mesh_type>
    struct simplex_traits<mesh_type, 0>
    {
        typedef typename mesh_type::node_type             simplex_type;
    };
    
    template<typename mesh_type>
    struct simplex_traits<mesh_type, 1>
    {
        typedef typename mesh_type::edge_type             simplex_type;
    };
    
    template<typename mesh_type>
    struct simplex_traits<mesh_type, 2>
    {
        typedef typename mesh_type::face_type             simplex_type;
    };
    
    template<typename mesh_type>
    struct simplex_traits<mesh_type, 3>
    {
        typedef typename mesh_type::tetrahedron_type      simplex_type;
    };
    
    /**
     * A data structure for managing a Simplicial Complex. Based on the work
     * by de Floriani and Hui, the Incidence Simplicial.
     * The complex is specialixed for 3-dimensional simplices, and can only
     * store 0-, 1-, 2-, nad 3-simplices.
     * Simplices are explicitly stored using a different memory kernel for
     * each type.h
     */
    template<
    typename NodeTraits
    , typename TetrahedronTraits
    , typename EdgeTraits
    , typename FaceTraits
    >
    class t4mesh
    {
    public:
        typedef t4mesh<NodeTraits, TetrahedronTraits, EdgeTraits, FaceTraits>        mesh_type;
        
        typedef NodeKey                                                              node_key_type;
        typedef EdgeKey                                                              edge_key_type;
        typedef FaceKey                                                              face_key_type;
        typedef TetrahedronKey                                                       tetrahedron_key_type;
        
        typedef NodeTraits                                                           node_traits;
        typedef TetrahedronTraits                                                    tetrahedron_traits;
        typedef EdgeTraits                                                           edge_traits;
        typedef FaceTraits                                                           face_traits;
        
        typedef Node<NodeTraits, mesh_type>                                         node_type;
        typedef Edge<EdgeTraits, mesh_type>                                         edge_type;
        typedef Face<FaceTraits, mesh_type>                                         face_type;
        typedef Tetrahedron<TetrahedronTraits, mesh_type>                           tetrahedron_type;
        
        typedef          simplex_set<node_key_type, edge_key_type
        , face_key_type, tetrahedron_key_type>            simplex_set_type;
        
    public:
        typedef          kernel<node_type, node_key_type>                            node_kernel_type;
        typedef          kernel<edge_type, edge_key_type>                            edge_kernel_type;
        typedef          kernel<face_type, face_key_type>                            face_kernel_type;
        typedef          kernel<tetrahedron_type, tetrahedron_key_type>              tetrahedron_kernel_type;
        
    public:
        
        typedef typename node_kernel_type::iterator                                  node_iterator;
        typedef typename edge_kernel_type::iterator                                  edge_iterator;
        typedef typename face_kernel_type::iterator                                  face_iterator;
        typedef typename tetrahedron_kernel_type::iterator                           tetrahedron_iterator;
        
        //expose some types from the kernel
        typedef typename node_kernel_type::size_type                                 size_type;
        
    private:
        node_kernel_type*                  m_node_kernel;
        edge_kernel_type*                  m_edge_kernel;
        face_kernel_type*                  m_face_kernel;
        tetrahedron_kernel_type*           m_tetrahedron_kernel;
        
        size_type                          m_uncompressed; //an estimate of the numbers of uncompressed simplices in the mesh
        
    public:
        
        /**
         *
         */
        node_type & lookup_simplex(node_key_type const & k)
        {
            return m_node_kernel->find(k);
        }
        
        /**
         *
         */
        edge_type & lookup_simplex(edge_key_type const & k)
        {
            return m_edge_kernel->find(k);
        }
        
        /**
         *
         */
        face_type & lookup_simplex(face_key_type const & k)
        {
            return m_face_kernel->find(k);
        }
        
        /**
         *
         */
        tetrahedron_type & lookup_simplex(tetrahedron_key_type const & k)
        {
            return m_tetrahedron_kernel->find(k);
        }
        
        /**
         * Helper function for labeling simplices in a connected component.
         * Is doubly recursive with traverse_co_boundary
         * Marek: Will it still work for the non-manifold case? I ain't sure.
         */
        template<typename seed_simplex_key
        , typename simplex_type>
        void label_co_co_bound(seed_simplex_key& ssk, simplex_type& s, int label)
        {
            assert(s.get_label() == 0 || !"traverse_co_co_bound called with simplex already labled");
            s.set_label(label);
            
            for(auto bound_it : *s.get_boundary())
            {
                auto& simplex = lookup_simplex(bound_it);
                if (simplex.get_label() == 0)
                {
                    simplex_set_type s_boundary;
                    boundary(bound_it, s_boundary);
                    if (s_boundary.contains(ssk))
                        label_co_bound(ssk, simplex, label);
                }
            }
        }
        
        /**
         * Helper function for labeling simplices in a connected component.
         * Is doubly recursive with traverse_co_co_boundary
         * Marek: Will it still work for the non-manifold case? I ain't sure.
         */
        template<typename seed_simplex_key
        , typename simplex_type>
        void  label_co_bound(seed_simplex_key& ssk, simplex_type& s, int label)
        {
            assert(s.get_label() == 0 || !"traverse_co_bound called with simplex already labled");
            s.set_label(label);
            
            for(auto co_bound_it : *s.get_co_boundary())
            {
                auto& simplex = lookup_simplex(co_bound_it);
                if (simplex.get_label() == 0)
                {
                    simplex_set_type s_boundary;
                    boundary(co_bound_it, s_boundary);
                    if (s_boundary.contains(ssk))
                        label_co_co_bound(ssk, simplex, label);
                }
            }
        }
        
        /**
         * Marek
         * Helper function for boundary and closure.
         */
        void boundary_helper(tetrahedron_key_type const & k, simplex_set_type & set)
        {
            typename tetrahedron_type::boundary_list s_boundary = lookup_simplex(k).get_boundary();
            typename tetrahedron_type::boundary_iterator it = s_boundary->begin();
            while (it != s_boundary->end())
            {
                set.insert(*it);
                boundary_helper(*it, set);
                ++it;
            }
        }
        
        void boundary_helper(face_key_type const & k, simplex_set_type & set)
        {
            typename face_type::boundary_list s_boundary = lookup_simplex(k).get_boundary();
            typename face_type::boundary_iterator it = s_boundary->begin();
            while (it != s_boundary->end())
            {
                set.insert(*it);
                boundary_helper(*it, set);
                ++it;
            }
        }
        
        void boundary_helper(edge_key_type const & k, simplex_set_type & set)
        {
            typename edge_type::boundary_list s_boundary = lookup_simplex(k).get_boundary();
            typename edge_type::boundary_iterator it = s_boundary->begin();
            while (it != s_boundary->end())
            {
                set.insert(*it);
                ++it;
            }
        }
        
        void boundary_helper(node_key_type const & k, simplex_set_type & set) {}
        
        
        void closure_helper(simplex_set_type & input_set, simplex_set_type & set)
        {
            typename simplex_set_type::tetrahedron_set_iterator tit = input_set.tetrahedra_begin();
            
            while (tit != input_set.tetrahedra_end())
            {
                set.insert(*tit);
                ++tit;
            }
            
            tit = set.tetrahedra_begin();
            while (tit != set.tetrahedra_end())
            {
                typename tetrahedron_type::boundary_list t_boundary = lookup_simplex(*tit).get_boundary();
                typename tetrahedron_type::boundary_iterator it = t_boundary->begin();
                while (it != t_boundary->end())
                {
                    set.insert(*it);
                    ++it;
                }
                ++tit;
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
                typename face_type::boundary_list f_boundary = lookup_simplex(*fit).get_boundary();
                typename face_type::boundary_iterator it = f_boundary->begin();
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
                typename edge_type::boundary_list e_boundary = lookup_simplex(*eit).get_boundary();
                typename edge_type::boundary_iterator it = e_boundary->begin();
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
            if (k.get_dim() == 0) return;
            
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
        void orient_face_helper(const TetrahedronKey& sk, const FaceKey& fk, bool consistently)
        {
            auto simplex_boundary = lookup_simplex(sk).get_boundary();
            auto face_boundary = lookup_simplex(fk).get_boundary();
            
            auto sb_it = simplex_boundary->begin();
            auto fb_it = face_boundary->begin();
            
            std::vector<edge_key_type> new_face_boundary(face_boundary->size());
            
            unsigned char f_index = 0, i = 0;
            
            while (sb_it != simplex_boundary->end())
            {
                if (*sb_it == fk)
                {
                    f_index = i+1;
                }
                else
                {
                    edge_key_type ek;
                    bool res = get_intersection(fk, *sb_it, ek);
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
                invert_orientation(fk);
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
            
            std::vector<node_key_type> new_face_boundary(face_boundary->size());
            
            unsigned char f_index = 0, i = 0;
            
            while (sb_it != simplex_boundary->end())
            {
                if (*sb_it == eid)
                {
                    f_index = i+1;
                }
                else
                {
                    node_key_type ek;
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
            assert(sk.get_dim() == (cfk.get_dim() - 1) || !"sk is not a boundary face of cfk.");
            assert(sk.get_dim() < 3 || !"No simplices of dimension more than three.");
            assert(sk.get_dim() > 0 || !"Vertices are not oriented.");
            
            auto coface_boundary = lookup_simplex(cfk).get_boundary();
            
            std::map<edge_key_type, face_key_type> face_to_simplex;
            
            for (auto cfb_it : *coface_boundary)
            {
                if (cfb_it != sk)
                {
                    edge_key_type k;
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
        void orient_coface_helper(const edge_key_type& sk, const face_key_type& cfk, bool consistently)
        {
            assert(sk.get_dim() == (cfk.get_dim() - 1) || !"sk is not a boundary face of cfk.");
            assert(sk.get_dim() < 3 || !"No simplices of dimension more than three.");
            assert(sk.get_dim() > 0 || !"Vertices are not oriented.");
            
            auto coface_boundary = lookup_simplex(cfk).get_boundary();
            
            std::map<node_key_type, edge_key_type> face_to_simplex;
            
            for (auto cfb_it : *coface_boundary)
            {
                if (cfb_it != sk)
                {
                    node_key_type k;
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
         */
    protected:
        node_key_type split_tetrahedron_helper(const tetrahedron_key_type & t,
                                               std::map<tetrahedron_key_type, tetrahedron_key_type> & new_tets)
        {
            orient_faces_oppositely(t);
            simplex_set_type t_boundary;
            boundary(t, t_boundary);
            
            unsafe_remove(t);
            
            node_key_type n = insert_node();
            lookup_simplex(n).set_compact(true);
            std::map<node_key_type, edge_key_type> node_2_edge_map;
            simplex_set_type::node_set_iterator nit = t_boundary.nodes_begin();
            while (nit != t_boundary.nodes_end())
            {
                node_2_edge_map[*nit] = unsafe_insert_edge(*nit, n);
                lookup_simplex(node_2_edge_map[*nit]).set_compact(true);
                if (nit == t_boundary.nodes_begin())
                {
                    typename node_type::co_boundary_set n_coboundary = lookup_simplex(n).get_co_boundary();
                    n_coboundary->insert(node_2_edge_map[*nit]);
                }
                ++nit;
            }
            std::map<edge_key_type, face_key_type> edge_2_face_map;
            simplex_set_type::edge_set_iterator eit = t_boundary.edges_begin();
            while (eit != t_boundary.edges_end())
            {
                typename edge_type::boundary_list e_boundary = lookup_simplex(*eit).get_boundary();
                assert(e_boundary->size() == 2 || !"Edge boundary corrupted");
                typename edge_type::boundary_iterator ebit = e_boundary->begin();
                node_key_type n1 = *ebit; ++ebit;
                node_key_type n2 = *ebit;
                edge_key_type e1 = node_2_edge_map[n1];
                edge_key_type e2 = node_2_edge_map[n2];
                face_key_type f = unsafe_insert_face(*eit, e1, e2);
                edge_2_face_map[*eit] = f;
                typename edge_type::co_boundary_set e1_coboundary = lookup_simplex(e1).get_co_boundary();
                typename edge_type::co_boundary_set e2_coboundary = lookup_simplex(e2).get_co_boundary();
                if (e1_coboundary->empty()) e1_coboundary->insert(f);
                if (e2_coboundary->empty()) e2_coboundary->insert(f);
                ++eit;
            }
            simplex_set_type::face_set_iterator fit = t_boundary.faces_begin();
            while (fit != t_boundary.faces_end())
            {
                typename face_type::boundary_list f_boundary = lookup_simplex(*fit).get_boundary();
                assert(f_boundary->size() == 3 || !"Face boundary corrupted");
                typename face_type::boundary_iterator fbit = f_boundary->begin();
                edge_key_type e1 = *fbit; ++fbit;
                edge_key_type e2 = *fbit; ++fbit;
                edge_key_type e3 = *fbit;
                face_key_type f1 = edge_2_face_map[e1];
                face_key_type f2 = edge_2_face_map[e2];
                face_key_type f3 = edge_2_face_map[e3];
                tetrahedron_key_type tet = unsafe_insert_tetrahedron(*fit, f1, f2, f3);
                new_tets[tet] = t;
                orient_coface_oppositely(*fit, tet);
                ++fit;
            }
            
            return n;
        }
    public:
        /**
         * Marek
         */
        node_key_type split_face_helper(face_key_type const & f,
                                        std::map<tetrahedron_key_type, tetrahedron_key_type> & new_tets)
        {
            simplex_set_type region, shell;
            star(f, region);
            region.insert(f);
            closure(region, shell);
            shell.difference(region);
            repair_co_boundaries(region, shell);
            
            simplex_set_type::tetrahedron_set_iterator tit = region.tetrahedra_begin();
            std::map<face_key_type, tetrahedron_key_type> face_2_tet_map;
            while (tit != region.tetrahedra_end())
            {
                orient_faces_oppositely(*tit);
                typename tetrahedron_type::boundary_list tbnd = find_tetrahedron(*tit).get_boundary();
                typename tetrahedron_type::boundary_iterator tbit = tbnd->begin();
                while (tbit != tbnd->end())
                {
                    if (*tbit != f)
                        face_2_tet_map[*tbit] = *tit;
                    ++tbit;
                }
                ++tit;
            }
            
            simplex_set_type region_boundary, f_closure;
            closure(f, f_closure);
            boundary(region, region_boundary);
            
            tit = region.tetrahedra_begin();
            while (tit != region.tetrahedra_end())
            {
                unsafe_remove(*tit);
                ++tit;
            }
            unsafe_erase(f);
            
            node_key_type n = insert_node();
            lookup_simplex(n).set_compact(true);
            std::map<node_key_type, edge_key_type> node_2_edge_map;
            simplex_set_type::node_set_iterator nit = region_boundary.nodes_begin();
            while (nit != region_boundary.nodes_end())
            {
                node_2_edge_map[*nit] = unsafe_insert_edge(*nit, n);
                lookup_simplex(node_2_edge_map[*nit]).set_compact(true);
                if (nit == region_boundary.nodes_begin())
                {
                    typename node_type::co_boundary_set n_coboundary = lookup_simplex(n).get_co_boundary();
                    n_coboundary->insert(node_2_edge_map[*nit]);
                }
                ++nit;
            }
            
            std::map<edge_key_type, face_key_type> edge_2_face_map;
            simplex_set_type::edge_set_iterator eit = region_boundary.edges_begin();
            while (eit != region_boundary.edges_end())
            {
                typename edge_type::boundary_list e_boundary = lookup_simplex(*eit).get_boundary();
                assert(e_boundary->size() == 2 || !"Edge boundary corrupted");
                typename edge_type::boundary_iterator ebit = e_boundary->begin();
                node_key_type n1 = *ebit; ++ebit;
                node_key_type n2 = *ebit;
                edge_key_type e1 = node_2_edge_map[n1];
                edge_key_type e2 = node_2_edge_map[n2];
                face_key_type f = unsafe_insert_face(*eit, e1, e2);
                edge_2_face_map[*eit] = f;
                typename edge_type::co_boundary_set e1_coboundary = lookup_simplex(e1).get_co_boundary();
                typename edge_type::co_boundary_set e2_coboundary = lookup_simplex(e2).get_co_boundary();
                if (e1_coboundary->empty()) e1_coboundary->insert(f);
                if (e2_coboundary->empty()) e2_coboundary->insert(f);
                ++eit;
            }
            
            simplex_set_type::face_set_iterator fit = region_boundary.faces_begin();
            
            while (fit != region_boundary.faces_end())
            {
                if (*fit == f)
                {
                    ++fit;
                    continue;
                }
                typename face_type::boundary_list f_boundary = lookup_simplex(*fit).get_boundary();
                assert(f_boundary->size() == 3 || !"Face boundary corrupted");
                typename face_type::boundary_iterator fbit = f_boundary->begin();
                edge_key_type e1 = *fbit; ++fbit;
                edge_key_type e2 = *fbit; ++fbit;
                edge_key_type e3 = *fbit;
                face_key_type f1 = edge_2_face_map[e1];
                face_key_type f2 = edge_2_face_map[e2];
                face_key_type f3 = edge_2_face_map[e3];
                tetrahedron_key_type tet = unsafe_insert_tetrahedron(*fit, f1, f2, f3);
                new_tets[tet] = face_2_tet_map[*fit];
                orient_coface_oppositely(*fit, tet);
                ++fit;
            }
            
            simplex_set_type::edge_set_iterator ceit = f_closure.edges_begin();
            while (ceit != f_closure.edges_end())
            {
                typename edge_type::co_boundary_set e_co_boundary = lookup_simplex(*ceit).get_co_boundary();
                typename edge_type::co_boundary_iterator f_pos = e_co_boundary->find(f);
                if (f_pos != e_co_boundary->end())
                {
                    e_co_boundary->erase(f_pos);
                    e_co_boundary->insert(edge_2_face_map[*ceit]);
                }
                ++ceit;
            }
            
            return n;
        }
        
    public:
        /**
         * Marek
         */
        node_key_type split_edge_helper(const edge_key_type & edge,
                                        std::map<tetrahedron_key_type, tetrahedron_key_type> & new_tets)
        {
            edge_key_type e = edge;
            
            simplex_set_type st_e;
            star(e, st_e);
            st_e.insert(e);
            simplex_set_type shell;
            closure(st_e, shell);
            shell.difference(st_e);
            
            repair_co_boundaries(st_e, shell);
            
            std::map<face_key_type, tetrahedron_key_type> face_2_tet_map;
            simplex_set_type::tetrahedron_set_iterator tit = st_e.tetrahedra_begin();
            while (tit != st_e.tetrahedra_end())
            {
                typename tetrahedron_type::boundary_list tbnd = find_tetrahedron(*tit).get_boundary();
                typename tetrahedron_type::boundary_iterator tbit = tbnd->begin();
                while (tbit != tbnd->end())
                {
                    face_2_tet_map[*tbit] = *tit;
                    ++tbit;
                }
                ++tit;
            }
            
            std::map<edge_key_type, face_key_type> old_edge_2_face_map;
            std::map<edge_key_type, bool> non_link_edge;
            simplex_set_type::face_set_iterator fit = st_e.faces_begin();
            while (fit != st_e.faces_end())
            {
                typename face_type::boundary_list f_boundary = lookup_simplex(*fit).get_boundary();
                typename face_type::boundary_iterator fbit = f_boundary->begin();
                while (fbit != f_boundary->end())
                {
                    if (*fbit != e)
                    {
                        old_edge_2_face_map[*fbit] = *fit;
                        non_link_edge[*fbit] = true;
                    }
                    ++fbit;
                }
                orient_faces_oppositely(*fit);
                ++fit;
            }
            
            tit = st_e.tetrahedra_begin();
            while (tit != st_e.tetrahedra_end())
            {
                orient_faces_oppositely(*tit);
                ++tit;
            }
            
            typename edge_type::boundary_list e_boundary = lookup_simplex(e).get_boundary();
            typename edge_type::boundary_iterator ebit = e_boundary->begin();
            node_key_type n1 = *ebit;
            ++ebit;
            node_key_type n2 = *ebit;
            
            node_key_type n = insert_node();
            find_node(n).set_compact(false);
            std::map<node_key_type, edge_key_type> node_2_edge_map;
            typename node_type::co_boundary_set n_coboundary = lookup_simplex(n).get_co_boundary();
            
            edge_key_type e1 = unsafe_insert_edge(n, n1);
            find_edge(e1).set_compact(false);
            node_2_edge_map[n1] = e1;
            n_coboundary->insert(e1);
            
            edge_key_type e2 = unsafe_insert_edge(n2, n);
            find_edge(e2).set_compact(false);
            node_2_edge_map[n2] = e2;
            n_coboundary->insert(e2);
            
            simplex_set_type::node_set_iterator nit = shell.nodes_begin();
            
            while (nit != shell.nodes_end())
            {
                if ((*nit != n1) && (*nit != n2))
                {
                    edge_key_type new_edge = unsafe_insert_edge(*nit, n);
                    find_edge(new_edge).set_compact(false);
                    node_2_edge_map[*nit] = new_edge;
                    n_coboundary->insert(new_edge);
                }
                ++nit;
            }
            
            std::map<edge_key_type, face_key_type> edge_2_face_map;
            simplex_set_type::edge_set_iterator eit = shell.edges_begin();
            while (eit != shell.edges_end())
            {
                typename edge_type::boundary_list e_boundary = lookup_simplex(*eit).get_boundary();
                assert(e_boundary->size() == 2 || !"Edge boundary corrupted");
                typename edge_type::boundary_iterator ebit = e_boundary->begin();
                node_key_type n1 = *ebit; ++ebit;
                node_key_type n2 = *ebit;
                edge_key_type e1 = node_2_edge_map[n1];
                edge_key_type e2 = node_2_edge_map[n2];
                face_key_type f = unsafe_insert_face(*eit, e1, e2);
                edge_2_face_map[*eit] = f;
                typename edge_type::co_boundary_set e1_coboundary = lookup_simplex(e1).get_co_boundary();
                typename edge_type::co_boundary_set e2_coboundary = lookup_simplex(e2).get_co_boundary();
                e1_coboundary->insert(f);
                e2_coboundary->insert(f);
                typename edge_type::co_boundary_set e_coboundary = lookup_simplex(*eit).get_co_boundary();
                //assert (*(e_coboundary->begin()) != old_edge_2_face_map[*eit]);
                if (non_link_edge[*eit] &&
                    (e_coboundary->find(old_edge_2_face_map[*eit]) != e_coboundary->end()))
                {
                    e_coboundary->erase(old_edge_2_face_map[*eit]);
                    e_coboundary->insert(f);
                }
                orient_coface_oppositely(*eit, f);
                ++eit;
            }
            
            simplex_set_type::face_set_iterator sfit = shell.faces_begin();
            
            while (sfit != shell.faces_end())
            {
                typename face_type::boundary_list f_boundary = lookup_simplex(*sfit).get_boundary();
                assert(f_boundary->size() == 3 || !"Face boundary corrupted");
                typename face_type::boundary_iterator fbit = f_boundary->begin();
                edge_key_type e1 = *fbit; ++fbit;
                edge_key_type e2 = *fbit; ++fbit;
                edge_key_type e3 = *fbit;
                face_key_type f1 = edge_2_face_map[e1];
                face_key_type f2 = edge_2_face_map[e2];
                face_key_type f3 = edge_2_face_map[e3];
                tetrahedron_key_type tet = unsafe_insert_tetrahedron(*sfit, f1, f2, f3);
                new_tets[tet] = face_2_tet_map[*sfit];
                orient_coface_oppositely(*sfit, tet);
                ++sfit;
            }
            
            /******* unsafe_remove(e); *******/
            tit = st_e.tetrahedra_begin();
            while (tit != st_e.tetrahedra_end())
            {
                unsafe_remove(*tit);
                ++tit;
            }
            fit = st_e.faces_begin();
            while (fit != st_e.faces_end())
            {
                unsafe_erase(*fit);
                ++fit;
            }
            unsafe_erase(e);
            /**********************************/
            
            simplex_set_type starn;
            star(n, starn);
            starn.insert(n);
            compress(starn);
            
            return n;
        }
        
    private:
        
        /**
         * Marek
         */
        void repair_co_boundaries(simplex_set_type & interior, simplex_set_type & boundary)
        {
            std::map<node_key_type, char> node_repaired;
            std::map<edge_key_type, char> edge_repaired;
            simplex_set_type::face_set_iterator fit = boundary.faces_begin();
            while (fit != boundary.faces_end())
            {
                typename face_type::boundary_list f_boundary = lookup_simplex(*fit).get_boundary();
                typename face_type::boundary_iterator fbit = f_boundary->begin();
                while (fbit != f_boundary->end())
                {
                    if (edge_repaired[*fbit] == 0)
                    {
                        typename edge_type::co_boundary_set e_co_boundary = lookup_simplex(*fbit).get_co_boundary();
                        typename edge_type::co_boundary_iterator ecit = e_co_boundary->begin();
                        while (ecit != e_co_boundary->end())
                        {
                            if (interior.contains(*ecit))
                            {
                                e_co_boundary->insert(*fit);
                                e_co_boundary->erase(*ecit);
                                break;
                            }
                            ++ecit;
                        }
                        edge_repaired[*fbit] = 1;
                    }
                    ++fbit;
                }
                ++fit;
            }
            
            simplex_set_type::edge_set_iterator eit = boundary.edges_begin();
            while (eit != boundary.edges_end())
            {
                typename edge_type::boundary_list e_boundary = lookup_simplex(*eit).get_boundary();
                typename edge_type::boundary_iterator ebit = e_boundary->begin();
                while (ebit != e_boundary->end())
                {
                    if (node_repaired[*ebit] == 0)
                    {
                        typename node_type::co_boundary_set n_co_boundary = lookup_simplex(*ebit).get_co_boundary();
                        typename node_type::co_boundary_iterator ncit = n_co_boundary->begin();
                        while (ncit != n_co_boundary->end())
                        {
                            if (interior.contains(*ncit))
                            {
                                n_co_boundary->insert(*eit);
                                n_co_boundary->erase(*ncit);
                                break;
                            }
                            ++ncit;
                        }
                        node_repaired[*ebit] = 1;
                    }
                    ++ebit;
                }
                ++eit;
            }
        }
        
        /**
         * Marek
         * so far for the manifold edge ONLY!
         */
        bool edge_collapse_precond(edge_key_type & e,
                                   node_key_type const & n1,
                                   node_key_type const & n2)
        {
            simplex_set_type lk_e, lk1, lk12;
            link(e, lk_e);
            link(n1, lk1);
            link(n2, lk12);
            lk12.intersection(lk1);
            lk12.difference(lk_e);
            if (lk12.size_nodes() == 0 &&
                lk12.size_edges() == 0 &&
                lk12.size_faces() == 0)
                return true;
            return false;
        }
        
    public:
        /**
         * Marek
         * so far for the manifold edge ONLY!
         */
        node_key_type edge_collapse_helper(edge_key_type & e,
                                           node_key_type const & n1,
                                           node_key_type const & n2)
        {
            if (!edge_collapse_precond(e,n1,n2))
                return NodeKey();
            
            simplex_set_type st2, lk1, lk2, st_e, lk_e, st_e_boundary;
            star(n2, st2);
            link(n1, lk1);
            link(n2, lk2);
            star(e, st_e);
            st_e.insert(e);
            link(e, lk_e);
            boundary(st_e, st_e_boundary);
            st_e_boundary.difference(st_e);
            
            repair_co_boundaries(st_e, st_e_boundary);
            
            std::map<edge_key_type, edge_key_type> edge_2_edge_map;
            std::map<face_key_type, face_key_type> face_2_face_map;
            
            edge_collapse_clear_interior(n1, n2, e, st_e, st_e_boundary, lk1, lk2, edge_2_edge_map, face_2_face_map);
            edge_collapse_sew_hole_up(n1, n2, st2, st_e, lk1, edge_2_edge_map, face_2_face_map);
            
            return n1;
        }
        
    private:
        /**
         * Marek
         */
        void edge_collapse_clear_interior(node_key_type const & n1,
                                          node_key_type const & n2,
                                          edge_key_type & e,
                                          simplex_set_type & st_e,
                                          simplex_set_type & st_e_boundary,
                                          simplex_set_type & lk1,
                                          simplex_set_type & lk2,
                                          std::map<edge_key_type, edge_key_type> & edge_2_edge_map,
                                          std::map<face_key_type, face_key_type> & face_2_face_map)
        {
            simplex_set_type::tetrahedron_set_iterator tit = st_e.tetrahedra_begin();
            while (tit != st_e.tetrahedra_end())
            {
                face_key_type f1, f2;
                typename tetrahedron_type::boundary_list t_boundary = lookup_simplex(*tit).get_boundary();
                typename tetrahedron_type::boundary_iterator tbit = t_boundary->begin();
                while (tbit != t_boundary->end())
                {
                    if (!st_e.contains(*tbit))
                    {
                        if (lk1.contains(*tbit))
                        {
                            f2 = *tbit;
                        }
                        if (lk2.contains(*tbit))
                        {
                            f1 = *tbit;
                        }
                    }
                    ++tbit;
                }
                face_2_face_map[f2] = f1;
                unsafe_remove(*tit);
                ++tit;
            }
            
            simplex_set_type::face_set_iterator fit = st_e.faces_begin();
            while (fit != st_e.faces_end())
            {
                edge_key_type e1, e2;
                typename face_type::boundary_list f_boundary = lookup_simplex(*fit).get_boundary();
                typename face_type::boundary_iterator fbit = f_boundary->begin();
                while (fbit != f_boundary->end())
                {
                    if (*fbit != e)
                    {
                        if (lk1.contains(*fbit))
                        {
                            e2 = *fbit;
                        }
                        if (lk2.contains(*fbit))
                        {
                            e1 = *fbit;
                        }
                    }
                    ++fbit;
                }
                edge_2_edge_map[e2] = e1;
                unsafe_erase(*fit);
                ++fit;
            }
            
            unsafe_erase(e);
        }
        
        /**
         * Marek
         */
        void edge_collapse_sew_hole_up(node_key_type const & n1,
                                       node_key_type const & n2,
                                       simplex_set_type & st2,
                                       simplex_set_type & st_e,
                                       simplex_set_type & lk1,
                                       std::map<edge_key_type, edge_key_type> & edge_2_edge_map,
                                       std::map<face_key_type, face_key_type> & face_2_face_map)
        {
            simplex_set_type to_be_removed;
            std::map<edge_key_type, edge_key_type>::iterator meeit = edge_2_edge_map.begin();
            while (meeit != edge_2_edge_map.end())
            {
                to_be_removed.insert(meeit->first);
                ++meeit;
            }
            std::map<face_key_type, face_key_type>::iterator mffit = face_2_face_map.begin();
            while (mffit != face_2_face_map.end())
            {
                to_be_removed.insert(mffit->first);
                ++mffit;
            }
            
            simplex_set_type::edge_set_iterator eit = st2.edges_begin();
            while (eit != st2.edges_end())
            {
                if ((!st_e.contains(*eit)) && (!to_be_removed.contains(*eit)))
                {
                    typename edge_type::boundary_list e_boundary = lookup_simplex(*eit).get_boundary();
                    typename edge_type::boundary_iterator ebit = e_boundary->begin();
                    while (ebit != e_boundary->end())
                    {
                        if (*ebit == n2)
                        {
                            *ebit = n1;
                            break;
                        }
                        ++ebit;
                    }
                }
                ++eit;
            }
            
            simplex_set_type::face_set_iterator fit = st2.faces_begin();
            while (fit != st2.faces_end())
            {
                if ((!st_e.contains(*fit)) && (!to_be_removed.contains(*fit)))
                {
                    typename face_type::boundary_list f_boundary = lookup_simplex(*fit).get_boundary();
                    typename face_type::boundary_iterator fbit = f_boundary->begin();
                    while (fbit != f_boundary->end())
                    {
                        if (to_be_removed.contains(*fbit))
                        {
                            *fbit = edge_2_edge_map[*fbit];
                            //break;
                        }
                        ++fbit;
                    }
                }
                ++fit;
            }
            
            simplex_set_type::tetrahedron_set_iterator tit = st2.tetrahedra_begin();
            while (tit != st2.tetrahedra_end())
            {
                if (!st_e.contains(*tit))
                {
                    typename tetrahedron_type::boundary_list t_boundary = lookup_simplex(*tit).get_boundary();
                    typename tetrahedron_type::boundary_iterator tbit = t_boundary->begin();
                    while (tbit != t_boundary->end())
                    {
                        if (to_be_removed.contains(*tbit))
                        {
                            face_key_type f = face_2_face_map[*tbit];
                            *tbit = f;
                            lookup_simplex(f).get_co_boundary()->insert(*tit);
                            //break;
                        }
                        ++tbit;
                    }
                }
                ++tit;
            }
            
            std::map<face_key_type, face_key_type>::iterator ffit = face_2_face_map.begin();
            while (ffit != face_2_face_map.end())
            {
                typename face_type::boundary_list f_boundary = lookup_simplex(ffit->first).get_boundary();
                typename face_type::boundary_iterator fbit = f_boundary->begin();
                while (fbit != f_boundary->end())
                {
                    typename edge_type::co_boundary_set e_coboundary = lookup_simplex(*fbit).get_co_boundary();
                    assert (e_coboundary->size() == 1);
                    typename edge_type::co_boundary_iterator ecit = e_coboundary->begin();
                    typename edge_type::co_boundary_set new_coboundary = new typename edge_type::set_type;
                    while (ecit != e_coboundary->end())
                    {
                        if (*ecit == ffit->first)
                        {
                            new_coboundary->insert(ffit->second);
                            //*ecit = ffit->second;
                            //break;
                        }
                        else
                        {
                            new_coboundary->insert(*ecit);
                        }
                        ++ecit;
                    }
                    lookup_simplex(*fbit).get_co_boundary()->clear();
                    typename edge_type::co_boundary_iterator ncbit = new_coboundary->begin();
                    while (ncbit != new_coboundary->end())
                    {
                        lookup_simplex(*fbit).get_co_boundary()->insert(*ncbit);
                        ++ncbit;
                    }
                    delete new_coboundary;
                    ++fbit;
                }
                unsafe_erase(ffit->first);
                ++ffit;
            }
            
            std::map<edge_key_type, edge_key_type>::iterator eeit = edge_2_edge_map.begin();
            while (eeit != edge_2_edge_map.end())
            {
                typename edge_type::boundary_list e_boundary = lookup_simplex(eeit->first).get_boundary();
                typename edge_type::boundary_iterator ebit = e_boundary->begin();
                while (ebit != e_boundary->end())
                {
                    typename node_type::co_boundary_set n_coboundary = lookup_simplex(*ebit).get_co_boundary();
                    assert (n_coboundary->size() == 1);
                    typename node_type::co_boundary_iterator ncit = n_coboundary->begin();
                    typename node_type::co_boundary_set new_coboundary = new typename node_type::set_type;
                    while (ncit != n_coboundary->end())
                    {
                        if (*ncit == eeit->first)
                        {
                            new_coboundary->insert(eeit->second);
                            //*ncit = eeit->second;
                            //break;
                        }
                        else
                            new_coboundary->insert(*ncit);
                        ++ncit;
                    }
                    find_node(*ebit).get_co_boundary()->clear();
                    typename node_type::co_boundary_iterator ncbit = new_coboundary->begin();
                    while (ncbit != new_coboundary->end())
                    {
                        find_node(*ebit).get_co_boundary()->insert(*ncbit);
                        ++ncbit;
                    }
                    delete new_coboundary;
                    ++ebit;
                }
                unsafe_erase(eeit->first);
                ++eeit;
            }
            
            unsafe_erase(n2);
        }
        
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
            if (f.get_dim() >= s.get_dim())
            {
                //cannot be in boundary as dimensions are wrong
                return false;
            }
            auto b_list = lookup_simplex(s).get_boundary();
            if (f.get_dim()+1 < s.get_dim())
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
        bool in_boundary(key_type_face const & f, node_key_type const & s)
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
            assert(s.get_dim() < t.get_dim() || !"Star traversed a wrong dimension simplex");
            
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
        void star_helper_recurse_up(simplex_key_s const & s, tetrahedron_key_type const & t, simplex_set_type & set) {}
        
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
        void star_helper_recurse_down(edge_key_type const & s,
                                      tetrahedron_key_type const & t,
                                      simplex_set_type & set)
        {
            star_helper_recurse_down_(s, t, set);
        }
        
        /**
         *
         */
        void star_helper_recurse_down(node_key_type const & s,
                                      tetrahedron_key_type const & t,
                                      simplex_set_type & set)
        {
            star_helper_recurse_down_(s, t, set);
        }
        
        /**
         *
         */
        void star_helper_recurse_down(node_key_type const & s,
                                      face_key_type const & t,
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
        
        /**
         * Inserts an edge into the mesh without updating any of the boundary nodes co-boundary relation.
         * Might leave the mesh in an inconsistent state.
         */
        edge_key_type unsafe_insert_edge(node_key_type node1, node_key_type node2)
        {
            edge_iterator edge = m_edge_kernel->create();
            //edge->set_compact(true);
            //the nodes might be uncompressed
            edge->add_face(node1);
            edge->add_face(node2);
            return edge.key();
        }
        
        /**
         * Inserts a face into the mesh without updating any of the boundary edges co-boundary relation
         * Might leave the mesh in an inconsistent state.
         */
        face_key_type unsafe_insert_face(edge_key_type edge1, edge_key_type edge2, edge_key_type edge3)
        {
            face_iterator face = m_face_kernel->create();
            face->add_face(edge1);
            face->add_face(edge2);
            face->add_face(edge3);
            return face.key();
        }
        
        /**
         * Inserts a tetrahedron into the mesh.
         * Updates boundary relation both ways, but don't uncompresses anything.
         * Might leave the mesh in an inconsistent state.
         */
        tetrahedron_key_type unsafe_insert_tetrahedron(face_key_type face1, face_key_type face2, face_key_type face3, face_key_type face4)
        {
            tetrahedron_iterator tetrahedron = m_tetrahedron_kernel->create();
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
        void unsafe_erase(tetrahedron_key_type const & key)
        {
            m_tetrahedron_kernel->erase(key);
        }
        
        /**
         *
         */
        void unsafe_erase(face_key_type const & key)
        {
            m_face_kernel->erase(key);
        }
        
        /**
         *
         */
        void unsafe_erase(edge_key_type const & key)
        {
            m_edge_kernel->erase(key);
        }
        
        /**
         *
         */
        void unsafe_erase(node_key_type const & key)
        {
            m_node_kernel->erase(key);
        }
        
    public:
        /**
         *
         */
        void unsafe_remove(tetrahedron_key_type const & key)
        {
            tetrahedron_type& tet = lookup_simplex(key);
            typename tetrahedron_type::boundary_iterator itr = tet.get_boundary()->begin();
            for( ; itr != tet.get_boundary()->end(); ++itr)
            {
                lookup_simplex(*itr).remove_co_face(key);
            }
            m_tetrahedron_kernel->erase(key);
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
            m_uncompressed = 0;
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
        void clear()
        {
            m_node_kernel->clear();
            m_edge_kernel->clear();
            m_face_kernel->clear();
            m_tetrahedron_kernel->clear();
        }
        
        /**
         *
         */
        node_type & find_node(const node_key_type k) { return m_node_kernel->find(k); }
        
        /**
         *
         */
        edge_type & find_edge(const edge_key_type k) { return m_edge_kernel->find(k); }
        
        /**
         *
         */
        face_type & find_face(const face_key_type k) { return m_face_kernel->find(k); }
        
        /**
         *
         */
        tetrahedron_type & find_tetrahedron(const tetrahedron_key_type k) { return m_tetrahedron_kernel->find(k); }
        
        /**
         *
         */
        typename node_kernel_type::iterator nodes_begin() { return m_node_kernel->begin(); }
        
        /**
         *
         */
        typename node_kernel_type::iterator nodes_end() { return m_node_kernel->end(); }
        
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
        node_key_type insert_node(bool is_compact =false)
        {
            node_iterator node = m_node_kernel->create();
            node->set_compact(is_compact);
            return node.key();
        }
        
        /**
         * Inserts an edge into the mesh. Updates the co-boundary of the boundary nodes with the newly created edge.
         * Leaves the closure of the edge in an uncompressed state.
         */
        edge_key_type insert_edge(node_key_type node1, node_key_type node2, bool is_compact =false)
        {
            edge_iterator edge = m_edge_kernel->create();
            //first uncompress boundary
            uncompress(node1);
            uncompress(node2);
            //add the new simplex to the co-boundary relation of the boundary simplices
            m_node_kernel->find(node1).add_co_face(edge.key());
            m_node_kernel->find(node2).add_co_face(edge.key());
            //set the boundary relation
            edge->add_face(node1);
            edge->add_face(node2);
            edge->set_compact(is_compact);
            return edge.key();
        }
        
        /**
         * Inserts a face into the mesh. Updates the co-boundary of the boundary faces with the newly created face.
         * Leaves the closure of the face in an uncompressed state.
         */
        face_key_type insert_face(edge_key_type edge1, edge_key_type edge2, edge_key_type edge3)
        {
            face_iterator face = m_face_kernel->create();
            //first uncompress full boundary - fast if allready uncompressed, else might be heavy
            simplex_set_type set;
            closure(edge1, set);
            closure(edge2, set);
            closure(edge3, set);
            uncompress(set);
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
        tetrahedron_key_type insert_tetrahedron(face_key_type face1, face_key_type face2, face_key_type face3, face_key_type face4)
        {
            tetrahedron_iterator tetrahedron = m_tetrahedron_kernel->create();
            //first uncompress full boundary - fast if allready uncompressed, else might be heavy
            simplex_set_type set;
            closure(face1, set);
            closure(face2, set);
            closure(face3, set);
            closure(face4, set);
            uncompress(set);
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
        
        size_type size_nodes() { return m_node_kernel->size(); }
        size_type size_edges() { return m_edge_kernel->size(); }
        size_type size_faces() { return m_face_kernel->size(); }
        size_type size_tetrahedra() { return m_tetrahedron_kernel->size(); }
        size_type size() { return size_nodes() + size_edges() + size_faces() + size_tetrahedra(); }
        
        /**
         * Makes sure that the mesh is in compressed format, according to the
         * principles behind the Incidence Simplicial data structure, as described
         * by Hui and de Floriani.
         */
        void compress_all()
        {
            //call compress with a new simplex set containing the entire mesh
            //only compress if more than 10% is uncompressed. Only nodes and edges can be uncompressed
            if ( (m_uncompressed * 5) > (size_nodes() + size_edges())  )
            {
                std::cout << "Compressing all..." << std::endl;
                simplex_set_type set;
                for (node_iterator iter = m_node_kernel->begin(); iter != m_node_kernel->end(); ++iter)
                    set.insert(iter.key());
                for (edge_iterator iter = m_edge_kernel->begin(); iter != m_edge_kernel->end(); ++iter)
                    set.insert(iter.key());
                for (face_iterator iter = m_face_kernel->begin(); iter != m_face_kernel->end(); ++iter)
                    set.insert(iter.key());
                for (tetrahedron_iterator iter = m_tetrahedron_kernel->begin(); iter != m_tetrahedron_kernel->end(); ++iter)
                    set.insert(iter.key());
                compress(set);
            }
        }
        
        /**
         * Returns the restricted star of a simplex.
         */
        void star(node_key_type const & key, simplex_set_type & s_set)
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
        
        
        void star(edge_key_type const & key, simplex_set_type & s_set)
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
        
        void star(face_key_type const & f, simplex_set_type & s_set)
        {
            face_type& simplex = lookup_simplex(f);
            typename face_type::co_boundary_set co_boundary = simplex.get_co_boundary();
            typename face_type::co_boundary_iterator co_bit = co_boundary->begin();
            for( ; co_bit != co_boundary->end() ; ++co_bit )
            {
                s_set.insert(*co_bit);
            }
        }
        
        void star(tetrahedron_key_type const &, simplex_set_type &)
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
         * Initial compress - might need to be refactored, a lot...
         * Only compresses a node that is non-compact.
         */
        void compress(simplex_set_type & s)
        {
            //by definition tetrahedra and faces are already "compressed".. so only compress nodes and edges.
            //we handle nodes first to not destroy the full coboudary of edges while we handle nodes
            typename simplex_set_type::node_set_iterator node_itr = s.nodes_begin();
            for( ; node_itr != s.nodes_end() ; ++node_itr)
            {
                //first we need to label all edges for each connected component in the star of the node
                node_type& node = m_node_kernel->find(*node_itr);
                if (node.is_compact()) continue;
                node.set_compact(true);
                typename node_type::co_boundary_set cob_set = node.get_co_boundary();
                std::vector<edge_key_type> edge_vec(cob_set->size());
                std::copy(cob_set->begin(), cob_set->end(), edge_vec.begin());
                typename std::vector<edge_key_type>::iterator edge_itr = edge_vec.begin();
                //co_boundary_iterator edge_itr = cob_set->begin();
                int label = 1;
                for( ;edge_itr != edge_vec.end() ; ++edge_itr)
                {
                    //find an iterator to the actual edge
                    edge_type& edge = m_edge_kernel->find(*edge_itr);
                    //only proceed if we havn't assigned a label yet
                    if (edge.get_label() == 0)
                    {
                        //enter labeling routine - after this call the entire connected component has been labeled.
                        label_co_bound(*node_itr, edge, label);
                        ++label; //used that label, increase to next...
                    }
                }
                
                //with all the components now labeled, we can proceed to compress the co-boundary relation
                edge_itr = edge_vec.begin(); //reset iterator
                label = 1;               //reset label
                node.get_co_boundary()->clear(); //reset the co-boundary
                //labels should come in increasing order starting from 1
                for( ;edge_itr != edge_vec.end() ; ++edge_itr)
                {
                    edge_type& edge = m_edge_kernel->find(*edge_itr);
                    assert(edge.get_label() <= label || !"Label ordering is wacked while compressing nodes");
                    if (edge.get_label() == label)
                    {
                        node.add_co_face(*edge_itr);
                        ++label;
                    }
                    //reset labels
                    reset_label(edge);
                }
                //done compressing that node's co_boundary.. next node...
                m_uncompressed = 0; //reset counting label
            }
            
            //now for edges... more or less same procedure - could be generalized into function, but...
            typename simplex_set_type::edge_set_iterator edge_itr = s.edges_begin();
            for( ; edge_itr != s.edges_end() ; ++edge_itr)
            {
                //first we need to label all edges for each connected component in the star of the node
                edge_type& edge = m_edge_kernel->find(*edge_itr);
                if (edge.is_compact()) continue;
                edge.set_compact(true);
                typename edge_type::co_boundary_set cob_set = edge.get_co_boundary();
                std::vector<face_key_type> face_vec(cob_set->size());
                std::copy(cob_set->begin(), cob_set->end(), face_vec.begin());
                typename std::vector<face_key_type>::iterator face_itr = face_vec.begin();
                //co_boundary_iterator edge_itr = cob_set->begin();
                int label = 1;
                for( ;face_itr != face_vec.end() ; ++face_itr)
                {
                    //find an iterator to the actual edge
                    face_type& face = m_face_kernel->find(*face_itr);
                    //only proceed if we havn't assigned a label yet
                    if (face.get_label() == 0)
                    {
                        //enter labeling routine - after this call the entire connected component has been labeled.
                        label_co_bound(*edge_itr, face, label);
                        ++label; //used that label, increase to next...
                    }
                }
                
                //with all the components now labeled, we can proceed to compress the co-boundary relation
                face_itr = face_vec.begin(); //reset iterator
                label = 1;               //reset label
                edge.get_co_boundary()->clear(); //reset the co-boundary
                //labels should come in increasing order starting from 1
                for( ;face_itr != face_vec.end() ; ++face_itr)
                {
                    face_type& face = m_face_kernel->find(*face_itr);
                    assert(face.get_label() <= label || !"Label ordering is wacked while compressing nodes");
                    if (face.get_label() == label)
                    {
                        edge.add_co_face(*face_itr);
                        ++label;
                    }
                    reset_label(face);
                }
                //done compressing that node's co_boundary.. next node...
            }
        } //compress(simplex_set_type)
        
        
        void uncompress(const tetrahedron_key_type & t) {}
        void uncompress(const face_key_type & f) {}
        
        //// NOT TESTED!!!!!
        void uncompress(const edge_key_type & edge_k)
        {
            //assert(0);
            ++m_uncompressed;
            edge_type& edge = lookup_simplex(edge_k);
            if (!edge.is_compact()) return; //only uncompress compressed nodes.
            edge.set_compact(false);
            simplex_set_type star_set;
            star(edge_k, star_set);
            typename simplex_set_type::face_set_iterator face_itr = star_set.faces_begin();
            for( ; face_itr != star_set.faces_end(); ++face_itr)
            {
                if(in_boundary(edge_k, *face_itr))
                    edge.add_co_face(*face_itr);
            }
        }
        
        //// NOT TESTED!!!!!!
        void uncompress(const node_key_type & node_k)
        {
            //assert(0);
            ++m_uncompressed;
            node_type& node = lookup_simplex(node_k);
            if (!node.is_compact()) return; //only uncompress compressed nodes.
            node.set_compact(false);
            simplex_set_type star_set;
            star(node_k, star_set);
            typename simplex_set_type::edge_set_iterator edge_itr = star_set.edges_begin();
            for( ; edge_itr != star_set.edges_end(); ++edge_itr)
            {
                if(in_boundary(node_k, *edge_itr))
                    node.add_co_face(*edge_itr);
            }
        }
        
        
        void uncompress(simplex_set_type & s)
        {
            //first uncompress edges
            for (typename simplex_set_type::edge_set_iterator edge_itr = s.edges_begin(); edge_itr != s.edges_end() ; ++edge_itr )
            {
                uncompress(*edge_itr);
            }
            //uncompress nodes
            for (typename simplex_set_type::node_set_iterator node_itr = s.nodes_begin(); node_itr != s.nodes_end() ; ++node_itr )
            {
                uncompress(*node_itr);
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
            std::map<face_key_type, char> face_occurrences;
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
        void orient_faces_consistently(const tetrahedron_key_type& tid)
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
            assert(k1.get_dim() > 0 || !"Cannot intersect two vertices");
            
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
        
        /**
         * Marek
         */
        node_key_type split_tetrahedron(const tetrahedron_key_type & t)
        {
            std::map<tetrahedron_key_type, tetrahedron_key_type> new_tets;
            return split_tetrahedron_helper(t, new_tets);
        }
        
        /**
         * Marek
         */
        node_key_type split_face(face_key_type & f)
        {
            std::map<tetrahedron_key_type, tetrahedron_key_type> new_tets;
            return split_face_helper(f, new_tets);
        }
        
        /**
         * Marek
         */
        node_key_type split_edge(edge_key_type & e)
        {
            std::map<tetrahedron_key_type, tetrahedron_key_type> new_tets;
            return split_edge_helper(e, new_tets);
        }
        
        /**
         * Marek
         */
        node_key_type edge_collapse(edge_key_type & e)
        {
            typename edge_type::boundary_list e_boundary = lookup_simplex(e).get_boundary();
            typename edge_type::boundary_iterator ebit = e_boundary->begin();
            node_key_type n1 = *ebit;
            ++ebit;
            node_key_type n2 = *ebit;
            return edge_collapse_helper(e, n1, n2);
        }
        
        void link(tetrahedron_key_type const & k, simplex_set_type & result){}
        
        void link(face_key_type const & f, simplex_set_type & result)
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
                
        void link(edge_key_type const & e, simplex_set_type & result)
        {
            simplex_set_type st_e, cl_e, temp;
            star(e, st_e);
            closure(st_e, temp);
            closure(e, cl_e);
            node_key_type n1, n2;
            simplex_set_type::node_set_iterator nit = cl_e.nodes_begin();
            n1 = *nit;  ++nit;  n2 = *nit;
            temp.difference(st_e);
            temp.difference(cl_e);
            simplex_set_type::edge_set_iterator eit = temp.edges_begin();
            while (eit != temp.edges_end())
            {
                typename edge_type::boundary_list ebnd = find_edge(*eit).get_boundary();
                typename edge_type::boundary_iterator ebit = ebnd->begin();
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
        
        void link(node_key_type const & n, simplex_set_type & result)
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
        bool is_boundary(tetrahedron_key_type & t) { return false; }
        
        /**
         *
         */
        bool is_boundary(face_key_type & f)
        {
            typename face_type::co_boundary_set f_coboundary = lookup_simplex(f).get_co_boundary();
            if (f_coboundary->size() == 2) return false;
            return true;
        }
        
        /**
         *
         */
        bool is_boundary(edge_key_type & e)
        {
            simplex_set_type ste;
            star(e, ste);
            simplex_set_type::face_set_iterator fit = ste.faces_begin();
            if (ste.faces_begin() == ste.faces_end()) return true;
            while (fit != ste.faces_end())
            {
                if (is_boundary(*fit)) return true;
                ++fit;
            }
            return false;
        }
        
        /**
         *
         */
        bool is_boundary(node_key_type & n)
        {
            simplex_set_type stn;
            star(n, stn);
            simplex_set_type::face_set_iterator fit = stn.faces_begin();
            if (stn.faces_begin() == stn.faces_end()) return true;
            while (fit != stn.faces_end())
            {
                if (is_boundary(*fit)) return true;
                ++fit;
            }
            return false;
        }
        
        /**
         *
         */
        bool exists(tetrahedron_key_type const & t)
        {
            return m_tetrahedron_kernel->is_valid(t);
        }
        
        /**
         *
         */
        bool exists(face_key_type const & f)
        {
            return m_face_kernel->is_valid(f);
        }
        
        /**
         *
         */
        bool exists(edge_key_type const & e)
        {
            return m_edge_kernel->is_valid(e);
        }
        
        /**
         *
         */
        bool exists(node_key_type const & n)
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
