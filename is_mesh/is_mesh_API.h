//
//  is_mesh_API.h
//  DSC
//
//  Created by Asger Nyman Christiansen on 6/14/13.
//  Copyright (c) 2013 Asger Nyman Christiansen. All rights reserved.
//

#ifndef DSC_is_mesh_API_h
#define DSC_is_mesh_API_h

#include "simplicial_complex.h"
#include <is_mesh/is_mesh.h>
#include <is_mesh/io/is_mesh_lists_read.h>

#include "attributes.h"

template <typename MT>
class ISMesh
{
    typedef typename OpenTissue::is_mesh::t4mesh< default_node_traits<MT>, default_tetrahedron_traits<MT>, default_edge_traits<MT>, default_face_traits<MT> > Mesh;
    
public:
    typedef typename Mesh::node_key_type    node_key_type;
    typedef typename Mesh::edge_key_type    edge_key_type;
    typedef typename Mesh::face_key_type    face_key_type;
    typedef typename Mesh::tetrahedron_key_type tetrahedron_key_type;
    
    typedef typename Mesh::simplex_set_type simplex_set;
    
    
    const node_key_type NULL_NODE;
    const edge_key_type NULL_EDGE;
    const face_key_type NULL_FACE;
    const tetrahedron_key_type NULL_TETRAHEDRON;
    
private:
    Mesh mesh;
    
public:
    ISMesh(std::vector<typename MT::real_type> & points, std::vector<int> & tets): NULL_NODE(-1), NULL_EDGE(-1), NULL_FACE(-1), NULL_TETRAHEDRON(-1)
    {
        vectors_read(points, tets, mesh);
        
        check_validity();
        check_validity();
    }
    
    ///////////////
    // ITERATIRS //
    ///////////////
    
    typename Mesh::node_iterator nodes_begin()
    {
        return mesh.nodes_begin();
    }
    
    typename Mesh::node_iterator nodes_end()
    {
        return mesh.nodes_end();
    }
    
    typename Mesh::edge_iterator edges_begin()
    {
        return mesh.edges_begin();
    }
    
    typename Mesh::edge_iterator edges_end()
    {
        return mesh.edges_end();
    }
    
    typename Mesh::face_iterator faces_begin()
    {
        return mesh.faces_begin();
    }
    
    typename Mesh::face_iterator faces_end()
    {
        return mesh.faces_end();
    }
    
    typename Mesh::tetrahedron_iterator tetrahedra_begin()
    {
        return mesh.tetrahedra_begin();
    }
    
    typename Mesh::tetrahedron_iterator tetrahedra_end()
    {
        return mesh.tetrahedra_end();
    }
    
    
    //////////////////////
    // GETTER FUNCTIONS //
    //////////////////////
    
    typename Mesh::node_type & get(const node_key_type& k)
    {
        return mesh.find_node(k);
    }
    
    typename Mesh::edge_type & get(const edge_key_type& k)
    {
        return mesh.find_edge(k);
    }
    
    typename Mesh::face_type & get(const face_key_type& k)
    {
        return mesh.find_face(k);
    }
    
    typename Mesh::tetrahedron_type & get(const tetrahedron_key_type& k)
    {
        return mesh.find_tetrahedron(k);
    }
    
    void get_nodes(const edge_key_type& e, std::vector<node_key_type>& nodes)
    {
        nodes = std::vector<node_key_type>(2);
        mesh.vertices(e, nodes);
    }
    
    void get_nodes(const face_key_type& f, std::vector<node_key_type>& nodes)
    {
        nodes = std::vector<node_key_type>(3);
        mesh.vertices(f, nodes);
    }
    
    void get_nodes(const tetrahedron_key_type& t, std::vector<node_key_type>& nodes)
    {
        nodes = std::vector<node_key_type>(4);
        mesh.vertices(t, nodes);
    }
    
    edge_key_type get_edge(const node_key_type& n1, const node_key_type& n2)
    {
        simplex_set st1, st2;
        star(n1, st1);
        star(n2, st2);
        st1.intersection(st2);
        
        if (st1.size_edges() != 1)
        {
            return NULL_EDGE;
        }
        return *(st1.edges_begin());
    }
    
    edge_key_type get_edge(const face_key_type& f1, const face_key_type& f2)
    {
        simplex_set cl1, cl2;
        closure(f1, cl1);
        closure(f2, cl2);
        cl1.intersection(cl2);
        
        if (cl1.size_edges() != 1)
        {
            return NULL_EDGE;
        }
        return *(cl1.edges_begin());
    }
    
    void get_apices(const face_key_type& f, std::vector<node_key_type>& apices)
    {
        apices = std::vector<node_key_type>(0);
        simplex_set lk_f;
        link(f, lk_f);
        for(auto nit = lk_f.nodes_begin(); nit != lk_f.nodes_end(); nit++)
        {
            apices.push_back(*nit);
        }
    }
    
    ////////////////////
    // MESH FUNCTIONS //
    ////////////////////
    
    template<typename Key>
    bool exists(const Key& k)
    {
        return mesh.exists(k);
    }
    
    void star(const node_key_type &n, simplex_set& set)
    {
        mesh.star(n, set);
    }
    
    void star(const edge_key_type &e, simplex_set& set)
    {
        mesh.star(e, set);
    }
    
    void star(const face_key_type &f, simplex_set& set)
    {
        mesh.star(f, set);
    }
    
    void star(const tetrahedron_key_type &t, simplex_set& set)
    {
        mesh.star(t, set);
    }
    
    void star(simplex_set &set_, simplex_set& set)
    {
        mesh.star(set_, set);
    }
    
    void closure(const node_key_type &n, simplex_set& set)
    {
        mesh.closure(n, set);
    }
    
    void closure(const edge_key_type &e, simplex_set& set)
    {
        mesh.closure(e, set);
    }
    
    void closure(const face_key_type &f, simplex_set& set)
    {
        mesh.closure(f, set);
    }
    
    void closure(const tetrahedron_key_type &t, simplex_set& set)
    {
        mesh.closure(t, set);
    }
    
    void closure(simplex_set &set_, simplex_set& set)
    {
        mesh.closure(set_, set);
    }
    
    template<typename Key>
    void link(const Key& k, simplex_set& set)
    {
        mesh.link(k, set);
    }
    
    /**
     * Ensures consistent orientation of all faces to the two tetrahedra which are in the star of f.
     */
    void orient_face(const face_key_type& f)
    {
        if (get(f).is_interface())
        {
            simplex_set st_f;
            star(f, st_f);
            int label = -100;
            for (auto tit = st_f.tetrahedra_begin(); tit != st_f.tetrahedra_end(); tit++)
            {
                int tl = get(*tit).label;
                
                if (tl > label)
                {
                    mesh.orient_faces_consistently(*tit);
                }
                label = tl;
            }
        }
    }
    
    node_key_type split(tetrahedron_key_type & t)
    {
        return mesh.split_tetrahedron(t);
    }
    
    node_key_type split(const edge_key_type & e)
    {
        std::map<tetrahedron_key_type, int> tt;
        simplex_set st_e;
        star(e, st_e);
        
        for (auto tit = st_e.tetrahedra_begin(); tit != st_e.tetrahedra_end(); tit++)
        {
            tt[*tit] = get(*tit).label;
        }
        
        std::map<tetrahedron_key_type, tetrahedron_key_type> new_tets;
        node_key_type n = mesh.split_edge_helper(e, new_tets);
        
        for (auto it = new_tets.begin(); it != new_tets.end(); it++)
        {
            get(it->first).label = tt[it->second];
        }
        return n;
    }
    
    node_key_type split(const face_key_type& f)
    {
        std::map<tetrahedron_key_type, int> tt;
        simplex_set st_f;
        star(f, st_f);
        
        for (auto tit = st_f.tetrahedra_begin(); tit != st_f.tetrahedra_end(); tit++)
        {
            tt[*tit] = get(*tit).label;
        }
        
        std::map<tetrahedron_key_type, tetrahedron_key_type> new_tets;
        node_key_type n = mesh.split_face_helper(f, new_tets);
        
        for(auto it = new_tets.begin(); it != new_tets.end(); it++)
        {
            get(it->first).label = tt[it->second];
        }
        return n;
    }
    
    node_key_type collapse(edge_key_type & e, node_key_type const & n1, node_key_type const & n2)
    {
        return mesh.edge_collapse_helper(e, n1, n2);
    }
    
    ///////////////////////////////////
    // SHOULD BE REMOVED OR REPLACED //
    ///////////////////////////////////
    
    void remove_edge(const edge_key_type& e, std::vector<node_key_type> & new_edges, simplex_set& new_simplices)
    {
        mesh.remove_edge(e, new_edges, new_simplices);
    }
    
    void find_min_multi_face(const face_key_type& f, simplex_set & multi_face, simplex_set & min_multi_face)
    {
        mesh.find_min_multi_face(f, multi_face, min_multi_face);
    }
    
    void multi_face_retriangulation(simplex_set & removed_faces, std::vector<node_key_type> & new_edges_desc,
                                    simplex_set & new_faces, simplex_set & new_simplices)
    {
        mesh.multi_face_retriangulation(removed_faces, new_edges_desc, new_faces, new_simplices);
    }
    
    void multi_face_remove(simplex_set & removed_faces, simplex_set & new_simplices)
    {
        mesh.multi_face_remove(removed_faces, new_simplices);
    }
    
    void boundary_2manifold(simplex_set & faces, simplex_set & result_set)
    {
        mesh.boundary_2manifold(faces, result_set);
    }
    
    bool is_connected(const face_key_type& f1, const face_key_type& f2)
    {
        edge_key_type e;
        return mesh.get_intersection(f1, f2, e);
    }
    
    ///////////////////////
    // UTILITY FUNCTIONS //
    ///////////////////////
        
    void garbage_collect()
    {
        mesh.garbage_collect();
    }
    
private:
    void check_validity()
    {
        std::cout << "Validity check" << std::endl;
        
        bool valid = true;
        
        for (auto fit = mesh.faces_begin(); fit != mesh.faces_end(); fit++)
        {
            simplex_set st_f;
            mesh.star(fit.key(), st_f);
            if (st_f.size_tetrahedra() > 2)
            {
                std::map<tetrahedron_key_type, simplex_set> t_bnds;
                
                std::cout << fit.key() << ": " << std::endl;
                auto tit = st_f.tetrahedra_begin();
                while (tit != st_f.tetrahedra_end())
                {
                    simplex_set bnd_t;
                    mesh.boundary(*tit, bnd_t);
                    t_bnds[*tit] = bnd_t;
                    
                    auto ftit = bnd_t.nodes_begin();
                    while (ftit != bnd_t.nodes_end())
                    {
                        std::cout << *ftit << " ";
                        ++ftit;
                    }
                    std::cout << std::endl;
                    ++tit;
                }
                
                while (!t_bnds.empty())
                {
                    typename std::map<tetrahedron_key_type, simplex_set>::iterator it;
                    it = t_bnds.begin();
                    tetrahedron_key_type t = it->first;
                    simplex_set bnd0;
                    bnd0.add(it->second);
                    ++it;
                    while (it != t_bnds.end())
                    {
                        simplex_set bnd1;
                        bnd1.add(it->second);
                        bnd1.difference(bnd0);
                        if (bnd1.size_faces() == 0)
                        {
                            mesh.unsafe_remove(t); // FUCKED
                            break;
                        }
                        ++it;
                    }
                    t_bnds.erase(t_bnds.begin());
                }
                
                valid = false;
            }
        }
        
        if (!valid)
            std::cout << "Input mesh invalid" << std::endl;
        else
            std::cout << "Input mesh valid" << std::endl;
    }
};

#endif
