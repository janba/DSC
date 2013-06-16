//
//  is_mesh_API.h
//  DSC
//
//  Created by Asger Nyman Christiansen on 6/14/13.
//  Copyright (c) 2013 Asger Nyman Christiansen. All rights reserved.
//

#ifndef DSC_is_mesh_API_h
#define DSC_is_mesh_API_h

#include <is_mesh/is_mesh.h>
#include <is_mesh/io/is_mesh_lists_read.h>

#include "attributes.h"

template <typename MT>
class ISMesh
{
    typedef typename OpenTissue::is_mesh::t4mesh< default_node_traits<MT>, default_tetrahedron_traits<MT>, default_edge_traits<MT>, default_face_traits<MT> > Mesh;
    
    typedef typename MT::real_type      T;
    typedef typename MT::vector3_type   V;
    
public:
    typedef typename Mesh::node_key_type            node_key;
    typedef typename Mesh::edge_key_type            edge_key;
    typedef typename Mesh::face_key_type            face_key;
    typedef typename Mesh::tetrahedron_key_type     tet_key;
    typedef typename Mesh::simplex_set_type         simplex_set;
    
    
    const node_key NULL_NODE;
    const edge_key NULL_EDGE;
    const face_key NULL_FACE;
    const tet_key NULL_TETRAHEDRON;
    
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
public:
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
protected:
    typename Mesh::node_type & get(const node_key& k)
    {
        return mesh.find_node(k);
    }
    
    typename Mesh::edge_type & get(const edge_key& k)
    {
        return mesh.find_edge(k);
    }
    
    typename Mesh::face_type & get(const face_key& k)
    {
        return mesh.find_face(k);
    }
    
    typename Mesh::tetrahedron_type & get(const tet_key& k)
    {
        return mesh.find_tetrahedron(k);
    }
    
public:
    void get_nodes(const edge_key& e, std::vector<node_key>& nodes)
    {
        nodes = std::vector<node_key>(2);
        mesh.vertices(e, nodes);
    }
    
    void get_nodes(const face_key& f, std::vector<node_key>& nodes)
    {
        nodes = std::vector<node_key>(3);
        mesh.vertices(f, nodes);
    }
    
    void get_nodes(const tet_key& t, std::vector<node_key>& nodes)
    {
        nodes = std::vector<node_key>(4);
        mesh.vertices(t, nodes);
    }
    
    edge_key get_edge(const node_key& n1, const node_key& n2)
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
    
    edge_key get_edge(const face_key& f1, const face_key& f2)
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
    
    void get_apices(const face_key& f, std::vector<node_key>& apices)
    {
        apices = std::vector<node_key>(0);
        simplex_set lk_f;
        link(f, lk_f);
        for(auto nit = lk_f.nodes_begin(); nit != lk_f.nodes_end(); nit++)
        {
            apices.push_back(*nit);
        }
    }
    
    /////////////////////////
    // ATTRIBUTE FUNCTIONS //
    /////////////////////////
public:
    /**
     * Returns the position of node n.
     */
    V get_pos(const node_key& n)
    {
        V p = get(n).get_pos();
        assert(!MT::is_nan(p[0]) && !MT::is_nan(p[1]) && !MT::is_nan(p[2]));
        return p;
    }
    
    /// Returns the positions of the nodes of edge e in verts.
    void get_pos(const edge_key & e, std::vector<V>& verts)
    {
        verts = std::vector<V>(2);
        std::vector<node_key> nodes(2);
        get_nodes(e, nodes);
        for (int k = 0; k < 2; ++k)
        {
            verts[k] = get_pos(nodes[k]);
        }
    }
    
    /// Returns the positions of the nodes of face f in verts.
    void get_pos(const face_key & f, std::vector<V>& verts)
    {
        verts = std::vector<V>(3);
        orient_face(f);
        std::vector<node_key> nodes(3);
        get_nodes(f, nodes);
        for (int k = 0; k < 3; ++k)
        {
            verts[k] = get_pos(nodes[k]);
        }
    }
    
    /// Returns the positions of the nodes of tetrahedron t in verts.
    void get_pos(const tet_key& t, std::vector<V>& verts)
    {
        verts = std::vector<V>(4);
        std::vector<node_key> nodes(4);
        get_nodes(t, nodes);
        for (int k = 0; k < 4; ++k)
        {
            verts[k] = get_pos(nodes[k]);
        }
    }
    
protected:
    /**
     * Sets the position of node n.
     */
    void set_pos(const node_key& n, V p)
    {
        get(n).set_pos(p);
    }
    
public:
    /**
     * Sets the destination where the node n is moved to when deform() is called.
     */
    void set_destination(const node_key& n, V p)
    {
        get(n).set_destination(p);
    }
    
    bool is_interface(const node_key& n)
    {
        return get(n).is_interface();
    }
    
    bool is_interface(const edge_key& e)
    {
        return get(e).is_interface();
    }
    
    bool is_interface(const face_key& f)
    {
        return get(f).is_interface();
    }
    
    ////////////////////
    // MESH FUNCTIONS //
    ////////////////////
    
    template<typename Key>
    bool exists(const Key& k)
    {
        return mesh.exists(k);
    }
    
    void star(const node_key &n, simplex_set& set)
    {
        mesh.star(n, set);
    }
    
    void star(const edge_key &e, simplex_set& set)
    {
        mesh.star(e, set);
    }
    
    void star(const face_key &f, simplex_set& set)
    {
        mesh.star(f, set);
    }
    
    void star(const tet_key &t, simplex_set& set)
    {
        mesh.star(t, set);
    }
    
    void star(simplex_set &set_, simplex_set& set)
    {
        mesh.star(set_, set);
    }
    
    void closure(const node_key &n, simplex_set& set)
    {
        mesh.closure(n, set);
    }
    
    void closure(const edge_key &e, simplex_set& set)
    {
        mesh.closure(e, set);
    }
    
    void closure(const face_key &f, simplex_set& set)
    {
        mesh.closure(f, set);
    }
    
    void closure(const tet_key &t, simplex_set& set)
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
    void orient_face(const face_key& f)
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
    
    node_key split(const edge_key & e)
    {
        std::map<tet_key, int> tt;
        simplex_set st_e;
        star(e, st_e);
        
        for (auto tit = st_e.tetrahedra_begin(); tit != st_e.tetrahedra_end(); tit++)
        {
            tt[*tit] = get(*tit).label;
        }
        
        std::map<tet_key, tet_key> new_tets;
        node_key n = mesh.split_edge_helper(e, new_tets);
        
        for (auto it = new_tets.begin(); it != new_tets.end(); it++)
        {
            get(it->first).label = tt[it->second];
        }
        return n;
    }
    
    node_key split(const face_key& f)
    {
        std::map<tet_key, int> tt;
        simplex_set st_f;
        star(f, st_f);
        
        for (auto tit = st_f.tetrahedra_begin(); tit != st_f.tetrahedra_end(); tit++)
        {
            tt[*tit] = get(*tit).label;
        }
        
        std::map<tet_key, tet_key> new_tets;
        node_key n = mesh.split_face_helper(f, new_tets);
        
        for(auto it = new_tets.begin(); it != new_tets.end(); it++)
        {
            get(it->first).label = tt[it->second];
        }
        return n;
    }
    
    node_key split(tet_key & t)
    {
        int label = get(t).label;
        
        node_key n = mesh.split_tetrahedron(t);
        
        simplex_set st_n;
        star(n, st_n);
        for (auto tit = st_n.tetrahedra_begin(); tit != st_n.tetrahedra_end(); tit++)
        {
            get(*tit).label = label;
        }
        return n;
    }
    
    node_key collapse(edge_key & e, node_key const & n1, node_key const & n2)
    {
        return mesh.edge_collapse_helper(e, n1, n2);
    }
    
    ///////////////////////////////////
    // SHOULD BE REMOVED OR REPLACED //
    ///////////////////////////////////
    
    void remove_edge(const edge_key& e, std::vector<node_key> & new_edges, simplex_set& new_simplices)
    {
        mesh.remove_edge(e, new_edges, new_simplices);
    }
    
    void find_min_multi_face(const face_key& f, simplex_set & multi_face, simplex_set & min_multi_face)
    {
        mesh.find_min_multi_face(f, multi_face, min_multi_face);
    }
    
    void multi_face_retriangulation(simplex_set & removed_faces, std::vector<node_key> & new_edges_desc,
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
                std::map<tet_key, simplex_set> t_bnds;
                
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
                    typename std::map<tet_key, simplex_set>::iterator it;
                    it = t_bnds.begin();
                    tet_key t = it->first;
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
