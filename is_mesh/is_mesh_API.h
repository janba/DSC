//
//  is_mesh_API.h
//  DSC
//
//  Created by Asger Nyman Christiansen on 6/14/13.
//  Copyright (c) 2013 Asger Nyman Christiansen. All rights reserved.
//

#pragma once

#include <is_mesh/is_mesh.h>
#include <is_mesh/is_mesh_lists_read.h>

namespace is_mesh {
    
    template <typename node_traits, typename edge_traits, typename face_traits, typename tet_traits>
    class ISMesh
    {
        typedef typename is_mesh::t4mesh< node_traits, tet_traits, edge_traits, face_traits> Mesh;
        
    public:
        typedef NodeKey            node_key;
        typedef EdgeKey            edge_key;
        typedef FaceKey            face_key;
        typedef TetrahedronKey     tet_key;
        typedef typename Mesh::simplex_set_type         simplex_set;
        
    private:
        Mesh mesh;
        
    public:
        template<typename real>
        ISMesh(std::vector<real> & points, std::vector<int> & tets)
        {
            vectors_read(points, tets, mesh);
            init();
//            validity_check();
        }
        
        ///////////////
        // ITERATORS //
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
        
        /////////////////////
        // LABEL FUNCTIONS //
        /////////////////////
    public:
        template<typename key>
        bool is_interface(const key& k)
        {
            return get(k).is_interface();
        }
        
        template<typename key>
        bool is_boundary(const key& k)
        {
            return get(k).is_boundary();
        }
        
        template<typename key>
        bool is_crossing(const key& k)
        {
            return get(k).is_crossing();
        }
        
        int get_label(const tet_key& t)
        {
            return get(t).label();
        }
        
    private:
        template<typename key>
        void set_interface(const key& k, bool b)
        {
            return get(k).set_interface(b);
        }
        
        template<typename key>
        void set_boundary(const key& k, bool b)
        {
            return get(k).set_boundary(b);
        }
        
        template<typename key>
        void set_crossing(const key& k, bool b)
        {
            return get(k).set_crossing(b);
        }
        
    public:
        void set_label(const tet_key& t, int label)
        {
            get(t).label(label);
            simplex_set cl_t;
            closure(t, cl_t);
            update(cl_t);
        }
        
    private:
        /**
         * Perform an initial update of flags for all nodes, edges and faces.
         */
        void init()
        {
            for (auto fit = faces_begin(); fit != faces_end(); fit++)
            {
                update_flag(fit.key());
            }
            
            for (auto eit = edges_begin(); eit != edges_end(); eit++)
            {
                update_flag(eit.key());
            }
            
            for (auto nit = nodes_begin(); nit != nodes_end(); nit++)
            {
                update_flag(nit.key());
            }
        }
        
        /**
         * Updates the flags (is interface, is boundary, is crossing) of simplices in set.
         */
        void update(simplex_set & set)
        {
            // Update faces
            for (auto fit = set.faces_begin(); fit != set.faces_end(); fit++)
            {
                if (exists(*fit))
                {
                    update_flag(*fit);
                }
            }
            
            // Update edges
            for (auto eit = set.edges_begin(); eit != set.edges_end(); eit++)
            {
                if (exists(*eit))
                {
                    update_flag(*eit);
                }
            }
            
            // Update nodes
            for (auto nit = set.nodes_begin(); nit != set.nodes_end(); nit++)
            {
                if (exists(*nit))
                {
                    update_flag(*nit);
                }
            }
        }
        
        void update_flag(const face_key & f)
        {
            set_interface(f, false);
            set_boundary(f, false);
            
            simplex_set st_f;
            star(f, st_f);
            if (st_f.size_tetrahedra() == 1)
            {
                // On the boundary
                set_boundary(f, true);
                if (get_label(*(st_f.tetrahedra_begin())) != 0)
                {
                    set_interface(f, true);
                }
            }
            else if(st_f.size_tetrahedra() == 2)
            {
                auto tit = st_f.tetrahedra_begin();
                int label0 = get_label(*tit);   ++tit;
                int label1 = get_label(*tit);
                if (label0 != label1)
                {
                    // On the interface
                    set_interface(f, true);
                }
            }
        }
        
        void update_flag(const edge_key & e)
        {
            set_boundary(e, false);
            set_interface(e, false);
            set_crossing(e, false);
            
            simplex_set st_e;
            star(e, st_e);
            
            int i = 0;
            for (auto fit = st_e.faces_begin(); fit != st_e.faces_end(); fit++)
            {
                if (exists(*fit))
                {
                    if (is_boundary(*fit))
                    {
                        set_boundary(e, true);
                    }
                    if (is_interface(*fit))
                    {
                        set_interface(e, true);
                        i++;
                    }
                }
            }
            if(i > 2)
            {
                set_crossing(e, true);
            }
        }
        
        void connected_component(simplex_set& st_n, const tet_key& t)
        {
            int label = get_label(t);
            st_n.erase(t);
            simplex_set cl_t;
            closure(t, cl_t);
            
            for(auto fit = cl_t.faces_begin(); fit != cl_t.faces_end(); fit++)
            {
                tet_key t2 = get_tet(t, *fit);
                if(st_n.contains(t2) && label == get_label(t2))
                {
                    connected_component(st_n, t2);
                }
            }
        }
        
        bool crossing(const node_key& n)
        {
            simplex_set st_n;
            star(n, st_n);
            
            int c = 0;
            while (st_n.size_tetrahedra() > 0)
            {
                if(c == 2)
                {
                    return true;
                }
                tet_key t = *st_n.tetrahedra_begin();
                connected_component(st_n, t);
                c++;
            }
            return false;
        }
        
        void update_flag(const node_key & n)
        {
            set_interface(n, false);
            set_boundary(n, false);
            set_crossing(n, false);
            
            simplex_set st_n;
            star(n, st_n);
            for (auto eit = st_n.edges_begin(); eit != st_n.edges_end(); eit++)
            {
                if (exists(*eit))
                {
                    if (is_interface(*eit))
                    {
                        set_interface(n, true);
                    }
                    if (is_boundary(*eit))
                    {
                        set_boundary(n, true);
                    }
                    if (is_crossing(*eit))
                    {
                        set_crossing(n, true);
                    }
                }
            }
            if(!is_crossing(n) && is_interface(n) && crossing(n))
            {
                set_crossing(n, true);
            }
        }
        
        //////////////////////
        // GETTER FUNCTIONS //
        //////////////////////
    public:
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
        
        SimplexSet<NodeKey> get_boundary(const EdgeKey& eid)
        {
            std::vector<NodeKey> temp = *mesh.lookup_simplex(eid).get_boundary();
            SimplexSet<NodeKey> res;
            for(auto k : temp)
            {
                res.push_back(k);
            }
            return res;
        }
        
        SimplexSet<EdgeKey> get_boundary(const FaceKey& fid)
        {
            std::vector<EdgeKey> temp = *mesh.lookup_simplex(fid).get_boundary();
            SimplexSet<EdgeKey> res;
            for(auto k : temp)
            {
                res.push_back(k);
            }
            return res;
        }
        
        SimplexSet<FaceKey> get_boundary(const TetrahedronKey& tid)
        {
            std::vector<FaceKey> temp =  *mesh.lookup_simplex(tid).get_boundary();
            SimplexSet<FaceKey> res;
            for(auto k : temp)
            {
                res.push_back(k);
            }
            return res;
        }
        
        
        SimplexSet<EdgeKey> get_co_boundary(const NodeKey& nid)
        {
            return *mesh.lookup_simplex(nid).get_co_boundary();
        }
        
        SimplexSet<FaceKey> get_co_boundary(const EdgeKey& eid)
        {
            return *mesh.lookup_simplex(eid).get_co_boundary();
        }
        
        SimplexSet<TetrahedronKey> get_co_boundary(const FaceKey& fid)
        {
            return *mesh.lookup_simplex(fid).get_co_boundary();
        }
        
        SimplexSet<NodeKey> get_boundary(const SimplexSet<EdgeKey>& set)
        {
            SimplexSet<NodeKey> res;
            for (auto &k : set)
            {
                res += get_boundary(k);
            }
            return res;
        }
        
        SimplexSet<EdgeKey> get_boundary(const SimplexSet<FaceKey>& set)
        {
            SimplexSet<EdgeKey> res;
            for (auto &k : set)
            {
                res += get_boundary(k);
            }
            return res;
        }
        
        SimplexSet<FaceKey> get_boundary(const SimplexSet<TetrahedronKey>& set)
        {
            SimplexSet<FaceKey> res;
            for (auto &k : set)
            {
                res += get_boundary(k);
            }
            return res;
        }
        
        SimplexSet<EdgeKey> get_co_boundary(const SimplexSet<NodeKey>& set)
        {
            SimplexSet<EdgeKey> res;
            for (auto &k : set)
            {
                res += get_co_boundary(k);
            }
            return res;
        }
        
        SimplexSet<FaceKey> get_co_boundary(const SimplexSet<EdgeKey>& set)
        {
            SimplexSet<FaceKey> res;
            for (auto &k : set)
            {
                res += get_co_boundary(k);
            }
            return res;
        }
        
        SimplexSet<TetrahedronKey> get_co_boundary(const SimplexSet<FaceKey>& set)
        {
            SimplexSet<TetrahedronKey> res;
            for (auto &k : set)
            {
                res += get_co_boundary(k);
            }
            return res;
        }
        
        std::vector<node_key> get_nodes(const edge_key& eid)
        {
            std::vector<node_key> nodes;
            for (auto nid : *mesh.lookup_simplex(eid).get_boundary()) {
                nodes.push_back(nid);
            }
            return nodes;
        }
        
        std::vector<node_key> get_nodes(const face_key& fid, bool sort = true)
        {
            assert(exists(fid));
            std::vector<node_key> nodes;
            for (auto eid : *mesh.lookup_simplex(fid).get_boundary()) {
                assert(exists(eid));
                if(sort)
                {
                    mesh.orient_face_helper(fid, eid, true);
                }
                nodes.push_back(get_nodes(eid)[0]);
            }
            return nodes;
        }
        
        std::vector<node_key> get_nodes(const tet_key& tid, bool sort = true)
        {
            assert(exists(tid));
            std::vector<node_key> nodes;
            auto bit = mesh.lookup_simplex(tid).get_boundary();
            face_key fid = *bit->begin();
            assert(exists(fid));
            if(sort)
            {
                mesh.orient_face_helper(tid, fid, true);
            }
            for (auto nid : get_nodes(fid)) {
                nodes.push_back(nid);
            }
            
            for (auto nid : get_nodes(*(bit->begin()+1)))
            {
                if(std::find(nodes.begin(), nodes.end(), nid) == nodes.end())
                {
                    nodes.push_back(nid);
                    break;
                }
            }
            return nodes;
        }
        
        std::vector<edge_key> get_edges(const node_key& nid)
        {
            auto coboundary = *mesh.lookup_simplex(nid).get_co_boundary();
            return std::vector<edge_key>(coboundary.begin(), coboundary.end());
        }
        
        std::vector<edge_key> get_edges(const face_key& fid)
        {
            std::vector<edge_key> edges;
            for (auto eid : *mesh.lookup_simplex(fid).get_boundary()) {
                edges.push_back(eid);
            }
            return edges;
        }
        
        std::vector<edge_key> get_edges(const tet_key& tid, bool sort = true)
        {
            std::vector<edge_key> edges;
            int j = 0;
            for(auto fid : *mesh.lookup_simplex(tid).get_boundary())
            {
                if(sort)
                {
                    mesh.orient_face_helper(tid, fid, true);
                }
                auto f_edges = get_edges(fid);
                if(edges.size() == 0)
                {
                    for (auto eid : f_edges) {
                        edges.push_back(eid);
                    }
                    edges.resize(6);
                }
                else {
                    for (int i = 0; i < 3; i++)
                    {
                        if(f_edges[i] == edges[0] || f_edges[i] == edges[1] || f_edges[i] == edges[2])
                        {
                            edges[3+j] = f_edges[(i+1)%3];
                            j++;
                            break;
                        }
                    }
                }
            }
            return edges;
        }
        
        edge_key get_edge(const node_key& n1, const node_key& n2)
        {
            simplex_set st1, st2;
            star(n1, st1);
            star(n2, st2);
            st1.intersection(st2);
            
            if (st1.size_edges() != 1)
            {
                return EdgeKey();
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
                return EdgeKey();
            }
            return *(cl1.edges_begin());
        }
        
        std::vector<edge_key> get_edges(const std::vector<tet_key>& tets)
        {
            std::vector<edge_key> edges;
            for (auto t : tets)
            {
                edges = uni(get_edges(t), edges);
            }
            return edges;
        }
        
        std::vector<face_key> get_faces(const node_key& nid)
        {
            std::vector<face_key> fids;
            for(auto e : *mesh.lookup_simplex(nid).get_co_boundary())
            {
                auto e_fids = get_faces(e);
                fids.insert(fids.end(), e_fids.begin(), e_fids.end());
            }
            return fids;
        }
        
        std::vector<face_key> get_faces(const edge_key& eid)
        {
            auto coboundary = *mesh.lookup_simplex(eid).get_co_boundary();
            return std::vector<face_key>(coboundary.begin(), coboundary.end());
        }
        
        std::vector<face_key> get_faces(const tet_key& tid)
        {
            auto boundary = *mesh.lookup_simplex(tid).get_boundary();
            return std::vector<face_key>(boundary.begin(), boundary.end());
        }
        
        std::vector<face_key> get_faces(const std::vector<tet_key>& tids)
        {
            std::vector<face_key> faces;
            for (auto tid : tids)
            {
                faces = uni(get_faces(tid), faces);
            }
            return faces;
        }
        
        face_key get_face(const node_key& n1, const node_key& n2, const node_key& n3)
        {
            simplex_set st1, st2, st3;
            star(n1, st1);
            star(n2, st2);
            star(n3, st3);
            
            st1.intersection(st2);
            st1.intersection(st3);
            
            if (st1.size_faces() != 1)
            {
                return FaceKey();
            }
            return *(st1.faces_begin());
        }
        
        face_key get_face(const tet_key& t1, const tet_key& t2)
        {
            simplex_set cl1, cl2;
            closure(t1, cl1);
            closure(t2, cl2);
            cl1.intersection(cl2);
            
            if (cl1.size_faces() != 1)
            {
                return FaceKey();
            }
            return *(cl1.faces_begin());
        }
        
        std::vector<tet_key> get_tets(const node_key& nid)
        {
            std::vector<tet_key> tets;
            for (auto eid : *mesh.lookup_simplex(nid).get_co_boundary()) {
                tets = uni(tets, get_tets(eid));
            }
            return tets;
        }
        
        std::vector<tet_key> get_tets(const edge_key& eid)
        {
            std::vector<tet_key> tets;
            for (auto fid : *mesh.lookup_simplex(eid).get_co_boundary()) {
                for (auto tid : *mesh.lookup_simplex(fid).get_co_boundary()) {
                    if(std::find(tets.begin(), tets.end(), tid) == tets.end())
                    {
                        tets.push_back(tid);
                    }
                }
            }
            return tets;
        }
        
        std::vector<tet_key> get_tets(const face_key& fid)
        {
            auto coboundary = *mesh.lookup_simplex(fid).get_co_boundary();
            return std::vector<tet_key>(coboundary.begin(), coboundary.end());
        }
        
        tet_key get_tet(const tet_key& t, const face_key& f)
        {
            simplex_set st_f;
            star(f, st_f);
            for(auto tit = st_f.tetrahedra_begin(); tit != st_f.tetrahedra_end(); tit++)
            {
                if(*tit != t)
                {
                    return *tit;
                }
            }
            return TetrahedronKey();
        }
        
        node_key get_apex(const tet_key& t, const face_key& f)
        {
            simplex_set cl_f, cl_t;
            closure(t, cl_t);
            closure(f, cl_f);
            cl_t.difference(cl_f);
            return *cl_t.nodes_begin();
        }
        
        node_key get_apex(const face_key& f, const edge_key& e)
        {
            simplex_set cl_f, cl_e;
            closure(f, cl_f);
            closure(e, cl_e);
            cl_f.difference(cl_e);
#ifdef DEBUG
            assert(cl_f.size_nodes() == 1);
#endif
            return *cl_f.nodes_begin();
        }
        
        std::vector<node_key> get_apices(const face_key& f)
        {
            std::vector<node_key> apices;
            simplex_set lk_f;
            link(f, lk_f);
            for(auto nit = lk_f.nodes_begin(); nit != lk_f.nodes_end(); nit++)
            {
                apices.push_back(*nit);
            }
            return apices;
        }
        
        ////////////////////
        // MESH FUNCTIONS //
        ////////////////////
    public:
        
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
        void orient_face(const face_key& fid)
        {
            if (is_interface(fid))
            {
                simplex_set st_f;
                star(fid, st_f);
                int label = -100;
                for (auto tit = st_f.tetrahedra_begin(); tit != st_f.tetrahedra_end(); tit++)
                {
                    int tl = get_label(*tit);
                    
                    if (tl > label)
                    {
                        mesh.orient_faces_consistently(*tit);
                    }
                    label = tl;
                }
            }
            else {
                simplex_set st_f;
                star(fid, st_f);
                mesh.orient_faces_consistently(*st_f.tetrahedra_begin());
            }
        }
        
        node_key split(const edge_key & e)
        {
            std::map<tet_key, int> tt;
            simplex_set st_e;
            star(e, st_e);
            
            for (auto tit = st_e.tetrahedra_begin(); tit != st_e.tetrahedra_end(); tit++)
            {
                tt[*tit] = get_label(*tit);
            }
            
            std::map<tet_key, tet_key> new_tets;
            node_key n = mesh.split_edge_helper(e, new_tets);
            
            for (auto it = new_tets.begin(); it != new_tets.end(); it++)
            {
                set_label(it->first, tt[it->second]);
            }
            
            simplex_set st_n;
            star(n, st_n);
            st_n.insert(n);
            update(st_n);
            return n;
        }
        
        node_key split(const face_key& f)
        {
            std::map<tet_key, int> tt;
            simplex_set st_f;
            star(f, st_f);
            
            for (auto tit = st_f.tetrahedra_begin(); tit != st_f.tetrahedra_end(); tit++)
            {
                tt[*tit] = get_label(*tit);
            }
            
            std::map<tet_key, tet_key> new_tets;
            node_key n = mesh.split_face_helper(f, new_tets);
            
            for(auto it = new_tets.begin(); it != new_tets.end(); it++)
            {
                set_label(it->first, tt[it->second]);
            }
            
            simplex_set st_n;
            star(n, st_n);
            st_n.insert(n);
            update(st_n);
            return n;
        }
        
        node_key split_new(const edge_key& eid)
        {
            std::cout << "SPLIT" << std::endl;
            auto nids = get_nodes(eid);
            auto fids = get_faces(eid);
            auto tids = get_tets(eid);
            auto nid = nids[1];
            
            // Split edge
            auto new_nid = mesh.insert_node();
            mesh.lookup_simplex(new_nid).set_pos(0.5*(mesh.lookup_simplex(nids[0]).get_pos() + mesh.lookup_simplex(nids[1]).get_pos()));
            mesh.lookup_simplex(new_nid).set_destination(0.5*(mesh.lookup_simplex(nids[0]).get_destination() + mesh.lookup_simplex(nids[1]).get_destination()));
            
            mesh.lookup_simplex(nid).remove_co_face(eid);
            mesh.lookup_simplex(new_nid).add_co_face(eid);
            auto& edge = mesh.lookup_simplex(eid);
            edge.remove_face(nid);
            edge.add_face(new_nid);
            
            auto new_eid = mesh.insert_edge(new_nid, nids[1]);
            
            // Update faces, create faces
            std::vector<edge_key> new_f_eids;
            std::vector<face_key> new_fids;
            for (auto f : fids)
            {
                auto eids = get_edges(f);
                assert(eids.size() == 3);
                auto temp = uni(uni(get_nodes(eids[0]), get_nodes(eids[1])), get_nodes(eids[2]));
                assert(temp.size() == 4);
                auto e_nids = difference(temp, nids);
                assert(e_nids.size() == 2);
                new_f_eids.push_back(mesh.insert_edge(e_nids[0], e_nids[1]));
                
                EdgeKey f_eid;
                assert(get_edges(f).size() == 3);
                for (auto e : get_edges(f)) {
                    if(get_nodes(e)[0] == nid || get_nodes(e)[1] == nid)
                    {
                        f_eid = e;
                        break;
                    }
                }
                
                mesh.lookup_simplex(f_eid).remove_co_face(f);
                mesh.lookup_simplex(new_f_eids.back()).add_co_face(f);
                auto& face = mesh.lookup_simplex(f);
                face.remove_face(f_eid);
                face.add_face(new_f_eids.back());
                new_fids.push_back(mesh.insert_face(f_eid, new_f_eids.back(), new_eid));
            }
            
            // Update tetrahedra, create tetrahedra
            std::vector<tet_key> new_tids;
            for (auto t : tids)
            {
                std::vector<edge_key> t_eids;
                assert(get_faces(t).size() == 4);
                for (auto f : get_faces(t))
                {
                    t_eids = uni(t_eids, get_edges(f));
                }
                
                assert(t_eids.size() == 8);
                auto f_eids = intersection(t_eids, new_f_eids);
                assert(f_eids.size() == 2);
                for(auto e : t_eids)
                {
                    if(is_neighbour(e, f_eids))
                    {
                        f_eids = uni(f_eids, {e});
                    }
                }
                assert(f_eids.size() == 3);
                auto new_t_fid = mesh.insert_face(f_eids[0], f_eids[1], f_eids[2]);
                
                FaceKey t_fid;
                for (auto f : get_faces(t)) {
                    if(get_nodes(f)[0] == nid || get_nodes(f)[1] == nid || get_nodes(f)[2] == nid)
                    {
                        t_fid = f;
                        break;
                    }
                }
                assert(t_fid.is_valid());
                mesh.lookup_simplex(t_fid).remove_co_face(t);
                mesh.lookup_simplex(new_t_fid).add_co_face(t);
                auto& tet = mesh.lookup_simplex(t);
                tet.remove_face(t_fid);
                tet.add_face(new_t_fid);
                
                std::vector<face_key> t_fids = {new_t_fid, t_fid};
                for (auto f : new_fids)
                {
                    if(is_neighbour(f, t_fids))
                    {
                        t_fids.push_back(f);
                    }
                }
                assert(t_fids.size() == 4);
                new_tids.push_back(mesh.insert_tetrahedron(t_fids[0], t_fids[1], t_fids[2], t_fids[3]));
            }
            
            // Update flags
            for (int i = 0; i < tids.size(); i++)
            {
                set_label(new_tids[i], get_label(tids[i]));
            }
            
            for(auto t : uni(new_tids, tids))
            {
                if(is_inverted(t))
                {
                    mesh.lookup_simplex(t).invert_orientation();
                }
            }
            
            return new_nid;
        }
        
        
        node_key collapse_new(const edge_key& eid)
        {
            auto nids = get_nodes(eid);
            auto fids = get_faces(eid);
            auto tids = get_tets(eid);
            node_key n = nids[1];
            
            // Remove edge
            mesh.remove(eid);
            
            // Remove faces
            std::vector<std::vector<edge_key>> merge_edges;
            for(auto f : fids)
            {
                auto eids = get_edges(f);
                assert(eids.size() == 2);
                auto nodes = get_nodes(eids[0]);
                assert(nodes.size() == 2);
                
                if(nodes[0] != n && nodes[1] != n)
                {
                    assert(get_nodes(eids[1])[0] == n || get_nodes(eids[1])[1] == n);
                    std::swap(eids[0], eids[1]);
                }
                merge_edges.push_back(eids);
                mesh.remove(f);
            }
            
            // Remove tetrahedra
            std::vector<std::vector<face_key>> merge_faces;
            for(auto t : tids)
            {
                auto fids = get_faces(t);
                assert(fids.size() == 2);
                auto nodes = get_nodes(fids[0]);
                assert(nodes.size() == 3);
                
                if(nodes[0] != n && nodes[1] != n && nodes[2] != n)
                {
                    assert(get_nodes(fids[1])[0] == n || get_nodes(fids[1])[1] == n || get_nodes(fids[1])[2] == n);
                    std::swap(fids[0], fids[1]);
                }
                merge_faces.push_back(fids);
                mesh.remove(t);
            }
            
            // Merge nodes
            mesh.merge(n, nids[0]);
            
            for (auto &eids : merge_edges)
            {
                assert(eids.size() == 2);
                mesh.merge(eids[0], eids[1]);
            }
            
            for (auto &fids : merge_faces)
            {
                assert(fids.size() == 2);
                mesh.merge(fids[0], fids[1]);
            }
            
            // Update flags
            simplex_set st_n, cl_st_n;
            star(n, st_n);
            closure(st_n, cl_st_n);
            update(cl_st_n);
            
            for(auto tit = st_n.tetrahedra_begin(); tit != st_n.tetrahedra_end(); tit++)
            {
                if(is_inverted(*tit))
                {
                    mesh.lookup_simplex(*tit).invert_orientation();
                }
            }
            return n;
        }
        
        node_key collapse(edge_key& e)
        {
            auto nodes = get_nodes(e);
#ifdef DEBUG
            assert(nodes[0].is_valid());
            assert(nodes[1].is_valid());
#endif
            node_key n = mesh.edge_collapse_helper(e, nodes[0], nodes[1]);
            if (!n.is_valid()) {
                return n;
            }
            simplex_set st_n, cl_st_n;
            star(n, st_n);
            closure(st_n, cl_st_n);
            update(cl_st_n);
            return n;
        }
        
        node_key flip_32(const edge_key& e)
        {
#ifdef DEBUG
            assert(!is_interface(e) && !is_boundary(e));
#endif
            simplex_set lk_e;
            link(e, lk_e);
#ifdef DEBUG
            assert(lk_e.size_nodes() == 3);
#endif
            node_key n1 = *lk_e.nodes_begin();
            node_key n2 = split(e);
            edge_key e2 = get_edge(n1, n2);
#ifdef DEBUG
            assert(e2.is_valid());
#endif
            node_key n3 = collapse(e2);
#ifdef DEBUG
            assert(n3.is_valid());
            assert(n1 == n3);
#endif
            return n3;
        }
        
        template<typename key_type>
        std::vector<key_type> difference(const std::vector<key_type>& keys1, const std::vector<key_type>& keys2)
        {
            std::vector<key_type> keys;
            for (auto &k1 : keys1) {
                if(std::find(keys2.begin(), keys2.end(), k1) == keys2.end())
                {
                    keys.push_back(k1);
                }
            }
            
            for (auto &k2 : keys2) {
                if(std::find(keys1.begin(), keys1.end(), k2) == keys1.end())
                {
                    keys.push_back(k2);
                }
            }
            return keys;
        }
        
        template<typename key_type>
        std::vector<key_type> uni(const std::vector<key_type>& keys1, const std::vector<key_type>& keys2)
        {
            std::vector<key_type> keys;
            for (auto &k : keys1) {
                if(std::find(keys.begin(), keys.end(), k) == keys.end())
                {
                    keys.push_back(k);
                }
            }
            for (auto &k : keys2) {
                if(std::find(keys.begin(), keys.end(), k) == keys.end())
                {
                    keys.push_back(k);
                }
            }
            return keys;
        }
        
        template<typename key_type>
        std::vector<key_type> intersection(const std::vector<key_type>& keys1, const std::vector<key_type>& keys2)
        {
            std::vector<key_type> keys;
            for (auto &k1 : keys1) {
                if(std::find(keys2.begin(), keys2.end(), k1) != keys2.end())
                {
                    keys.push_back(k1);
                }
            }
            return keys;
        }
        
        bool is_neighbour(const NodeKey& nid1, const NodeKey& nid2)
        {
            auto coboundary = get_co_boundary(nid1);
            for (auto n : get_co_boundary(nid2))
            {
                if(coboundary.contains(n))
                {
                    return true;
                }
            }
            return false;
        }
        
        template<typename key_type>
        bool is_neighbour(const key_type& key1, const key_type& key2)
        {
            auto boundary = get_boundary(key1);
            for (auto k : get_boundary(key2))
            {
                if(boundary.contains(k))
                {
                    return true;
                }
            }
            return false;
        }
        
        template<typename key_type>
        bool is_neighbour(const key_type& key, const std::vector<key_type>& keys)
        {
            int i = 0;
            for (auto k1 : *mesh.lookup_simplex(key).get_boundary())
            {
                for (auto k2 : keys)
                {
                    auto boundary = *mesh.lookup_simplex(k2).get_boundary();
                    if(std::find(boundary.begin(), boundary.end(), k1) != boundary.end())
                    {
                        i++;
                        break;
                    }
                }
            }
            return i == keys.size();
        }
        
        std::vector<face_key> create_faces(const edge_key& interior_edge, const std::vector<edge_key>& exterior_edges)
        {
            assert(exterior_edges.size()%3 == 0);
            std::vector<std::vector<edge_key>> faces_edges(exterior_edges.size()/3);
            for(auto& face_edges : faces_edges)
            {
                face_edges.push_back(interior_edge);
            }
            for (auto e : exterior_edges)
            {
                for(auto& face_edges : faces_edges)
                {
                    if(is_neighbour(e, face_edges))
                    {
                        face_edges.push_back(e);
                        break;
                    }
                }
            }
            
            std::vector<face_key> new_faces;
            for(auto& face_edges : faces_edges)
            {
                assert(face_edges.size() == 3);
                new_faces.push_back(mesh.insert_face(face_edges[0], face_edges[1], face_edges[2]));
            }
            return new_faces;
        }
        
        bool is_inverted(const tet_key& tid)
        {
            auto nodes = get_nodes(tid);
            std::vector<typename node_traits::vec3> verts;
            for(auto &n : nodes)
            {
                assert(exists(n));
                verts.push_back(mesh.lookup_simplex(n).get_pos());
            }
            return dot(verts[0]-verts[3], cross(verts[1]-verts[3], verts[2]-verts[3])) < 0.;
        }
        
        /**
         * Inserts a tetrahedron into the mesh. Updates the co-boundary of the boundary edges with the newly created tetrahedron.
         * Leaves the closure of the tetrahedron in an uncompressed state.
         */
        TetrahedronKey insert_tetrahedron(FaceKey face1, FaceKey face2, FaceKey face3, FaceKey face4)
        {
            auto tetrahedron = mesh.m_tetrahedron_kernel->create();
            //update relations
            mesh.m_face_kernel->find(face1).add_co_face(tetrahedron.key());
            mesh.m_face_kernel->find(face2).add_co_face(tetrahedron.key());
            mesh.m_face_kernel->find(face3).add_co_face(tetrahedron.key());
            mesh.m_face_kernel->find(face4).add_co_face(tetrahedron.key());
            tetrahedron->add_face(face1);
            tetrahedron->add_face(face2);
            tetrahedron->add_face(face3);
            tetrahedron->add_face(face4);
            
            if(is_inverted(tetrahedron.key()))
            {
                tetrahedron->invert_orientation();
            }
            
            return tetrahedron.key();
        }
        
        tet_key create_tetrahedron(const std::vector<face_key>& interior_faces, std::vector<face_key>& exterior_faces)
        {
            std::vector<face_key> tet_faces {exterior_faces.back()};
            exterior_faces.pop_back();
            for (auto& f : interior_faces)
            {
                if(std::find(tet_faces.begin(), tet_faces.end(), f) == tet_faces.end() && is_neighbour(f, tet_faces))
                {
                    tet_faces.push_back(f);
                }
            }
            
            for (auto fit = exterior_faces.begin(); fit != exterior_faces.end(); fit++)
            {
                if(is_neighbour(*fit, tet_faces))
                {
                    tet_faces.push_back(*fit);
                    fit--;
                    exterior_faces.erase(fit+1);
                }
            }
            assert(tet_faces.size() == 4);
            return insert_tetrahedron(tet_faces[0], tet_faces[1], tet_faces[2], tet_faces[3]);
        }
        
        std::vector<tet_key> create_tetrahedra(const std::vector<face_key>& interior_faces, const std::vector<face_key>& exterior_faces)
        {
            int N_tets = (2*static_cast<int>(interior_faces.size()) + static_cast<int>(exterior_faces.size()))/4;
            assert((2*interior_faces.size()+exterior_faces.size())%4 == 0);
            
            std::vector<face_key> exterior_faces_(exterior_faces.begin(), exterior_faces.end());
            std::vector<tet_key> new_tets;
            for(int i = 0; i < N_tets; i++)
            {
                new_tets.push_back(create_tetrahedron(interior_faces, exterior_faces_));
            }
            assert(exterior_faces_.size() == 0);
            assert(new_tets.size() == N_tets);
            return new_tets;
        }
        
        void flip_32_new(const edge_key& eid)
        {
            auto faces = get_faces(eid);
            auto tets = get_tets(eid);
            assert(tets.size() == 3);
            int label = get_label(tets[0]);
            assert(label == get_label(tets[1]));
            assert(label == get_label(tets[1]));
            
            // Find edges for the new face
            auto edges = get_edges(tets);
            assert(edges.size() == 10);
            std::vector<edge_key> face_edges;
            for(auto& e : edges)
            {
                if(e != eid && !is_neighbour(e, eid))
                {
                    face_edges.push_back(e);
                }
            }
            assert(face_edges.size() == 3);
            
            // Remove edge
            mesh.remove(eid);
            
            // Remove faces
            for(auto f : faces)
            {
                mesh.remove(f);
            }
            
            // Create face
            auto new_face = mesh.insert_face(face_edges[0], face_edges[1], face_edges[2]);
            
            // Create tets
            auto exterior_faces = get_faces(tets);
            assert(exterior_faces.size() == 6);
            auto new_tets = create_tetrahedra({new_face}, exterior_faces);
            assert(new_tets.size() == 2);
            
            // Remove tets
            for(auto t : tets)
            {
                mesh.remove(t);
            }
            
            // Update flags
            for (auto t : new_tets) {
                set_label(t, label);
            }
        }
        
        void flip_23_new(const face_key& fid)
        {
            auto nodes = get_apices(fid);
            auto edges = get_edges(fid);
            auto tets = get_tets(fid);
            int label = get_label(tets[0]);
            assert(label == get_label(tets[1]));
            assert(label == get_label(tets[1]));
            
            // Create edge
            edge_key new_edge = mesh.insert_edge(nodes[0], nodes[1]);
            
            // Create faces
            std::vector<edge_key> exterior_edges = get_edges(tets);
            assert(exterior_edges.size() == 9);
            
            std::vector<face_key> new_faces = create_faces(new_edge, exterior_edges);
            assert(new_faces.size() == 3);
            
            // Remove face
            mesh.remove(fid);
            
            // Create tetrahedra
            std::vector<face_key> exterior_faces = get_faces(tets);
            assert(exterior_faces.size() == 6);
            auto new_tets = create_tetrahedra(new_faces, exterior_faces);
            assert(new_tets.size() == 3);
            
            // Remove tetrahedra
            for (auto t : tets) {
                mesh.remove(t);
            }
            
            // Update flags
            for (auto t : new_tets) {
                set_label(t, label);
            }
        }
        
        node_key flip_23(const face_key& f)
        {
#ifdef DEBUG
            assert(!is_interface(f) && !is_boundary(f));
#endif
            simplex_set lk_f;
            link(f, lk_f);
            node_key n1 = *lk_f.nodes_begin();
#ifdef DEBUG
            assert(lk_f.size_nodes() == 2);
#endif
            node_key n2 = split(f);
            edge_key e = get_edge(n1, n2);
#ifdef DEBUG
            assert(e.is_valid());
#endif
            node_key n3 = collapse(e);
#ifdef DEBUG
            assert(n3.is_valid());
            assert(n1 == n3);
#endif
            return n3;
        }
        
        template<typename child_key, typename parent_key>
        void connect(const child_key& ck, const parent_key& pk)
        {
            mesh.lookup_simplex(ck).add_co_face(pk);
            mesh.lookup_simplex(pk).add_face(ck);
        }
        
        template<typename child_key, typename parent_key>
        void disconnect(const child_key& ck, const parent_key& pk)
        {
            mesh.lookup_simplex(ck).remove_co_face(pk);
            mesh.lookup_simplex(pk).remove_face(ck);
        }
        
        template<typename child_key, typename parent_key>
        void swap(const child_key& ck1, const parent_key& pk1, const child_key& ck2, const parent_key& pk2)
        {
            if(!contains(*mesh.lookup_simplex(pk1).get_boundary(), ck1))
            {
                assert(contains(*mesh.lookup_simplex(pk1).get_boundary(), ck2));
                assert(contains(*mesh.lookup_simplex(pk2).get_boundary(), ck1));
                
                disconnect(ck1, pk2);
                disconnect(ck2, pk1);
                connect(ck1, pk1);
                connect(ck2, pk2);
            }
            else {
                assert(contains(*mesh.lookup_simplex(pk1).get_boundary(), ck1));
                assert(contains(*mesh.lookup_simplex(pk2).get_boundary(), ck2));
                
                disconnect(ck1, pk1);
                disconnect(ck2, pk2);
                connect(ck1, pk2);
                connect(ck2, pk1);
            }
        }
        
        template<typename key_type>
        bool contains(const std::vector<key_type>& keys, const key_type& key)
        {
            return std::find(keys.begin(), keys.end(), key) != keys.end();
        }
        
        void flip_22_new(const face_key& fid1, const face_key& fid2)
        {
            assert(false || !"2-2 flip is not working yet");
            edge_key eid = intersection(get_edges(fid1), get_edges(fid2)).front();
            node_key nid1 = difference(get_nodes(eid), get_nodes(fid1)).front();
            node_key nid2 = difference(get_nodes(eid), get_nodes(fid2)).front();
            
            auto nids = get_nodes(eid);
            auto fids = get_faces(eid);
            auto tids = get_tets(eid);
            
            // Find apices
            
            std::vector<node_key> apices;
            for(auto f : difference(fids, {fid1, fid2}))
            {
                apices = uni(difference( get_nodes(f), nids), apices);
            }
            assert(apices.size() == 2);
            
            // Reconnect edge
            disconnect(nids[0], eid);
            disconnect(nids[1], eid);
            connect(nid1, eid);
            connect(nid2, eid);
            
            // Reconnect faces
            edge_key eid1 = intersection(get_edges(nids[0]), get_edges(nid1)).front();
            edge_key eid2 = intersection(get_edges(nids[1]), get_edges(nid2)).front();
            swap(eid1, fid1, eid2, fid2);
            
            for(auto f : difference(fids, {fid1, fid2}))
            {
                auto rm_eids = difference(get_edges(f), {eid});
                assert(rm_eids.size() == 2);
                auto f_tids = get_tets(f);
                assert(f_tids.size() == 2);
                auto t_eids = uni(get_edges(f_tids[0], false), get_edges(f_tids[1], false));
                assert(t_eids.size() == 9);
                
                auto temp = uni(uni(get_edges(f), get_edges(fid1)), get_edges(fid2));
                assert(temp.size() == 7);
                
                auto add_eids = difference(t_eids, temp);
                assert(add_eids.size() == 2);
                disconnect(rm_eids[0], f);
                disconnect(rm_eids[1], f);
                connect(add_eids[0], f);
                connect(add_eids[1], f);
            }
            
            // Reconnect tetrahedra
            for (auto n : apices) {
                auto fid1 = intersection(intersection(get_faces(n), get_faces(nids[0])), get_faces(nid1));
                auto fid2 = intersection(intersection(get_faces(n), get_faces(nids[1])), get_faces(nid2));
                assert(fid1.size() == 1);
                assert(fid2.size() == 1);
                auto tids_n = intersection(tids, get_tets(n));
                assert(tids_n.size() == 2);
                swap(fid1[0], tids_n[0], fid2[0], tids_n[1]);
            }
            validity_check();
            
            // Ensure correct orientation
            for(auto t : tids)
            {
                if(is_inverted(t))
                {
                    mesh.lookup_simplex(t).invert_orientation();
                }
            }
            
            // Update flags
            for (auto t : tids) {
                simplex_set cl_t;
                closure(t, cl_t);
                update(cl_t);
            }
        }
        
        void flip_44_new(const face_key& fid1, const face_key& fid2)
        {
            edge_key eid = intersection(get_edges(fid1), get_edges(fid2)).front();
            node_key nid1 = difference(get_nodes(eid), get_nodes(fid1)).front();
            node_key nid2 = difference(get_nodes(eid), get_nodes(fid2)).front();
            auto faces = get_faces(eid);
            auto tets = get_tets(eid);
            assert(tets.size() == 4);
            std::vector<int> labels = {get_label(tets[0]), get_label(tets[1]), get_label(tets[2]), get_label(tets[3])};
            
            // Find the edges for creating the faces
            std::vector<edge_key> exterior_edges = get_edges(tets);
            exterior_edges = difference(exterior_edges, {eid});
            assert(exterior_edges.size() == 12);
            
            // Create edge
            auto new_edge = mesh.insert_edge(nid1, nid2);
            
            // Remove edge
            mesh.remove(eid);
            
            // Create faces
            std::vector<face_key> new_faces = create_faces(new_edge, exterior_edges);
            assert(new_faces.size() == 4);
            
            for (auto f1 : new_faces) {
                for (auto f2 : new_faces) {
                    if(f1 != f2)
                    {
                        assert(intersection(get_edges((f1)), get_edges(f2)).size() == 1);
                    }
                }
            }
            
            // Remove faces
            for (auto f : faces)
            {
                mesh.remove(f);
            }
            
            // Create tetrahedra
            std::vector<face_key> exterior_faces = get_faces(tets);
            assert(exterior_faces.size() == 8);
            auto new_tets = create_tetrahedra(new_faces, exterior_faces);
            assert(new_tets.size() == 4);
            
            // Remove tetrahedra
            for(auto t : tets)
            {
                mesh.remove(t);
            }
            
            // Update flags
            for(int i = 0; i < 4; i++)
            {
                set_label(new_tets[i], labels[i]);
            }
        }
        
        node_key flip_22(const face_key& f1, const face_key& f2)
        {
            return flip_44(f1, f2);
        }
        
        node_key flip_44(const face_key& f1, const face_key& f2)
        {
#ifdef DEBUG
            assert((is_interface(f1) && is_interface(f2)) || (!is_interface(f1) && !is_interface(f2)));
            assert((is_boundary(f1) && is_boundary(f2)) || (!is_boundary(f1) && !is_boundary(f2)));
#endif
            edge_key e1 = get_edge(f1, f2);
            node_key n1 = get_apex(f1, e1);
            node_key n2 = split(e1);
            edge_key e2 = get_edge(n1, n2);
#ifdef DEBUG
            assert(e2.is_valid());
#endif
            node_key n3 = collapse(e2);
#ifdef DEBUG
            assert(n3.is_valid());
            assert(n1 == n3);
#endif
            return n3;
        }
        
        ///////////////////////
        // UTILITY FUNCTIONS //
        ///////////////////////
        
        void garbage_collect()
        {
            mesh.garbage_collect();
        }
                
        void validity_check()
        {
            std::cout << "Validity check" << std::endl;
            for(auto tit = tetrahedra_begin(); tit != tetrahedra_end(); tit++)
            {
                assert(exists(tit.key()));
                // Check faces:
                auto faces = get_faces(tit.key());
                assert(faces.size() == 4);
                for (auto f : faces) {
                    assert(exists(f));
                    auto cotets = get_tets(f);
                    assert((is_boundary(f) && cotets.size() == 1) || (!is_boundary(f) && cotets.size() == 2));
                    assert(std::find(cotets.begin(), cotets.end(), tit.key()) != cotets.end());
                    for (auto f2 : faces) {
                        assert(f == f2 || is_neighbour(f, f2));
                    }
                    
                    // Check edges:
                    auto edges = get_edges(f);
                    assert(edges.size() == 3);
                    for (auto e : edges)
                    {
                        assert(exists(e));
                        auto cofaces = get_faces(e);
                        assert(std::find(cofaces.begin(), cofaces.end(), f) != cofaces.end());
                        for (auto e2 : edges) {
                            assert(e == e2 || is_neighbour(e, e2));
                        }
                        
                        // Check nodes:
                        auto nodes = get_nodes(e);
                        assert(nodes.size() == 2);
                        for (auto n : nodes)
                        {
                            assert(exists(n));
                            auto coedges = get_edges(n);
                            assert(std::find(coedges.begin(), coedges.end(), e) != coedges.end());
                        }
                    }
                    
                }
                
                assert(get_edges(tit.key()).size() == 6);
                assert(get_nodes(tit.key()).size() == 4);
            }
            
            std::cout << "Input mesh valid" << std::endl;
        }
    };
    
}
