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
            validity_check();
            simplex_set_test();
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
        void set_label(const tet_key& tid, int label)
        {
            get(tid).label(label);
            SimplexSet<tet_key> tids = {tid};
            update(tids);
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
        void update(const SimplexSet<tet_key>& tids)
        {
            // Update faces
            auto fids = get_boundary(tids);
            for (auto f : fids)
            {
                if (exists(f))
                {
                    update_flag(f);
                }
            }
            
            // Update edges
            auto eids = get_boundary(fids);
            for (auto e : eids)
            {
                if (exists(e))
                {
                    update_flag(e);
                }
            }
            
            // Update nodes
            auto nids = get_boundary(eids);
            for (auto n : nids)
            {
                if (exists(n))
                {
                    update_flag(n);
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
            return mesh.lookup_simplex(k);
        }
        
        typename Mesh::edge_type & get(const edge_key& k)
        {
            return mesh.lookup_simplex(k);
        }
        
        typename Mesh::face_type & get(const face_key& k)
        {
            return mesh.lookup_simplex(k);
        }
        
        typename Mesh::tetrahedron_type & get(const tet_key& k)
        {
            return mesh.lookup_simplex(k);
        }
        
        const SimplexSet<NodeKey>& get_boundary(const EdgeKey& eid)
        {
            return *get(eid).get_boundary();
        }
        
        const SimplexSet<EdgeKey>& get_boundary(const FaceKey& fid)
        {
            return *get(fid).get_boundary();
        }
        
        const SimplexSet<FaceKey>& get_boundary(const TetrahedronKey& tid)
        {
            return *get(tid).get_boundary();
        }
        
        const SimplexSet<EdgeKey>& get_co_boundary(const NodeKey& nid)
        {
            return *get(nid).get_co_boundary();
        }
        
        const SimplexSet<FaceKey>& get_co_boundary(const EdgeKey& eid)
        {
            return *get(eid).get_co_boundary();
        }
        
        const SimplexSet<TetrahedronKey>& get_co_boundary(const FaceKey& fid)
        {
            return *get(fid).get_co_boundary();
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
        
        const SimplexSet<NodeKey>& get_nodes(const edge_key& eid)
        {
            return *get(eid).get_boundary();
        }
        
        SimplexSet<NodeKey> get_nodes(const face_key& fid, bool sort = true)
        {
            SimplexSet<NodeKey> nids;
            for (auto e : get_edges(fid)) {
                if(sort)
                {
                    mesh.orient_face_helper(fid, e, true);
                }
                nids += get_nodes(e)[0];
            }
            return nids;
        }
        
        SimplexSet<NodeKey> get_nodes(const tet_key& tid, bool sort = true)
        {
            SimplexSet<NodeKey> nids;
            auto fids = get_faces(tid);
            if(sort)
            {
                mesh.orient_face_helper(tid, fids[0], true);
            }
            nids += get_nodes(fids[0]);
            nids += get_nodes(fids[1]);
            assert(nids.size() == 4);
            return nids;
        }
        
        SimplexSet<NodeKey> get_nodes(const SimplexSet<tet_key>& tids)
        {
            SimplexSet<NodeKey> nids;
            for(auto t : tids)
            {
                nids += get_nodes(t);
            }
            return nids;
        }
        
        const SimplexSet<EdgeKey>& get_edges(const node_key& nid)
        {
            return *get(nid).get_co_boundary();
        }
        
        const SimplexSet<EdgeKey>& get_edges(const face_key& fid)
        {
            return *get(fid).get_boundary();
        }
        
        SimplexSet<EdgeKey> get_edges(const tet_key& tid, bool sort = true)
        {
            SimplexSet<EdgeKey> eids;
            for(auto f : get_faces(tid))
            {
                if(sort)
                {
                    mesh.orient_face_helper(tid, f, true);
                }
                SimplexSet<EdgeKey> f_eids = get_edges(f);
                if(eids.size() == 0)
                {
                    eids += f_eids;
                }
                else {
                    for(int i = 0; i < 3; i++)
                    {
                        if(eids.contains(f_eids[i]))
                        {
                            eids += f_eids[(i+1)%3];
                        }
                    }
                }
            }
            return eids;
        }
        
        EdgeKey get_edge(const NodeKey& nid1, const NodeKey& nid2)
        {
            SimplexSet<EdgeKey> eid = get_edges(nid1) & get_edges(nid2);
            assert(eid.size() == 1);
            return eid.front();
        }
        
        EdgeKey get_edge(const FaceKey& fid1, const FaceKey& fid2)
        {
            SimplexSet<EdgeKey> eid = get_edges(fid1) & get_edges(fid2);
            assert(eid.size() == 1);
            return eid.front();
        }
        
        SimplexSet<FaceKey> get_faces(const NodeKey& nid)
        {
            SimplexSet<FaceKey> fids;
            for(auto e : get_edges(nid))
            {
                fids += get_faces(e);
            }
            return fids;
        }
        
        const SimplexSet<FaceKey>& get_faces(const EdgeKey& eid)
        {
            return *get(eid).get_co_boundary();
        }
        
        const SimplexSet<FaceKey>& get_faces(const TetrahedronKey& tid)
        {
            return *get(tid).get_boundary();
        }
        
        FaceKey get_face(const NodeKey& nid1, const NodeKey& nid2, const NodeKey& nid3)
        {
            SimplexSet<FaceKey> fid = (get_faces(nid1) & get_faces(nid2)) & get_faces(nid3);
            assert(fid.size() == 1);
            return fid.front();
        }
        
        FaceKey get_face(const TetrahedronKey& tid1, const TetrahedronKey& tid2)
        {
            SimplexSet<FaceKey> fid = get_faces(tid1) & get_faces(tid2);
            assert(fid.size() == 1);
            return fid.front();
        }
        
        SimplexSet<TetrahedronKey> get_tets(const NodeKey& nid)
        {
            SimplexSet<TetrahedronKey> tids;
            for (auto e : get_edges(nid)) {
                tids += get_tets(e);
            }
            return tids;
        }
        
        SimplexSet<TetrahedronKey> get_tets(const EdgeKey& eid)
        {
            SimplexSet<TetrahedronKey> tids;
            for (auto f : get_faces(eid)) {
                tids += get_tets(f);
            }
            return tids;
        }
        
        const SimplexSet<TetrahedronKey>& get_tets(const FaceKey& fid)
        {
            return *get(fid).get_co_boundary();
        }
        
        TetrahedronKey get_tet(const TetrahedronKey& tid, const FaceKey& fid)
        {
            SimplexSet<TetrahedronKey> t = get_tets(fid) - tid;
            assert(t.size() == 1);
            return t.front();
        }
        
        NodeKey get_apex(const FaceKey& fid, const EdgeKey& e)
        {
            SimplexSet<NodeKey> n = get_nodes(fid) - get_nodes(e);
            assert(n.size() == 1);
            return n.front();
        }
        
        SimplexSet<NodeKey> get_apices(const FaceKey& fid)
        {
            return get_nodes(get_tets(fid)) - get_nodes(fid);
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
        
        node_key split(const edge_key& eid)
        {
            auto nids = get_nodes(eid);
            auto fids = get_faces(eid);
            auto tids = get_tets(eid);
            auto nid = nids[1];
            
            // Split edge
            auto new_nid = mesh.insert_node();
            get(new_nid).set_pos(0.5*(get(nids[0]).get_pos() + get(nids[1]).get_pos()));
            get(new_nid).set_destination(0.5*(get(nids[0]).get_destination() + get(nids[1]).get_destination()));
            
            get(nid).remove_co_face(eid);
            get(new_nid).add_co_face(eid);
            auto& edge = get(eid);
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
                
                get(f_eid).remove_co_face(f);
                get(new_f_eids.back()).add_co_face(f);
                auto& face = get(f);
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
                get(t_fid).remove_co_face(t);
                get(new_t_fid).add_co_face(t);
                auto& tet = get(t);
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
                    get(t).invert_orientation();
                }
            }
            
            return new_nid;
        }
        
        
        node_key collapse(const edge_key& eid)
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
            
            // Update flags and ensure no inverted tetrahedra.
            SimplexSet<tet_key> changed_tids = get_co_boundary(get_co_boundary(get_co_boundary(n)));
            for(auto t : changed_tids)
            {
                if(is_inverted(t))
                {
                    get(t).invert_orientation();
                }
            }
            update(changed_tids);
            return n;
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
            for (auto k1 : *get(key).get_boundary())
            {
                for (auto k2 : keys)
                {
                    auto boundary = *get(k2).get_boundary();
                    if(std::find(boundary.begin(), boundary.end(), k1) != boundary.end())
                    {
                        i++;
                        break;
                    }
                }
            }
            return i == keys.size();
        }
        
        bool is_inverted(const tet_key& tid)
        {
            auto nodes = get_nodes(tid);
            std::vector<typename node_traits::vec3> verts;
            for(auto &n : nodes)
            {
                assert(exists(n));
                verts.push_back(get(n).get_pos());
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
        
        face_key flip_32(const edge_key& eid)
        {
            auto e_nids = get_boundary(eid);
            auto e_fids = get_co_boundary(eid);
            assert(e_fids.size() == 3);
            auto e_tids = get_co_boundary(e_fids);
            assert(e_tids.size() == 3);
            int label = get_label(e_tids[0]);
            assert(label == get_label(e_tids[1]));
            assert(label == get_label(e_tids[2]));
            
            // Remove edge
            mesh.remove(eid);
            
            // Create face
            auto f_eids = get_boundary(get_boundary(e_tids)) - get_boundary(e_fids);
            assert(f_eids.size() == 3);
            auto new_fid = mesh.insert_face(f_eids[0], f_eids[1], f_eids[2]);
            
            // Remove faces
            for(face_key& f : e_fids)
            {
                mesh.remove(f);
            }
            
            // Create tetrahedra
            auto exterior_fids = get_boundary(e_tids);
            for (node_key& n : e_nids)
            {
                auto t_fids = exterior_fids & get_co_boundary(get_co_boundary(n));
                assert(t_fids.size() == 3);
                insert_tetrahedron(t_fids[0], t_fids[1], t_fids[2], new_fid);
            }
            
            // Remove tetrahedra
            for(tet_key& t : e_tids)
            {
                mesh.remove(t);
            }
            
            // Update flags
            assert(get_co_boundary(new_fid).size() == 2);
            for (auto t : get_co_boundary(new_fid)) {
                set_label(t, label);
            }
            return new_fid;
        }
        
        edge_key flip_23(const face_key& fid)
        {
            auto f_tids = get_co_boundary(fid);
            assert(f_tids.size() == 2);
            auto f_eids = get_boundary(fid);
            assert(f_eids.size() == 3);
            auto f_nids = get_boundary(f_eids);
            int label = get_label(f_tids[0]);
            assert(label == get_label(f_tids[1]));
            
            // Create edge
            auto e_nids = get_boundary(get_boundary(get_boundary(f_tids))) - f_nids;
            assert(e_nids.size() == 2);
            edge_key new_eid = mesh.insert_edge(e_nids[0], e_nids[1]);
            
            // Create faces
            auto new_fs_eids = get_boundary(get_boundary(f_tids)) - f_eids;
            assert(new_fs_eids.size() == 6);
            for (node_key& n : f_nids)
            {
                auto new_f_eids = new_fs_eids & get_co_boundary(n);
                assert(new_f_eids.size() == 2);
                mesh.insert_face(new_f_eids[0], new_f_eids[1], new_eid);
            }
            
            // Remove face
            mesh.remove(fid);
            
            // Create tetrahedra
            auto new_ts_fids1 = get_boundary(f_tids);
            auto new_ts_fids2 = get_co_boundary(new_eid);
            assert(new_ts_fids1.size() == 6);
            assert(new_ts_fids2.size() == 3);
            for (edge_key& e : f_eids)
            {
                auto new_t_fids = (new_ts_fids1 & get_co_boundary(e)) + (new_ts_fids2 & get_co_boundary(get_co_boundary(get_boundary(e))));
                assert(new_t_fids.size() == 4);
                insert_tetrahedron(new_t_fids[0], new_t_fids[1], new_t_fids[2], new_t_fids[3]);
            }
            
            // Remove tetrahedra
            for (auto t : f_tids) {
                mesh.remove(t);
            }
            
            // Update flags
            assert(get_co_boundary(get_co_boundary(new_eid)).size() == 3);
            for (auto t : get_co_boundary(get_co_boundary(new_eid))) {
                set_label(t, label);
            }
            return new_eid;
        }
        
        template<typename child_key, typename parent_key>
        void connect(const child_key& ck, const parent_key& pk)
        {
            get(ck).add_co_face(pk);
            get(pk).add_face(ck);
        }
        
        template<typename child_key, typename parent_key>
        void disconnect(const child_key& ck, const parent_key& pk)
        {
            get(ck).remove_co_face(pk);
            get(pk).remove_face(ck);
        }
        
        template<typename child_key, typename parent_key>
        void swap(const child_key& ck1, const parent_key& pk1, const child_key& ck2, const parent_key& pk2)
        {
            if(!get_boundary(pk1).contains(ck1))
            {
                assert(get_boundary(pk1).contains(ck2));
                assert(get_boundary(pk2).contains(ck1));
                
                disconnect(ck1, pk2);
                disconnect(ck2, pk1);
                connect(ck1, pk1);
                connect(ck2, pk2);
            }
            else {
                assert(get_boundary(pk1).contains(ck1));
                assert(get_boundary(pk2).contains(ck2));
                
                disconnect(ck1, pk1);
                disconnect(ck2, pk2);
                connect(ck1, pk2);
                connect(ck2, pk1);
            }
        }
        
        void flip(const edge_key& eid, const face_key& fid1, const face_key& fid2)
        {
            SimplexSet<face_key> fids = {fid1, fid2};
            auto e_nids = get_boundary(eid);
            auto e_fids = get_co_boundary(eid);
            auto e_tids = get_co_boundary(e_fids);
            
            // Reconnect edge
            auto new_e_nids = get_boundary(get_boundary(fids)) - e_nids;
            assert(new_e_nids.size() == 2);
            
            disconnect(e_nids[0], eid);
            disconnect(e_nids[1], eid);
            connect(new_e_nids[0], eid);
            connect(new_e_nids[1], eid);
            
            // Reconnect faces
            auto swap_eids = (get_co_boundary(e_nids[0]) & get_co_boundary(new_e_nids[0])) + (get_co_boundary(e_nids[1]) & get_co_boundary(new_e_nids[1]));
            assert(swap_eids.size() == 2);
            swap(swap_eids[0], fids[0], swap_eids[1], fids[1]);
            
            assert((e_fids - fids).size() <= 2);
            for(auto f : (e_fids - fids))
            {
                auto rm_eids = get_boundary(f) - eid;
                assert(rm_eids.size() == 2);
                auto apex = get_boundary(get_boundary(f)) - (new_e_nids + e_nids);
                assert(apex.size() == 1);
                auto add_eids = get_co_boundary(apex) & get_co_boundary(new_e_nids);
                assert(add_eids.size() == 2);
                
                disconnect(rm_eids[0], f);
                disconnect(rm_eids[1], f);
                connect(add_eids[0], f);
                connect(add_eids[1], f);
                
                // Reconnect tetrahedra
                auto tids = get_co_boundary(f);
                assert(tids.size() == 2);
                auto swap_fids = (get_boundary(tids) - e_fids) & get_co_boundary(swap_eids);
                assert(swap_eids.size() == 2);
                swap(swap_fids[0], tids[0], swap_fids[1], tids[1]);
            }
            
            // Ensure correct orientation
            for(auto t : e_tids)
            {
                if(is_inverted(t))
                {
                    get(t).invert_orientation();
                }
            }
            
            // Update flags
            update(e_tids);
        }
        
        
        void flip_22(const face_key& fid1, const face_key& fid2)
        {
            SimplexSet<edge_key> eid = get_boundary(fid1) & get_boundary(fid2);
            assert(eid.size() == 1);
            assert(get_co_boundary(eid).size() == 3);
            assert(get_co_boundary(get_co_boundary(eid)).size() == 2);
            
            flip(eid[0], fid1, fid2);
        }
        
        void flip_44(const face_key& fid1, const face_key& fid2)
        {
            SimplexSet<edge_key> eid = get_boundary(fid1) & get_boundary(fid2);
            assert(eid.size() == 1);
            assert(get_co_boundary(eid).size() == 4);
            assert(get_co_boundary(get_co_boundary(eid)).size() == 4);
            
            flip(eid[0], fid1, fid2);
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
            std::cout << "Testing validity of simplicial complex: ";
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
            
            std::cout << "PASSED" << std::endl;
        }
    };
    
}
