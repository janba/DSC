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
        
    private:
        Mesh mesh;
        
    public:
        template<typename real>
        ISMesh(std::vector<real> & points, std::vector<int> & tets)
        {
            vectors_read(points, tets, *this);
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
            SimplexSet<FaceKey> fids = get_faces(tids);
            for (auto f : fids)
            {
                if (exists(f))
                {
                    update_flag(f);
                }
            }
            
            // Update edges
            SimplexSet<EdgeKey> eids = get_edges(fids);
            for (auto e : eids)
            {
                if (exists(e))
                {
                    update_flag(e);
                }
            }
            
            // Update nodes
            SimplexSet<NodeKey> nids = get_nodes(eids);
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
            
            SimplexSet<TetrahedronKey> tids = get_tets(f);
            if (tids.size() == 1)
            {
                // On the boundary
                set_boundary(f, true);
                if (get_label(tids.front()) != 0)
                {
                    set_interface(f, true);
                }
            }
            else if(tids.size() == 2)
            {
                if (get_label(tids.front()) != get_label(tids.back()))
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
            
            int i = 0;
            for (auto f : get_faces(e))
            {
                if (exists(f))
                {
                    if (is_boundary(f))
                    {
                        set_boundary(e, true);
                    }
                    if (is_interface(f))
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
        
        void connected_component(SimplexSet<TetrahedronKey>& tids, const tet_key& tid)
        {
            int label = get_label(tid);
            tids -= tid;
            
            for(auto f : get_faces(tid))
            {
                tet_key tid2 = (get_tets(f) - tid).front();
                if(tids.contains(tid2) && label == get_label(tid2))
                {
                    connected_component(tids, tid2);
                }
            }
        }
        
        bool crossing(const node_key& n)
        {
            SimplexSet<TetrahedronKey> tids = get_tets(n);
            
            int c = 0;
            while (tids.size() > 0)
            {
                if(c == 2)
                {
                    return true;
                }
                tet_key tid = tids.front();
                connected_component(tids, tid);
                c++;
            }
            return false;
        }
        
        
        
        void update_flag(const node_key & n)
        {
            set_interface(n, false);
            set_boundary(n, false);
            set_crossing(n, false);
            
            for (auto e : get_edges(n))
            {
                if (exists(e))
                {
                    if (is_interface(e))
                    {
                        set_interface(n, true);
                    }
                    if (is_boundary(e))
                    {
                        set_boundary(n, true);
                    }
                    if (is_crossing(e))
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
        
        // Getters for getting the boundary/coboundary of a simplex:
        const SimplexSet<NodeKey>& get_nodes(const edge_key& eid)
        {
            return *get(eid).get_boundary();
        }
        
        const SimplexSet<EdgeKey>& get_edges(const node_key& nid)
        {
            return *get(nid).get_co_boundary();
        }
        
        const SimplexSet<EdgeKey>& get_edges(const face_key& fid)
        {
            return *get(fid).get_boundary();
        }
        
        const SimplexSet<FaceKey>& get_faces(const EdgeKey& eid)
        {
            return *get(eid).get_co_boundary();
        }
        
        const SimplexSet<FaceKey>& get_faces(const TetrahedronKey& tid)
        {
            return *get(tid).get_boundary();
        }
        
        const SimplexSet<TetrahedronKey>& get_tets(const FaceKey& fid)
        {
            return *get(fid).get_co_boundary();
        }
        
        // Getters for getting the boundary of a boundary etc.
        SimplexSet<NodeKey> get_sorted_nodes(const face_key& fid)
        {
            SimplexSet<NodeKey> nids;
            for (auto e : get_edges(fid)) {
                mesh.orient_face_helper(fid, e);
                nids += get_nodes(e)[0];
            }
            return nids;
        }
        
        SimplexSet<NodeKey> get_sorted_nodes(const tet_key& tid)
        {
            SimplexSet<NodeKey> nids;
            auto fids = get_faces(tid);
            mesh.orient_face_helper(tid, fids[0]);
            nids += get_sorted_nodes(fids[0]);
            nids += get_sorted_nodes(fids[1]);
            assert(nids.size() == 4);
            return nids;
        }
        
        SimplexSet<EdgeKey> get_sorted_edges(const tet_key& tid)
        {
            SimplexSet<EdgeKey> eids;
            for(auto f : get_faces(tid))
            {
                mesh.orient_face_helper(tid, f, true);
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
        
        SimplexSet<NodeKey> get_nodes(const FaceKey& fid)
        {
            return get_nodes(get_edges(fid));
        }
        
        SimplexSet<NodeKey> get_nodes(const tet_key& tid)
        {
            return get_nodes(get_edges(get_faces(tid)));
        }
        
        SimplexSet<EdgeKey> get_edges(const tet_key& tid)
        {
            return get_edges(get_faces(tid));
        }
        
        SimplexSet<FaceKey> get_faces(const NodeKey& nid)
        {
            return get_faces(get_edges(nid));
        }
        
        SimplexSet<TetrahedronKey> get_tets(const NodeKey& nid)
        {
            return get_tets(get_faces(get_edges(nid)));
        }
        
        SimplexSet<TetrahedronKey> get_tets(const EdgeKey& eid)
        {
            return get_tets(get_faces(eid));
        }
        
        // Getters which have a SimplexSet as input
        template<typename key_type>
        SimplexSet<NodeKey> get_nodes(const SimplexSet<key_type>& keys)
        {
            SimplexSet<NodeKey> nids;
            for(auto k : keys)
            {
                nids += get_nodes(k);
            }
            return nids;
        }
        
        template<typename key_type>
        SimplexSet<EdgeKey> get_edges(const SimplexSet<key_type>& keys)
        {
            SimplexSet<EdgeKey> eids;
            for(auto k : keys)
            {
                eids += get_edges(k);
            }
            return eids;
        }
        
        template<typename key_type>
        SimplexSet<FaceKey> get_faces(const SimplexSet<key_type>& keys)
        {
            SimplexSet<FaceKey> fids;
            for(auto k : keys)
            {
                fids += get_faces(k);
            }
            return fids;
        }
        
        template<typename key_type>
        SimplexSet<TetrahedronKey> get_tets(const SimplexSet<key_type>& keys)
        {
            SimplexSet<TetrahedronKey> tids;
            for(auto k : keys)
            {
                tids += get_tets(k);
            }
            return tids;
        }
        
        // Other getter functions
        
        NodeKey get_node(const EdgeKey& eid1, const EdgeKey& eid2)
        {
            SimplexSet<NodeKey> nid = get_nodes(eid1) & get_nodes(eid2);
            assert(nid.size() == 1);
            return nid.front();
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
        
        bool is_clockwise_order(const node_key& nid, SimplexSet<node_key>& nids)
        {
            auto x = get(nid).get_pos() - get(nids[0]).get_pos();
            auto y = get(nids[1]).get_pos() - get(nids[0]).get_pos();
            auto z = get(nids[2]).get_pos() - get(nids[0]).get_pos();
            auto val = dot(x, cross(y,z));
            
#ifdef DEBUG
            assert(val != 0.);
#endif
            return val > 0.;
        }
        
        /*
         * Orient the nodes in a counter clockwise order seen from the node a.
         */
        void orient_cc(const node_key& nid, SimplexSet<node_key>& nids)
        {
            if(is_clockwise_order(nid, nids))
            {
                std::swap(nids[0], nids[2]);
            }
        }
        
        /**
         * Ensures consistent orientation of the face fid if fid is an interface or boundary face.
         */
        void orient_face(const face_key& fid)
        {
            if (is_interface(fid))
            {
                int label = -100;
                tet_key tid;
                for (auto t : get_tets(fid))
                {
                    int tl = get_label(t);
                    if (tl > label)
                    {
                        label = tl;
                        tid = t;
                    }
                }
                mesh.orient_face_helper(tid, fid);
            }
            else if (is_boundary(fid))
            {
                mesh.orient_face_helper(get_tets(fid).front(), fid);
            }
        }
        
        node_key split(const edge_key& eid, const typename node_traits::vec3& pos, const typename node_traits::vec3& destination)
        {
            auto nids = get_nodes(eid);
            auto fids = get_faces(eid);
            auto tids = get_tets(eid);
            
            // Split edge
            auto new_nid = insert_node();
            get(new_nid).set_pos(pos);
            get(new_nid).set_destination(destination);
            
            disconnect(nids[1], eid);
            connect(new_nid, eid);
            
            auto new_eid = insert_edge(new_nid, nids[1]);
            
            // Update faces, create faces
            for (auto f : fids)
            {
                EdgeKey f_eid = (get_edges(f) & get_edges(nids[1])).front();
                disconnect(f_eid, f);
                
                SimplexSet<NodeKey> new_e_nids = get_nodes(get_edges(f)) - nids[0];
                assert(new_e_nids.size() == 2);
                EdgeKey new_f_eid = insert_edge(new_e_nids[0], new_e_nids[1]);
                connect(new_f_eid, f);
                
                insert_face(f_eid, new_f_eid, new_eid);
            }
            
            // Update tetrahedra, create tetrahedra
            SimplexSet<TetrahedronKey> new_tids;
            for (auto t : tids)
            {
                FaceKey t_fid = (get_faces(t) - get_faces(nids[0])).front();
                disconnect(t_fid, t);
                
                SimplexSet<EdgeKey> new_f_eids = get_edges(get_faces(t)) - get_edges(nids[0]);
                assert(new_f_eids.size() == 3);
                FaceKey new_t_fid = insert_face(new_f_eids[0], new_f_eids[1], new_f_eids[2]);
                connect(new_t_fid, t);
                
                SimplexSet<FaceKey> t_fids = get_faces(new_eid) & get_faces(get_edges(t_fid));
                assert(t_fids.size() == 2);
                new_tids += insert_tetrahedron(t_fids[0], t_fids[1], new_t_fid, t_fid);
            }
            
            // Update flags
            for (int i = 0; i < tids.size(); i++)
            {
                set_label(new_tids[i], get_label(tids[i]));
            }
            
            for(auto t : tids)
            {
                if(is_inverted(t))
                {
                    get(t).invert_orientation();
                }
            }
            
            return new_nid;
        }
        
        
        node_key collapse(const edge_key& eid, const typename node_traits::vec3& pos, const typename node_traits::vec3& destination)
        {
            auto nids = get_nodes(eid);
            auto fids = get_faces(eid);
            auto tids = get_tets(eid);
            
            // Remove edge
            mesh.remove(eid);
            
            // Merge nodes
            node_key nid = merge(nids[1], nids[0]);
            get(nid).set_pos(pos);
            get(nid).set_destination(destination);
            
            // Remove faces
            for(auto f : fids)
            {
                auto eids = get_edges(f);
                mesh.remove(f);
                merge(eids[0], eids[1]);
            }
            
            // Remove tetrahedra
            for(auto t : tids)
            {
                auto fids = get_faces(t);
                mesh.remove(t);
                merge(fids[0], fids[1]);
            }
            
            // Update flags and ensure no inverted tetrahedra.
            SimplexSet<tet_key> changed_tids = get_tets(nid);
            for(auto t : changed_tids)
            {
                if(is_inverted(t))
                {
                    get(t).invert_orientation();
                }
            }
            update(changed_tids);
            return nid;
        }
        
        /**
         * Returns whether the tetrahedron with ID tid is inverted.
         */
        bool is_inverted(const tet_key& tid)
        {
            auto nodes = get_sorted_nodes(tid);
            std::vector<typename node_traits::vec3> verts;
            for(auto &n : nodes)
            {
                verts.push_back(get(n).get_pos());
            }
            return dot(verts[0]-verts[3], cross(verts[1]-verts[3], verts[2]-verts[3])) < 0.;
        }
        
    private:
        
        NodeKey merge(const NodeKey& key1, const NodeKey& key2)
        {
            for(auto e : get_edges(key2))
            {
                connect(key1, e);
            }
            mesh.remove(key2);
            return key1;
        }
        
        template<typename key_type>
        key_type merge(const key_type& key1, const key_type& key2)
        {
            auto& simplex = get(key2);
            for(auto k : *simplex.get_co_boundary())
            {
                connect(key1, k);
            }
            for(auto k : *simplex.get_boundary())
            {
                connect(k, key1);
            }
            mesh.remove(key2);
            return key1;
        }
        
    public:
        /**
         * Inserts a node into the mesh. Trivial.
         */
        NodeKey insert_node()
        {
            auto node = mesh.m_node_kernel->create();
            return node.key();
        }
        
        /**
         * Inserts an edge into the mesh. Updates the co-boundary of the boundary nodes with the newly created edge.
         * Leaves the closure of the edge in an uncompressed state.
         */
        EdgeKey insert_edge(NodeKey node1, NodeKey node2)
        {
            auto edge = mesh.m_edge_kernel->create();
            //add the new simplex to the co-boundary relation of the boundary simplices
            get(node1).add_co_face(edge.key());
            get(node2).add_co_face(edge.key());
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
            auto face = mesh.m_face_kernel->create();
            //update relations
            get(edge1).add_co_face(face.key());
            get(edge2).add_co_face(face.key());
            get(edge3).add_co_face(face.key());
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
            auto tetrahedron = mesh.m_tetrahedron_kernel->create();
            //update relations
            get(face1).add_co_face(tetrahedron.key());
            get(face2).add_co_face(tetrahedron.key());
            get(face3).add_co_face(tetrahedron.key());
            get(face4).add_co_face(tetrahedron.key());
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
        
    public:
        
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
            auto new_fid = insert_face(f_eids[0], f_eids[1], f_eids[2]);
            
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
            edge_key new_eid = insert_edge(e_nids[0], e_nids[1]);
            
            // Create faces
            auto new_fs_eids = get_boundary(get_boundary(f_tids)) - f_eids;
            assert(new_fs_eids.size() == 6);
            for (node_key& n : f_nids)
            {
                auto new_f_eids = new_fs_eids & get_co_boundary(n);
                assert(new_f_eids.size() == 2);
                insert_face(new_f_eids[0], new_f_eids[1], new_eid);
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
                        assert(f == f2 || get_edge(f, f2).is_valid());
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
                            assert(e == e2 || get_node(e, e2).is_valid());
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
