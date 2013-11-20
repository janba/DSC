//
//  Deformabel Simplicial Complex (DSC) method
//  Copyright (C) 2013  Technical University of Denmark
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  See licence.txt for a copy of the GNU General Public License.

#pragma once

#include <is_mesh/kernel.h>
#include <is_mesh/simplex.h>
#include <is_mesh/simplex_set.h>
#include <is_mesh/is_mesh_lists_read.h>

namespace is_mesh {
    
    template <typename node_traits, typename edge_traits, typename face_traits, typename tet_traits>
    class ISMesh
    {
        typedef Node<node_traits>                                         node_type;
        typedef Edge<edge_traits>                                         edge_type;
        typedef Face<face_traits>                                         face_type;
        typedef Tetrahedron<tet_traits>                           tetrahedron_type;
        
        kernel<node_type, NodeKey>* m_node_kernel;
        kernel<edge_type, EdgeKey>*                  m_edge_kernel;
        kernel<face_type, FaceKey>*                  m_face_kernel;
        kernel<tetrahedron_type, TetrahedronKey>*           m_tetrahedron_kernel;
        
    public:
        template<typename real>
        ISMesh(std::vector<real> & points, std::vector<int> & tets, const std::vector<int>& tet_labels = std::vector<int>())
        {
            m_node_kernel = new kernel<node_type, NodeKey>();
            m_edge_kernel = new kernel<edge_type, EdgeKey>();
            m_face_kernel = new kernel<face_type, FaceKey>();
            m_tetrahedron_kernel = new kernel<tetrahedron_type, TetrahedronKey>();
            
            vectors_read(points, tets, *this);
            init_flags(tet_labels);
            validity_check();
            simplex_set_test();
        }
        
        ~ISMesh()
        {
            delete m_tetrahedron_kernel;
            delete m_face_kernel;
            delete m_edge_kernel;
            delete m_node_kernel;
        }
        
        ///////////////
        // ITERATORS //
        ///////////////
    public:
        typename kernel<node_type, NodeKey>::iterator nodes_begin()
        {
            return m_node_kernel->begin();
        }
        
        typename kernel<node_type, NodeKey>::iterator nodes_end()
        {
            return m_node_kernel->end();
        }
        
        typename kernel<edge_type, EdgeKey>::iterator edges_begin()
        {
            return m_edge_kernel->begin();
        }
        
        typename kernel<edge_type, EdgeKey>::iterator edges_end()
        {
            return m_edge_kernel->end();
        }
        
        typename kernel<face_type, FaceKey>::iterator faces_begin()
        {
            return m_face_kernel->begin();
        }
        
        typename kernel<face_type, FaceKey>::iterator faces_end()
        {
            return m_face_kernel->end();
        }
        
        typename kernel<tetrahedron_type, TetrahedronKey>::iterator tetrahedra_begin()
        {
            return m_tetrahedron_kernel->begin();
        }
        
        typename kernel<tetrahedron_type, TetrahedronKey>::iterator tetrahedra_end()
        {
            return m_tetrahedron_kernel->end();
        }
        
        /////////////////////
        // LABEL FUNCTIONS //
        /////////////////////
    public:
        
        int get_label(const TetrahedronKey& t)
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
        void set_label(const TetrahedronKey& tid, int label)
        {
            get(tid).label(label);
            SimplexSet<TetrahedronKey> tids = {tid};
            update(tids);
        }
        
    private:
        /**
         * Perform an initial update of flags for all nodes, edges and faces.
         */
        void init_flags(const std::vector<int>& tet_labels)
        {
            if(tet_labels.size() > 0)
            {
                for (auto tit = tetrahedra_begin(); tit != tetrahedra_end(); tit++)
                {
                    tit->label(tet_labels[tit.key()]);
                }
            }
            
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
        void update(const SimplexSet<TetrahedronKey>& tids)
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
        
        void update_flag(const FaceKey & f)
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
        
        void update_flag(const EdgeKey & e)
        {
            set_boundary(e, false);
            set_interface(e, false);
            set_crossing(e, false);
            
            int i = 0;
            for (auto f : get_faces(e))
            {
                if (exists(f))
                {
                    if (get(f).is_boundary())
                    {
                        set_boundary(e, true);
                    }
                    if (get(f).is_interface())
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
        
        void connected_component(SimplexSet<TetrahedronKey>& tids, const TetrahedronKey& tid)
        {
            int label = get_label(tid);
            tids -= tid;
            
            for(auto f : get_faces(tid))
            {
                if(get_tets(f).size() == 2)
                {
                    TetrahedronKey tid2 = get_tet(tid, f);
                    if(tids.contains(tid2) && label == get_label(tid2))
                    {
                        connected_component(tids, tid2);
                    }
                }
            }
        }
        
        bool crossing(const NodeKey& n)
        {
            SimplexSet<TetrahedronKey> tids = get_tets(n);
            
            int c = 0;
            while (tids.size() > 0)
            {
                if(c == 2)
                {
                    return true;
                }
                TetrahedronKey tid = tids.front();
                connected_component(tids, tid);
                c++;
            }
            return false;
        }
        
        
        
        void update_flag(const NodeKey & n)
        {
            set_interface(n, false);
            set_boundary(n, false);
            set_crossing(n, false);
            
            for (auto e : get_edges(n))
            {
                if (exists(e))
                {
                    if (get(e).is_interface())
                    {
                        set_interface(n, true);
                    }
                    if (get(e).is_boundary())
                    {
                        set_boundary(n, true);
                    }
                    if (get(e).is_crossing())
                    {
                        set_crossing(n, true);
                    }
                }
            }
            if(!get(n).is_crossing() && get(n).is_interface() && crossing(n))
            {
                set_crossing(n, true);
            }
        }
        
        //////////////////////
        // GETTER FUNCTIONS //
        //////////////////////
    public:
        node_type & get(const NodeKey& nid)
        {
            return m_node_kernel->find(nid);
        }
        
        edge_type & get(const EdgeKey& eid)
        {
            return m_edge_kernel->find(eid);
        }
        
        face_type & get(const FaceKey& fid)
        {
            return m_face_kernel->find(fid);
        }
        
        tetrahedron_type & get(const TetrahedronKey& tid)
        {
            return m_tetrahedron_kernel->find(tid);
        }
        
        // Getters for getting the boundary/coboundary of a simplex:
        const SimplexSet<NodeKey>& get_nodes(const EdgeKey& eid)
        {
            return get(eid).get_boundary();
        }
        
        const SimplexSet<EdgeKey>& get_edges(const NodeKey& nid)
        {
            return get(nid).get_co_boundary();
        }
        
        const SimplexSet<EdgeKey>& get_edges(const FaceKey& fid)
        {
            return get(fid).get_boundary();
        }
        
        const SimplexSet<FaceKey>& get_faces(const EdgeKey& eid)
        {
            return get(eid).get_co_boundary();
        }
        
        const SimplexSet<FaceKey>& get_faces(const TetrahedronKey& tid)
        {
            return get(tid).get_boundary();
        }
        
        const SimplexSet<TetrahedronKey>& get_tets(const FaceKey& fid)
        {
            return get(fid).get_co_boundary();
        }
        
        // Getters for getting the boundary of a boundary etc.
        
        SimplexSet<NodeKey> get_sorted_nodes(const FaceKey& fid, const TetrahedronKey& tid)
        {
            SimplexSet<NodeKey> nids = get_nodes(fid);
            NodeKey apex = (get_nodes(tid) - nids).front();
            orient_cc(apex, nids);
            return nids;
        }
        
        SimplexSet<NodeKey> get_sorted_nodes(const FaceKey& fid)
        {
            if (get(fid).is_interface())
            {
                int label = -100;
                TetrahedronKey tid;
                for (auto t : get_tets(fid))
                {
                    int tl = get_label(t);
                    if (tl > label)
                    {
                        label = tl;
                        tid = t;
                    }
                }
                return get_sorted_nodes(fid, tid);
            }
            else if (get(fid).is_boundary())
            {
                TetrahedronKey tid = get_tets(fid).front();
                return get_sorted_nodes(fid, tid);
            }
            return get_nodes(fid);
        }
        
        SimplexSet<NodeKey> get_nodes(const FaceKey& fid)
        {
            const SimplexSet<EdgeKey>& eids = get_edges(fid);
            SimplexSet<NodeKey> nids = get_nodes(eids[0]);
            nids += get_nodes(eids[1]);
            return nids;
        }
        
        SimplexSet<NodeKey> get_nodes(const TetrahedronKey& tid)
        {
            const SimplexSet<FaceKey>& fids = get_faces(tid);
            SimplexSet<NodeKey> nids = get_nodes(fids[0]);
            nids += get_nodes(fids[1]);
            return nids;
        }
        
        SimplexSet<EdgeKey> get_edges(const TetrahedronKey& tid)
        {
            SimplexSet<EdgeKey> eids;
            for(const FaceKey& f : get_faces(tid))
            {
                eids += get_edges(f);
            }
            return eids;
        }
        
        SimplexSet<FaceKey> get_faces(const NodeKey& nid)
        {
            SimplexSet<FaceKey> fids;
            for(const EdgeKey& e : get_edges(nid))
            {
                fids += get_faces(e);
            }
            return fids;
        }
        
        SimplexSet<TetrahedronKey> get_tets(const NodeKey& nid)
        {
            SimplexSet<TetrahedronKey> tids;
            for(const EdgeKey& e : get_edges(nid))
            {
                for(const FaceKey& f : get_faces(e))
                {
                    tids += get_tets(f);
                }
            }
            return tids;
        }
        
        SimplexSet<TetrahedronKey> get_tets(const EdgeKey& eid)
        {
            SimplexSet<TetrahedronKey> tids;
            for(const FaceKey& f : get_faces(eid))
            {
                tids += get_tets(f);
            }
            return tids;
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
        
        /**
         *  Returns the node shared by the edges eid1 and eid2.
         */
        NodeKey get_node(const EdgeKey& eid1, const EdgeKey& eid2)
        {
            const SimplexSet<NodeKey>& nids = get_nodes(eid2);
            for (auto& n : get_nodes(eid1)) {
                if(nids.contains(n))
                {
                    return n;
                }
            }
            return NodeKey();
        }
        
        /**
         *  Returns the edge between the nodes nid1 and nid2.
         */
        EdgeKey get_edge(const NodeKey& nid1, const NodeKey& nid2)
        {
            const SimplexSet<EdgeKey>& eids = get_edges(nid2);
            for (auto& e : get_edges(nid1)) {
                if(eids.contains(e))
                {
                    return e;
                }
            }
            return EdgeKey();
        }
        
        /**
         *  Returns the edge shared by the faces fid1 and fid2.
         */
        EdgeKey get_edge(const FaceKey& fid1, const FaceKey& fid2)
        {
            const SimplexSet<EdgeKey>& eids = get_edges(fid2);
            for (auto& e : get_edges(fid1)) {
                if(eids.contains(e))
                {
                    return e;
                }
            }
            return EdgeKey();
        }
        
        /**
         *  Returns the face between the nodes nid1, nid2 and nid3.
         */
        FaceKey get_face(const NodeKey& nid1, const NodeKey& nid2, const NodeKey& nid3)
        {
            SimplexSet<FaceKey> fid = (get_faces(nid1) & get_faces(nid2)) & get_faces(nid3);
            assert(fid.size() == 1);
            return fid.front();
        }
        
        /**
         *  Returns the face shared by the tetrahedra tid1 and tid2.
         */
        FaceKey get_face(const TetrahedronKey& tid1, const TetrahedronKey& tid2)
        {
            const SimplexSet<EdgeKey>& fids = get_faces(tid2);
            for (auto& f : get_faces(tid1)) {
                if(fids.contains(f))
                {
                    return f;
                }
            }
            assert(false);
            return FaceKey();
        }
        
        /**
         *  Returns the tetrahedron which shares the face fid with tid, i.e. the neighbour to tid.
         */
        TetrahedronKey get_tet(const TetrahedronKey& tid, const FaceKey& fid)
        {
            for(TetrahedronKey t : get_tets(fid))
            {
                if(t != tid)
                {
                    return t;
                }
            }
            return TetrahedronKey();
        }
        
        /**
         * Returns the position of node nid.
         */
        vec3 get_pos(const NodeKey& nid)
        {
            return get(nid).get_pos();
        }
        
        /**
         * Returns the positions of nodes nids.
         */
        std::vector<vec3> get_pos(const SimplexSet<NodeKey>& nids)
        {
            std::vector<vec3> verts(nids.size());
            for (int i = 0; i < nids.size(); i++)
            {
                verts[i] = get_pos(nids[i]);
            }
            return verts;
        }
        
        //////////////////////
        // EXISTS FUNCTIONS //
        //////////////////////
    public:
        
        /**
         *
         */
        bool exists(const TetrahedronKey& t)
        {
            return m_tetrahedron_kernel->is_valid(t);
        }
        
        /**
         *
         */
        bool exists(const FaceKey& f)
        {
            return m_face_kernel->is_valid(f);
        }
        
        /**
         *
         */
        bool exists(const EdgeKey& e)
        {
            return m_edge_kernel->is_valid(e);
        }
        
        /**
         *
         */
        bool exists(const NodeKey& n)
        {
            return m_node_kernel->is_valid(n);
        }
        
        
        ///////////////////////////
        // ORIENTATION FUNCTIONS //
        ///////////////////////////
        
    public:
        bool is_clockwise_order(const NodeKey& nid, SimplexSet<NodeKey>& nids)
        {
            auto x = get(nid).get_pos() - get(nids[0]).get_pos();
            auto y = get(nids[1]).get_pos() - get(nids[0]).get_pos();
            auto z = get(nids[2]).get_pos() - get(nids[0]).get_pos();
            auto val = dot(x, cross(y,z));
            
            return val > 0.;
        }
        
        /*
         * Orient the nodes in a counter clockwise order seen from the node a.
         */
        void orient_cc(const NodeKey& nid, SimplexSet<NodeKey>& nids)
        {
            if(is_clockwise_order(nid, nids))
            {
                nids.swap();
            }
        }
        
        /**
         * Returns whether the tetrahedron with ID tid is inverted.
         */
        bool is_inverted(const TetrahedronKey& tid)
        {
            for(const FaceKey& f : get_faces(tid))
            {
                const SimplexSet<TetrahedronKey>& tids = get_tets(f);
                if(tids.size() == 2) // Check that f is not a boundary face.
                {
                    SimplexSet<NodeKey> nids = get_nodes(f);
                    SimplexSet<NodeKey> apices = get_nodes(tids) - nids;
                    auto normal = cross(get_pos(nids[0]) - get_pos(nids[2]), get_pos(nids[1]) - get_pos(nids[2]));
                    auto d1 = dot(get_pos(apices[0]) - get_pos(nids[2]), normal);
                    auto d2 = dot(get_pos(apices[1]) - get_pos(nids[2]), normal);
                    if((d1 < 0. && d2 < 0) || (d1 > 0. && d2 > 0.))
                    {
                        return true;
                    }
                }
            }
            return false;
        }
        
        ////////////////////
        // MESH FUNCTIONS //
        ////////////////////
    
    public:
        /**
         * Inserts a node into the mesh. Trivial.
         */
        template<typename vec3>
        NodeKey insert_node(const vec3& p)
        {
            auto node = m_node_kernel->create(node_traits(p));
            return node.key();
        }
        
        /**
         * Inserts an edge into the mesh. Updates the co-boundary of the boundary nodes with the newly created edge.
         * Leaves the closure of the edge in an uncompressed state.
         */
        EdgeKey insert_edge(NodeKey node1, NodeKey node2)
        {
            auto edge = m_edge_kernel->create(edge_traits());
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
            auto face = m_face_kernel->create(face_traits());
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
            auto tetrahedron = m_tetrahedron_kernel->create(tet_traits());
            //update relations
            get(face1).add_co_face(tetrahedron.key());
            get(face2).add_co_face(tetrahedron.key());
            get(face3).add_co_face(tetrahedron.key());
            get(face4).add_co_face(tetrahedron.key());
            tetrahedron->add_face(face1);
            tetrahedron->add_face(face2);
            tetrahedron->add_face(face3);
            tetrahedron->add_face(face4);
            
            return tetrahedron.key();
        }
        
    private:
        /**
         *
         */
        void remove(const NodeKey& nid)
        {
            for(auto e : get_edges(nid))
            {
                get(e).remove_face(nid);
            }
            m_node_kernel->erase(nid);
        }
        
        /**
         *
         */
        void remove(const EdgeKey& eid)
        {
            for(auto f : get_faces(eid))
            {
                get(f).remove_face(eid);
            }
            for(auto n : get_nodes(eid))
            {
                get(n).remove_co_face(eid);
            }
            m_edge_kernel->erase(eid);
        }
        
        /**
         *
         */
        void remove(const FaceKey& fid)
        {
            for(auto t : get_tets(fid))
            {
                get(t).remove_face(fid);
            }
            for(auto e : get_edges(fid))
            {
                get(e).remove_co_face(fid);
            }
            m_face_kernel->erase(fid);
        }
        
        /**
         *
         */
        void remove(const TetrahedronKey& tid)
        {
            for(auto f : get_faces(tid))
            {
                get(f).remove_co_face(tid);
            }
            m_tetrahedron_kernel->erase(tid);
        }
        
        NodeKey merge(const NodeKey& key1, const NodeKey& key2)
        {
            for(auto e : get_edges(key2))
            {
                connect(key1, e);
            }
            remove(key2);
            return key1;
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
            if(!get(pk1).get_boundary().contains(ck1))
            {
                assert(get(pk1).get_boundary().contains(ck2));
                assert(get(pk2).get_boundary().contains(ck1));
                
                disconnect(ck1, pk2);
                disconnect(ck2, pk1);
                connect(ck1, pk1);
                connect(ck2, pk2);
            }
            else {
                assert(get(pk1).get_boundary().contains(ck1));
                assert(get(pk2).get_boundary().contains(ck2));
                
                disconnect(ck1, pk1);
                disconnect(ck2, pk2);
                connect(ck1, pk2);
                connect(ck2, pk1);
            }
        }
        
        template<typename key_type>
        key_type merge(const key_type& key1, const key_type& key2)
        {
            auto& simplex = get(key2);
            for(auto k : simplex.get_co_boundary())
            {
                connect(key1, k);
            }
            for(auto k : simplex.get_boundary())
            {
                connect(k, key1);
            }
            remove(key2);
            return key1;
        }
        
    protected:
        
        virtual void split(const NodeKey& nid_new, const NodeKey& nid1, const NodeKey& nid2)
        {
            
        }
        
        NodeKey split(const EdgeKey& eid, const vec3& pos, const vec3& destination)
        {
            auto nids = get_nodes(eid);
            auto fids = get_faces(eid);
            auto tids = get_tets(eid);
            
            // Split edge
            auto new_nid = insert_node(pos);
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
            
            split(new_nid, nids[0], nids[1]);
            
            return new_nid;
        }
        
        virtual void collapse(const NodeKey& nid_new)
        {
            
        }
        
        NodeKey collapse(const EdgeKey& eid, const vec3& pos, const vec3& destination)
        {
            auto nids = get_nodes(eid);
            auto fids = get_faces(eid);
            auto tids = get_tets(eid);
            
            // Remove edge
            remove(eid);
            
            // Merge nodes
            NodeKey nid = merge(nids[1], nids[0]);
            get(nid).set_pos(pos);
            get(nid).set_destination(destination);
            
            // Remove faces
            for(auto f : fids)
            {
                auto eids = get_edges(f);
                remove(f);
                merge(eids[0], eids[1]);
            }
            
            // Remove tetrahedra
            for(auto t : tids)
            {
                auto fids = get_faces(t);
                remove(t);
                merge(fids[0], fids[1]);
            }
            
            // Update flags.
            update(get_tets(nid));
            
            collapse(nid);
            
            return nid;
        }
        
        FaceKey flip_32(const EdgeKey& eid)
        {
            SimplexSet<NodeKey> e_nids = get_nodes(eid);
            SimplexSet<FaceKey> e_fids = get_faces(eid);
            assert(e_fids.size() == 3);
            SimplexSet<TetrahedronKey> e_tids = get_tets(e_fids);
            assert(e_tids.size() == 3);
            int label = get_label(e_tids[0]);
            assert(label == get_label(e_tids[1]));
            assert(label == get_label(e_tids[2]));
            
            // Remove edge
            remove(eid);
            
            // Create face
            SimplexSet<EdgeKey> f_eids = get_edges(e_tids) - get_edges(e_fids);
            assert(f_eids.size() == 3);
            FaceKey new_fid = insert_face(f_eids[0], f_eids[1], f_eids[2]);
            
            // Remove faces
            for(const FaceKey& f : e_fids)
            {
                remove(f);
            }
            
            // Create tetrahedra
            SimplexSet<FaceKey> exterior_fids = get_faces(e_tids);
            for (const NodeKey& n : e_nids)
            {
                SimplexSet<FaceKey> t_fids = exterior_fids & get_faces(n);
                assert(t_fids.size() == 3);
                insert_tetrahedron(t_fids[0], t_fids[1], t_fids[2], new_fid);
            }
            
            // Remove tetrahedra
            for(const TetrahedronKey& t : e_tids)
            {
                remove(t);
            }
            
            // Update flags
            assert(get_tets(new_fid).size() == 2);
            for (auto t : get_tets(new_fid)) {
                set_label(t, label);
            }
            return new_fid;
        }
        
        EdgeKey flip_23(const FaceKey& fid)
        {
            SimplexSet<TetrahedronKey> f_tids = get_tets(fid);
            assert(f_tids.size() == 2);
            SimplexSet<EdgeKey> f_eids = get_edges(fid);
            assert(f_eids.size() == 3);
            SimplexSet<NodeKey> f_nids = get_nodes(f_eids);
            int label = get_label(f_tids[0]);
            assert(label == get_label(f_tids[1]));
            
            // Create edge
            SimplexSet<NodeKey> e_nids = get_nodes(f_tids) - f_nids;
            assert(e_nids.size() == 2);
            EdgeKey new_eid = insert_edge(e_nids[0], e_nids[1]);
            
            // Create faces
            SimplexSet<EdgeKey> new_fs_eids = get_edges(f_tids) - f_eids;
            assert(new_fs_eids.size() == 6);
            for (const NodeKey& n : f_nids)
            {
                auto new_f_eids = new_fs_eids & get_edges(n);
                assert(new_f_eids.size() == 2);
                insert_face(new_f_eids[0], new_f_eids[1], new_eid);
            }
            
            // Remove face
            remove(fid);
            
            // Create tetrahedra
            SimplexSet<FaceKey> new_ts_fids1 = get_faces(f_tids);
            SimplexSet<FaceKey> new_ts_fids2 = get_faces(new_eid);
            assert(new_ts_fids1.size() == 6);
            assert(new_ts_fids2.size() == 3);
            for (const EdgeKey& e : f_eids)
            {
                SimplexSet<FaceKey> new_t_fids = (new_ts_fids1 & get_faces(e)) + (new_ts_fids2 & get_faces(get_nodes(e)));
                assert(new_t_fids.size() == 4);
                insert_tetrahedron(new_t_fids[0], new_t_fids[1], new_t_fids[2], new_t_fids[3]);
            }
            
            // Remove tetrahedra
            for (auto t : f_tids) {
                remove(t);
            }
            
            // Update flags
            assert(get_tets(new_eid).size() == 3);
            for (auto t : get_tets(new_eid)) {
                set_label(t, label);
            }
            return new_eid;
        }
        
        void flip(const EdgeKey& eid, const FaceKey& fid1, const FaceKey& fid2)
        {
            SimplexSet<FaceKey> fids = {fid1, fid2};
            SimplexSet<NodeKey> e_nids = get_nodes(eid);
            SimplexSet<FaceKey> e_fids = get_faces(eid);
            SimplexSet<TetrahedronKey> e_tids = get_tets(e_fids);
            
            // Reconnect edge
            SimplexSet<NodeKey> new_e_nids = get_nodes(fids) - e_nids;
            assert(new_e_nids.size() == 2);
            
            disconnect(e_nids[0], eid);
            disconnect(e_nids[1], eid);
            connect(new_e_nids[0], eid);
            connect(new_e_nids[1], eid);
            
            // Reconnect faces
            SimplexSet<EdgeKey> swap_eids = {get_edge(e_nids[0], new_e_nids[0]), get_edge(e_nids[1], new_e_nids[1])};
            swap(swap_eids[0], fid1, swap_eids[1], fid2);
            
            assert((e_fids - fids).size() <= 2);
            for(FaceKey f : (e_fids - fids))
            {
                SimplexSet<EdgeKey> rm_eids = get_edges(f) - eid;
                assert(rm_eids.size() == 2);
                SimplexSet<NodeKey> apex = get_nodes(f) - (new_e_nids + e_nids);
                assert(apex.size() == 1);
                
                SimplexSet<EdgeKey> add_eids = {get_edge(apex[0], new_e_nids[0]), get_edge(apex[0], new_e_nids[1])};
                
                disconnect(rm_eids[0], f);
                disconnect(rm_eids[1], f);
                connect(add_eids[0], f);
                connect(add_eids[1], f);
                
                // Reconnect tetrahedra
                SimplexSet<TetrahedronKey> tids = get_tets(f);
                assert(tids.size() == 2);
                SimplexSet<FaceKey> swap_fids = (get_faces(tids) - e_fids) & get_faces(swap_eids);
                assert(swap_fids.size() == 2);
                swap(swap_fids[0], tids[0], swap_fids[1], tids[1]);
            }
            
            // Update flags
            update(e_tids);
        }
        
        
        void flip_22(const FaceKey& fid1, const FaceKey& fid2)
        {
            SimplexSet<EdgeKey> eid = get_edges(fid1) & get_edges(fid2);
            assert(eid.size() == 1);
            assert(get_faces(eid).size() == 3);
            assert(get_tets(eid).size() == 2);
            
            flip(eid[0], fid1, fid2);
        }
        
        void flip_44(const FaceKey& fid1, const FaceKey& fid2)
        {
            SimplexSet<EdgeKey> eid = get_edges(fid1) & get_edges(fid2);
            assert(eid.size() == 1);
            assert(get_faces(eid).size() == 4);
            assert(get_tets(eid).size() == 4);
            
            flip(eid[0], fid1, fid2);
        }
        
        ///////////////////////
        // UTILITY FUNCTIONS //
        ///////////////////////
    public:
        
        void garbage_collect()
        {
            m_node_kernel->garbage_collect();
            m_edge_kernel->garbage_collect();
            m_face_kernel->garbage_collect();
            m_tetrahedron_kernel->garbage_collect();
        }
        
        void extract_surface_mesh(std::vector<vec3>& verts, std::vector<int>& indices)
        {
            garbage_collect();
            
            std::map<NodeKey, int> vert_index;
            
            // Extract vertices
            for (auto nit = nodes_begin(); nit != nodes_end(); nit++)
            {
                if (nit->is_interface())
                {
                    verts.push_back(get_pos(nit.key()));
                    vert_index[nit.key()] = static_cast<int>(verts.size());
                }
            }
            
            // Extract faces
            for (auto fit = faces_begin(); fit != faces_end(); fit++)
            {
                if (fit->is_interface())
                {
                    SimplexSet<NodeKey> nodes = get_sorted_nodes(fit.key());
                    
                    indices.push_back(vert_index[nodes[0]]);
                    indices.push_back(vert_index[nodes[1]]);
                    indices.push_back(vert_index[nodes[2]]);
                }
            }
        }
        
        void extract_tet_mesh(std::vector<vec3>& points, std::vector< std::vector<int> >& tets)
        {
            garbage_collect();
            
            std::map<NodeKey, int> indices;
            int counter = 0;
            for (auto nit = nodes_begin(); nit != nodes_end(); nit++)
            {
                indices[nit.key()] = counter;
                points.push_back(nit->get_pos());
                ++counter;
            }
            
            for (auto tit = tetrahedra_begin(); tit != tetrahedra_end(); tit++)
            {
                std::vector<int> tet;
                auto nodes = get_nodes(tit.key());
                
                for (auto &n : nodes)
                {
                    tet.push_back(indices[n]);
                }
                tet.push_back(get_label(tit.key()));
                tets.push_back(tet);
            }
        }
                
        void validity_check()
        {
            std::cout << "Testing connectivity of simplicial complex: ";
            for(auto tit = tetrahedra_begin(); tit != tetrahedra_end(); tit++)
            {
                assert(exists(tit.key()));
                // Check faces:
                auto faces = get_faces(tit.key());
                assert(faces.size() == 4);
                for (auto f : faces) {
                    assert(exists(f));
                    auto cotets = get_tets(f);
                    assert((get(f).is_boundary() && cotets.size() == 1) || (!get(f).is_boundary() && cotets.size() == 2));
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
            
            std::cout << "Testing for inverted tetrahedra: ";
            for (auto tit = tetrahedra_begin(); tit != tetrahedra_end(); tit++)
            {
                assert(!is_inverted(tit.key()));
            }
            std::cout << "PASSED" << std::endl;
            
            std::cout << "Testing for corrupted interface or boundary: ";
            for (auto eit = edges_begin(); eit != edges_end(); eit++)
            {
                int boundary = 0;
                int interface = 0;
                for (auto f : get_faces(eit.key())) {
                    if(get(f).is_boundary())
                    {
                        boundary++;
                    }
                    if(get(f).is_interface())
                    {
                        interface++;
                    }
                }
                assert((eit->is_interface() && interface >= 2) || (!eit->is_interface() && interface == 0)); // Check that the interface is not corrupted
                assert((eit->is_boundary() && boundary == 2) || (!eit->is_boundary() && boundary == 0)); // Check that the boundary is not corrupted
            }
            std::cout << "PASSED" << std::endl;
        }
    };
    
}
