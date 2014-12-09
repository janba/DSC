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

#include "util.h"
#include "kernel.h"
#include "simplex.h"
#include "simplex_set.h"

namespace is_mesh {

    template <typename node_traits>
    class NodeIterator {
        const kernel<Node, NodeKey> *m_node_kernel;
    public:
        NodeIterator(const kernel<Node, NodeKey> *m_node_kernel) : m_node_kernel(m_node_kernel) {
        }

        typename kernel<Node, NodeKey>::iterator begin() const {
            return m_node_kernel->begin();
        }

        typename kernel<Node, NodeKey>::iterator end() const {
            return m_node_kernel->end();
        }
    };


    template <typename edge_traits>
    class EdgeIterator {
        const kernel<Edge, EdgeKey>*                  m_edge_kernel;
    public:
        EdgeIterator(const kernel<Edge, EdgeKey> *m_edge_kernel) : m_edge_kernel(m_edge_kernel) {
        }

        typename kernel<Edge, EdgeKey>::iterator begin() const {
            return m_edge_kernel->begin();
        }

        typename kernel<Edge, EdgeKey>::iterator end() const {
            return m_edge_kernel->end();
        }
    };

    template <typename face_traits>
    class FaceIterator {
        const kernel<Face, FaceKey>*                  m_face_kernel;
    public:
        FaceIterator(const kernel<Face, FaceKey> *m_face_kernel) : m_face_kernel(m_face_kernel) {
        }

        typename kernel<Face, FaceKey>::iterator begin() const {
            return m_face_kernel->begin();
        }

        typename kernel<Face, FaceKey>::iterator end() const {
            return m_face_kernel->end();
        }

    };

    template <typename tet_traits>
    class TetrahedronIterator {
        const kernel<Tetrahedron, TetrahedronKey>*           m_tetrahedron_kernel;
    public:
        TetrahedronIterator(const kernel<Tetrahedron, TetrahedronKey> *m_tetrahedron_kernel)
                : m_tetrahedron_kernel(m_tetrahedron_kernel) {
        }

        typename kernel<Tetrahedron, TetrahedronKey>::iterator begin() const{
            return m_tetrahedron_kernel->begin();
        }

        typename kernel<Tetrahedron, TetrahedronKey>::iterator end() const {
            return m_tetrahedron_kernel->end();
        }

    };



    class ISMesh
    {


        kernel<Node, NodeKey>* m_node_kernel;
        kernel<Edge, EdgeKey>*                  m_edge_kernel;
        kernel<Face, FaceKey>*                  m_face_kernel;
        kernel<Tetrahedron, TetrahedronKey>*           m_tetrahedron_kernel;
        
    public:
        ISMesh(std::vector<vec3> & points, std::vector<int> & tets, const std::vector<int>& tet_labels);

        ~ISMesh();

        unsigned int get_no_nodes() const;

        unsigned int get_no_edges() const;

        unsigned int get_no_faces() const;

        unsigned int get_no_tets() const;
        
        ///////////////
        // ITERATORS //
        ///////////////
    public:
        NodeIterator<NodeAttributes> nodes() const;

        typename kernel<Node, NodeKey>::iterator nodes_begin();

        typename kernel<Node, NodeKey>::iterator nodes_end();

        EdgeIterator<EdgeAttributes> edges() const;

        typename kernel<Edge, EdgeKey>::iterator edges_begin();

        typename kernel<Edge, EdgeKey>::iterator edges_end();

        FaceIterator<FaceAttributes> faces() const;

        typename kernel<Face, FaceKey>::iterator faces_begin();

        typename kernel<Face, FaceKey>::iterator faces_end();

        TetrahedronIterator<TetAttributes> tetrahedra() const;

        typename kernel<Tetrahedron, TetrahedronKey>::iterator tetrahedra_begin();

        typename kernel<Tetrahedron, TetrahedronKey>::iterator tetrahedra_end();
        
        /////////////////////
        // LABEL FUNCTIONS //
        /////////////////////
    public:

        virtual int get_label(const TetrahedronKey& t);
        
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
        virtual void set_label(const TetrahedronKey& tid, int label);
        
    private:

        struct edge_key {
            int k1, k2;
            edge_key(int i, int j) : k1(i), k2(j) {}
            bool operator<(const edge_key& k) const
            {
                return k1 < k.k1 || (k1 == k.k1 && k2 < k.k2);
            }
        };
        struct face_key
        {
            int k1, k2, k3;
            face_key(int i, int j, int k) : k1(i), k2(j), k3(k){}
            bool operator<(const face_key& k) const
            {
                //        return k1 < k.k1 || (k1 == k.k1 && k2 < k.k2 || (k2 == k.k2 && k3 < k.k3)) ;
                if (k1 < k.k1) return true;
                if (k1 == k.k1 && k2 < k.k2) return true;
                if (k1 == k.k1 && k2 == k.k2 && k3 < k.k3) return true;
                return false;
            }
        };

        template<typename map_type, typename mesh_type>
        inline int create_edge(int i, int j, map_type& edge_map, mesh_type& mesh)
        {
            int a,b;
            if (i <= j) { a = i; b = j; }
            else { a = j; b = i; }
            edge_key key(a, b);
            typename map_type::iterator it = edge_map.find(key);
            if (it == edge_map.end())
            {
                int n = mesh.insert_edge(i, j); //non-sorted
                it = edge_map.insert(std::pair<edge_key,int>(key, n)).first;
            }
            return it->second;
        }
        
        template<typename map_type, typename mesh_type>
        inline int create_face(int i, int j, int k, map_type& face_map, mesh_type& mesh)
        {
            int a[3] = {i, j, k};
            std::sort(a, a+3);
            face_key key(a[0], a[1], a[2]);
            typename map_type::iterator it = face_map.find(key); //lookup in sorted order
            if (it == face_map.end())
            {
                int a = mesh.insert_face(i, j, k); //create in supplied order
                it = face_map.insert(std::pair<face_key,int>(key, a)).first;
            }
            return it->second;
        }

        bool create(const std::vector<vec3>& points, const std::vector<int>& tets);

        void init_flags(const std::vector<int>& tet_labels);

        void update(const SimplexSet<TetrahedronKey>& tids);

        void update_flag(const FaceKey & f);

        void update_flag(const EdgeKey & e);
        
        void connected_component(SimplexSet<TetrahedronKey>& tids, const TetrahedronKey& tid);

        bool crossing(const NodeKey& n);


        void update_flag(const NodeKey & n);
        
        //////////////////////
        // GETTER FUNCTIONS //
        //////////////////////
    public:
        Node & get(const NodeKey& nid);

        Edge & get(const EdgeKey& eid);

        Face & get(const FaceKey& fid);

        Tetrahedron & get(const TetrahedronKey& tid);

        const SimplexSet<NodeKey> & get_nodes(const EdgeKey& eid);

        const SimplexSet<EdgeKey> & get_edges(const NodeKey& nid);

        const SimplexSet<EdgeKey> & get_edges(const FaceKey& fid);

        const SimplexSet<FaceKey> & get_faces(const EdgeKey& eid);

        const SimplexSet<FaceKey> & get_faces(const TetrahedronKey& tid);

        const SimplexSet<TetrahedronKey> & get_tets(const FaceKey& fid);
        
        // Getters for getting the boundary of a boundary etc.

        SimplexSet<NodeKey> get_sorted_nodes(const FaceKey& fid, const TetrahedronKey& tid);

        SimplexSet<NodeKey> get_sorted_nodes(const FaceKey& fid);

        SimplexSet<NodeKey> get_nodes(const FaceKey& fid);

        SimplexSet<NodeKey> get_nodes(const TetrahedronKey& tid);

        SimplexSet<EdgeKey> get_edges(const TetrahedronKey& tid);

        SimplexSet<FaceKey> get_faces(const NodeKey& nid);

        SimplexSet<TetrahedronKey> get_tets(const NodeKey& nid);

        SimplexSet<TetrahedronKey> get_tets(const EdgeKey& eid);
        
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

        NodeKey get_node(const EdgeKey& eid1, const EdgeKey& eid2);

        NodeKey get_node(const EdgeKey& eid, const NodeKey& nid);

        EdgeKey get_edge(const NodeKey& nid1, const NodeKey& nid2);

        EdgeKey get_edge(const FaceKey& fid1, const FaceKey& fid2);

        FaceKey get_face(const NodeKey& nid1, const NodeKey& nid2, const NodeKey& nid3);

        FaceKey get_face(const TetrahedronKey& tid1, const TetrahedronKey& tid2);

        TetrahedronKey get_tet(const TetrahedronKey& tid, const FaceKey& fid);

        vec3 get_pos(const NodeKey& nid);

        std::vector<vec3> get_pos(const SimplexSet<NodeKey>& nids);
        
        //////////////////////
        // EXISTS FUNCTIONS //
        //////////////////////
    public:

        bool exists(const TetrahedronKey& t);

        bool exists(const FaceKey& f);

        bool exists(const EdgeKey& e);

        bool exists(const NodeKey& n);
        
        
        ///////////////////////////
        // ORIENTATION FUNCTIONS //
        ///////////////////////////
        
    public:
        bool is_clockwise_order(const NodeKey& nid, SimplexSet<NodeKey>& nids);

        void orient_cc(const NodeKey& nid, SimplexSet<NodeKey>& nids);

        bool is_inverted(const TetrahedronKey& tid);

        bool is_inverted_destination(const TetrahedronKey& tid);
        
        ////////////////////
        // MESH FUNCTIONS //
        ////////////////////
    
    public:
        NodeKey insert_node(const vec3& p);

        EdgeKey insert_edge(NodeKey node1, NodeKey node2);

        FaceKey insert_face(EdgeKey edge1, EdgeKey edge2, EdgeKey edge3);

        TetrahedronKey insert_tetrahedron(FaceKey face1, FaceKey face2, FaceKey face3, FaceKey face4);
        
    private:
        void remove(const NodeKey& nid);

        void remove(const EdgeKey& eid);

        void remove(const FaceKey& fid);

        void remove(const TetrahedronKey& tid);

        NodeKey merge(const NodeKey& key1, const NodeKey& key2);
        
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

        virtual void update_split(const NodeKey& nid_new, const NodeKey& nid1, const NodeKey& nid2);

        NodeKey split(const EdgeKey& eid, const vec3& pos, const vec3& destination);

        virtual void update_collapse(const NodeKey& nid, const NodeKey& nid_removed, double weight);

        void collapse(const EdgeKey& eid, const NodeKey& nid, double weight = 0.5);

        FaceKey flip_32(const EdgeKey& eid);

        EdgeKey flip_23(const FaceKey& fid);

        void flip(const EdgeKey& eid, const FaceKey& fid1, const FaceKey& fid2);


        void flip_22(const FaceKey& fid1, const FaceKey& fid2);

        void flip_44(const FaceKey& fid1, const FaceKey& fid2);
        
        ///////////////////////
        // UTILITY FUNCTIONS //
        ///////////////////////
    public:

        void garbage_collect();

        virtual void scale(const vec3& s);

        void extract_surface_mesh(std::vector<vec3>& points, std::vector<int>& faces);

        void extract_tet_mesh(std::vector<vec3>& points, std::vector<int>& tets, std::vector<int>& tet_labels);

        void validity_check();
    };
    
}
