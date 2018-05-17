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

#include <functional>
#include <thread>
#include "util.h"
#include "kernel.h"
#include "simplex.h"
#include "simplex_set.h"
#include "is_mesh_iterator.h"
#include "geometry.h"
#include "attribute_vector.h"
#include "cache.h"

namespace is_mesh {

    struct GarbageCollectDeletions {
        std::vector<NodeKey> nodeKeys;
        std::vector<EdgeKey> edgeKeys;
        std::vector<FaceKey> faceKeys;
        std::vector<TetrahedronKey> tetrahedronKeys;
    };

    class Geometry;

    class ISMesh
    {
        kernel<NodeKey,Node> m_node_kernel;
        kernel<EdgeKey,Edge> m_edge_kernel;
        kernel<FaceKey,Face> m_face_kernel;
        kernel<TetrahedronKey,Tetrahedron> m_tetrahedron_kernel;

        std::shared_ptr<Geometry> subdomain;

        std::map<long,std::function<void(const GarbageCollectDeletions&)>> m_gc_listeners;
        std::map<long,std::function<void(const TetrahedronKey& tid, unsigned int oldValue)>> m_set_label_listeners;
        std::map<long,std::function<void(const NodeKey& nid_new, const NodeKey& nid1, const NodeKey& nid2)>> m_split_listeners;
        std::map<long,std::function<void(const NodeKey& nid, const NodeKey& nid_removed, double weight)>> m_collapse_listeners;

        unsigned int m_number_of_threads = std::thread::hardware_concurrency();
        
        //////////////////////
        // Caching
        //////////////////////
#ifdef DSC_CACHE
    public:
        void invalidate_cache(SimplexSet<TetrahedronKey> tets);
        void update_cache_size();
        is_mesh::SimplexSet<is_mesh::TetrahedronKey> * get_tets_cache(is_mesh::NodeKey nid);
        is_mesh::SimplexSet<is_mesh::TetrahedronKey>get_tets_cache(is_mesh::SimplexSet<is_mesh::NodeKey> nids);
        is_mesh::SimplexSet<is_mesh::FaceKey> * get_faces_cache(is_mesh::NodeKey nid);
        is_mesh::SimplexSet<is_mesh::FaceKey> * get_link(is_mesh::NodeKey nk);
        is_mesh::SimplexSet<is_mesh::NodeKey> * get_nodes_cache(is_mesh::NodeKey nk);
        is_mesh::SimplexSet<is_mesh::NodeKey> * get_nodes_cache(is_mesh::TetrahedronKey tk);
        is_mesh::SimplexSet<is_mesh::NodeKey> * get_nodes_cache(is_mesh::FaceKey fk);
        
        is_mesh::SimplexSet<is_mesh::NodeKey> get_nodes_cache(const is_mesh::SimplexSet<TetrahedronKey> & keys);
#endif //DSC_CACHE
    public:
        ISMesh(std::vector<vec3> & points, std::vector<int> & tets, const std::vector<int>& tet_labels);

        // copy of ISMesh is rare. marked explicit to avoid copying the object by mistake (such as pass by value)
        // if mesh has subdomain is is not copied
        explicit ISMesh(const ISMesh& mesh);

        ~ISMesh();

        unsigned int get_no_nodes() const;

        unsigned int get_no_edges() const;

        unsigned int get_no_faces() const;

        unsigned int get_no_tets() const;

        unsigned int get_max_node_key() const;

        unsigned int get_max_edge_key() const;

        unsigned int get_max_face_key() const;

        unsigned int get_max_tet_key() const;

        std::shared_ptr<Geometry> get_subdomain();

        void clear_subdomain();

        void set_subdomain(std::shared_ptr<Geometry> subdomain);

        unsigned int get_number_of_threads() const;

        void set_number_of_threads(unsigned int m_number_of_threads);

        ///////////////
        // ITERATORS //
        ///////////////
    public:
        std::vector<TetrahedronKey> find_par_tet(std::function<bool(Tetrahedron&)> include){ return find_par<TetrahedronKey,Tetrahedron>(include); }
        std::vector<FaceKey> find_par_face(std::function<bool(Face&)> include){ return find_par<FaceKey,Face>(include); }
        std::vector<EdgeKey> find_par_edge(std::function<bool(Edge&)> include){ return find_par<EdgeKey,Edge>(include); }
        std::vector<NodeKey> find_par_node(std::function<bool(Node&)> include){ return find_par<NodeKey,Node>(include); }

        template<typename key_type, typename value_type>
        std::vector<key_type> find_par(std::function<bool(value_type&)> include);

        // Runs the function fn on each node simultaneously on many threads
        // Number of threads used is std::thread::hardware_concurrency()
        template<typename value_type>
        void for_each_par(std::function<void(value_type&,int)> fn);

        // Map each element to a return type using map_fn and then reduce the return value into a single value using reduce_fn
        template<typename simplex_type, typename return_type>
        return_type map_reduce_par(std::function<return_type(simplex_type&)> map_fn, std::function<return_type(return_type,return_type)> reduce_fn, return_type default_value = {});

        // Space partitioned parallel for each
        // Runs the function fn on each node simultaneously on many threads
        // Number of threads used is std::thread::hardware_concurrency()
        // partitionsize must be larger than twice the maximum node size
        // the function is evaluated once for each simplex
        // dimension is on which axis the space is partitioned (x,y or z)
        template<typename value_type>
        void for_each_par_sp(double partitionsize,  int dimension, std::function<void(value_type& node, int threadid)> fn);

        NodeIterator nodes() const;

        EdgeIterator edges() const;

        FaceIterator faces() const;

        TetrahedronIterator tetrahedra() const;

    public:
        void set_label(const TetrahedronKey& tid, int label);

    private:
        template<typename key_type, typename value_type>
        kernel<key_type,value_type>& get_kernel();

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
            auto it = edge_map.find(key);
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
            auto it = face_map.find(key); //lookup in sorted order
            if (it == face_map.end())
            {
                int index = mesh.insert_face(i, j, k); //create in supplied order
                it = face_map.insert(std::pair<face_key,int>(key, index)).first;
            }
            return it->second;
        }

        bool create(const std::vector<vec3>& points, const std::vector<int>& tets);
        /**
        * Perform an initial update of flags for all nodes, edges and faces.
        */
        void init_flags(const std::vector<int>& tet_labels);

        /**
        * Updates the flags (is interface, is boundary, is crossing) of all
        */
        void update_flag();

        /**
        * Updates the flags (is interface, is boundary, is crossing) of simplices in set.
        */
        void update_flag(const SimplexSet<TetrahedronKey>& tids);

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

        const Node & get(const NodeKey& nid) const ;

        const Edge & get(const EdgeKey& eid) const ;

        const Face & get(const FaceKey& fid) const ;

        const Tetrahedron & get(const TetrahedronKey& tid) const ;

        // Getters for getting the boundary of a boundary etc.

        SimplexSet<NodeKey> get_sorted_nodes(const FaceKey& fid, const TetrahedronKey& tid);

        SimplexSet<NodeKey> get_sorted_nodes(const FaceKey& fid);

        // Getters which have a SimplexSet as input
        template<typename key_type>
        SimplexSet<NodeKey> get_nodes(const SimplexSet<key_type>& keys)
        {
            SimplexSet<NodeKey> nids;
            for(auto k : keys)
            {
                nids += get(k).node_keys();
            }
            return nids;
        }

        template<typename key_type>
        SimplexSet<EdgeKey> get_edges(const SimplexSet<key_type>& keys)
        {
            SimplexSet<EdgeKey> eids;
            for(auto k : keys)
            {
                eids += get(k).edge_keys();
            }
            return eids;
        }


        template<typename key_type>
        SimplexSet<FaceKey> get_faces(const SimplexSet<key_type>& keys)
        {
            SimplexSet<FaceKey> fids;
            for(auto k : keys)
            {
                fids += get(k).face_keys();
            }
            return fids;
        }

        template<typename key_type>
        SimplexSet<TetrahedronKey> get_tets(const SimplexSet<key_type>& keys)
        {
            SimplexSet<TetrahedronKey> tids;
            for(auto k : keys)
            {
                tids += get(k).tet_keys();
            }
            return tids;
        }
        
        // Other getter functions
        /**
        *  Returns the node shared by the edges eid1 and eid2.
        */
        NodeKey get_node(const EdgeKey& eid1, const EdgeKey& eid2);
        /**
        *  Returns the node adjacent to edge eid which is not nid.
        */
        NodeKey get_node(const EdgeKey& eid, const NodeKey& nid);
        /**
        *  Returns the edge between the nodes nid1 and nid2.
        */
        EdgeKey get_edge(const NodeKey& nid1, const NodeKey& nid2);
        /**
        *  Returns the edge shared by the faces fid1 and fid2.
        */
        EdgeKey get_edge(const FaceKey& fid1, const FaceKey& fid2);
        /**
        *  Returns the face between the nodes nid1, nid2 and nid3.
        */
        FaceKey get_face(const NodeKey& nid1, const NodeKey& nid2, const NodeKey& nid3);
        /**
        *  Returns the face shared by the tetrahedra tid1 and tid2.
        */
        FaceKey get_face(const TetrahedronKey& tid1, const TetrahedronKey& tid2);
        /**
        *  Returns the tetrahedron which shares the face fid with tid, i.e. the neighbour to tid.
        */
        TetrahedronKey get_tet(const TetrahedronKey& tid, const FaceKey& fid);
        /**
        * Returns the positions of nodes nids.
        */
        std::vector<vec3> get_pos(const SimplexSet<NodeKey>& nids);
        

        
        //////////////////////
        // EXISTS FUNCTIONS //
        //////////////////////
    public:

        bool exists(const TetrahedronKey& t);

        bool exists(const FaceKey& f);

        bool exists(const EdgeKey& e);

        bool exists(const NodeKey& n);


        bool excluded(const TetrahedronKey& t);

        bool excluded(const FaceKey& f);

        bool excluded(const EdgeKey& e);

        bool excluded(const NodeKey& n);


        ///////////////////////////
        // ORIENTATION FUNCTIONS //
        ///////////////////////////
        
    public:
        bool is_clockwise_order(const NodeKey& nid, SimplexSet<NodeKey>& nids);
        /*
         * Orient the nodes in a counter clockwise order seen from the node a.
         */
        void orient_cc(const NodeKey& nid, SimplexSet<NodeKey>& nids);
        /**
        * Returns whether the tetrahedron with ID tid is inverted.
        */
        bool is_inverted(const TetrahedronKey& tid);

        bool is_inverted_destination(const TetrahedronKey& tid);
        
        ////////////////////
        // MESH FUNCTIONS //
        ////////////////////
    
    public:
        /**
        * Inserts a node into the mesh. Trivial.
        */
        NodeKey insert_node(const vec3& p);
        /**
        * Inserts an edge into the mesh. Updates the co-boundary of the boundary nodes with the newly created edge.
        * Leaves the closure of the edge in an uncompressed state.
        */
        EdgeKey insert_edge(NodeKey node1, NodeKey node2);
        /**
        * Inserts a face into the mesh. Updates the co-boundary of the boundary faces with the newly created face.
        * Leaves the closure of the face in an uncompressed state.
        */
        FaceKey insert_face(EdgeKey edge1, EdgeKey edge2, EdgeKey edge3);
        /**
        * Inserts a tetrahedron into the mesh. Updates the co-boundary of the boundary edges with the newly created tetrahedron.
        * Leaves the closure of the tetrahedron in an uncompressed state.
        */
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
        
    public:
        /**
        * Calculates the average position of the nodes in the simplex set nids.
        * If interface is true, the average position is only calculated among the nodes which are interface.
        */
        vec3 get_barycenter(const SimplexSet<NodeKey>& nids, bool interface = false);

        NodeKey split(const EdgeKey& eid, const vec3& pos, const vec3& destination);


        /**
        *  Collapses the edge eid. The node nid must be adjacent to eid before the collapse. The node nid survives, while the other is removed. The weight parameter specifies how the attributes of the old nodes are weighted in the surviving node. For example the position of the surviving node is given by (1.-weight)*get(nid).get_pos() + weight*get(nid_remove).get_pos(). This means that if weight is 0, the surviving node retain its attributes.
        */
        void collapse(const EdgeKey& eid, const NodeKey& nid, double weight = 0.5);

        FaceKey flip_32(const EdgeKey& eid);

        EdgeKey flip_23(const FaceKey& fid);

        void flip(const EdgeKey& eid, const FaceKey& fid1, const FaceKey& fid2);


        void flip_22(const FaceKey& fid1, const FaceKey& fid2);

        void flip_44(const FaceKey& fid1, const FaceKey& fid2);

    private:
        void update_split(const NodeKey& nid_new, const NodeKey& nid1, const NodeKey& nid2);

        void update_collapse(const NodeKey& nid, const NodeKey& nid_removed, double weight);

        ///////////////////////
        // UTILITY FUNCTIONS //
        ///////////////////////
    public:

        double volume_destination(const is_mesh::SimplexSet<is_mesh::NodeKey>& nids);

        double signed_volume_destination(const is_mesh::SimplexSet<is_mesh::NodeKey>& nids);

        void garbage_collect();

        virtual void scale(const vec3& s);

        void extract_surface_mesh(std::vector<vec3>& points, std::vector<int>& faces);

        // extract all tests as tet surfaces
        void extract_surface_mesh_debug(std::vector<vec3>& points, std::vector<int>& faces);

        void extract_tet_mesh(std::vector<vec3>& points, std::vector<int>& tets, std::vector<int>& tet_labels);

        void validity_check(bool skip_boundary_check = false);

        // subscribe for gc events - needed if you are working with persistent attribute_vectors
        // returns listener id
        long add_gc_listener(std::function<void(const GarbageCollectDeletions&)> fn);

        // remove listener by id
        bool remove_gc_listener(long id);

        long add_label_listener(std::function<void(const TetrahedronKey& tid, unsigned int oldValue)> fn);

        bool remove_label_listener(long id);

        long add_split_listener(std::function<void(const NodeKey& nid_new, const NodeKey& nid1, const NodeKey& nid2)> fn);

        bool remove_split_listener(long id);

        long add_collapse_listener(std::function<void(const NodeKey& nid, const NodeKey& nid_removed, double weight)> fn);

        bool remove_collapse_listener(long id);

        friend class Node;
        friend class Edge;
        friend class Face;
        friend class Tetrahedron;
    };

    template<typename key_type, typename value_type>
    inline void run_for_each_par(std::function<void(value_type&,int)> fn, kernel<key_type,value_type> *kernel, int from, int to, int threadid){
        auto begin = kernel->find_valid_iterator(key_type(std::min(from, (int)kernel->capacity())));
        auto end = kernel->find_valid_iterator(key_type(std::min(to, (int)kernel->capacity())));
        for (auto iter = begin;iter!=end;iter++){
            auto& n = *iter;
            fn(n,threadid);
        }
    }

    template<typename value_type>
    inline void ISMesh::for_each_par(std::function<void(value_type&,int)> fn) {
        using KeyType = decltype(std::declval<value_type>().key()); // the type of value().key()
        auto kernel = &get_kernel<KeyType, value_type>();
        kernel->readonly = true;

        int thread_count = m_number_of_threads;
        if (thread_count<=1){
            for (auto &n : *kernel){
                fn(n,0);
            }
        } else {
            std::vector<std::thread *> threads;
            int chunk_size = (int) ceil(kernel->capacity() / (float) (thread_count));

            for (int i = 0; i < thread_count; i++) {
                threads.push_back(new std::thread(run_for_each_par<KeyType, value_type>, fn, kernel, i * chunk_size, (1 + i) * chunk_size, i));
            }

            for (int i = 0; i < thread_count; i++) {
                threads[i]->join();
                delete threads[i];
            }
        }
        kernel->readonly = false;
    }

    template<typename simplex_type, typename return_type>
    inline return_type ISMesh::map_reduce_par(std::function<return_type(simplex_type&)> map_fn, std::function<return_type(return_type,return_type)> reduce_fn, return_type default_value){
        std::vector<return_type> res(m_number_of_threads, default_value);
        for_each_par<simplex_type>([&](simplex_type& v, int index){
            auto value = map_fn(v);
            res[index] = reduce_fn(res[index],value);
        });
        for (int i=1;i<res.size();i++){
            res[i] = reduce_fn(res[i-1], res[i]);
        }
        return res[res.size()-1];
    }

    template<typename key_type, typename value_type>
    inline std::vector<key_type> ISMesh::find_par(std::function<bool(value_type&)> include){
        std::vector<std::vector<key_type>> res_array(m_number_of_threads, {});
        for_each_par<value_type>([&](value_type &value, int threadid){
            if (include(value))
            {
                res_array[threadid].push_back(value.key());
            }
        });
        std::vector<key_type> res;
        for (auto & t:res_array){
            res.insert(res.end(), t.begin(), t.end());
        }
        return res;
    }

    template<typename value_type, typename kernel_type, typename key_type>
    inline void run_for_each_par_sp(std::function<void(value_type&,int)> fn, kernel_type* kernel,  int threadid, int actualthread, AttributeVector<key_type, int> *attributeVector){
        for (int i=0;i< attributeVector->size();i++){
            key_type key(i);
            int run_in_thread = (*attributeVector)[key];
            if (run_in_thread == threadid){
                fn(kernel->get(key), actualthread);
            }
        }
    }

    template<typename value_type>
    inline void ISMesh::for_each_par_sp(double partitionsize, int dimension, std::function<void(value_type&,int)> fn) {
        using KeyType = decltype(std::declval<value_type>().key()); // the type of value().key()
        using KernelType = kernel<KeyType, value_type>;
        auto kernel = &get_kernel<KeyType, value_type>();
        kernel->readonly = true;

        if (m_number_of_threads <=1){
            for (auto &n : *kernel){
                fn(n,0);
            }
        } else {
            AttributeVector<KeyType, int> attributeVector(kernel->capacity(), -1);
            for_each_par<value_type>([&](value_type& n, int t){
                double p = n.get_center()[dimension];
                double concur_partitionsize = partitionsize * m_number_of_threads * 2;
                if (p < 0) {
                    p += concur_partitionsize * (int) (ceil(-p / concur_partitionsize));
                }
                p = fmod(p, concur_partitionsize);
                int run_in_thread = std::min((int) floor(p / partitionsize), ((int) m_number_of_threads * 2) - 1);
                attributeVector[n.key()] = run_in_thread;
            });

            std::vector<std::thread*> threads;
            for (int i=0;i< m_number_of_threads;i++){
                threads.push_back(new std::thread(run_for_each_par_sp<value_type, KernelType, KeyType>, fn, kernel, i*2,i, &attributeVector));
            }
            for (int i=0;i< m_number_of_threads;i++){
                threads[i]->join();
            }
            for (int i=0;i< m_number_of_threads;i++){
                delete threads[i];
                threads[i] = new std::thread(run_for_each_par_sp<value_type, KernelType, KeyType>, fn, kernel, 1+i*2,i, &attributeVector);
            }
            for (int i=0;i< m_number_of_threads;i++){
                threads[i]->join();
                delete threads[i];
            }
        }
        kernel->readonly = false;
    }

    template<>
    inline kernel<NodeKey,Node>& ISMesh::get_kernel(){
        return m_node_kernel;
    }

    template<>
    inline kernel<EdgeKey,Edge>& ISMesh::get_kernel(){
        return m_edge_kernel;
    }

    template<>
    inline kernel<FaceKey,Face>& ISMesh::get_kernel(){
        return m_face_kernel;
    }

    template<>
    inline kernel<TetrahedronKey,Tetrahedron>& ISMesh::get_kernel(){
        return m_tetrahedron_kernel;
    }
}
