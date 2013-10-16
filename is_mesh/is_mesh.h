
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
        
        size_type size_nodes() { return m_node_kernel->size(); }
        size_type size_edges() { return m_edge_kernel->size(); }
        size_type size_faces() { return m_face_kernel->size(); }
        size_type size_tetrahedra() { return m_tetrahedron_kernel->size(); }
        size_type size() { return size_nodes() + size_edges() + size_faces() + size_tetrahedra(); }
        
        
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
