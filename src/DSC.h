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

#include "is_mesh.h"
#include "util.h"
#include "attributes.h"
#include "design_domain.h"

namespace DSC {
    
    template <typename node_att = NodeAttributes, typename edge_att = EdgeAttributes, typename face_att = FaceAttributes, typename tet_att = TetAttributes>
    class DeformableSimplicialComplex : public is_mesh::ISMesh<node_att, edge_att, face_att, tet_att>
    {
        friend class ObjectGenerator;
        typedef is_mesh::ISMesh<node_att, edge_att, face_att, tet_att> ISMesh;
    public:
        
        typedef is_mesh::NodeKey      node_key;
        typedef is_mesh::EdgeKey      edge_key;
        typedef is_mesh::FaceKey      face_key;
        typedef is_mesh::TetrahedronKey       tet_key;
        
    private:
        DesignDomain *design_domain;
        
        // Input parameter
        real AVG_EDGE_LENGTH;
        
        // Should be eliminated
        real FLIP_EDGE_INTERFACE_FLATNESS;
        
    protected:
        // Thresholds on the quality of edges
        real DEG_EDGE_QUALITY;
        real MIN_EDGE_QUALITY;
        
        // Thresholds on the quality of faces.
        real DEG_FACE_QUALITY;
        real MIN_FACE_QUALITY;
        
        // Thresholds on the quality of tetrahedra.
        real DEG_TET_QUALITY;
        real MIN_TET_QUALITY;
        
        // Thresholds on the length of edges.
        real MIN_LENGTH;
        real MAX_LENGTH;
        
        // Thresholds on the area of faces.
        real MIN_AREA;
        real MAX_AREA;
        
        // Thresholds on the volume of tetrahedra.
        real MIN_VOLUME;
        real MAX_VOLUME;
        
        // As close as a node can get to an opposite face before movement is stopped.
        real MIN_DEFORMATION;
        
        //////////////////////////
        // INITIALIZE FUNCTIONS //
        //////////////////////////
        
    public:
        
        /// SimplicialComplex constructor.
        DeformableSimplicialComplex(real _AVG_EDGE_LENGTH, std::vector<real> & points, std::vector<int> & tets, DesignDomain *domain = nullptr):
            ISMesh(points, tets), design_domain(domain)
        {
            AVG_EDGE_LENGTH = _AVG_EDGE_LENGTH;
            MIN_DEFORMATION = 0.0001 * AVG_EDGE_LENGTH;
            
            DEG_EDGE_QUALITY = 0.1;
            MIN_EDGE_QUALITY = 0.5;
            
            DEG_FACE_QUALITY = 1. - cos(2.*M_PI/180.);
            MIN_FACE_QUALITY = 1. - cos(10.*M_PI/180.);
            
            DEG_TET_QUALITY = 0.02;
            MIN_TET_QUALITY = 0.3;
            
            FLIP_EDGE_INTERFACE_FLATNESS = 0.995;
            
            MIN_LENGTH = 0.;
            MAX_LENGTH = 2. * AVG_EDGE_LENGTH;
            
            real area_avg = AVG_EDGE_LENGTH*AVG_EDGE_LENGTH*0.5;
            MIN_AREA = 0.2*area_avg;
            MAX_AREA = 5.*area_avg;
            
            real vol_avg = AVG_EDGE_LENGTH*AVG_EDGE_LENGTH*AVG_EDGE_LENGTH*sqrt(2.)/12.;
            MIN_VOLUME = 0.2*vol_avg;
            MAX_VOLUME = INFINITY;
            
            //        fix_complex();
            //        resize_complex();
        }
        
        DeformableSimplicialComplex()
        {
            delete design_domain;
        }
        
        using ISMesh::get;
        
        using ISMesh::get_nodes;
        using ISMesh::get_edges;
        using ISMesh::get_faces;
        using ISMesh::get_tets;
        
        using ISMesh::get_pos;
    private:
        
        // For debugging!
        void print(const node_key& n)
        {
            std::cout << "Node: " << n << std::endl;
            vec3 p = get_pos(n);
            vec3 d = get(n).get_destination();
            std::cout << "P = " << p[0] << ", " << p[1] << ", " << p[2] << std::endl;
            std::cout << "D = " << d[0] << ", " << d[1] << ", " << d[2]  << std::endl;
            
            std::cout << "\nStar_edges = [";
            for(auto e : ISMesh::get_edges(n))
            {
                auto verts = get_pos(get_nodes(e));
                vec3 p1 = verts[0];
                vec3 p2 = verts[1];
                std::cout << p1[0] << ", " << p1[1] << ", " << p1[2] << "; " << std::endl;
                std::cout << p2[0] << ", " << p2[1] << ", " << p2[2] << "; " << std::endl;
            }
            std::cout << "];" << std::endl;
            
            std::cout << "\nStar_Iedges = [";
            for(auto e : ISMesh::get_edges(n))
            {
                if(get(e).is_interface())
                {
                    auto verts = get_pos(get_nodes(e));
                    vec3 p1 = verts[0];
                    vec3 p2 = verts[1];
                    std::cout << p1[0] << ", " << p1[1] << ", " << p1[2] << "; " << std::endl;
                    std::cout << p2[0] << ", " << p2[1] << ", " << p2[2] << "; " << std::endl;
                }
            }
            std::cout << "];" << std::endl;
            
            std::cout << "\nedges = [";
            auto eids = get_edges(get_tets(n)) - get_edges(n);
            for(auto e : eids)
            {
                auto verts = get_pos(get_nodes(e));
                vec3 p1 = verts[0];
                vec3 p2 = verts[1];
                std::cout << p1[0] << ", " << p1[1] << ", " << p1[2] << "; " << std::endl;
                std::cout << p2[0] << ", " << p2[1] << ", " << p2[2] << "; " << std::endl;
            }
            std::cout << "];" << std::endl;
            
            std::cout << "\nIedges = [";
            for(auto e : eids)
            {
                if(get(e).is_interface())
                {
                    auto verts = ISMesh::get_pos(ISMesh::get_nodes(e));
                    vec3 p1 = verts[0];
                    vec3 p2 = verts[1];
                    std::cout << p1[0] << ", " << p1[1] << ", " << p1[2] << "; " << std::endl;
                    std::cout << p2[0] << ", " << p2[1] << ", " << p2[2] << "; " << std::endl;
                }
            }
            std::cout << "];" << std::endl;
        }
        
        /////////////////////////
        // ATTRIBUTE FUNCTIONS //
        /////////////////////////
    public:
        virtual void update_attributes()
        {
            
        }
        
    protected:
        
        virtual bool is_unsafe_editable(const node_key& nid)
        {
            return ISMesh::exists(nid) && !get(nid).is_boundary();
        }
        
        virtual bool is_unsafe_editable(const edge_key& eid)
        {
            return ISMesh::exists(eid) && !get(eid).is_boundary();
        }
        
        virtual bool is_unsafe_editable(const face_key& fid)
        {
            return ISMesh::exists(fid) && !get(fid).is_boundary();
        }
        
        virtual bool is_unsafe_editable(const tet_key& tid)
        {
            return ISMesh::exists(tid);
        }
        
        virtual bool is_safe_editable(const node_key& nid)
        {
            return is_unsafe_editable(nid) && !get(nid).is_interface();
        }
        
        virtual bool is_safe_editable(const edge_key& eid)
        {
            return is_unsafe_editable(eid) && !get(eid).is_interface();
        }
        
        virtual bool is_safe_editable(const face_key& fid)
        {
            return is_unsafe_editable(fid) && !get(fid).is_interface();
        }
        
        virtual bool is_safe_editable(const tet_key& tid)
        {
            return is_unsafe_editable(tid);
        }
        
    public:
        virtual bool is_movable(const node_key& nid)
        {
            return is_unsafe_editable(nid) && get(nid).is_interface() && !get(nid).is_crossing();
        }
        
        int get_label(const tet_key& t)
        {
            return ISMesh::get_label(t);
        }
        
    protected:
        /**
         * Sets the position of node n.
         */
        void set_pos(const node_key& nid, const vec3& p)
        {
            ISMesh::get(nid).set_pos(p);
            if(!is_movable(nid))
            {
                ISMesh::get(nid).set_destination(p);
            }
        }
        
    public:
        /// Returns the destination of the node with id nid.
        vec3 get_dest(const node_key& nid)
        {
            if(is_movable(nid))
            {
                return get(nid).get_destination();
            }
            return get(nid).get_pos();
        }
        
        /// Returns the destinations of the nodes in the simplex set.
        std::vector<vec3> get_dest(const is_mesh::SimplexSet<node_key>& nids)
        {
            std::vector<vec3> verts;
            for (auto n : nids) {
                verts.push_back(get_dest(n));
            }
            return verts;
        }
        
        /**
         * Sets the destination where the node n is moved to when deform() is called.
         */
        void set_destination(const node_key& nid, vec3 dest)
        {
            if(is_movable(nid))
            {
                if(design_domain)
                {
                    vec3 p = get_pos(nid);
                    vec3 vec = dest - p;
                    design_domain->clamp_vector(p, vec);
                    get(nid).set_destination(p + vec);
                }
                else {
                    get(nid).set_destination(dest);
                }
            }
        }
        
        /////////////
        // GETTERS //
        /////////////
    public:
        
        vec3 get_center() const
        {
            return vec3(0.);
        }
        
        real get_min_tet_quality() const
        {
            return MIN_TET_QUALITY;
        }
        
        real get_deg_tet_quality() const
        {
            return DEG_TET_QUALITY;
        }
        
        real get_deg_face_quality() const
        {
            return DEG_FACE_QUALITY;
        }
        
        real get_min_face_quality() const
        {
            return MIN_FACE_QUALITY;
        }
        
        real get_avg_edge_length() const
        {
            return AVG_EDGE_LENGTH;
        }
        
        const DesignDomain* get_design_domain() const
        {
            return design_domain;
        }
        
        ////////////////////////
        // FIX MESH FUNCTIONS //
        ////////////////////////
    private:
        
        //////////////////////////////
        // TOPOLOGICAL EDGE REMOVAL //
        //////////////////////////////
        
        /**
         * Build a table K for the dynamic programming method by Klincsek (see Shewchuk "Two Discrete Optimization Algorithms
         * for the Topological Improvement of Tetrahedral Meshes" article for details).
         */
        real build_table(const edge_key& e, const is_mesh::SimplexSet<node_key>& polygon, std::vector<std::vector<int>>& K)
        {
            is_mesh::SimplexSet<node_key> nids = get_nodes(e);
            
            const int m = (int) polygon.size();
            
            auto Q = std::vector<std::vector<real>>(m-1, std::vector<real>(m, 0.));
            K = std::vector<std::vector<int>>(m-1, std::vector<int>(m, 0) );
            
            for(int i = 0; i < m-1; i++)
            {
                Q[i][i+1] = INFINITY;
            }
            
            for (int i = m-3; i >= 0; i--)
            {
                for (int j = i+2; j < m; j++)
                {
                    for (int k = i+1; k < j; k++)
                    {
                        real q2 = Util::quality<real>(get_pos(polygon[i]), get_pos(polygon[k]), get_pos(nids[0]), get_pos(polygon[j]));
                        real q1 = Util::quality<real>(get_pos(polygon[k]), get_pos(polygon[i]), get_pos(nids[1]), get_pos(polygon[j]));
                        real q = Util::min(q1, q2);
                        if (k < j-1)
                        {
                            q = Util::min(q, Q[k][j]);
                        }
                        if (k > i+1)
                        {
                            q = Util::min(q, Q[i][k]);
                        }
                        
                        if (k == i+1 || q > Q[i][j])
                        {
                            Q[i][j] = q;
                            K[i][j] = k;
                        }
                    }
                }
            }
            
            return Q[0][m-1];
        }
        
        
        node_key get_next(const node_key& nid, is_mesh::SimplexSet<edge_key>& eids)
        {
            for (auto e : eids)
            {
                auto n = get_nodes(e) - nid;
                if(n.size() == 1)
                {
                    eids -= e;
                    return n.front();
                }
            }
            return is_mesh::NodeKey();
        }
        
        is_mesh::SimplexSet<node_key> get_polygon(is_mesh::SimplexSet<edge_key>& eids)
        {
            is_mesh::SimplexSet<node_key> polygon = {get_nodes(eids[0]).front()};
            node_key next_nid;
            do {
                next_nid = get_next(polygon.back(), eids);
                if(next_nid.is_valid() && !polygon.contains(next_nid))
                {
                    polygon.push_back(next_nid);
                }
            } while(next_nid.is_valid());
            
            do {
                next_nid = get_next(polygon.front(), eids);
                if(next_nid.is_valid() && !polygon.contains(next_nid))
                {
                    polygon.push_front(next_nid);
                }
            } while(next_nid.is_valid());
            return polygon;
        }
        
        std::vector<is_mesh::SimplexSet<node_key>> get_polygons(const edge_key& eid)
        {
            std::vector<is_mesh::SimplexSet<tet_key>> tid_groups;
            for (auto t : get_tets(eid))
            {
                bool found = false;
                for(auto& tids : tid_groups)
                {
                    if(get_label(t) == get_label(tids.front()))
                    {
                        tids += t;
                        found = true;
                        break;
                    }
                }
                if(!found)
                {
                    tid_groups.push_back({t});
                }
            }
            
            std::vector<is_mesh::SimplexSet<node_key>> polygons;
            is_mesh::SimplexSet<edge_key> m_eids = get_edges(get_faces(eid));
            for(auto& tids : tid_groups)
            {
                is_mesh::SimplexSet<edge_key> eids = get_edges(tids) - m_eids;
                is_mesh::SimplexSet<node_key> polygon = get_polygon(eids);
                check_consistency(get_nodes(eid), polygon);
                polygons.push_back(polygon);
            }
            
            struct {
                bool operator()(const is_mesh::SimplexSet<node_key>& a, const is_mesh::SimplexSet<node_key>& b)
                {
                    return a.size() > b.size();
                }
            } compare;
            std::sort(polygons.begin(), polygons.end(), compare);
            
            return polygons;
        }
        
        void flip_23_recursively(const is_mesh::SimplexSet<node_key>& polygon, const node_key& n1, const node_key& n2, std::vector<std::vector<int>>& K, int i, int j)
        {
            if(j >= i+2)
            {
                int k = K[i][j];
                flip_23_recursively(polygon, n1, n2, K, i, k);
                flip_23_recursively(polygon, n1, n2, K, k, j);
                ISMesh::flip_23(ISMesh::get_face(n1, n2, polygon[k]));
            }
        }
        
        void topological_edge_removal(const is_mesh::SimplexSet<node_key>& polygon, const node_key& n1, const node_key& n2, std::vector<std::vector<int>>& K)
        {
            const int m = static_cast<int>(polygon.size());
            int k = K[0][m-1];
            flip_23_recursively(polygon, n1, n2, K, 0, k);
            flip_23_recursively(polygon, n1, n2, K, k, m-1);
            ISMesh::flip_32(ISMesh::get_edge(n1, n2));
        }
        
        /**
         * Attempt to remove edge e by mesh reconnection using the dynamic programming method by Klincsek (see Shewchuk "Two Discrete Optimization Algorithms
         * for the Topological Improvement of Tetrahedral Meshes" article for details).
         */
        bool topological_edge_removal(const edge_key& eid)
        {
            std::vector<is_mesh::SimplexSet<node_key>> polygon = get_polygons(eid);
            assert(polygon.size() == 1 && polygon.front().size() > 2);
            
            std::vector<std::vector<int>> K;
            real q_new = build_table(eid, polygon.front(), K);
            
            if (q_new > min_quality(get_tets(eid)))
            {
                const is_mesh::SimplexSet<node_key>& nodes = get_nodes(eid);
                topological_edge_removal(polygon.front(), nodes[0], nodes[1], K);
                return true;
            }
            return false;
        }
        
        void topological_boundary_edge_removal(const is_mesh::SimplexSet<node_key>& polygon1, const is_mesh::SimplexSet<node_key>& polygon2, const edge_key& eid, std::vector<std::vector<int>>& K1, std::vector<std::vector<int>>& K2)
        {
            auto nids = get_nodes(eid);
            const int m1 = static_cast<int>(polygon1.size());
            const int m2 = static_cast<int>(polygon2.size());
            int k = K1[0][m1-1];
            flip_23_recursively(polygon1, nids[0], nids[1], K1, 0, k);
            flip_23_recursively(polygon1, nids[0], nids[1], K1, k, m1-1);
            
            if(m2 <= 2) {
                // Find the faces to flip about.
                face_key f1 = ISMesh::get_face(nids[0], nids[1], polygon1.front());
                face_key f2 = ISMesh::get_face(nids[0], nids[1], polygon1.back());
                assert(get(f1).is_boundary() && get(f2).is_boundary());
                
                if(precond_flip_edge(ISMesh::get_edge(f1, f2), f1, f2))
                {
                    ISMesh::flip_22(f1, f2);
                }
            }
            else {
                k = K2[0][m2-1];
                flip_23_recursively(polygon2, nids[0], nids[1], K2, 0, k);
                flip_23_recursively(polygon2, nids[0], nids[1], K2, k, m2-1);
                
                // Find the faces to flip about.
                face_key f1 = ISMesh::get_face(nids[0], nids[1], polygon1.front());
                face_key f2 = ISMesh::get_face(nids[0], nids[1], polygon1.back());
                
                if(precond_flip_edge(ISMesh::get_edge(f1, f2), f1, f2))
                {
                    ISMesh::flip_44(f1, f2);
                }
            }
        }
        
        bool topological_boundary_edge_removal(const edge_key& eid)
        {
            std::vector<is_mesh::SimplexSet<node_key>> polygons = get_polygons(eid);
            
            if(polygons.size() > 2 || polygons[0].size() <= 2)
            {
                return false;
            }
            
            std::vector<std::vector<int>> K1, K2;
            real q_new = build_table(eid, polygons[0], K1);
            
            if(polygons.size() == 2 && polygons[1].size() > 2)
            {
                q_new = Util::min(q_new, build_table(eid, polygons[1], K2));
            }
            else
            {
                polygons.push_back({polygons[0].front(), polygons[0].back()});
            }
            
            if (q_new > min_quality(get_tets(eid)))
            {
                topological_boundary_edge_removal(polygons[0], polygons[1], eid, K1, K2);
                return true;
            }
            return false;
        }
        
        /**
         * Improve tetrahedra quality by the topological operation (re-connection) edge removal. It do so only for tetrahedra of quality lower than MIN_TET_QUALITY.
         * For details see "Two Discrete Optimization Algorithms for the Topological Improvement of Tetrahedral Meshes" by Shewchuk.
         */
        void topological_edge_removal()
        {
            std::vector<tet_key> tets;
            for (auto tit = ISMesh::tetrahedra_begin(); tit != ISMesh::tetrahedra_end(); tit++)
            {
                if (quality(tit.key()) < MIN_TET_QUALITY)
                {
                    tets.push_back(tit.key());
                }
            }
            
            // Attempt to remove each edge of each tetrahedron in tets. Accept if it increases the minimum quality locally.
            int i = 0, j = 0, k = 0;
            for (auto &t : tets)
            {
                if (ISMesh::exists(t) && quality(t) < MIN_TET_QUALITY)
                {
                    for (auto e : ISMesh::get_edges(t))
                    {
                        if (ISMesh::exists(e))
                        {
                            if(is_safe_editable(e))
                            {
                                if(topological_edge_removal(e))
                                {
                                    i++;
                                }
                            }
                            else if((get(e).is_interface() || get(e).is_boundary()) && is_flippable(e))
                            {
                                if(topological_boundary_edge_removal(e))
                                {
                                    k++;
                                }
                            }
                            j++;
                        }
                    }
                }
            }
            std::cout << "Topological edge removals: " << i + k << "/" << j << " (" << k << " at interface)" << std::endl;
            ISMesh::garbage_collect();
        }
        
        //////////////////////////////
        // TOPOLOGICAL FACE REMOVAL //
        //////////////////////////////
        
        is_mesh::SimplexSet<edge_key> test_neighbour(const face_key& f, const node_key& a, const node_key& b, const node_key& u, const node_key& w, real& q_old, real& q_new)
        {
            edge_key e = ISMesh::get_edge(u,w);
            is_mesh::SimplexSet<face_key> g_set = get_faces(e) - get_faces(get_tets(f));
            real q = Util::quality<real>(get_pos(a), get_pos(b), get_pos(w), get_pos(u));
            
            if(g_set.size() == 1 && is_safe_editable(e))
            {
                face_key g = g_set.front();
                node_key v = (get_nodes(g) - get_nodes(e)).front();
                real V_uv = Util::signed_volume<real>(get_pos(a), get_pos(b), get_pos(v), get_pos(u));
                real V_vw = Util::signed_volume<real>(get_pos(a), get_pos(b), get_pos(w), get_pos(v));
                real V_wu = Util::signed_volume<real>(get_pos(a), get_pos(b), get_pos(u), get_pos(w));
                
                if((V_uv > 0. && V_vw > 0.) || (V_vw > 0. && V_wu > 0.) || (V_wu > 0. && V_uv > 0.))
                {
                    q_old = Util::min(Util::quality<real>(get_pos(a), get_pos(u), get_pos(w), get_pos(v)),
                                     Util::quality<real>(get_pos(u), get_pos(v), get_pos(b), get_pos(w)));
                    
                    real q_uv_old, q_uv_new, q_vw_old, q_vw_new;
                    is_mesh::SimplexSet<edge_key> uv_edges = test_neighbour(g, a, b, u, v, q_uv_old, q_uv_new);
                    is_mesh::SimplexSet<edge_key> vw_edges = test_neighbour(g, a, b, v, w, q_vw_old, q_vw_new);
                    
                    q_old = Util::min(Util::min(q_old, q_uv_old), q_vw_old);
                    q_new = Util::min(q_uv_new, q_vw_new);
                    
                    if(q_new > q_old || q_new > q)
                    {
                        is_mesh::SimplexSet<edge_key> edges = {ISMesh::get_edge(f, g)};
                        edges += uv_edges;
                        edges += vw_edges;
                        return edges;
                    }
                }
            }
            q_old = INFINITY;
            q_new = q;
            return {};
        }
        
        /**
         * Attempt to remove the faces sandwiched between the apices of f using multi-face removal. The face f is used as a starting point.
         */
        bool topological_face_removal(const face_key& f)
        {
            is_mesh::SimplexSet<node_key> nids = get_nodes(f);
            is_mesh::SimplexSet<node_key> apices = get_nodes(get_tets(f)) - nids;
            ISMesh::orient_cc(apices[0], nids);
            
            real q_01_new, q_01_old, q_12_new, q_12_old, q_20_new, q_20_old;
            is_mesh::SimplexSet<edge_key> e01 = test_neighbour(f, apices[0], apices[1], nids[0], nids[1], q_01_old, q_01_new);
            is_mesh::SimplexSet<edge_key> e12 = test_neighbour(f, apices[0], apices[1], nids[1], nids[2], q_12_old, q_12_new);
            is_mesh::SimplexSet<edge_key> e20 = test_neighbour(f, apices[0], apices[1], nids[2], nids[0], q_20_old, q_20_new);
            
            real q_old = Util::min(Util::min(Util::min(min_quality(get_tets(f)), q_01_old), q_12_old), q_20_old);
            real q_new = Util::min(Util::min(q_01_new, q_12_new), q_20_new);
            
            if(q_new > q_old)
            {
                ISMesh::flip_23(f);
                for(auto &e : e01)
                {
                    ISMesh::flip_32(e);
                }
                for(auto &e : e12)
                {
                    ISMesh::flip_32(e);
                }
                for(auto &e : e20)
                {
                    ISMesh::flip_32(e);
                }
                return true;
            }
            return false;
        }
        
        /**
         * Attempt to remove the faces sandwiched between the nodes apex1 and apex2 using multi-face removal.
         * The face which intersects with the line segment |apex1 apex2| is used as a starting point.
         */
        bool topological_face_removal(const node_key& apex1, const node_key& apex2)
        {
            is_mesh::SimplexSet<face_key> fids = get_faces(get_tets(apex1)) & get_faces(get_tets(apex2));
            vec3 p = get_pos(apex1);
            vec3 ray = get_pos(apex2) - p;
            for(auto f : fids)
            {
                if(is_safe_editable(f))
                {
                    auto nids = get_nodes(f);
                    this->orient_cc(apex2, nids);
                    
                    real t = Util::intersection_ray_triangle<real>(p, ray, get_pos(nids[0]), get_pos(nids[1]), get_pos(nids[2]));
                    if(0. < t && t < 1.)
                    {
                        if(topological_face_removal(f))
                        {
                            return true;
                        }
                    }
                }
            }
            return false;
        }
        
        /**
         * Improve tetrahedra quality by the topological operation (re-connection) multi-face removal. It do so only for tetrahedra of quality lower than MIN_TET_QUALITY.
         * For details see "Two Discrete Optimization Algorithms for the Topological Improvement of Tetrahedral Meshes" by Shewchuk.
         */
        void topological_face_removal()
        {
            std::vector<tet_key> tets;
            for (auto tit = ISMesh::tetrahedra_begin(); tit != ISMesh::tetrahedra_end(); tit++)
            {
                if (quality(tit.key()) < MIN_TET_QUALITY)
                {
                    tets.push_back(tit.key());
                }
            }
            
            // Attempt to remove each face of each remaining tetrahedron in tets using multi-face removal.
            // Accept if it increases the minimum quality locally.
            int i = 0, j = 0;
            for (auto &t : tets)
            {
                if (ISMesh::exists(t) && quality(t) < MIN_TET_QUALITY)
                {
                    for (auto f : ISMesh::get_faces(t))
                    {
                        if (is_safe_editable(f))
                        {
                            auto apices = ISMesh::get_nodes(ISMesh::get_tets(f)) - ISMesh::get_nodes(f);
                            if(topological_face_removal(apices[0], apices[1]))
                            {
                                i++;
                            }
                            j++;
                        }
                    }
                }
            }
            std::cout << "Topological face removals: " << i << "/" << j << std::endl;
            
            ISMesh::garbage_collect();
        }
        
        ////////////////
        // THICKENING //
        ////////////////
        
        /**
         * Splits all interface edges with a volume greater than MAX_EDGE_LENGTH by inserting a vertex.
         */
        void thickening_interface()
        {
            if(MAX_LENGTH == INFINITY)
            {
                return;
            }
            
            std::vector<edge_key> edges;
            for (auto eit = ISMesh::edges_begin(); eit != ISMesh::edges_end(); eit++)
            {
                if (is_unsafe_editable(eit.key()) && eit->is_interface() && length(eit.key()) > MAX_LENGTH)
                {
                    edges.push_back(eit.key());
                }
            }
            int i = 0, j = 0;
            for(auto &e : edges)
            {
                if (is_unsafe_editable(e) && get(e).is_interface() && length(e) > MAX_LENGTH)
                {
                    node_key nid = split(e);
                    if(nid.is_valid())
                    {
                        i++;
                    }
                    j++;
                }
            }
            std::cout << "Thickening interface splits: " << i << "/" << j << std::endl;
        }
        
        /**
         * Splits all tetrahedra with a volume greater than MAX_TET_VOLUME by inserting a vertex.
         */
        void thickening()
        {
            if(MAX_VOLUME == INFINITY)
            {
                return;
            }
            
            std::vector<tet_key> tetrahedra;
            for (auto tit = ISMesh::tetrahedra_begin(); tit != ISMesh::tetrahedra_end(); tit++)
            {
                if (volume(tit.key()) > MAX_VOLUME)
                {
                    tetrahedra.push_back(tit.key());
                }
            }
            int i = 0, j = 0;
            for(auto &t : tetrahedra)
            {
                if (ISMesh::exists(t) && volume(t) > MAX_VOLUME)
                {
                    if(split(t).is_valid())
                    {
                        i++;
                    }
                    j++;
                }
            }
            std::cout << "Thickening splits: " << i << "/" << j << std::endl;
        }
        
        //////////////
        // THINNING //
        //////////////
        
        void thinning_interface()
        {
            if(MIN_LENGTH <= 0.)
            {
                return;
            }
            
            std::vector<edge_key> edges;
            for (auto eit = ISMesh::edges_begin(); eit != ISMesh::edges_end(); eit++)
            {
                if (is_unsafe_editable(eit.key()) && eit->is_interface() && length(eit.key()) < MIN_LENGTH)
                {
                    edges.push_back(eit.key());
                }
            }
            int i = 0, j = 0;
            for(auto &e : edges)
            {
                if (is_unsafe_editable(e) && get(e).is_interface() && length(e) < MIN_LENGTH)
                {
                    if(collapse(e, false))
                    {
                        i++;
                    }
                    j++;
                }
            }
            std::cout << "Thinning interface splits: " << i << "/" << j << std::endl;
        }
        
        /**
         * Collapses all edges which is shorter than MIN_EDGE_LENGTH. However, the collapse is only performed if it does not alter the interface and it improves the minimum quality of the tetrahedra in the star of the edge.
         */
        void thinning()
        {
            if(MIN_VOLUME <= 0.)
            {
                return;
            }
            
            std::vector<tet_key> tetrahedra;
            for (auto tit = ISMesh::tetrahedra_begin(); tit != ISMesh::tetrahedra_end(); tit++)
            {
                if (volume(tit.key()) < MIN_VOLUME)
                {
                    tetrahedra.push_back(tit.key());
                }
            }
            int i = 0, j = 0;
            for(auto &t : tetrahedra)
            {
                if (ISMesh::exists(t) && volume(t) < MIN_VOLUME)
                {
                    if(collapse(t))
                    {
                        i++;
                    }
                    j++;
                }
            }
            std::cout << "Thinning collapses: " << i << "/" << j << std::endl;
        }
        
        /////////////////////////
        // REMOVE DEGENERACIES //
        /////////////////////////
        /**
         * Attempt to remove edges with lower quality than DEG_EDGE_QUALITY by collapsing them.
         */
        void remove_degenerate_edges()
        {
            std::list<edge_key> edges;
            for (auto eit = ISMesh::edges_begin(); eit != ISMesh::edges_end(); eit++)
            {
                if (quality(eit.key()) < DEG_EDGE_QUALITY)
                {
                    edges.push_back(eit.key());
                }
            }
            int i = 0, j = 0;
            for(auto e : edges)
            {
                if(ISMesh::exists(e) && quality(e) < DEG_EDGE_QUALITY)
                {
                    if(collapse(e, false))
                    {
                        i++;
                    }
                    j++;
                }
            }
            std::cout << "Removed " << i <<"/"<< j << " degenerate edges" << std::endl;
            ISMesh::garbage_collect();
        }
        
        void remove_degenerate_faces()
        {
            std::list<face_key> faces;
            
            for (auto fit = ISMesh::faces_begin(); fit != ISMesh::faces_end(); fit++)
            {
                if(quality(fit.key()) < DEG_FACE_QUALITY)
                {
                    faces.push_back(fit.key());
                }
            }
            
            int i = 0, j = 0;
            for (auto &f : faces)
            {
                if (ISMesh::exists(f) && quality(f) < DEG_FACE_QUALITY)
                {
                    if(collapse(f, false))
                    {
                        i++;
                    }
                    else {
                        split(longest_edge(get_edges(f)));
                    }
                    j++;
                }
            }
            std::cout << "Removed " << i <<"/"<< j << " degenerate faces" << std::endl;
            ISMesh::garbage_collect();
        }
        
        void remove_degenerate_tets()
        {
            std::vector<tet_key> tets;
            
            for (auto tit = ISMesh::tetrahedra_begin(); tit != ISMesh::tetrahedra_end(); tit++)
            {
                if (quality(tit.key()) < DEG_TET_QUALITY)
                {
                    tets.push_back(tit.key());
                }
            }
            int i = 0, j = 0;
            for (auto &t : tets)
            {
                if (ISMesh::exists(t) && quality(t) < DEG_TET_QUALITY)
                {
                    if(collapse(t, false))
                    {
                        i++;
                    }
                    else {
                        split(longest_edge(get_edges(t)));
                    }
                    j++;
                }
            }
            std::cout << "Removed " << i <<"/"<< j << " degenerate tets" << std::endl;
            ISMesh::garbage_collect();
        }
        
        //////////////////////////////////
        // REMOVE LOW QUALITY SIMPLICES //
        //////////////////////////////////
        
        /**
         * Attempt to remove edges with worse quality than MIN_EDGE_QUALITY by safely collapsing them.
         */
        void remove_edges()
        {
            std::list<edge_key> edges;
            for (auto eit = ISMesh::edges_begin(); eit != ISMesh::edges_end(); eit++)
            {
                if (quality(eit.key()) < MIN_EDGE_QUALITY)
                {
                    edges.push_back(eit.key());
                }
            }
            int i = 0, j = 0;
            for(auto e : edges)
            {
                if(ISMesh::exists(e) && quality(e) < MIN_EDGE_QUALITY)
                {
                    if(collapse(e))
                    {
                        i++;
                    }
                    j++;
                }
            }
            std::cout << "Removed " << i <<"/"<< j << " low quality edges" << std::endl;
            ISMesh::garbage_collect();
        }
        
        /**
         * Attempt to remove the cap f by splitting the longest edge and collapsing it with cap's apex.
         */
        bool remove_cap(const face_key& fid)
        {
            // Find longest edge
            edge_key eid = longest_edge(ISMesh::get_edges(fid));
            
            // Find apex
            node_key apex = (ISMesh::get_nodes(fid) - ISMesh::get_nodes(eid)).front();
            // Find the projected position of the apex
            auto verts = get_pos(get_nodes(eid));
            vec3 p = Util::project(get_pos(apex), verts[0], verts[1]);
            
            // Split longest edge
            node_key n = ISMesh::split(eid, p, p);
            
            // Collapse new edge
            edge_key e_rem = ISMesh::get_edge(apex, n);
            return collapse(e_rem);
        }
        
        /**
         * Attempt to remove the cap f by splitting the longest edge and collapsing it with cap's apex.
         */
        bool remove_needle(const face_key& fid)
        {
            // Find shortest edge
            edge_key e = shortest_edge(ISMesh::get_edges(fid));
            
            // Remove edge
            return collapse(e);
        }
        
        /**
         * Attempt to remove the face f by first determining whether it's a cap or a needle.
         */
        bool remove_face(const face_key& f)
        {
            if(max_angle(f) > 0.9*M_PI)
            {
                return remove_cap(f);
            }
            return remove_needle(f);
        }
        
        /**
         * Attempts to remove degenerate faces (faces with a minimum angle smaller than MIN_ANGLE).
         */
        void remove_faces()
        {
            std::list<face_key> faces;
            
            for (auto fit = ISMesh::faces_begin(); fit != ISMesh::faces_end(); fit++)
            {
                if(quality(fit.key()) < MIN_FACE_QUALITY)
                {
                    faces.push_back(fit.key());
                }
            }
            
            int i = 0, j = 0;
            for (auto &f : faces)
            {
                if (ISMesh::exists(f) && quality(f) < MIN_FACE_QUALITY)
                {
                    if(remove_face(f))
                    {
                        i++;
                    }
                    j++;
                }
            }
            std::cout << "Removed " << i <<"/"<< j << " low quality faces" << std::endl;
            ISMesh::garbage_collect();
        }
        
        /**
         * Remove a degenerate tetrahedron of a type "sliver" by splitting the two longest edges
         * and collapsing the newly created vertices together. Return true if successful.
         */
        bool remove_sliver(const tet_key& tid)
        {
            is_mesh::SimplexSet<edge_key> eids = get_edges(tid);
            edge_key e1 = longest_edge(eids);
            eids -= e1;
            edge_key e2 = longest_edge(eids);
            
            node_key n1 = split(e1);
            node_key n2 = split(e2);
            
            edge_key e = ISMesh::get_edge(n1, n2);
            return collapse(e);
        }
        
        /**
         * Remove a degenerate tetrahedron of a type "cap" by splitting the face opposite cap's apex and collapsing cap's apex with the newly created vertex.
         * Return true if successful.
         */
        bool remove_cap(const tet_key& tid)
        {
            // Find the largest face
            face_key fid = largest_face(ISMesh::get_faces(tid));
            
            // Find the apex
            node_key apex = (get_nodes(tid) - get_nodes(fid)).front();
            
            // Project the apex
            auto verts = ISMesh::get_pos(get_nodes(fid));
            vec3 p = Util::project(ISMesh::get_pos(apex), verts[0], verts[1], verts[2]);
            
            // Split the face
            node_key n = ISMesh::split(fid, p, p);
            
            // Collapse edge
            edge_key e = ISMesh::get_edge(n, apex);
            return collapse(e);
        }
        
        /**
         * Remove a degenerate tetrahedron of a type "wedge" or "needle" by collapsing the shortest edge.
         * Return true if successful.
         */
        bool remove_wedge(const tet_key& tid)
        {
            is_mesh::SimplexSet<edge_key> eids = ISMesh::get_edges(tid);
            while(eids.size() > 2)
            {
                edge_key e = shortest_edge(eids);
                if(collapse(e))
                {
                    return true;
                }
                eids -= e;
            }
            return false;
            
            //        simplex_set cl_t;
            //        ISMesh::closure(t, cl_t);
            //        edge_key e1 = longest_edge(cl_t);
            //        cl_t.erase(e1);
            //        edge_key e2 = longest_edge(cl_t);
            //
            //        node_key n1 = split(e1);
            //        node_key n2 = split(e2);
            //
            //        edge_key e = ISMesh::get_edge(n1, n2);
            //        return collapse(e);
        }
        
        /**
         * Remove a tetrahedron of a type "needle" by splitting the tetrahedron.
         * Return true if successful.
         */
        bool remove_needle(const tet_key& tid)
        {
            is_mesh::SimplexSet<edge_key> eids = ISMesh::get_edges(tid);
            while(eids.size() > 1)
            {
                edge_key e = shortest_edge(eids);
                if(collapse(e))
                {
                    return true;
                }
                eids -= e;
            }
            return false;
            //        split(t);
            //        return true;
        }
        
        /**
         * Destroy degenerate (nearly flat) tetrahedron t by splits and collapses.
         * This function detects what type of degeneracy tetrahedron t is (sliver, cap, needle or wedge)
         * and selects appropriate degeneracy removal routine.
         */
        bool remove_tet(const tet_key& tid)
        {
            // Find the largest face
            is_mesh::SimplexSet<face_key> fids = get_faces(tid);
            face_key fid = largest_face(fids);
            is_mesh::SimplexSet<node_key> nids = get_nodes(fid);
            
            // Find the apex
            node_key apex = (get_nodes(tid) - nids).front();
            
            // Project the apex
            auto verts = ISMesh::get_pos(nids);
            vec3 proj_apex = Util::project(ISMesh::get_pos(apex), verts[0], verts[1], verts[2]);
            
            // Find barycentric coordinates
            std::vector<real> barycentric_coords = Util::barycentric_coords<real>(proj_apex, verts[0], verts[1], verts[2]);
            
            if(barycentric_coords[0] > 0.2 && barycentric_coords[1] > 0.2 && barycentric_coords[2] > 0.2) // The tetrahedron is a cap
            {
                return remove_cap(tid);
            }
            else if(barycentric_coords[0] < -0.2 || barycentric_coords[1] < -0.2 || barycentric_coords[2] < -0.2) // The tetrahedron is a sliver
            {
                return remove_sliver(tid);
            }
            
            real mean_dist = 0.;
            for(vec3 &p : verts)
            {
                mean_dist += (p-proj_apex).length()/3.;
            }
            int close = 0;
            for(vec3 &p : verts)
            {
                if((p-proj_apex).length() < mean_dist)
                {
                    close++;
                }
            }
            
            if(close == 2) // The tetrahedron is a needle
            {
                return remove_needle(tid);
            }
            else if(close == 1) // The tetrahedron is a wedge
            {
                return remove_wedge(tid);
            }
            return false;
        }
        
        /**
         * Attempt to remove tetrahedra with quality lower than MIN_TET_QUALITY.
         */
        void remove_tets()
        {
            std::vector<tet_key> tets;
            
            for (auto tit = ISMesh::tetrahedra_begin(); tit != ISMesh::tetrahedra_end(); tit++)
            {
                if (quality(tit.key()) < MIN_TET_QUALITY)
                {
                    tets.push_back(tit.key());
                }
            }
            int i = 0, j=0;
            for (auto &tet : tets)
            {
                if (ISMesh::exists(tet) && quality(tet) < MIN_TET_QUALITY)
                {
                    if(remove_tet(tet))
                    {
                        i++;
                    }
                    j++;
                }
            }
            std::cout << "Removed " << i <<"/"<< j << " low quality tets" << std::endl;
            ISMesh::garbage_collect();
        }
        
        ///////////////
        // SMOOTHING //
        ///////////////
    private:
        /**
         * Performs Laplacian smoothing if it improves the minimum tetrahedron quality locally.
         */
        bool smart_laplacian(const node_key& nid, real alpha = 1.)
        {
            is_mesh::SimplexSet<tet_key> tids = get_tets(nid);
            is_mesh::SimplexSet<face_key> fids = get_faces(tids) - get_faces(nid);
            
            vec3 old_pos = get_pos(nid);
            vec3 avg_pos = get_barycenter(get_nodes(fids));
            vec3 new_pos = old_pos + alpha * (avg_pos - old_pos);
            
            if(min_quality_improvement(fids, old_pos, new_pos) > 0.)
            {
                set_pos(nid, new_pos);
                return true;
            }
            return false;
        }
        
        void smooth()
        {
            int i = 0, j = 0;
            for (auto nit = ISMesh::nodes_begin(); nit != ISMesh::nodes_end(); nit++)
            {
                if (is_safe_editable(nit.key()))
                {
                    if (smart_laplacian(nit.key()))
                    {
                        i++;
                    }
                    j++;
                }
            }
            std::cout << "Smoothed: " << i << "/" << j << std::endl;
        }
        
        ///////////////////
        // FIX FUNCTIONS //
        ///////////////////
        
        void fix_complex()
        {
            smooth();
            
            topological_edge_removal();
            topological_face_removal();
            
//            remove_tets();
//            remove_faces();
//            remove_edges();
            
            remove_degenerate_tets();
            remove_degenerate_faces();
            remove_degenerate_edges();
        }
        
        void resize_complex()
        {
            thickening_interface();
            
            thinning_interface();
            
            thickening();
            
            thinning();
            
            fix_complex();
        }
        
        ////////////////////
        // MOVE FUNCTIONS //
        ////////////////////
    public:
        /**
         * Moves all the vertices to their destination which can be set by the set_destination() function.
         */
        void deform(int num_steps = 10)
        {
            std::cout << std::endl << "********************************" << std::endl;
            int missing;
            int step = 0;
            do {
                std::cout << "\nMove vertices step " << step << std::endl;
                missing = 0;
                int movable = 0;
                for (auto nit = ISMesh::nodes_begin(); nit != ISMesh::nodes_end(); nit++)
                {
                    if (is_movable(nit.key()))
                    {
                        if(!move_vertex(nit.key()))
                        {
                            missing++;
                        }
                        movable++;
                    }
                }
                std::cout << "Vertices missing to be moved: " << missing <<"/" << movable << std::endl;
                
                fix_complex();
                
                ++step;
            } while (missing > 0 && step < num_steps);
            
            resize_complex();
            
            ISMesh::garbage_collect();
            update_attributes();
            ISMesh::validity_check();
        }
        
    private:
        
        /**
         * Tries moving the node n to the new position new_pos. Returns true if it succeeds.
         */
        bool move_vertex(const node_key & n)
        {
            vec3 pos = get_pos(n);
            vec3 destination = get(n).get_destination();
            real l = Util::length(destination - pos);
            
            if (l < 1e-4*AVG_EDGE_LENGTH) // The vertex is not moved
            {
                return true;
            }
            
            real max_l = l*intersection_with_link(n, destination) - MIN_DEFORMATION;
            l = Util::max(Util::min(0.5*max_l, l), 0.);
            set_pos(n, pos + l*Util::normalize(destination - pos));
            
            if (Util::length(destination - get_pos(n)) < 1e-4*AVG_EDGE_LENGTH)
            {
                return true;
            }
            return false;
        }
        
        
    public:
        /**
         * Returns the intersection point (= pos + t*(new_pos-pos)) with the link of the node n and
         * when moving the node n to the new position new_pos.
         */
        real intersection_with_link(const node_key & n, const vec3& destination)
        {
            vec3 pos = get_pos(n);
            vec3 ray = destination - pos;

            real min_t = INFINITY;
            auto fids = ISMesh::get_faces(ISMesh::get_tets(n)) - ISMesh::get_faces(n);
            for(auto f : fids)
            {
                auto face_pos = ISMesh::get_pos(ISMesh::get_nodes(f));
                real t = Util::intersection_ray_plane<real>(pos, ray, face_pos[0], face_pos[1], face_pos[2]);
                if (0. <= t)
                {
                    min_t = Util::min(t, min_t);
                }
            }
            assert(min_t < INFINITY);
            return min_t;
        }
        
        ///////////
        // FLIPS //
        ///////////
    private:
        
        
        /**
         * Returns whether it is possible to flip the edge e or not, i.e. whether the edge is a feature edge
         * (it is not a feature edge if its neighborhood is sufficiently flat).
         */
        bool is_flippable(const edge_key & eid)
        {
            is_mesh::SimplexSet<face_key> fids;
            for(auto f : get_faces(eid))
            {
                if (get(f).is_interface() || get(f).is_boundary())
                {
                    fids += f;
                }
            }
            if(fids.size() != 2)
            {
                return false;
            }
            
            is_mesh::SimplexSet<node_key> e_nids = get_nodes(eid);
            is_mesh::SimplexSet<node_key> new_e_nids = (get_nodes(fids[0]) + get_nodes(fids[1])) - e_nids;
            assert(new_e_nids.size() == 2);
            
            // Check that there does not already exist an edge.
            if(ISMesh::get_edge(new_e_nids[0], new_e_nids[1]).is_valid())
            {
                return false;
            }
            
            // Check that the edge is not a feature edge if it is a part of the interface.
            if(get(eid).is_interface())
            {
                real angle = cos_dihedral_angle(fids[0], fids[1]);
                if(angle < FLIP_EDGE_INTERFACE_FLATNESS)
                {
                    return false;
                }
            }
            
            return true;
        }
        
        bool precond_flip_edge(const edge_key& eid, const face_key& f1, const face_key& f2)
        {
            is_mesh::SimplexSet<node_key> e_nids = get_nodes(eid);
            is_mesh::SimplexSet<node_key> new_e_nids = (get_nodes(f1) + get_nodes(f2)) - e_nids;
            is_mesh::SimplexSet<node_key> apices = (get_nodes(get_faces(eid)) - e_nids) - new_e_nids;
            assert(e_nids.size() == 2);
            assert(new_e_nids.size() == 2);
            
            // Check that there does not already exist an edge.
            if(ISMesh::get_edge(new_e_nids[0], new_e_nids[1]).is_valid())
            {
                return false;
            }
            
            vec3 p = get_pos(new_e_nids[0]);
            vec3 r = get_pos(new_e_nids[1]) - p;
            vec3 a = get_pos(e_nids[0]);
            vec3 b = get_pos(e_nids[1]);
            
            for (node_key n : apices) {
                vec3 c = get_pos(n);
                real t = Util::intersection_ray_plane<real>(p, r, a, b, c);
                if(t > 0. && t < 1.)
                {
                    std::vector<real> coords = Util::barycentric_coords<real>(p + t*r, c, a, b);
                    if(coords[0] > MIN_DEFORMATION && coords[1] > MIN_DEFORMATION && coords[2] >= 0.)
                    {
                        return true;
                    }
                }
            }
            
            return false;
        }
        
        ////////////
        // SPLITS //
        ////////////
    public:
        /**
         * Split a tetrahedron t and returns the new node which is positioned at the barycenter of the vertices of t.
         */
        node_key split(const tet_key& tid)
        {
            is_mesh::SimplexSet<edge_key> eids = ISMesh::get_edges(tid);
            edge_key eid = longest_edge(eids);
            return split(eid);
        }
        
        /**
         * Split a face f and returns the new node which is positioned at the barycenter of the vertices of f.
         */
        node_key split(const face_key& fid)
        {
            is_mesh::SimplexSet<edge_key> eids = ISMesh::get_edges(fid);
            edge_key eid = longest_edge(eids);
            return split(eid);
        }
        
        /**
         * Split an edge e and returns the new node which is placed at the middle of e.
         */
        node_key split(const edge_key& eid)
        {
            auto verts = ISMesh::get_pos(ISMesh::get_nodes(eid));
            vec3 pos = Util::barycenter(verts[0], verts[1]);
            vec3 destination = pos;
            if(get(eid).is_interface())
            {
                auto dests = get_dest(get_nodes(eid));
                destination = Util::barycenter(dests[0], dests[1]);
            }
            
            return ISMesh::split(eid, pos, destination);
        }
        
        ///////////////
        // COLLAPSES //
        ///////////////
    private:
        /**
         * Collapses the edge e and places the new node at the most optimal position of the position of either end node or their barycenter.
         * If the parameter safe is true, the method if the nodes of edge e are editable, i.e. not a part of the interface, and will therefore not change the interface.
         * Returns whether the collapse was successful.
         */
        bool collapse(edge_key& eid, bool safe = true)
        {
            if (!ISMesh::exists(eid) || !eid.is_valid())
            {
                return node_key();
            }
            is_mesh::SimplexSet<node_key> nids = get_nodes(eid);
            is_mesh::SimplexSet<tet_key> e_tids = get_tets(eid);
            is_mesh::SimplexSet<face_key> fids0 = get_faces(get_tets(nids[0]) - e_tids) - get_faces(nids[0]);
            is_mesh::SimplexSet<face_key> fids1 = get_faces(get_tets(nids[1]) - e_tids) - get_faces(nids[1]);
            
            real q0 = min_quality(fids0, get_pos(nids[0]));
            real q1 = min_quality(fids1, get_pos(nids[1]));
            
            vec3 pos_opt, destination_opt;
            real q_max = -INFINITY;
            
            if ((safe && is_safe_editable(nids[0]) && is_safe_editable(nids[1])) || (!safe && is_unsafe_editable(nids[0]) && is_unsafe_editable(nids[1])))
            {
                vec3 p = Util::barycenter(get_pos(nids[0]), get_pos(nids[1]));
                real q = Util::min(min_quality(fids0, get_pos(nids[0]), p), min_quality(fids1, get_pos(nids[1]), p));
                
                if (q > q_max)
                {
                    destination_opt = Util::barycenter(get_dest(nids[0]), get_dest(nids[1]));
                    pos_opt = p;
                    q_max = q;
                }
            }
            
            if ((safe && is_safe_editable(nids[0])) || (!safe && is_unsafe_editable(nids[0])))
            {
                vec3 p = get_pos(nids[1]);
                real q = Util::min(min_quality(fids0, get_pos(nids[0]), p), q1);
                
                if (q > q_max)
                {
                    destination_opt = get_dest(nids[1]);
                    pos_opt = p;
                    q_max = q;
                }
            }
            
            if ((safe && is_safe_editable(nids[1])) || (!safe && is_unsafe_editable(nids[1])))
            {
                vec3 p = ISMesh::get_pos(nids[0]);
                real q = Util::min(q0, min_quality(fids1, get_pos(nids[1]), p));
                
                if (q > q_max)
                {
                    destination_opt = get_dest(nids[0]);
                    pos_opt = p;
                    q_max = q;
                }
            }
            
            if(!safe && q_max > EPSILON)
            {
                ISMesh::collapse(eid, pos_opt, destination_opt);
                return true;
            }
            real q_old = Util::min(Util::min(min_quality(e_tids), q1), q0);
            if(q_max > Util::min(q_old, MIN_TET_QUALITY) + EPSILON)
            {
                ISMesh::collapse(eid, pos_opt, destination_opt);
                return true;
            }
            return false;
        }
        
        bool collapse(is_mesh::SimplexSet<edge_key>& eids, bool safe)
        {
            while(eids.size() > 0)
            {
                edge_key e = shortest_edge(eids);
                if(collapse(e, safe))
                {
                    return true;
                }
                eids -= e;
            }
            return false;
        }
        
        bool collapse(const face_key& fid, bool safe = true)
        {
            is_mesh::SimplexSet<edge_key> eids = get_edges(fid);
            return collapse(eids, safe);
        }
        
        bool collapse(const tet_key& tid, bool safe = true)
        {
            is_mesh::SimplexSet<edge_key> eids = get_edges(tid);
            return collapse(eids, safe);
        }
        
        //////////////////////
        // GETTER FUNCTIONS //
        //////////////////////
    public:
        /**
         Returns the normal to interface face fid.
         */
        vec3 get_normal(const face_key& fid)
        {
            auto pos = ISMesh::get_pos(ISMesh::get_sorted_nodes(fid));
            return Util::normal_direction(pos[0], pos[1], pos[2]);
        }
        
        /**
         Returns the normal to interface node n.
         */
        vec3 get_normal(const node_key& nid)
        {
            vec3 result(0.);
            for (auto f : ISMesh::get_faces(nid))
            {
                if (get(f).is_interface())
                {
                    result += get_normal(f);
                }
            }
            if (Util::length(result) < EPSILON) {
                return vec3(0.);
            }
#ifdef DEBUG
            assert(!Util::isnan(result[0]) && !Util::isnan(result[1]) && !Util::isnan(result[2]));
#endif
            return Util::normalize(result);
        }
        
        /**
         * Calculates the average position of the nodes in the simplex set nids.
         * If interface is true, the average position is only calculated among the nodes which are interface.
         */
        vec3 get_barycenter(const is_mesh::SimplexSet<node_key>& nids, bool interface = false)
        {
            vec3 avg_pos(0.);
            int i = 0;
            for (auto n : nids)
            {
                if (!interface || get(n).is_interface())
                {
                    avg_pos += get_pos(n);
                    i++;
                }
            }
#ifdef DEBUG
            assert(i != 0);
#endif
            return avg_pos / static_cast<real>(i);
        }
        
        /**
         * Calculates the average position of the neighbouring nodes to node n.
         * If interface is true, the average position is only calculated among the neighbouring nodes which are interface.
         */
        vec3 get_barycenter(const node_key& nid, bool interface = false)
        {
            if(interface && !get(nid).is_interface())
            {
                return get_pos(nid);
            }
            
            is_mesh::SimplexSet<node_key> nids = get_nodes(get_tets(nid)) - nid;
            return get_barycenter(nids, interface);
        }
        
        ///////////////////////
        // UTILITY FUNCTIONS //
        ///////////////////////
        
    public:
        
        real length(const edge_key& eid)
        {
            is_mesh::SimplexSet<node_key> nids = get_nodes(eid);
            return Util::length(get_pos(nids[0]) - get_pos(nids[1]));
        }
        
        real length_destination(const edge_key& eid)
        {
            is_mesh::SimplexSet<node_key> nids = get_nodes(eid);
            return Util::length(get_dest(nids[0]) - get_dest(nids[1]));
        }
        
        real area(const face_key& fid)
        {
            is_mesh::SimplexSet<node_key> nids = get_nodes(fid);
            return Util::area<real>(get_pos(nids[0]), get_pos(nids[1]), get_pos(nids[2]));
        }
        
        real area_destination(const face_key& fid)
        {
            is_mesh::SimplexSet<node_key> nids = get_nodes(fid);
            return Util::area<real>(get_dest(nids[0]), get_dest(nids[1]), get_dest(nids[2]));
        }
        
        real volume(const tet_key& tid)
        {
            is_mesh::SimplexSet<node_key> nids = get_nodes(tid);
            return Util::volume<real>(get_pos(nids[0]), get_pos(nids[1]), get_pos(nids[2]), get_pos(nids[3]));
        }
        
        real volume_destination(const tet_key& tid)
        {
            is_mesh::SimplexSet<node_key> nids = get_nodes(tid);
            return Util::volume<real>(get_dest(nids[0]), get_dest(nids[1]), get_dest(nids[2]), get_dest(nids[3]));
        }
        
        real quality(const tet_key& tid)
        {
            is_mesh::SimplexSet<node_key> nids = get_nodes(tid);
            return std::abs(Util::quality<real>(get_pos(nids[0]), get_pos(nids[1]), get_pos(nids[2]), get_pos(nids[3])));
        }
        
        real min_angle(const face_key& fid)
        {
            is_mesh::SimplexSet<node_key> nids = get_nodes(fid);
            return Util::min_angle<real>(get_pos(nids[0]), get_pos(nids[1]), get_pos(nids[2]));
        }
        
        real max_angle(const face_key& fid)
        {
            is_mesh::SimplexSet<node_key> nids = get_nodes(fid);
            return Util::max_angle<real>(get_pos(nids[0]), get_pos(nids[1]), get_pos(nids[2]));
        }
        
        real quality(const face_key& fid)
        {
            is_mesh::SimplexSet<node_key> nids = get_nodes(fid);
            auto angles = Util::cos_angles<real>(get_pos(nids[0]), get_pos(nids[1]), get_pos(nids[2]));
            real worst_a;
            for(auto a : angles)
            {
                worst_a = std::max(worst_a, std::abs(a));
            }
            return 1. - worst_a;
        }
        
        real quality(const edge_key& eid)
        {
            return length(eid)/AVG_EDGE_LENGTH;
        }
        
        /**
         * Returns the largest face in the simplex set.
         */
        face_key largest_face(const is_mesh::SimplexSet<face_key>& fids)
        {
            real max_a = -INFINITY;
            face_key max_f;
            for(auto f : fids)
            {
                real a = area(f);
                if(a > max_a)
                {
                    max_a = a;
                    max_f = f;
                }
            }
            return max_f;
        }
        
        /**
         * Returns the shortest edge in the simplex set.
         */
        edge_key shortest_edge(const is_mesh::SimplexSet<edge_key>& eids)
        {
            real min_l = INFINITY;
            edge_key min_e;
            for(auto e : eids)
            {
                real l = length(e);
                if(l < min_l)
                {
                    min_l = l;
                    min_e = e;
                }
            }
            return min_e;
        }
        
        /**
         * Returns the longest edge in the simplex set.
         */
        edge_key longest_edge(const is_mesh::SimplexSet<edge_key>& eids)
        {
            real max_l = -INFINITY;
            edge_key max_e;
            for(auto e : eids)
            {
                real l = length(e);
                if(l > max_l)
                {
                    max_l = l;
                    max_e = e;
                }
            }
            return max_e;
        }
        
        /**
         * Returns the minimum quality of the tetrahedra in simplex set s.
         */
        real min_quality(const is_mesh::SimplexSet<tet_key>& tids)
        {
            real q_min = INFINITY;
            for (auto t : tids)
            {
                q_min = Util::min(quality(t), q_min);
            }
            return q_min;
        }
        
        /**
         * Returns the minimum tetrahedral quality of a node with position pos. The faces in the link of the node should be passed in fids.
         */
        real min_quality(const is_mesh::SimplexSet<face_key>& fids, const vec3& pos)
        {
            real min_q = INFINITY;
            for (auto f : fids)
            {
                is_mesh::SimplexSet<node_key> nids = get_nodes(f);
                min_q = Util::min(min_q, std::abs(Util::quality<real>(get_pos(nids[0]), get_pos(nids[1]), get_pos(nids[2]), pos)));
            }
            return min_q;
        }
        
        /**
         * Returns the new minimum tetrahedral quality when moving a node from old_pos to new_pos. The faces in the link of the node should be passed in fids.
         */
        real min_quality(const is_mesh::SimplexSet<face_key>& fids, const vec3& pos_old, const vec3& pos_new)
        {
            real min_q = INFINITY;
            for (auto f : fids)
            {
                is_mesh::SimplexSet<node_key> nids = get_nodes(f);
                if(Util::sign(Util::signed_volume<real>(get_pos(nids[0]), get_pos(nids[1]), get_pos(nids[2]), pos_old)) !=
                   Util::sign(Util::signed_volume<real>(get_pos(nids[0]), get_pos(nids[1]), get_pos(nids[2]), pos_new)))
                {
                    return -INFINITY;
                }
                min_q = Util::min(min_q, std::abs(Util::quality<real>(get_pos(nids[0]), get_pos(nids[1]), get_pos(nids[2]), pos_new)));
            }
            return min_q;
        }
        
        /**
         * Returns the improvement in minimum tetrahedral quality of moving a node from old_pos to new_pos. The faces in the link of the node should be passed in fids.
         */
        real min_quality_improvement(const is_mesh::SimplexSet<face_key>& fids, const vec3& pos_old, const vec3& pos_new)
        {
            real min_q_old = INFINITY;
            real min_q_new = INFINITY;
            for (auto f : fids)
            {
                is_mesh::SimplexSet<node_key> nids = get_nodes(f);
                if(Util::sign(Util::signed_volume<real>(get_pos(nids[0]), get_pos(nids[1]), get_pos(nids[2]), pos_old)) !=
                   Util::sign(Util::signed_volume<real>(get_pos(nids[0]), get_pos(nids[1]), get_pos(nids[2]), pos_new)))
                {
                    return -INFINITY;
                }
                min_q_old = Util::min(min_q_old, std::abs(Util::quality<real>(get_pos(nids[0]), get_pos(nids[1]), get_pos(nids[2]), pos_old)));
                min_q_new = Util::min(min_q_new, std::abs(Util::quality<real>(get_pos(nids[0]), get_pos(nids[1]), get_pos(nids[2]), pos_new)));
            }
            return min_q_new - min_q_old;
        }
        
        
    private:
        
        /**
         * Check if the sequence of vertices in polygon is consistent with positive orientation of tetrahedra in the mesh
         * with respect to the ordered pair of vertices in vv. If not, reverse the order of vertices in polygon.
         */
        void check_consistency(const is_mesh::SimplexSet<node_key>& nids, is_mesh::SimplexSet<node_key>& polygon)
        {
            unsigned int n = static_cast<unsigned int>(polygon.size());
            
            real sum = 0;
            for (unsigned int i = 0; i < n; ++i)
            {
                sum += Util::signed_volume<real>(get_pos(nids[0]), get_pos(nids[1]), get_pos(polygon[i]), get_pos(polygon[(i+1)%n]));
            }
            
            if (sum < 0.)
            {
                for (unsigned int i = 0; i < n/2; ++i)
                {
                    polygon.swap(i, n-1-i);
                }
            }
        }
        
        ////////////////////////
        // DOCUMENT FUNCTIONS //
        ////////////////////////
    public:
        ///
        void extract_interface(std::vector<vec3>& verts, std::vector<int>& indices)
        {
            std::map<node_key, int> vert_index;
            
            // Extract vertices
            for (auto nit = ISMesh::nodes_begin(); nit != ISMesh::nodes_end(); nit++)
            {
                if (nit->is_interface())
                {
                    verts.push_back(ISMesh::get_pos(nit.key()));
                    vert_index[nit.key()] = static_cast<int>(verts.size());
                }
            }
            
            // Extract faces
            for (auto fit = ISMesh::faces_begin(); fit != ISMesh::faces_end(); fit++)
            {
                if (fit->is_interface())
                {
                    auto nodes = get_sorted_nodes(fit.key());
                    
                    indices.push_back(vert_index[nodes[0]]);
                    indices.push_back(vert_index[nodes[1]]);
                    indices.push_back(vert_index[nodes[2]]);
                }
            }
        }
        
        void extract_tet_mesh(std::vector<vec3>& points, std::vector< std::vector<int> >& tets)
        {
            ISMesh::garbage_collect();
            
            std::map<node_key, int> indices;
            int counter = 0;
            for (auto nit = ISMesh::nodes_begin(); nit != ISMesh::nodes_end(); nit++)
            {
                indices[nit.key()] = counter;
                points.push_back(nit->get_pos());
                ++counter;
            }
            
            for (auto tit = ISMesh::tetrahedra_begin(); tit != ISMesh::tetrahedra_end(); tit++)
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
        
        /**
         * Returns the cosine to the dihedral angle between face f1 and face f2.
         */
        real cos_dihedral_angle(const face_key& f1, const face_key& f2)
        {
            auto nids1 = ISMesh::get_nodes(f1);
            auto nids2 = ISMesh::get_nodes(f2);
            is_mesh::SimplexSet<node_key> nids = nids1 & nids2;
            is_mesh::SimplexSet<node_key> apices = (nids1 + nids2) - nids;
            
            return Util::cos_dihedral_angle<real>(get_pos(nids[0]), get_pos(nids[1]), get_pos(apices[0]), get_pos(apices[1]));
        }
        
        /**
         * Returns the dihedral angle between face f1 and face f2.
         */
        real dihedral_angle(const face_key& f1, const face_key& f2)
        {
            return acos(cos_dihedral_angle(f1, f2));
        }
        
        std::vector<real> cos_dihedral_angles(const tet_key& tid)
        {
            auto verts = ISMesh::get_pos(ISMesh::get_nodes(tid));
            std::vector<real> angles;
            std::vector<int> apices;
            for (int i = 0; i < verts.size(); i++) {
                for (int j = 0; j < verts.size(); j++) {
                    if(i < j)
                    {
                        apices.clear();
                        for (int k = 0; k < verts.size(); k++) {
                            if(k != i && k != j)
                            {
                                apices.push_back(k);   
                            }
                        }
                        angles.push_back(Util::cos_dihedral_angle<real>(verts[i], verts[j], verts[apices[0]], verts[apices[1]]));
                    }
                }
            }
            return angles;
        }
        
        /**
         * Returns the cosine of the minimum dihedral angle between the faces of tetrahedron t.
         */
        real min_cos_dihedral_angle(const tet_key& t)
        {
            real min_angle = -1.;
            std::vector<real> angles = cos_dihedral_angles(t);
            for(auto a : angles)
            {
                min_angle = Util::max(min_angle, a);
            }
            return min_angle;
        }
        
        /**
         * Returns the minimum dihedral angle between the faces of tetrahedron t.
         */
        real min_dihedral_angle(const tet_key& t)
        {
            return acos(min_cos_dihedral_angle(t));
        }
        
        void get_qualities(std::vector<int>& histogram, real& min_quality)
        {
            min_quality = INFINITY;
            
            histogram = std::vector<int>(100);
            for (int i = 0; i < 100; ++i)
            {
                histogram[i] = 0;
            }
            
            for (auto tit = ISMesh::tetrahedra_begin(); tit != ISMesh::tetrahedra_end(); tit++)
            {
                real q = quality(tit.key());
                min_quality = Util::min(min_quality, q);
                int index = static_cast<int>(floor(q*100.));
#ifdef DEBUG
                assert(index < 100 && index >= 0);
#endif
                histogram[index] += 1;
            }
        }
        
        /**
         * Calculates the dihedral angles in the SimplicialComplex and returns these in a histogram,
         * along with the minimum and maximum dihedral angles.
         */
        void get_dihedral_angles(std::vector<int> & histogram, real & min_angle, real & max_angle)
        {
            max_angle = -INFINITY, min_angle = INFINITY;
            
            histogram = std::vector<int>(180);
            for (int i = 0; i < 180; ++i)
            {
                histogram[i] = 0;
            }
            
            for (auto tit = ISMesh::tetrahedra_begin(); tit != ISMesh::tetrahedra_end(); tit++)
            {
                std::vector<real> angles = cos_dihedral_angles(tit.key());
                for(auto cos_a : angles)
                {
                    real a = acos(cos_a)*180./M_PI;
                    min_angle = Util::min(min_angle, a);
                    max_angle = Util::max(max_angle, a);
                    histogram[(int)floor(a)] += 1;
                }
            }
        }
        
        real min_quality()
        {
            real min_q = INFINITY;
            for (auto tit = ISMesh::tetrahedra_begin(); tit != ISMesh::tetrahedra_end(); tit++)
            {
                min_q = Util::min(min_q, quality(tit.key()));
            }
            return min_q;
        }
        
        /// Counts the total number of nodes and the number of nodes on the interface(s).
        void count_nodes(int & total, int & object)
        {
            total = 0, object = 0;
            for (auto nit = ISMesh::nodes_begin(); nit != ISMesh::nodes_end(); nit++)
            {
                total++;
                if (nit->is_interface())
                {
                    object++;
                }
            }
        }
        
        /// Counts the total number of edges and the number of edges on the interface(s).
        void count_edges(int & total, int & object)
        {
            total = 0, object = 0;
            for (auto eit = ISMesh::edges_begin(); eit != ISMesh::edges_end(); eit++)
            {
                total++;
                if (eit->is_interface())
                {
                    object++;
                }
            }
        }
        
        /// Counts the total number of faces and the number of faces on the interface(s).
        void count_faces(int & total, int & object)
        {
            total = 0, object = 0;
            for (auto fit = ISMesh::faces_begin(); fit != ISMesh::faces_end(); fit++)
            {
                total++;
                if (fit->is_interface())
                {
                    object++;
                }
            }
        }
        
        /// Counts the total number of tetrahedra and the number of tetrahedra in the object(s).
        void count_tetrahedra(int & total, int & object)
        {
            total = 0, object = 0;
            for (auto tit = ISMesh::tetrahedra_begin(); tit != ISMesh::tetrahedra_end(); tit++)
            {
                total++;
                if (tit->label() != 0)
                {
                    object++;
                }
            }
        }
        
    public:
        void test_split_collapse()
        {
            is_mesh::SimplexSet<edge_key> eids;
            for (auto eit = ISMesh::edges_begin(); eit != ISMesh::edges_end(); eit++)
            {
                auto neighbours = ISMesh::get_edges(ISMesh::get_faces(eit.key()));
                bool ok = true;
                for(auto e : neighbours)
                {
                    if(eids.contains(e))
                    {
                        ok = false;
                    }
                }
                if (ok)
                {
                    eids += eit.key();
                }
            }
            
            int j = 0;
            std::cout << "Split test # = " << eids.size();
            is_mesh::SimplexSet<edge_key> new_eids;
            std::vector<vec3> verts;
            for (auto e : eids) {
                auto nids = ISMesh::get_nodes(e);
                auto new_nid = split(e);
                auto new_eid = (ISMesh::get_edges(nids) & ISMesh::get_edges(new_nid)) - e;
                assert(new_eid.size() == 1);
                new_eids += new_eid[0];
                auto old_nid = ISMesh::get_nodes(new_eid) - new_nid;
                assert(old_nid.size() == 1);
                verts.push_back(ISMesh::get_pos(old_nid[0]));
                j++;
                if(j%1000 == 0)
                {
                    std::cout << ".";
                }
            }
            std::cout << " DONE" << std::endl;
            ISMesh::garbage_collect();
            ISMesh::validity_check();
            
            std::cout << "Collapse test # = " << new_eids.size();
            j = 0;
            for (int i = 0; i < new_eids.size(); i++) {
                assert(ISMesh::exists(new_eids[i]));
                auto nid = ISMesh::collapse(new_eids[i], verts[i], verts[i]);
                assert(nid.is_valid());
                j++;
                if(j%1000 == 0)
                {
                    std::cout << ".";
                }
            }
            std::cout << " DONE" << std::endl;
            ISMesh::garbage_collect();
            ISMesh::validity_check();
        }
        
        void test_flip23_flip32()
        {
            is_mesh::SimplexSet<face_key> fids;
            for (auto fit = ISMesh::faces_begin(); fit != ISMesh::faces_end(); fit++)
            {
                if(is_safe_editable(fit.key()))
                {
                    auto nids = get_nodes(fit.key());
                    nids += get_nodes(get_tets(fit.key()));
                    real t = Util::intersection_ray_triangle<real>(get_pos(nids[3]), get_pos(nids[4]) - get_pos(nids[3]), get_pos(nids[0]), get_pos(nids[1]), get_pos(nids[2]));
                    
                    auto neighbours = get_faces(get_tets(fit.key()));
                    bool ok = true;
                    for(auto f : neighbours)
                    {
                        if(fids.contains(f))
                        {
                            ok = false;
                        }
                    }
                    if (ok && 0 < t && t < 1)
                    {
                        fids += fit.key();
                    }
                }
            }
            
            std::cout << "Flip 2-3 test # = " << fids.size();
            is_mesh::SimplexSet<edge_key> new_eids;
            int i = 0;
            for (auto f : fids) {
                assert(ISMesh::exists(f));
                auto new_eid = ISMesh::flip_23(f);
                assert(new_eid.is_valid());
                new_eids += new_eid;
                i++;
                if(i%1000 == 0)
                {
                    std::cout << ".";
                }
            }
            std::cout << " DONE" << std::endl;
            ISMesh::garbage_collect();
            ISMesh::validity_check();
            
            i=0;
            std::cout << "Flip 3-2 test # = " << new_eids.size();
            for (auto e : new_eids) {
                auto new_fid = ISMesh::flip_32(e);
                assert(new_fid.is_valid());
                i++;
                if(i%1000 == 0)
                {
                    std::cout << ".";
                }
            }
            std::cout << " DONE" << std::endl;
            ISMesh::garbage_collect();
            ISMesh::validity_check();
        }
        
        void test_flip44()
        {
            is_mesh::SimplexSet<edge_key> eids;
            for (auto eit = ISMesh::edges_begin(); eit != ISMesh::edges_end(); eit++)
            {
                if(is_unsafe_editable(eit.key()) && eit->is_interface() && get_faces(eit.key()).size() == 4)
                {
                    auto neighbours = get_edges(get_tets(eit.key()));
                    bool ok = true;
                    for(auto e : neighbours)
                    {
                        if(eids.contains(e))
                        {
                            ok = false;
                        }
                    }
                    is_mesh::SimplexSet<face_key> flip_fids;
                    for(auto f : get_faces(eit.key()))
                    {
                        if(get(f).is_interface())
                        {
                            flip_fids += f;
                        }
                    }
                    assert(flip_fids.size() == 2);
                    
                    if (ok && precond_flip_edge(eit.key(), flip_fids[0], flip_fids[1]))
                    {
                        eids += eit.key();
                    }
                }
            }
            
            for(int t = 0; t < 2; t++)
            {
                std::cout << "Flip 4-4 test # = " << eids.size();
                int i = 0;
                for (auto e : eids) {
                    is_mesh::SimplexSet<face_key> flip_fids;
                    for(auto f : get_faces(e))
                    {
                        if(get(f).is_interface())
                        {
                            flip_fids += f;
                        }
                    }
                    assert(flip_fids.size() == 2);
                    assert(get_faces(e).size() == 4);
                    ISMesh::flip_44(flip_fids[0], flip_fids[1]);
                    i++;
                    if(i%100 == 0)
                    {
                        std::cout << ".";
                    }
                }
                std::cout << " DONE" << std::endl;
                ISMesh::garbage_collect();
                ISMesh::validity_check();
            }
        }
        
        void test_flip22()
        {
            is_mesh::SimplexSet<edge_key> eids;
            for (auto eit = ISMesh::edges_begin(); eit != ISMesh::edges_end(); eit++)
            {
                if(eit->is_boundary() && get_faces(eit.key()).size() == 3)
                {
                    auto neighbours = get_edges(get_tets(eit.key()));
                    bool ok = true;
                    for(auto e : neighbours)
                    {
                        if(eids.contains(e))
                        {
                            ok = false;
                        }
                    }
                    is_mesh::SimplexSet<face_key> flip_fids;
                    for(auto f : get_faces(eit.key()))
                    {
                        if(get(f).is_boundary())
                        {
                            flip_fids += f;
                        }
                    }
                    assert(flip_fids.size() == 2);
                    
                    if (ok && precond_flip_edge(eit.key(), flip_fids[0], flip_fids[1]))
                    {
                        eids += eit.key();
                    }
                }
            }
            
            for(int t = 0; t < 2; t++)
            {
                std::cout << "Flip 2-2 test # = " << eids.size();
                int i = 0;
                for (auto e : eids) {
                    assert(ISMesh::exists(e));
                    auto fids = get_faces(e);
                    assert(fids.size() == 3);
                    is_mesh::SimplexSet<face_key> flip_fids;
                    for(auto f : fids)
                    {
                        if(get(f).is_boundary())
                        {
                            flip_fids += f;
                        }
                    }
                    
                    assert(flip_fids.size() == 2);
                    ISMesh::flip_22(flip_fids[0], flip_fids[1]);
                    i++;
                    if(i%10 == 0)
                    {
                        std::cout << ".";
                    }
                }
                std::cout << " DONE" << std::endl;
                ISMesh::garbage_collect();
                ISMesh::validity_check();
            }
        }
        
    };
    
}
