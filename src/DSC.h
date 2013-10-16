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

#include "is_mesh_API.h"
#include "util.h"
#include "attributes.h"
#include "design_domain.h"

namespace DSC {
    
    template <typename node_att = NodeAttributes, typename edge_att = EdgeAttributes, typename face_att = FaceAttributes, typename tet_att = TetAttributes>
    class DeformableSimplicialComplex : public is_mesh::ISMesh<node_att, edge_att, face_att, tet_att>
    {
        friend class ObjectGenerator;
        typedef is_mesh::ISMesh<node_att, edge_att, face_att, tet_att> Complex;
    public:
        
        typedef typename Complex::node_key      node_key;
        typedef typename Complex::edge_key      edge_key;
        typedef typename Complex::face_key      face_key;
        typedef typename Complex::tet_key       tet_key;
        typedef typename Complex::simplex_set   simplex_set;
        
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
            Complex(points, tets), design_domain(domain)
        {
            AVG_EDGE_LENGTH = _AVG_EDGE_LENGTH;
            MIN_DEFORMATION = 0.25 * AVG_EDGE_LENGTH;
            
            DEG_EDGE_QUALITY = 0.1;
            MIN_EDGE_QUALITY = 0.5;
            
            DEG_FACE_QUALITY = 1. - cos(5.*M_PI/180.);
            MIN_FACE_QUALITY = 1. - cos(10.*M_PI/180.);
            
            DEG_TET_QUALITY = 0.01;
            MIN_TET_QUALITY = 0.3;
            
            FLIP_EDGE_INTERFACE_FLATNESS = 0.995;
            
            MIN_LENGTH = 0.5 * AVG_EDGE_LENGTH;
            MAX_LENGTH = 2. * AVG_EDGE_LENGTH;
            
            real area_avg = AVG_EDGE_LENGTH*AVG_EDGE_LENGTH*0.5;
            MIN_AREA = 0.2*area_avg;
            MAX_AREA = 5.*area_avg;
            
            real vol_avg = AVG_EDGE_LENGTH*AVG_EDGE_LENGTH*AVG_EDGE_LENGTH*sqrt(2.)/12.;
            MIN_VOLUME = 0.5*vol_avg;
            MAX_VOLUME = 10.*vol_avg;
            
            //        fix_complex();
            //        resize_complex();
        }
        
        DeformableSimplicialComplex()
        {
            delete design_domain;
        }
        
    private:
        
        // For debugging!
        void print(const node_key& n)
        {
            std::cout << "Node: " << n << std::endl;
            vec3 p = get_pos(n);
            vec3 d = get_dest(n);
            std::cout << "P = " << p[0] << ", " << p[1] << ", " << p[2] << std::endl;
            std::cout << "D = " << d[0] << ", " << d[1] << ", " << d[2]  << std::endl;
            
            simplex_set st_n;
            Complex::star(n, st_n);
            std::cout << "\nStar_edges = [";
            for(auto eit = st_n.edges_begin(); eit != st_n.edges_end(); eit++)
            {
                auto verts = get_pos(*eit);
                vec3 p1 = verts[0];
                vec3 p2 = verts[1];
                std::cout << p1[0] << ", " << p1[1] << ", " << p1[2] << "; " << std::endl;
                std::cout << p2[0] << ", " << p2[1] << ", " << p2[2] << "; " << std::endl;
            }
            std::cout << "];" << std::endl;
            
            std::cout << "\nStar_Iedges = [";
            for(auto eit = st_n.edges_begin(); eit != st_n.edges_end(); eit++)
            {
                if(is_interface(*eit))
                {
                    auto verts = get_pos(*eit);
                    vec3 p1 = verts[0];
                    vec3 p2 = verts[1];
                    std::cout << p1[0] << ", " << p1[1] << ", " << p1[2] << "; " << std::endl;
                    std::cout << p2[0] << ", " << p2[1] << ", " << p2[2] << "; " << std::endl;
                }
            }
            std::cout << "];" << std::endl;
            
            std::cout << "\nedges = [";
            auto eids = Complex::get_edges(Complex::get_tets(n)) - Complex::get_edges(n);
            for(auto e : eids)
            {
                auto verts = get_pos(e);
                vec3 p1 = verts[0];
                vec3 p2 = verts[1];
                std::cout << p1[0] << ", " << p1[1] << ", " << p1[2] << "; " << std::endl;
                std::cout << p2[0] << ", " << p2[1] << ", " << p2[2] << "; " << std::endl;
            }
            std::cout << "];" << std::endl;
            
            std::cout << "\nIedges = [";
            for(auto e : eids)
            {
                if(is_interface(e))
                {
                    auto verts = get_pos(e);
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
    
        template<typename key>
        bool is_interface(const key& k)
        {
            return Complex::is_interface(k);
        }
        
        template<typename key>
        bool is_boundary(const key& k)
        {
            return Complex::is_boundary(k);
        }
        
        template<typename key>
        bool is_crossing(const key& k)
        {
            return Complex::is_crossing(k);
        }
        
    protected:
        
        virtual bool is_unsafe_editable(const node_key& nid)
        {
            return Complex::exists(nid) && !is_boundary(nid);
        }
        
        virtual bool is_unsafe_editable(const edge_key& eid)
        {
            return Complex::exists(eid) && !is_boundary(eid);
        }
        
        virtual bool is_unsafe_editable(const face_key& fid)
        {
            return Complex::exists(fid) && !is_boundary(fid);
        }
        
        virtual bool is_unsafe_editable(const tet_key& tid)
        {
            return Complex::exists(tid);
        }
        
        virtual bool is_safe_editable(const node_key& nid)
        {
            return is_unsafe_editable(nid) && !is_interface(nid);
        }
        
        virtual bool is_safe_editable(const edge_key& eid)
        {
            return is_unsafe_editable(eid) && !is_interface(eid);
        }
        
        virtual bool is_safe_editable(const face_key& fid)
        {
            return is_unsafe_editable(fid) && !is_interface(fid);
        }
        
        virtual bool is_safe_editable(const tet_key& tid)
        {
            return is_unsafe_editable(tid);
        }
        
    public:
        virtual bool is_movable(const node_key& nid)
        {
            return is_unsafe_editable(nid) && is_interface(nid) && !is_crossing(nid);
        }
        
        int get_label(const tet_key& t)
        {
            return Complex::get_label(t);
        }
        
        /**
         * Returns the position of node n.
         */
        vec3 get_pos(const node_key& n)
        {
            vec3 p = Complex::get(n).get_pos();
#ifdef DEBUG
            assert(!Util::isnan(p[0]) && !Util::isnan(p[1]) && !Util::isnan(p[2]));
#endif
            return p;
        }
        
        /// Returns the positions of the nodes of edge e.
        std::vector<vec3> get_pos(const edge_key & e)
        {
            std::vector<vec3> verts(2);
            auto nodes = Complex::get_nodes(e);
            for (int k = 0; k < 2; ++k)
            {
                verts[k] = get_pos(nodes[k]);
            }
            return verts;
        }
        
        /// Returns the positions of the nodes of face f.
        std::vector<vec3> get_pos(const face_key & f)
        {
            std::vector<vec3> verts(3);
            auto nodes = Complex::get_sorted_nodes(f);
            for (int k = 0; k < 3; ++k)
            {
                verts[k] = get_pos(nodes[k]);
            }
            return verts;
        }
        
        /// Returns the positions of the nodes of tetrahedron t.
        std::vector<vec3> get_pos(const tet_key& t)
        {
            std::vector<vec3> verts(4);
            auto nodes = Complex::get_sorted_nodes(t);
            for (int k = 0; k < 4; ++k)
            {
                verts[k] = get_pos(nodes[k]);
            }
            return verts;
        }
        
    protected:
        /**
         * Sets the position of node n.
         */
        void set_pos(const node_key& n, const vec3& p)
        {
            Complex::get(n).set_pos(p);
        }
        
    public:
        vec3 get_dest(const node_key& n)
        {
            if(is_movable(n))
            {
                return Complex::get(n).get_destination();
            }
            return get_pos(n);
        }
        
        /// Returns the destinations of the nodes of edge e.
        std::vector<vec3> get_dest(const edge_key & e)
        {
            std::vector<vec3> verts(2);
            auto nodes = Complex::get_nodes(e);
            for (int k = 0; k < 2; ++k)
            {
                verts[k] = get_dest(nodes[k]);
            }
            return verts;
        }
        
        /// Returns the destinations of the nodes of face f.
        std::vector<vec3> get_dest(const face_key & f)
        {
            std::vector<vec3> verts(3);
            auto nodes = Complex::get_nodes(f);
            for (int k = 0; k < 3; ++k)
            {
                verts[k] = get_dest(nodes[k]);
            }
            return verts;
        }
        
        /// Returns the destinations of the nodes of tetrahedron t.
        std::vector<vec3> get_dest(const tet_key& t)
        {
            std::vector<vec3> verts(4);
            auto nodes = Complex::get_nodes(t);
            for (int k = 0; k < 4; ++k)
            {
                verts[k] = get_dest(nodes[k]);
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
                    Complex::get(nid).set_destination(p + vec);
                }
                else {
                    Complex::get(nid).set_destination(dest);
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
        real build_table(const edge_key& e, const std::vector<node_key>& polygon, std::vector<std::vector<int>>& K)
        {
            auto verts = get_pos(e);
            
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
                        real q2 = Util::quality<real>(get_pos(polygon[i]), get_pos(polygon[k]), verts[0], get_pos(polygon[j]));
                        real q1 = Util::quality<real>(get_pos(polygon[k]), get_pos(polygon[i]), verts[1], get_pos(polygon[j]));
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
        
        std::vector<node_key> get_polygon(edge_key const & e)
        {
            simplex_set lk_e;
            Complex::link(e, lk_e);
            std::vector<node_key> polygon;
            
            sort_vertices(lk_e, polygon);
            check_consistency(get_pos(e), polygon);
            return polygon;
        }
        
        void get_half_polygon(const edge_key& e, std::vector<node_key>& polygon1, std::vector<node_key>& polygon2)
        {
            auto nodes = Complex::get_nodes(e);
            std::vector<node_key> polygon = get_polygon(e);
            int i1 = -1, i2 = -1;
            for(int i = 0; i < polygon.size(); i++)
            {
                face_key f = Complex::get_face(nodes[0], nodes[1], polygon[i]);
#ifdef DEBUG
                assert(f.is_valid());
#endif
                if(is_interface(f) || is_boundary(f))
                {
                    if(i2 != -1) // More than one boundary meets at edge e.
                    {
                        return;
                    }
                    
                    if(i1 == -1) {
                        i1 = i;
                    }
                    else {
                        i2 = i;
                    }
                }
            }
#ifdef DEBUG
            assert(i2 != -1);
#endif
            
            for(int i = i1; i != i2; i = (i+1)%polygon.size())
            {
                polygon1.push_back(polygon[i]);
            }
            polygon1.push_back(polygon[i2]);
            
            for(int i = i2; i != i1; i = (i+1)%polygon.size())
            {
                polygon2.push_back(polygon[i]);
            }
            polygon2.push_back(polygon[i1]);
            
            if(polygon1.size() <= 2)
            {
                swap(polygon1, polygon2);
            }
#ifdef DEBUG
            assert(polygon1.size() > 2);
#endif
        }
        
        
        void flip_23_recursively(const std::vector<node_key>& polygon, const node_key& n1, const node_key& n2, std::vector<std::vector<int>>& K, int i, int j)
        {
            if(j >= i+2)
            {
                int k = K[i][j];
                flip_23_recursively(polygon, n1, n2, K, i, k);
                flip_23_recursively(polygon, n1, n2, K, k, j);
                Complex::flip_23(Complex::get_face(n1, n2, polygon[k]));
            }
        }
        
        void topological_edge_removal(const std::vector<node_key>& polygon, const node_key& n1, const node_key& n2, std::vector<std::vector<int>>& K)
        {
            const int m = static_cast<int>(polygon.size());
            int k = K[0][m-1];
            flip_23_recursively(polygon, n1, n2, K, 0, k);
            flip_23_recursively(polygon, n1, n2, K, k, m-1);
            Complex::flip_32(Complex::get_edge(n1, n2));
        }
        
        /**
         * Attempt to remove edge e by mesh reconnection using the dynamic programming method by Klincsek (see Shewchuk "Two Discrete Optimization Algorithms
         * for the Topological Improvement of Tetrahedral Meshes" article for details).
         */
        bool topological_edge_removal(const edge_key& e)
        {
            std::vector<node_key> polygon = get_polygon(e);
            std::vector<std::vector<int>> K;
            real q_new = build_table(e, polygon, K);
            
            if (q_new > min_quality(e))
            {
                auto nodes = Complex::get_nodes(e);
                topological_edge_removal(polygon, nodes[0], nodes[1], K);
                return true;
            }
            return false;
        }
        
        void topological_boundary_edge_removal(const std::vector<node_key>& polygon1, const std::vector<node_key>& polygon2, const node_key& n1, const node_key& n2, std::vector<std::vector<int>>& K1, std::vector<std::vector<int>>& K2)
        {
            const int m1 = static_cast<int>(polygon1.size());
            const int m2 = static_cast<int>(polygon2.size());
            int k = K1[0][m1-1];
            flip_23_recursively(polygon1, n1, n2, K1, 0, k);
            flip_23_recursively(polygon1, n1, n2, K1, k, m1-1);
            
            if(m2 <= 2) {
                // Find the faces to flip about.
                face_key f1 = Complex::get_face(n1, n2, polygon1.front());
                face_key f2 = Complex::get_face(n1, n2, polygon1.back());
                
                if(is_boundary(f1) && is_boundary(f2))
                {
                    Complex::flip_22(f1, f2);
                }
            }
            else {
                k = K2[0][m2-1];
                flip_23_recursively(polygon2, n1, n2, K2, 0, k);
                flip_23_recursively(polygon2, n1, n2, K2, k, m2-1);
                
                // Find the faces to flip about.
                face_key f1 = Complex::get_face(n1, n2, polygon1.front());
                face_key f2 = Complex::get_face(n1, n2, polygon1.back());
                
                Complex::flip_44(f1, f2);
            }
        }
        
        bool topological_boundary_edge_removal(const edge_key& e)
        {
            std::vector<node_key> polygon1, polygon2;
            get_half_polygon(e, polygon1, polygon2);
            
            if(polygon1.size() <= 2)
            {
                return false;
            }
            
            std::vector<std::vector<int>> K1, K2;
            real q_new = build_table(e, polygon1, K1);
            
            if(polygon2.size() > 2)
            {
                q_new = Util::min(q_new, build_table(e, polygon2, K2));
            }
            
            if (q_new > min_quality(e))
            {
                auto nodes = Complex::get_nodes(e);
                topological_boundary_edge_removal(polygon1, polygon2, nodes[0], nodes[1], K1, K2);
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
            for (auto tit = Complex::tetrahedra_begin(); tit != Complex::tetrahedra_end(); tit++)
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
                if (Complex::exists(t) && quality(t) < MIN_TET_QUALITY)
                {
                    simplex_set cl_t;
                    Complex::closure(t, cl_t);
                    
                    for (auto eit = cl_t.edges_begin(); eit != cl_t.edges_end(); eit++)
                    {
                        if (Complex::exists(*eit))
                        {
                            if(is_safe_editable(*eit))
                            {
                                if(topological_edge_removal(*eit))
                                {
                                    i++;
                                }
                            }
                            else if(!is_boundary(*eit) && is_flippable(*eit))
                            {
                                if(topological_boundary_edge_removal(*eit))
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
            Complex::garbage_collect();
        }
        
        //////////////////////////////
        // TOPOLOGICAL FACE REMOVAL //
        //////////////////////////////
        
        face_key get_neighbour(const face_key& f, const edge_key& e)
        {
            simplex_set st_e;
            Complex::star(e, st_e);
            if(st_e.size_faces() != 4)
            {
                return face_key();
            }
            
            simplex_set st_f, cl_st_f;
            Complex::star(f, st_f);
            Complex::closure(st_f, cl_st_f);
            
            st_e.difference(cl_st_f);
#ifdef DEBUG
            assert(st_e.size_faces() == 1);
#endif
            return *st_e.faces_begin();
        }
        
        std::vector<edge_key> test_neighbour(const face_key& f, const node_key& a, const node_key& b, node_key& u, node_key& w, real& q_old, real& q_new)
        {
            edge_key e = Complex::get_edge(u,w);
            face_key g = get_neighbour(f, e);
            real q = Util::quality<real>(get_pos(a), get_pos(b), get_pos(w), get_pos(u));
            
            if(g.is_valid() && is_safe_editable(e))
            {
                node_key v = Complex::get_apex(g, e);
                real V_uv = Util::signed_volume<real>(get_pos(a), get_pos(b), get_pos(v), get_pos(u));
                real V_vw = Util::signed_volume<real>(get_pos(a), get_pos(b), get_pos(w), get_pos(v));
                real V_wu = Util::signed_volume<real>(get_pos(a), get_pos(b), get_pos(u), get_pos(w));
                
                if((V_uv > 0. && V_vw > 0.) || (V_vw > 0. && V_wu > 0.) || (V_wu > 0. && V_uv > 0.))
                {
                    q_old = Util::min(Util::quality<real>(get_pos(a), get_pos(u), get_pos(w), get_pos(v)),
                                     Util::quality<real>(get_pos(u), get_pos(v), get_pos(b), get_pos(w)));
                    
                    real q_uv_old, q_uv_new, q_vw_old, q_vw_new;
                    auto uv_edges = test_neighbour(g, a, b, u, v, q_uv_old, q_uv_new);
                    auto vw_edges = test_neighbour(g, a, b, v, w, q_vw_old, q_vw_new);
                    
                    q_old = Util::min(Util::min(q_old, q_uv_old), q_vw_old);
                    q_new = Util::min(q_uv_new, q_vw_new);
                    
                    if(q_new > q_old || q_new > q)
                    {
                        std::vector<edge_key> edges = {Complex::get_edge(f, g)};
                        edges.insert(edges.end(), uv_edges.begin(), uv_edges.end());
                        edges.insert(edges.end(), vw_edges.begin(), vw_edges.end());
                        return edges;
                    }
                }
            }
            q_old = INFINITY;
            q_new = q;
            return std::vector<edge_key>();
        }
        
        /*
         * Orient the nodes in a counter clockwise order seen from the node a.
         */
        void orient_cc(const node_key& a, std::vector<node_key>& nodes)
        {
            vec3 x = get_pos(a) - get_pos(nodes[0]);
            vec3 y = get_pos(nodes[1]) - get_pos(nodes[0]);
            vec3 z = get_pos(nodes[2]) - get_pos(nodes[0]);
            real val = Util::dot(x, Util::cross(y,z));
#ifdef DEBUG
            assert(val != 0.);
#endif
            if(val > 0.)
            {
                node_key t = nodes[0];
                nodes[0] = nodes[2];
                nodes[2] = t;
            }
        }
        
        /**
         * Attempt to remove the faces sandwiched between the apices of f using multi-face removal. The face f is used as a starting point.
         */
        bool topological_face_removal(const face_key& f)
        {
            auto apices = Complex::get_apices(f);
            auto nodes = Complex::get_nodes(f);
            orient_cc(apices[0], nodes);
            
            real q_01_new, q_01_old, q_12_new, q_12_old, q_20_new, q_20_old;
            auto e01 = test_neighbour(f, apices[0], apices[1], nodes[0], nodes[1], q_01_old, q_01_new);
            auto e12 = test_neighbour(f, apices[0], apices[1], nodes[1], nodes[2], q_12_old, q_12_new);
            auto e20 = test_neighbour(f, apices[0], apices[1], nodes[2], nodes[0], q_20_old, q_20_new);
            
            real q_old = Util::min(Util::min(Util::min(min_quality(f), q_01_old), q_12_old), q_20_old);
            real q_new = Util::min(Util::min(q_01_new, q_12_new), q_20_new);
            
            if(q_new > q_old)
            {
                Complex::flip_23(f);
                for(auto &e : e01)
                {
                    Complex::flip_32(e);
                }
                for(auto &e : e12)
                {
                    Complex::flip_32(e);
                }
                for(auto &e : e20)
                {
                    Complex::flip_32(e);
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
            is_mesh::SimplexSet<face_key> fids = Complex::get_faces(Complex::get_tets(apex1)) & Complex::get_faces(Complex::get_tets(apex2));
            for(auto f : fids)
            {
                if(is_safe_editable(f))
                {
                    auto nodes = Complex::get_nodes(f);
                    orient_cc(apex2, nodes);
                    
                    vec3 ray = get_pos(apex2) - get_pos(apex1);
                    real t = Util::intersection_ray_triangle<real>(get_pos(apex1), ray, get_pos(nodes[0]), get_pos(nodes[1]), get_pos(nodes[2]));
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
            for (auto tit = Complex::tetrahedra_begin(); tit != Complex::tetrahedra_end(); tit++)
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
                if (Complex::exists(t) && quality(t) < MIN_TET_QUALITY)
                {
                    simplex_set cl_t;
                    Complex::closure(t, cl_t);
                    
                    for (auto fit = cl_t.faces_begin(); fit != cl_t.faces_end(); fit++)
                    {
                        if (is_safe_editable(*fit))
                        {
                            auto apices = Complex::get_apices(*fit);
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
            
            Complex::garbage_collect();
        }
        
        ////////////////
        // THICKENING //
        ////////////////
        
        /**
         * Splits all interface edges with a volume greater than MAX_EDGE_LENGTH by inserting a vertex.
         */
        void thickening_interface()
        {
            std::vector<edge_key> edges;
            for (auto eit = Complex::edges_begin(); eit != Complex::edges_end(); eit++)
            {
                if (is_unsafe_editable(eit.key()) && is_interface(eit.key()) && length(eit.key()) > MAX_LENGTH)
                {
                    edges.push_back(eit.key());
                }
            }
            int i = 0, j = 0;
            for(auto &e : edges)
            {
                if (is_unsafe_editable(e) && is_interface(e) && length(e) > MAX_LENGTH)
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
            std::vector<tet_key> tetrahedra;
            for (auto tit = Complex::tetrahedra_begin(); tit != Complex::tetrahedra_end(); tit++)
            {
                if (volume(tit.key()) > MAX_VOLUME)
                {
                    tetrahedra.push_back(tit.key());
                }
            }
            int i = 0, j = 0;
            for(auto &t : tetrahedra)
            {
                if (Complex::exists(t) && volume(t) > MAX_VOLUME)
                {
                    if(split(t) != Complex::NULL_NODE)
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
        /**
         * Collapses all edges which is shorter than MIN_EDGE_LENGTH. However, the collapse is only performed if it does not alter the interface and it improves the minimum quality of the tetrahedra in the star of the edge.
         */
        void thinning()
        {
            std::vector<tet_key> tetrahedra;
            for (auto tit = Complex::tetrahedra_begin(); tit != Complex::tetrahedra_end(); tit++)
            {
                if (volume(tit.key()) < MIN_VOLUME)
                {
                    tetrahedra.push_back(tit.key());
                }
            }
            int i = 0, j = 0;
            for(auto &t : tetrahedra)
            {
                if (Complex::exists(t) && volume(t) < MIN_VOLUME)
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
            for (auto eit = Complex::edges_begin(); eit != Complex::edges_end(); eit++)
            {
                if (quality(eit.key()) < DEG_EDGE_QUALITY)
                {
                    edges.push_back(eit.key());
                }
            }
            int i = 0, j = 0;
            for(auto e : edges)
            {
                if(Complex::exists(e) && quality(e) < DEG_EDGE_QUALITY)
                {
                    if(collapse(e, false))
                    {
                        i++;
                    }
                    j++;
                }
            }
            std::cout << "Removed " << i <<"/"<< j << " degenerate edges" << std::endl;
            Complex::garbage_collect();
        }
        
        void remove_degenerate_faces()
        {
            std::list<face_key> faces;
            
            for (auto fit = Complex::faces_begin(); fit != Complex::faces_end(); fit++)
            {
                if(quality(fit.key()) < DEG_FACE_QUALITY)
                {
                    faces.push_back(fit.key());
                }
            }
            
            int i = 0, j = 0;
            for (auto &f : faces)
            {
                if (Complex::exists(f) && quality(f) < DEG_FACE_QUALITY)
                {
                    if(collapse(f, false))
                    {
                        i++;
                    }
                    else {
                        simplex_set cl_f;
                        Complex::closure(f, cl_f);
                        edge_key e = longest_edge(cl_f);
                        split(e);
                    }
                    j++;
                }
            }
            std::cout << "Removed " << i <<"/"<< j << " degenerate faces" << std::endl;
            Complex::garbage_collect();
        }
        
        void remove_degenerate_tets()
        {
            std::vector<tet_key> tets;
            
            for (auto tit = Complex::tetrahedra_begin(); tit != Complex::tetrahedra_end(); tit++)
            {
                if (quality(tit.key()) < DEG_TET_QUALITY)
                {
                    tets.push_back(tit.key());
                }
            }
            int i = 0, j = 0;
            for (auto &t : tets)
            {
                if (Complex::exists(t) && quality(t) < DEG_TET_QUALITY)
                {
                    if(collapse(t, false))
                    {
                        i++;
                    }
                    else {
                        simplex_set cl_t;
                        Complex::closure(t, cl_t);
                        edge_key e = longest_edge(cl_t);
                        split(e);
                    }
                    j++;
                }
            }
            std::cout << "Removed " << i <<"/"<< j << " degenerate tets" << std::endl;
            Complex::garbage_collect();
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
            for (auto eit = Complex::edges_begin(); eit != Complex::edges_end(); eit++)
            {
                if (quality(eit.key()) < MIN_EDGE_QUALITY)
                {
                    edges.push_back(eit.key());
                }
            }
            int i = 0, j = 0;
            for(auto e : edges)
            {
                if(Complex::exists(e) && quality(e) < MIN_EDGE_QUALITY)
                {
                    if(collapse(e) != Complex::NULL_NODE)
                    {
                        i++;
                    }
                    j++;
                }
            }
            std::cout << "Removed " << i <<"/"<< j << " low quality edges" << std::endl;
            Complex::garbage_collect();
        }
        
        /**
         * Attempt to remove the cap f by splitting the longest edge and collapsing it with cap's apex.
         */
        bool remove_cap(const face_key& f)
        {
            // Find longest edge
            simplex_set cl_f;
            Complex::closure(f, cl_f);
            edge_key e = longest_edge(cl_f);
            
            // Find apex
            simplex_set cl_e;
            Complex::closure(e, cl_e);
            cl_f.difference(cl_e);
            node_key apex = *cl_f.nodes_begin();
            
            // Find the projected position of the apex
            auto verts = get_pos(e);
            vec3 p = Util::project(get_pos(apex), verts[0], verts[1]);
            
            // Split longest edge
            node_key n = split(e);
            set_pos(n, p);
            
            // Collapse new edge
            edge_key e_rem = Complex::get_edge(apex, n);
            return collapse(e_rem) != Complex::NULL_NODE;
        }
        
        /**
         * Attempt to remove the cap f by splitting the longest edge and collapsing it with cap's apex.
         */
        bool remove_needle(const face_key& f)
        {
            // Find shortest edge
            simplex_set cl_f;
            Complex::closure(f, cl_f);
            edge_key e = shortest_edge(cl_f);
            
            // Remove edge
            return collapse(e) != Complex::NULL_NODE;
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
            
            for (auto fit = Complex::faces_begin(); fit != Complex::faces_end(); fit++)
            {
                if(quality(fit.key()) < MIN_FACE_QUALITY)
                {
                    faces.push_back(fit.key());
                }
            }
            
            int i = 0, j = 0;
            for (auto &f : faces)
            {
                if (Complex::exists(f) && quality(f) < MIN_FACE_QUALITY)
                {
                    if(remove_face(f))
                    {
                        i++;
                    }
                    j++;
                }
            }
            std::cout << "Removed " << i <<"/"<< j << " low quality faces" << std::endl;
            Complex::garbage_collect();
        }
        
        /**
         * Remove a degenerate tetrahedron of a type "sliver" by splitting the two longest edges
         * and collapsing the newly created vertices together. Return true if successful.
         */
        bool remove_sliver(const tet_key & t)
        {
            simplex_set cl_t;
            Complex::closure(t, cl_t);
            edge_key e1 = longest_edge(cl_t);
            cl_t.erase(e1);
            edge_key e2 = longest_edge(cl_t);
            
            node_key n1 = split(e1);
            node_key n2 = split(e2);
            
            edge_key e = Complex::get_edge(n1, n2);
            return collapse(e) != Complex::NULL_NODE;
        }
        
        /**
         * Remove a degenerate tetrahedron of a type "cap" by splitting the face opposite cap's apex and collapsing cap's apex with the newly created vertex.
         * Return true if successful.
         */
        bool remove_cap(const tet_key & t)
        {
            // Find the largest face
            simplex_set cl_t;
            Complex::closure(t, cl_t);
            face_key f = largest_face(cl_t);
            
            // Find the apex
            node_key apex = Complex::get_apex(t, f);
            
            // Project the apex
            auto verts = get_pos(f);
            vec3 p = Util::project(get_pos(apex), verts);
            
            // Split the face
            node_key n = split(f);
            set_pos(n, p);
            
            // Collapse edge
            edge_key e = Complex::get_edge(n, apex);
            return collapse(e) != Complex::NULL_NODE;
        }
        
        /**
         * Remove a degenerate tetrahedron of a type "wedge" or "needle" by collapsing the shortest edge.
         * Return true if successful.
         */
        bool remove_wedge(const tet_key & t)
        {
            simplex_set cl_t;
            Complex::closure(t, cl_t);
            while(cl_t.size_edges() > 2)
            {
                edge_key e = shortest_edge(cl_t);
                if(collapse(e) != Complex::NULL_NODE)
                {
                    return true;
                }
                cl_t.erase(e);
            }
            return false;
            
            //        simplex_set cl_t;
            //        Complex::closure(t, cl_t);
            //        edge_key e1 = longest_edge(cl_t);
            //        cl_t.erase(e1);
            //        edge_key e2 = longest_edge(cl_t);
            //
            //        node_key n1 = split(e1);
            //        node_key n2 = split(e2);
            //
            //        edge_key e = Complex::get_edge(n1, n2);
            //        return collapse(e) != Complex::NULL_NODE;
        }
        
        /**
         * Remove a tetrahedron of a type "needle" by splitting the tetrahedron.
         * Return true if successful.
         */
        bool remove_needle(const tet_key & t)
        {
            simplex_set cl_t;
            Complex::closure(t, cl_t);
            while(cl_t.size_edges() > 1)
            {
                edge_key e = shortest_edge(cl_t);
                if(collapse(e) != Complex::NULL_NODE)
                {
                    return true;
                }
                cl_t.erase(e);
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
        bool remove_tet(const tet_key & t)
        {
            // Find the largest face
            simplex_set cl_t;
            Complex::closure(t, cl_t);
            face_key f = largest_face(cl_t);
            
            // Find the apex
            node_key apex = Complex::get_apex(t, f);
            
            // Project the apex
            auto verts = get_pos(f);
            vec3 proj_apex = Util::project(get_pos(apex), verts);
            
            // Find barycentric coordinates
            std::vector<real> barycentric_coords = Util::barycentric_coords<real>(proj_apex, verts[0], verts[1], verts[2]);
            
            if(barycentric_coords[0] > 0.2 && barycentric_coords[1] > 0.2 && barycentric_coords[2] > 0.2) // The tetrahedron is a cap
            {
                return remove_cap(t);
            }
            else if(barycentric_coords[0] < -0.2 || barycentric_coords[1] < -0.2 || barycentric_coords[2] < -0.2) // The tetrahedron is a sliver
            {
                return remove_sliver(t);
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
                return remove_needle(t);
            }
            else if(close == 1) // The tetrahedron is a wedge
            {
                return remove_wedge(t);
            }
            return false;
        }
        
        /**
         * Attempt to remove tetrahedra with quality lower than MIN_TET_QUALITY.
         */
        void remove_tets()
        {
            std::vector<tet_key> tets;
            
            for (auto tit = Complex::tetrahedra_begin(); tit != Complex::tetrahedra_end(); tit++)
            {
                if (quality(tit.key()) < MIN_TET_QUALITY)
                {
                    tets.push_back(tit.key());
                }
            }
            int i = 0, j=0;
            for (auto &tet : tets)
            {
                if (Complex::exists(tet) && quality(tet) < MIN_TET_QUALITY)
                {
                    if(remove_tet(tet))
                    {
                        i++;
                    }
                    j++;
                }
            }
            std::cout << "Removed " << i <<"/"<< j << " low quality tets" << std::endl;
            Complex::garbage_collect();
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
            is_mesh::SimplexSet<tet_key> tids = Complex::get_tets(nid);
            real q_old = min_quality(tids);

            vec3 old_pos = get_pos(nid);
            vec3 avg_pos = get_barycenter(nid);
            set_pos(nid, old_pos + alpha * (avg_pos - old_pos));
            
            if (is_inverted(tids) || min_quality(tids) < q_old)
            {
                set_pos(nid, old_pos);
                return false;
            }
            return true;
        }
        
        void smooth()
        {
            int i = 0, j = 0;
            for (auto nit = Complex::nodes_begin(); nit != Complex::nodes_end(); nit++)
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
            std::cout << "Smooth." << std::endl;
            smooth();
            
            std::cout << "Topological removals." << std::endl;
            topological_edge_removal();
            topological_face_removal();
            
//            std::cout << "Low quality removal." << std::endl;
//            remove_tets();
//            remove_faces();
//            remove_edges();
            
            std::cout << "Degeneracy removal." << std::endl;
            remove_degenerate_tets();
            remove_degenerate_faces();
            remove_degenerate_edges();
            
            validity_check();
        }
        
        void resize_complex()
        {
            std::cout << "Thickening interface pass." << std::endl;
            thickening_interface();
            
//            std::cout << "Thinning interface pass." << std::endl;
//            thinning_interface();
            
//            std::cout << "Thickening pass." << std::endl;
//            thickening();
            
//            std::cout << "Thinning pass." << std::endl;
//            thinning();
            
            validity_check();
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
            int missing;
            int step = 0;
            do {
                std::cout << "Move vertices step " << step << std::endl;
                missing = 0;
                int movable = 0;
                for (auto nit = Complex::nodes_begin(); nit != Complex::nodes_end(); nit++)
                {
                    if (is_movable(nit.key()))
                    {
                        if(!move_vertex(nit.key()))
                        {
                            if(step == 9)
                            {
                                print(nit.key());
                            }
                            missing++;
                        }
                        movable++;
                    }
                }
                std::cout << "Vertices missing to be moved: " << missing <<"/" << movable << std::endl;
                validity_check();
                
                fix_complex();
                
                ++step;
            } while (missing > 0 && step < num_steps);
            
            resize_complex();
            
            Complex::garbage_collect();
            update_attributes();
        }
        
    private:
        
        /**
         * Tries moving the node n to the new position new_pos. Returns true if it succeeds.
         */
        bool move_vertex(const node_key & n)
        {
            vec3 pos = get_pos(n);
            vec3 destination = get_dest(n);
            real l = Util::length(destination - pos);
            
            if (l < 1e-4*AVG_EDGE_LENGTH) // The vertex is not moved
            {
                return true;
            }
            
            real t = intersection_with_link(n, destination);
            
            l = Util::max(Util::min(l*t - l*MIN_DEFORMATION, l), 0.);
            set_pos(n, pos + l*Util::normalize(destination - pos));
            
            if (Util::length(destination - get_pos(n)) < EPSILON)
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
            auto fids = Complex::get_faces(Complex::get_tets(n)) - Complex::get_faces(n);
            for(auto f : fids)
            {
                auto face_pos = get_pos(f);
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
        bool is_flippable(const edge_key & e)
        {
            if(!is_interface(e) && !is_boundary(e))
            {
                return false;
            }
            simplex_set st_e;
            Complex::star(e, st_e);
            std::vector<face_key> faces;
            for(auto fit = st_e.faces_begin(); fit != st_e.faces_end(); fit++)
            {
                if (is_interface(*fit) || is_boundary(*fit))
                {
                    faces.push_back(*fit);
                }
            }
            if(faces.size() > 2)
            {
                return false;
            }
#ifdef DEBUG
            assert(faces.size() == 2);
#endif
            
            real angle = cos_dihedral_angle(faces[0], faces[1]);
            if(angle > FLIP_EDGE_INTERFACE_FLATNESS)
            {
                return true;
            }
            return false;
        }
        
        /**
         * Flips the edge e (which is a special case of the edge remove operation in the embedding mesh).
         * Relabels the tetrahedra accordingly so that the interface mesh geometry does not change.
         */
        bool flip(const edge_key & e)
        {
            simplex_set st_e;
            Complex::star(e, st_e);
            std::vector<face_key> faces;
            for(auto fit = st_e.faces_begin(); fit != st_e.faces_end(); fit++)
            {
                if (is_interface(*fit) || is_boundary(*fit))
                {
                    faces.push_back(*fit);
                }
            }
            if(faces.size() > 2)
            {
                return false;
            }
#ifdef DEBUG
            assert(faces.size() == 2);
#endif
            
            node_key n = Complex::flip_44(faces[0], faces[1]);
            return n != Complex::NULL_NODE;
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
            auto verts = get_pos(tid);
            vec3 p = Util::barycenter(verts[0], verts[1], verts[2], verts[3]);
            
            node_key n = Complex::split(tid);
            set_pos(n, p);
            set_destination(n, p);
            return n;
        }
        
        /**
         * Split a face f and returns the new node which is positioned at the barycenter of the vertices of f.
         */
        node_key split(const face_key & f)
        {
            std::vector<node_key> nodes;
            Complex::get_nodes(f, nodes);
            vec3 p = Util::barycenter(get_pos(nodes[0]), get_pos(nodes[1]), get_pos(nodes[2]));
            vec3 p_new = p;
            if(is_interface(f))
            {
                p_new = Util::barycenter(get_dest(nodes[0]), get_dest(nodes[1]), get_dest(nodes[2]));
            }
            
            node_key n = Complex::split(f);
            set_pos(n, p);
            set_destination(n, p_new);
            
            return n;
        }
        
        /**
         * Split an edge e and returns the new node which is placed at the middle of e.
         */
        node_key split(const edge_key & e)
        {
            auto verts = get_pos(e);
            vec3 p = Util::barycenter(verts[0], verts[1]);
            vec3 p_new = p;
            if(is_interface(e))
            {
                auto dests = get_dest(e);
                p_new = Util::barycenter(dests[0], dests[1]);
            }
            
            node_key n = Complex::split(e);
            set_pos(n, p);
            set_destination(n, p_new);
            return n;
        }
        
        ///////////////
        // COLLAPSES //
        ///////////////
    private:
        /**
         * Collapses the edge e and moves the resulting node to the position p. Returns the new node if successful, otherwise NULL_NODE.
         */
        node_key collapse(edge_key& e, const vec3& p, const vec3& p_new)
        {
            node_key n_new = Complex::collapse(e);
            
            if (Complex::exists(n_new) && n_new.is_valid())
            {
                set_pos(n_new, p);
                set_destination(n_new, p_new);
            }
            return n_new;
        }
        
        /**
         * Returns true if the collapse of the edge e does not result in any inverted tetrahedra.
         * The merged nodes are assumed moved to p after the collapse.
         */
        bool precond_collapse(const edge_key& e, const vec3& p)
        {
            auto nodes = Complex::get_nodes(e);
            
            simplex_set lk_n0, lk_n1, lk_e;
            Complex::link(nodes[0], lk_n0);
            Complex::link(nodes[1], lk_n1);
            
            simplex_set st_e, cl_st_e;
            Complex::star(e, st_e);
            Complex::closure(st_e, cl_st_e);
            lk_n0.add(lk_n1);
            lk_n0.difference(cl_st_e);
            
            if(will_invert(nodes[0], p, lk_n0) || will_invert(nodes[1], p, lk_n0))
            {
                return false;
            }
            return true;
        }
        
        /**
         * Returns the minimum quality of neighbouring tetrahedra if the edge e is collapsed and the resulting node is moved to p_new.
         */
        real min_quality(const edge_key& eid, const vec3& p_new)
        {
            is_mesh::SimplexSet<node_key> nids = Complex::get_nodes(eid);
            is_mesh::SimplexSet<face_key> fids = Complex::get_faces(Complex::get_tets(nids)) - Complex::get_faces(nids);
            return min_quality(fids, p_new);
        }
        
        /**
         * Collapses the edge e and places the new node at the most optimal position of the position of either end node or their barycenter.
         * If the parameter safe is true, the method if the nodes of edge e are editable, i.e. not a part of the interface, and will therefore not change the interface.
         * If non of the nodes are editable or precond_collapse returns false, the method returns NULL_NODE.
         */
        node_key collapse(edge_key& e, bool safe = true)
        {
            if (!Complex::exists(e) || !e.is_valid())
            {
                return node_key();
            }
            auto nodes = Complex::get_nodes(e);
            
            vec3 p_opt, p_new_opt;
            real q_max = 0.;
            
            if (!is_boundary(nodes[0]) && !is_boundary(nodes[1]) && (!safe || (!is_interface(nodes[0]) && !is_interface(nodes[1]))))
            {
                vec3 p = Util::barycenter(get_pos(nodes[0]), get_pos(nodes[1]));
                real q = min_quality(e, p);
                if (precond_collapse(e, p) && q > q_max)
                {
                    p_new_opt = Util::barycenter(get_dest(nodes[0]), get_dest(nodes[1]));
                    p_opt = p;
                    q_max = q;
                }
            }
            
            if (!is_boundary(nodes[0]) && (!safe || !is_interface(nodes[0])))
            {
                vec3 p = get_pos(nodes[1]);
                real q = min_quality(e, p);
                
                if (precond_collapse(e, p) && q > q_max)
                {
                    p_new_opt = get_dest(nodes[1]);
                    p_opt = p;
                    q_max = q;
                }
            }
            
            if (!is_boundary(nodes[1]) && (!safe || !is_interface(nodes[1])))
            {
                vec3 p = get_pos(nodes[0]);
                real q = min_quality(e, p);
                
                if (precond_collapse(e, p) && q > q_max)
                {
                    p_new_opt = get_dest(nodes[0]);
                    p_opt = p;
                    q_max = q;
                }
            }
            real q = Util::min(min_quality(nodes[0]), min_quality(nodes[1]));
            if((!safe && q_max > EPSILON) || q_max > Util::min(q, MIN_TET_QUALITY))
            {
                return collapse(e, p_opt, p_new_opt);
            }
            return node_key();
        }
        
        bool collapse(const face_key& f, bool safe = true)
        {
            simplex_set cl_f;
            Complex::closure(f, cl_f);
            while(cl_f.size_edges() > 0)
            {
                edge_key e = shortest_edge(cl_f);
                node_key nid = collapse(e, safe);
                if(nid.is_valid())
                {
                    return true;
                }
                cl_f.erase(e);
            }
            return false;
        }
        
        bool collapse(const tet_key& t, bool safe = true)
        {
            simplex_set cl_t;
            Complex::closure(t, cl_t);
            while(cl_t.size_edges() > 0)
            {
                edge_key e = shortest_edge(cl_t);
                node_key nid = collapse(e, safe);
                if(nid.is_valid())
                {
                    return true;
                }
                cl_t.erase(e);
            }
            return false;
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
            Complex::orient_face(fid);
            auto pos = get_pos(fid);
            return Util::normal_direction(pos[0], pos[1], pos[2]);
        }
        
        /**
         Returns the normal to interface node n.
         */
        vec3 get_normal(const node_key& nid)
        {
            vec3 result(0.);
            for (auto f : Complex::get_faces(nid))
            {
                if (is_interface(f))
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
         * Calculates the average position of the neighbouring nodes to node n.
         * If interface is true, the average position is only calculated among the neighbouring nodes which are interface.
         */
        vec3 get_barycenter(const node_key& nid, bool interface = false)
        {
            if(interface && !is_interface(nid))
            {
                return get_pos(nid);
            }
            
            is_mesh::SimplexSet<node_key> nids = Complex::get_nodes(Complex::get_tets(nid)) - nid;
            vec3 avg_pos(0.);
            int i = 0;
            for (auto n : nids)
            {
                if (!interface || is_interface(n))
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
        
        ///////////////////////
        // UTILITY FUNCTIONS //
        ///////////////////////
        
    public:
        
        real length(const edge_key& eid)
        {
            auto verts = get_pos(eid);
            return Util::length(verts[0] - verts[1]);
        }
        
        real length_destination(const edge_key& eid)
        {
            auto dests = get_dest(eid);
            return Util::length(dests[0] - dests[1]);
        }
        
        real area(const face_key& fid)
        {
            auto verts = get_pos(fid);
            return Util::area<real>(verts[0], verts[1], verts[2]);
        }
        
        real area_destination(const face_key& fid)
        {
            auto dests = get_dest(fid);
            return Util::area<real>(dests[0], dests[1], dests[2]);
        }
        
        real volume(const tet_key& tid)
        {
            auto verts = get_pos(tid);
            return Util::volume<real>(verts[0], verts[1], verts[2], verts[3]);
        }
        
        real volume_destination(const tet_key& tid)
        {
            auto dests = get_dest(tid);
            return Util::volume<real>(dests[0], dests[1], dests[2], dests[3]);
        }
        
        real quality(const tet_key& tid)
        {
            is_mesh::SimplexSet<node_key> nids = Complex::get_nodes(tid);
            return std::abs(Util::quality<real>(get_pos(nids[0]), get_pos(nids[1]), get_pos(nids[2]), get_pos(nids[3])));
        }
        
        real min_angle(const face_key& f)
        {
            auto verts = get_pos(f);
            return Util::min_angle<real>(verts[0], verts[1], verts[2]);
        }
        
        real max_angle(const face_key& f)
        {
            auto verts = get_pos(f);
            return Util::max_angle<real>(verts[0], verts[1], verts[2]);
        }
        
        real quality(const face_key& fid)
        {
            auto verts = get_pos(fid);
            auto angles = Util::cos_angles<real>(verts[0], verts[1], verts[2]);
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
        face_key largest_face(simplex_set& set)
        {
            real max_a = -INFINITY;
            face_key max_f;
            for(auto f = set.faces_begin(); f != set.faces_end(); f++)
            {
                real a = area(*f);
                if(a > max_a)
                {
                    max_a = a;
                    max_f = *f;
                }
            }
            return max_f;
        }
        
        /**
         * Returns the shortest edge in the simplex set.
         */
        edge_key shortest_edge(simplex_set& set)
        {
            real min_l = INFINITY;
            edge_key min_e;
            for(auto e = set.edges_begin(); e != set.edges_end(); e++)
            {
                real l = length(*e);
                if(l < min_l)
                {
                    min_l = l;
                    min_e = *e;
                }
            }
            return min_e;
        }
        
        /**
         * Returns the shortest edge in the simplex set.
         */
        edge_key longest_edge(simplex_set& set)
        {
            real max_l = -INFINITY;
            edge_key max_e;
            for(auto e = set.edges_begin(); e != set.edges_end(); e++)
            {
                real l = length(*e);
                if(l > max_l)
                {
                    max_l = l;
                    max_e = *e;
                }
            }
            return max_e;
        }
        /**
         * Returns the minimum quality of the tetrahedra in simplex set s.
         */
        real min_quality(simplex_set& set)
        {
            real q_min = INFINITY;
            for (auto tit = set.tetrahedra_begin(); tit != set.tetrahedra_end(); tit++)
            {
                q_min = Util::min(quality(*tit), q_min);
            }
            return q_min;
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
         * Returns the minimum quality of the tetrahedra spanned by the vertices of the faces in s and pos.
         */
        real min_quality(const is_mesh::SimplexSet<face_key>& fids, const vec3& pos)
        {
            real min_q = INFINITY;
            for(auto f : fids)
            {
                auto verts = get_pos(f);
                min_q = Util::min(min_q, std::abs(Util::quality<real>(verts[0], verts[1], verts[2], pos)));
            }
            return min_q;
        }
        
        /**
         * Returns the minimum quality among tetrahedra in the star of the node nid.
         */
        real min_quality(const node_key& nid)
        {
            return min_quality(Complex::get_tets(nid));
        }
        
        /**
         * Returns the minimum quality among tetrahedra in the star of the edge eid.
         */
        real min_quality(const edge_key& eid)
        {
            return min_quality(Complex::get_tets(eid));
        }
        
        /**
         * Returns minimum quality among the tetrahedra adjacent to face fid.
         */
        real min_quality(const face_key& fid)
        {
            return min_quality(Complex::get_tets(fid));
        }
        
    private:
        
        /**
         * Check if the sequence of vertices in polygon is consistent with positive orientation of tetrahedra in the mesh
         * with respect to the ordered pair of vertices in vv. If not, reverse the order of vertices in polygon.
         */
        void check_consistency(const std::vector<vec3> & vv, std::vector<node_key> & polygon)
        {
            unsigned int n = static_cast<unsigned int>(polygon.size());
            std::vector<vec3> vp(n);
            
            for (unsigned int i = 0; i < n; ++i)
            {
                vp[i] = get_pos(polygon[i]);
            }
            
            real sum = 0;
            for (unsigned int i = 0; i < n; ++i)
            {
                sum += Util::signed_volume<real>(vv[0], vv[1], vp[i], vp[(i+1)%n]);
            }
            
            if (sum < 0.)
            {
                for (unsigned int i = 0; i < n/2; ++i)
                {
                    node_key temp = polygon[i];
                    polygon[i] = polygon[n-1-i];
                    polygon[n-1-i] = temp;
                }
            }
        }
        
        /// Check whether they are any inverted tetrahedra in the simplicial complex.
        bool simplicial_complex_criterion_check()
        {
            for (auto tit = Complex::tetrahedra_begin(); tit != Complex::tetrahedra_end(); tit++)
            {
                if (Complex::is_inverted(tit.key()))
                {
                    return false;
                }
            }
            return true;
        }
        
        /**
         * Returns whether any of the tetrahedra in the simplex set is inverted.
         */
        bool inverted(simplex_set& set)
        {
            for (auto tit = set.tetrahedra_begin(); tit != set.tetrahedra_end(); tit++)
            {
                if (Complex::is_inverted(*tit))
                {
                    return true;
                }
            }
            return false;
        }
        
        /**
         * Returns whether any of the tetrahedra in the simplex set is inverted.
         */
        bool is_inverted(const is_mesh::SimplexSet<tet_key>& tids)
        {
            for (auto t : tids)
            {
                if (Complex::is_inverted(t))
                {
                    return true;
                }
            }
            return false;
        }
        
        /**
         * Returns whether any of the tetrahedra in the simplex set will invert if node n is moved to p_new.
         */
        bool will_invert(const node_key& n, const vec3 p_new, simplex_set& set)
        {
            vec3 p = get_pos(n);
            for(auto fit = set.faces_begin(); fit != set.faces_end(); fit++)
            {
                auto verts = get_pos(*fit);
                real vol1 = Util::signed_volume<real>(verts[0], verts[1], verts[2], p);
                real vol2 = Util::signed_volume<real>(verts[0], verts[1], verts[2], p_new);
                if(Util::sign(vol1) !=  Util::sign(vol2))
                {
                    return true;
                }
            }
            return false;
        }
        
        void validity_check()
        {
            bool valid = simplicial_complex_criterion_check();
            assert(valid);
            for (auto tit = Complex::tetrahedra_begin(); tit != Complex::tetrahedra_end(); tit++)
            {
                valid = valid & Complex::exists(tit.key());
            }
            assert(valid);
            
            for (auto fit = Complex::faces_begin(); fit != Complex::faces_end(); fit++)
            {
                valid = valid & Complex::exists(fit.key());
            }
            assert(valid);
            
            for (auto nit = Complex::nodes_begin(); nit != Complex::nodes_end(); nit++)
            {
                valid = valid & Complex::exists(nit.key());
                valid = valid & !(is_interface(nit.key()) && is_boundary(nit.key())); // Check that the interface has not reached the boundary
            }
            assert(valid);
            
            for (auto eit = Complex::edges_begin(); eit != Complex::edges_end(); eit++)
            {
                valid = valid & Complex::exists(eit.key());
                simplex_set st_e;
                Complex::star(eit.key(), st_e);
                int boundary = 0;
                int interface = 0;
                for (auto fit = st_e.faces_begin(); fit != st_e.faces_end(); fit++) {
                    if(is_boundary(*fit))
                    {
                        boundary++;
                    }
                    if(is_interface(*fit))
                    {
                        interface++;
                    }
                }
                valid = valid & ((is_interface(eit.key()) && interface >= 2) || (!is_interface(eit.key()) && interface == 0)); // Check that the interface is not corrupted
                assert(valid);
                valid = valid & ((is_boundary(eit.key()) && boundary == 2) || (!is_boundary(eit.key()) && boundary == 0)); // Check that the boundary is not corrupted
                assert(valid);
            }
            assert(valid);
        }
        
        /**
         * Sort the vertices from set according to their connectivity and returned the sorted vertices in sorted_vertices.
         */
        void sort_vertices(simplex_set& set, std::vector<node_key>& sorted_vertices)
        {
            sorted_vertices = std::vector<node_key>();
            sorted_vertices.push_back(*(set.nodes_begin()));
            
            std::map<edge_key, bool> edge_used;
            for(auto eit = set.edges_begin(); eit != set.edges_end(); eit++)
            {
                edge_used[*eit] = false;
            }
            
            while(sorted_vertices.size() < set.size_nodes())
            {
                for (auto eit = set.edges_begin(); eit != set.edges_end(); eit++)
                {
                    if(!edge_used[*eit])
                    {
                        auto nodes = Complex::get_nodes(*eit);
                        
                        if (nodes[0] == sorted_vertices.back())
                        {
                            sorted_vertices.push_back(nodes[1]);
                            edge_used[*eit] = true;
                            break;
                        }
                        else if (nodes[1] == sorted_vertices.back())
                        {
                            sorted_vertices.push_back(nodes[0]);
                            edge_used[*eit] = true;
                            break;
                        }
                    }
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
            for (auto nit = Complex::nodes_begin(); nit != Complex::nodes_end(); nit++)
            {
                if (nit->is_interface())
                {
                    verts.push_back(get_pos(nit.key()));
                    vert_index[nit.key()] = static_cast<int>(verts.size());
                }
            }
            
            // Extract faces
            for (auto fit = Complex::faces_begin(); fit != Complex::faces_end(); fit++)
            {
                if (fit->is_interface())
                {
                    Complex::orient_face(fit.key());
                    auto nodes = Complex::get_nodes(fit.key());
                    
                    indices.push_back(vert_index[nodes[0]]);
                    indices.push_back(vert_index[nodes[1]]);
                    indices.push_back(vert_index[nodes[2]]);
                }
            }
        } // extract_interface
        
        void extract_tet_mesh(std::vector<vec3>& points, std::vector< std::vector<int> >& tets)
        {
            Complex::garbage_collect();
            
            std::map<node_key, int> indices;
            int counter = 0;
            for (auto nit = Complex::nodes_begin(); nit != Complex::nodes_end(); nit++)
            {
                indices[nit.key()] = counter;
                points.push_back(get_pos(nit.key()));
                ++counter;
            }
            
            for (auto tit = Complex::tetrahedra_begin(); tit != Complex::tetrahedra_end(); tit++)
            {
                std::vector<int> tet;
                auto nodes = Complex::get_nodes(tit.key());
                
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
            auto nodes = Complex::get_nodes(f1);
            auto temp = Complex::get_nodes(f2);
            nodes.insert(nodes.end(), temp.begin(), temp.end());
            
            std::vector<vec3> verts, apices;
            for(int i = 0; i < nodes.size(); i++)
            {
                bool found = false;
                for (int j = 0; j < nodes.size(); j++)
                {
                    if(i != j && nodes[i] == nodes[j])
                    {
                        if(i < j)
                        {
                            verts.push_back(get_pos(nodes[i]));
                        }
                        found = true;
                    }
                }
                if(!found)
                {
                    apices.push_back(get_pos(nodes[i]));
                }
            }
            
            return Util::cos_dihedral_angle<real>(verts[0], verts[1], apices[0], apices[1]);
        }
        
        /**
         * Returns the dihedral angle between face f1 and face f2.
         */
        real dihedral_angle(const face_key& f1, const face_key& f2)
        {
            return acos(cos_dihedral_angle(f1, f2));
        }
        
        std::vector<real> cos_dihedral_angles(const tet_key& t)
        {
            auto verts = get_pos(t);
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
            
            for (auto tit = Complex::tetrahedra_begin(); tit != Complex::tetrahedra_end(); tit++)
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
            
            for (auto tit = Complex::tetrahedra_begin(); tit != Complex::tetrahedra_end(); tit++)
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
            for (auto tit = Complex::tetrahedra_begin(); tit != Complex::tetrahedra_end(); tit++)
            {
                min_q = Util::min(min_q, quality(tit.key()));
            }
            return min_q;
        }
        
        /// Counts the total number of nodes and the number of nodes on the interface(s).
        void count_nodes(int & total, int & object)
        {
            total = 0, object = 0;
            for (auto nit = Complex::nodes_begin(); nit != Complex::nodes_end(); nit++)
            {
                total++;
                if (nit->is_interface())
                    object++;
            }
        }
        
        /// Counts the total number of edges and the number of edges on the interface(s).
        void count_edges(int & total, int & object)
        {
            total = 0, object = 0;
            for (auto eit = Complex::edges_begin(); eit != Complex::edges_end(); eit++)
            {
                total++;
                if (eit->is_interface())
                    object++;
            }
        }
        
        /// Counts the total number of faces and the number of faces on the interface(s).
        void count_faces(int & total, int & object)
        {
            total = 0, object = 0;
            for (auto fit = Complex::faces_begin(); fit != Complex::faces_end(); fit++)
            {
                total++;
                if (fit->is_interface())
                    object++;
            }
        }
        
        /// Counts the total number of tetrahedra and the number of tetrahedra in the object(s).
        void count_tetrahedra(int & total, int & object)
        {
            total = 0, object = 0;
            for (auto tit = Complex::tetrahedra_begin(); tit != Complex::tetrahedra_end(); tit++)
            {
                total++;
                if (get_label(tit.key()) != 0)
                    object++;
            }
        }
        
    public:
        void test_split_collapse()
        {
            is_mesh::SimplexSet<edge_key> eids;
            for (auto eit = Complex::edges_begin(); eit != Complex::edges_end(); eit++)
            {
                auto neighbours = Complex::get_boundary(Complex::get_co_boundary(eit.key()));
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
                auto nids = Complex::get_boundary(e);
                auto nid = Complex::split(e);
                auto new_eid = (Complex::get_co_boundary(nids) & Complex::get_co_boundary(nid)) - e;
                assert(new_eid.size() == 1);
                new_eids += new_eid[0];
                verts.push_back(get_pos(nids[1]));
                j++;
                if(j%100 == 0)
                {
                    std::cout << ".";
                }
            }
            std::cout << " DONE" << std::endl;
            
            std::cout << "Collapse test # = " << new_eids.size();
            j = 0;
            for (int i = 0; i < new_eids.size(); i++) {
                assert(Complex::exists(new_eids[i]));
                auto nid = Complex::collapse(new_eids[i]);
                assert(nid.is_valid());
                set_pos(nid, verts[i]);
                for(auto t : Complex::get_tets(nid))
                {
                    if(Complex::is_inverted(t))
                    {
                        Complex::get(t).invert_orientation();
                    }
                }
                j++;
                if(j%100 == 0)
                {
                    std::cout << ".";
                }
            }
            std::cout << " DONE" << std::endl;
            Complex::garbage_collect();
            Complex::validity_check();
            validity_check();
        }
        
        void test_flip23_flip32()
        {
            is_mesh::SimplexSet<face_key> fids;
            for (auto fit = Complex::faces_begin(); fit != Complex::faces_end(); fit++)
            {
                if(is_safe_editable(fit.key()))
                {
                    auto nids = Complex::get_boundary(Complex::get_boundary(fit.key()));
                    nids += Complex::get_boundary(Complex::get_boundary(Complex::get_boundary(Complex::get_co_boundary(fit.key()))));
                    real t = Util::intersection_ray_triangle<real>(get_pos(nids[3]), get_pos(nids[4]) - get_pos(nids[3]), get_pos(nids[0]), get_pos(nids[1]), get_pos(nids[2]));
                    
                    auto neighbours = Complex::get_boundary(Complex::get_co_boundary(fit.key()));
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
                assert(Complex::exists(f));
                auto new_eid = Complex::flip_23(f);
                assert(new_eid.is_valid());
                new_eids += new_eid;
                i++;
                if(i%1000 == 0)
                {
                    std::cout << ".";
                }
            }
            std::cout << " DONE" << std::endl;
            
            i=0;
            std::cout << "Flip 3-2 test # = " << new_eids.size();
            for (auto e : new_eids) {
                auto new_fid = Complex::flip_32(e);
                assert(new_fid.is_valid());
                i++;
                if(i%1000 == 0)
                {
                    std::cout << ".";
                }
            }
            std::cout << " DONE" << std::endl;
            Complex::garbage_collect();
            Complex::validity_check();
            validity_check();
        }
        
        void test_flip44()
        {
            is_mesh::SimplexSet<edge_key> eids;
            for (auto eit = Complex::edges_begin(); eit != Complex::edges_end(); eit++)
            {
                if(is_safe_editable(eit.key()) && Complex::get_co_boundary(eit.key()).size() == 4)
                {
                    auto neighbours = Complex::get_boundary(Complex::get_boundary(Complex::get_co_boundary(Complex::get_co_boundary(eit.key()))));
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
            }
            
            std::cout << "Flip 4-4 test # = " << eids.size();
            is_mesh::SimplexSet<face_key> flip_fids;
            int i = 0;
            for (auto e : eids) {
                assert(Complex::exists(e));
                auto fids = Complex::get_co_boundary(e);
                assert(fids.size() == 4);
                auto fid = fids - Complex::get_boundary(Complex::get_co_boundary(fids[0]));
                assert(fid.size() == 1);
                flip_fids += fid;
                flip_fids += fids[0];
                Complex::flip_44(fids[0], fid[0]);
                i++;
                if(i%1000 == 0)
                {
                    std::cout << ".";
                }
            }
            std::cout << " DONE" << std::endl;
            
            i=0;
            std::cout << "Flip 4-4 test # = " << eids.size();
            for (int j = 0; j < flip_fids.size(); j+=2)
            {
                Complex::flip_44(flip_fids[j], flip_fids[j+1]);
                i++;
                if(i%1000 == 0)
                {
                    std::cout << ".";
                }
            }
            std::cout << " DONE" << std::endl;
            Complex::garbage_collect();
            Complex::validity_check();
            validity_check();
        }
        
        void test_flip22()
        {
            is_mesh::SimplexSet<edge_key> eids;
            for (auto eit = Complex::edges_begin(); eit != Complex::edges_end(); eit++)
            {
                if(is_boundary(eit.key()) && Complex::get_co_boundary(eit.key()).size() == 3)
                {
                    auto neighbours = Complex::get_boundary(Complex::get_boundary(Complex::get_co_boundary(Complex::get_co_boundary(eit.key()))));
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
            }
            
            std::cout << "Flip 2-2 test # = " << eids.size();
            int i = 0;
            for (auto e : eids) {
                assert(Complex::exists(e));
                auto fids = Complex::get_co_boundary(e);
                assert(fids.size() == 3);
                is_mesh::SimplexSet<face_key> flip_fids;
                for(auto f : fids)
                {
                    if(is_boundary(f))
                    {
                        flip_fids += f;
                    }
                }
                
                assert(flip_fids.size() == 2);
                Complex::flip_22(flip_fids[0], flip_fids[1]);
                i++;
                if(i%1000 == 0)
                {
                    std::cout << ".";
                }
            }
            std::cout << " DONE" << std::endl;
            Complex::garbage_collect();
            Complex::validity_check();
            validity_check();
        }
        
    };
    
}
