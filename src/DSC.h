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
#include "attributes.h"
#include "geometry.h"
#include "mesh_io.h"

struct parameters {
    
    // Thresholds on the quality of edges
    double DEG_EDGE_QUALITY;
    double MIN_EDGE_QUALITY;
    
    // Thresholds on the quality of faces.
    double DEG_FACE_QUALITY;
    double MIN_FACE_QUALITY;
    
    // Thresholds on the quality of tetrahedra.
    double DEG_TET_QUALITY;
    double MIN_TET_QUALITY;
    
    // Thresholds on the length of edges.
    double MIN_LENGTH;
    double MAX_LENGTH;
    
    // Thresholds on the area of faces.
    double MIN_AREA;
    double MAX_AREA;
    
    // Thresholds on the volume of tetrahedra.
    double MIN_VOLUME;
    double MAX_VOLUME;
};

namespace DSC {
    

    class DeformableSimplicialComplex : public is_mesh::ISMesh
    {
    public:
        
        typedef is_mesh::NodeKey      node_key;
        typedef is_mesh::EdgeKey      edge_key;
        typedef is_mesh::FaceKey      face_key;
        typedef is_mesh::TetrahedronKey       tet_key;
        
    protected:
        is_mesh::MultipleGeometry design_domain;
        
        // Input parameter
        double AVG_LENGTH;
        double AVG_AREA;
        double AVG_VOLUME;
        
        // Should be eliminated
        double FLIP_EDGE_INTERFACE_FLATNESS = 0.995;
        
        parameters pars;
        
        //////////////////////////
        // INITIALIZE FUNCTIONS //
        //////////////////////////
        
    public:
        
        /// SimplicialComplex constructor.
        DeformableSimplicialComplex(std::vector<vec3> & points, std::vector<int> & tets, const std::vector<int>& tet_labels);

        ~DeformableSimplicialComplex();
        /*
        using is_mesh::ISMesh::get;
        using is_mesh::ISMesh::get_label;

        using is_mesh::ISMesh::nodes_begin;
        using is_mesh::ISMesh::nodes_end;
        using is_mesh::ISMesh::edges_begin;
        using is_mesh::ISMesh::edges_end;
        using is_mesh::ISMesh::faces_begin;
        using is_mesh::ISMesh::faces_end;
        using is_mesh::ISMesh::tetrahedra_begin;
        using is_mesh::ISMesh::tetrahedra_end;

        using is_mesh::ISMesh::get_pos;

        using is_mesh::ISMesh::get_nodes;
        using is_mesh::ISMesh::get_edges;
        using is_mesh::ISMesh::get_faces;
        using is_mesh::ISMesh::get_tets;

        using is_mesh::ISMesh<tet_att>::get_edge;
        using is_mesh::ISMesh<tet_att>::get_face;

        using is_mesh::ISMesh<tet_att>::validity_check;

    protected:
        using is_mesh::ISMesh<tet_att>::set_label;

    private:

        using is_mesh::ISMesh<tet_att>::flip_22;
        using is_mesh::ISMesh<tet_att>::flip_23;
        using is_mesh::ISMesh<tet_att>::flip_32;
        using is_mesh::ISMesh<tet_att>::flip_44;

        using is_mesh::ISMesh<tet_att>::split;
        using is_mesh::ISMesh<tet_att>::collapse;

        using is_mesh::ISMesh<tet_att>::garbage_collect;
        using is_mesh::ISMesh<tet_att>::exists;
             */
    public:

        virtual void set_avg_edge_length(double avg_edge_length = 0.);

        void set_parameters(parameters pars_);

        void set_design_domain(is_mesh::Geometry *geometry);

        virtual void set_labels(const is_mesh::Geometry& geometry, int label);

    private:

        void print(const node_key& n);
        
        /////////////////////////
        // ATTRIBUTE FUNCTIONS //
        /////////////////////////
        
    protected:

        virtual bool is_unsafe_editable(const node_key& nid);

        virtual bool is_unsafe_editable(const edge_key& eid);

        virtual bool is_unsafe_editable(const face_key& fid);

        virtual bool is_unsafe_editable(const tet_key& tid);

        virtual bool is_safe_editable(const node_key& nid);

        virtual bool is_safe_editable(const edge_key& eid);

        virtual bool is_safe_editable(const face_key& fid);

        virtual bool is_safe_editable(const tet_key& tid);

    public:
        virtual bool is_movable(const node_key& nid);

    protected:
        void set_pos(const node_key& nid, const vec3& p);
        
    public:
        virtual void set_destination(const node_key& nid, const vec3& dest);

        /////////////
        // GETTERS //
        /////////////
    public:

        vec3 get_center() const;

        double get_min_tet_quality() const;

        double get_deg_tet_quality() const;

        double get_deg_face_quality() const;

        double get_min_face_quality() const;

        double get_avg_edge_length() const;

        const is_mesh::MultipleGeometry & get_design_domain() const;
        
        ////////////////////////
        // FIX MESH FUNCTIONS //
        ////////////////////////
    private:
        
        //////////////////////////////
        // TOPOLOGICAL EDGE REMOVAL //
        //////////////////////////////

        double build_table(const edge_key& e, const is_mesh::SimplexSet<node_key>& polygon, std::vector<std::vector<int>>& K);


        node_key get_next(const node_key& nid, is_mesh::SimplexSet<edge_key>& eids);

        is_mesh::SimplexSet<node_key> get_polygon(is_mesh::SimplexSet<edge_key>& eids);
        
        std::vector<is_mesh::SimplexSet<node_key>> get_polygons(const edge_key& eid);

        void flip_23_recursively(const is_mesh::SimplexSet<node_key>& polygon, const node_key& n1, const node_key& n2, std::vector<std::vector<int>>& K, int i, int j);

        void topological_edge_removal(const is_mesh::SimplexSet<node_key>& polygon, const node_key& n1, const node_key& n2, std::vector<std::vector<int>>& K);

        bool topological_edge_removal(const edge_key& eid);

        void topological_boundary_edge_removal(const is_mesh::SimplexSet<node_key>& polygon1, const is_mesh::SimplexSet<node_key>& polygon2, const edge_key& eid, std::vector<std::vector<int>>& K1, std::vector<std::vector<int>>& K2);

        bool topological_boundary_edge_removal(const edge_key& eid);

        void topological_edge_removal();
        
        //////////////////////////////
        // TOPOLOGICAL FACE REMOVAL //
        //////////////////////////////

        is_mesh::SimplexSet<edge_key> test_neighbour(const face_key& f, const node_key& a, const node_key& b, const node_key& u, const node_key& w, double& q_old, double& q_new);

        bool topological_face_removal(const face_key& f);

        bool topological_face_removal(const node_key& apex1, const node_key& apex2);

        void topological_face_removal();
        
        ////////////////
        // THICKENING //
        ////////////////

        void thickening_interface();

        void thickening();
        
        //////////////
        // THINNING //
        //////////////

        void thinning_interface();

        void thinning();
        
        /////////////////////////
        // REMOVE DEGENERACIES //
        /////////////////////////
        void remove_degenerate_edges();

        void remove_degenerate_faces();

        void remove_degenerate_tets();
        
        //////////////////////////////////
        // REMOVE LOW QUALITY SIMPLICES //
        //////////////////////////////////

        void remove_edges();

        bool remove_cap(const face_key& fid);

        bool remove_needle(const face_key& fid);

        bool remove_face(const face_key& f);

        void remove_faces();

        bool remove_sliver(const tet_key& tid);

        bool remove_cap(const tet_key& tid);

        bool remove_wedge(const tet_key& tid);

        bool remove_needle(const tet_key& tid);

        bool remove_tet(const tet_key& tid);

        void remove_tets();
        
        ///////////////
        // SMOOTHING //
        ///////////////
    private:
        bool smart_laplacian(const node_key& nid, double alpha = 1.);

        void smooth();
        
        ///////////////////
        // FIX FUNCTIONS //
        ///////////////////

        void fix_complex();

        void resize_complex();
        
        ////////////////////
        // MOVE FUNCTIONS //
        ////////////////////
    public:
        void deform(int num_steps = 10);
        
    private:

        bool move_vertex(const node_key & n);
        
        
    public:
        double intersection_with_link(const node_key & n, const vec3& destination);
        
        ///////////
        // FLIPS //
        ///////////
    private:

        bool is_flat(const is_mesh::SimplexSet<face_key>& fids);

        bool is_flippable(const edge_key & eid);

        bool precond_flip_edge(const edge_key& eid, const face_key& f1, const face_key& f2);
        
        ////////////
        // SPLITS //
        ////////////
    public:
        void split(const tet_key& tid);

        void split(const face_key& fid);

        void split(const edge_key& eid);
        
        ///////////////
        // COLLAPSES //
        ///////////////
    private:
        bool is_collapsable(const edge_key& eid, const node_key& nid, bool safe);


        bool collapse(const edge_key& eid, bool safe = true);

        bool collapse(is_mesh::SimplexSet<edge_key>& eids, bool safe);

        bool collapse(const face_key& fid, bool safe = true);

        bool collapse(const tet_key& tid, bool safe = true);
        
        //////////////////////
        // GETTER FUNCTIONS //
        //////////////////////
    public:

        std::vector<vec3> get_interface_face_positions();

        vec3 get_normal(const face_key& fid);

        vec3 get_normal(const node_key& nid);

        vec3 get_barycenter(const is_mesh::SimplexSet<node_key>& nids, bool interface = false);

        vec3 get_barycenter(const node_key& nid, bool interface = false);
        
        ///////////////////////
        // UTILITY FUNCTIONS //
        ///////////////////////
        
    public:

        double length(const edge_key& eid);

        double length_destination(const edge_key& eid);

        double area(const face_key& fid);

        double area_destination(const face_key& fid);

        double volume(const tet_key& tid);

        double volume_destination(const tet_key& tid);

        double volume_destination(const is_mesh::SimplexSet<node_key>& nids);

        double signed_volume_destination(const is_mesh::SimplexSet<node_key>& nids);

        vec3 barycenter(const tet_key& tid);

        vec3 barycenter_destination(const tet_key& tid);

        double quality(const tet_key& tid);

        double min_angle(const face_key& fid);

        double max_angle(const face_key& fid);

        double quality(const face_key& fid);

        double quality(const edge_key& eid);

        face_key largest_face(const is_mesh::SimplexSet<face_key>& fids);

        edge_key shortest_edge(const is_mesh::SimplexSet<edge_key>& eids);

        edge_key longest_edge(const is_mesh::SimplexSet<edge_key>& eids);

        double min_quality(const is_mesh::SimplexSet<tet_key>& tids);

        double min_quality(const is_mesh::SimplexSet<face_key>& fids, const vec3& pos);

        double min_quality(const is_mesh::SimplexSet<face_key>& fids, const vec3& pos_old, const vec3& pos_new);

        void min_quality(const is_mesh::SimplexSet<face_key>& fids, const vec3& pos_old, const vec3& pos_new, double& min_q_old, double& min_q_new);
        
        
    private:

        void check_consistency(const is_mesh::SimplexSet<node_key>& nids, is_mesh::SimplexSet<node_key>& polygon);
        
        ////////////////////////
        // DOCUMENT FUNCTIONS //
        ////////////////////////
    public:

        double compute_avg_edge_length();

        double cos_dihedral_angle(const face_key& f1, const face_key& f2);

        double dihedral_angle(const face_key& f1, const face_key& f2);

        std::vector<double> cos_dihedral_angles(const tet_key& tid);

        double min_cos_dihedral_angle(const tet_key& t);

        double min_dihedral_angle(const tet_key& t);

        void get_qualities(std::vector<int>& histogram, double& min_quality);

        void get_dihedral_angles(std::vector<int> & histogram, double & min_angle, double & max_angle);

        double min_quality();

        void count_nodes(int & total, int & object);

        void count_edges(int & total, int & object);

        void count_faces(int & total, int & object);

        void count_tetrahedra(int & total, int & object);
        
    public:
        void test_split_collapse();

        void test_flip23_flip32();

        void test_flip44();

        void test_flip22();
        
    };
    
}
