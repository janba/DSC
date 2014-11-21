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
#include "util.h"

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

        DEPRECATED
        typedef is_mesh::NodeKey      node_key;
        DEPRECATED
        typedef is_mesh::EdgeKey      edge_key;
        DEPRECATED
        typedef is_mesh::FaceKey      face_key;
        DEPRECATED
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

    public:

        virtual void set_avg_edge_length(double avg_edge_length = 0.);

        void set_parameters(parameters pars_);

        void set_design_domain(is_mesh::Geometry *geometry);

        virtual void set_labels(const is_mesh::Geometry& geometry, int label);

    private:

        void print(const is_mesh::NodeKey& n);
        
        /////////////////////////
        // ATTRIBUTE FUNCTIONS //
        /////////////////////////
        
    protected:

        virtual bool is_unsafe_editable(const is_mesh::NodeKey& nid);

        virtual bool is_unsafe_editable(const is_mesh::EdgeKey& eid);

        virtual bool is_unsafe_editable(const is_mesh::FaceKey& fid);

        virtual bool is_unsafe_editable(const is_mesh::TetrahedronKey& tid);

        virtual bool is_safe_editable(const is_mesh::NodeKey& nid);

        virtual bool is_safe_editable(const is_mesh::EdgeKey& eid);

        virtual bool is_safe_editable(const is_mesh::FaceKey& fid);

        virtual bool is_safe_editable(const is_mesh::TetrahedronKey& tid);

    public:
        virtual bool is_movable(const is_mesh::NodeKey& nid);

    protected:
        void set_pos(const is_mesh::NodeKey& nid, const vec3& p);
        
    public:
        virtual void set_destination(const is_mesh::NodeKey& nid, const vec3& dest);

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

        double build_table(const is_mesh::EdgeKey& e, const is_mesh::SimplexSet<is_mesh::NodeKey>& polygon, std::vector<std::vector<int>>& K);


        is_mesh::NodeKey get_next(const is_mesh::NodeKey& nid, is_mesh::SimplexSet<is_mesh::EdgeKey>& eids);

        is_mesh::SimplexSet<is_mesh::NodeKey> get_polygon(is_mesh::SimplexSet<is_mesh::EdgeKey>& eids);
        
        std::vector<is_mesh::SimplexSet<is_mesh::NodeKey>> get_polygons(const is_mesh::EdgeKey& eid);

        void flip_23_recursively(const is_mesh::SimplexSet<is_mesh::NodeKey>& polygon, const is_mesh::NodeKey& n1, const is_mesh::NodeKey& n2, std::vector<std::vector<int>>& K, int i, int j);

        void topological_edge_removal(const is_mesh::SimplexSet<is_mesh::NodeKey>& polygon, const is_mesh::NodeKey& n1, const is_mesh::NodeKey& n2, std::vector<std::vector<int>>& K);

        bool topological_edge_removal(const is_mesh::EdgeKey& eid);

        void topological_boundary_edge_removal(const is_mesh::SimplexSet<is_mesh::NodeKey>& polygon1, const is_mesh::SimplexSet<is_mesh::NodeKey>& polygon2, const is_mesh::EdgeKey& eid, std::vector<std::vector<int>>& K1, std::vector<std::vector<int>>& K2);

        bool topological_boundary_edge_removal(const is_mesh::EdgeKey& eid);

        void topological_edge_removal();
        
        //////////////////////////////
        // TOPOLOGICAL FACE REMOVAL //
        //////////////////////////////

        is_mesh::SimplexSet<is_mesh::EdgeKey> test_neighbour(const is_mesh::FaceKey& f, const is_mesh::NodeKey& a, const is_mesh::NodeKey& b, const is_mesh::NodeKey& u, const is_mesh::NodeKey& w, double& q_old, double& q_new);

        bool topological_face_removal(const is_mesh::FaceKey& f);

        bool topological_face_removal(const is_mesh::NodeKey& apex1, const is_mesh::NodeKey& apex2);

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

        bool remove_cap(const is_mesh::FaceKey& fid);

        bool remove_needle(const is_mesh::FaceKey& fid);

        bool remove_face(const is_mesh::FaceKey& f);

        void remove_faces();

        bool remove_sliver(const is_mesh::TetrahedronKey& tid);

        bool remove_cap(const is_mesh::TetrahedronKey& tid);

        bool remove_wedge(const is_mesh::TetrahedronKey& tid);

        bool remove_needle(const is_mesh::TetrahedronKey& tid);

        bool remove_tet(const is_mesh::TetrahedronKey& tid);

        void remove_tets();
        
        ///////////////
        // SMOOTHING //
        ///////////////
    private:
        bool smart_laplacian(const is_mesh::NodeKey& nid, double alpha = 1.);

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

        bool move_vertex(const is_mesh::NodeKey & n);
        
        
    public:
        double intersection_with_link(const is_mesh::NodeKey & n, const vec3& destination);
        
        ///////////
        // FLIPS //
        ///////////
    private:

        bool is_flat(const is_mesh::SimplexSet<is_mesh::FaceKey>& fids);

        bool is_flippable(const is_mesh::EdgeKey & eid);

        bool precond_flip_edge(const is_mesh::EdgeKey& eid, const is_mesh::FaceKey& f1, const is_mesh::FaceKey& f2);
        
        ////////////
        // SPLITS //
        ////////////
    public:
        is_mesh::NodeKey split(const is_mesh::TetrahedronKey& tid);

        is_mesh::NodeKey split(const is_mesh::FaceKey& fid);

        is_mesh::NodeKey split(const is_mesh::EdgeKey& eid);

        ///////////////
        // COLLAPSES //
        ///////////////
    private:
        bool is_collapsable(const is_mesh::EdgeKey& eid, const is_mesh::NodeKey& nid, bool safe);

        bool collapse(const is_mesh::EdgeKey& eid, bool safe = true);

        bool collapse(is_mesh::SimplexSet<is_mesh::EdgeKey>& eids, bool safe);

        bool collapse(const is_mesh::FaceKey& fid, bool safe = true);

        bool collapse(const is_mesh::TetrahedronKey& tid, bool safe = true);
        
        //////////////////////
        // GETTER FUNCTIONS //
        //////////////////////
    public:

        std::vector<vec3> get_interface_face_positions();

        vec3 get_normal(const is_mesh::FaceKey& fid);

        vec3 get_normal(const is_mesh::NodeKey& nid);

        vec3 get_barycenter(const is_mesh::SimplexSet<is_mesh::NodeKey>& nids, bool interface = false);

        vec3 get_barycenter(const is_mesh::NodeKey& nid, bool interface = false);
        
        ///////////////////////
        // UTILITY FUNCTIONS //
        ///////////////////////
        
    public:

        double length(const is_mesh::EdgeKey& eid);

        double length_destination(const is_mesh::EdgeKey& eid);

        double area(const is_mesh::FaceKey& fid);

        double area_destination(const is_mesh::FaceKey& fid);

        double volume(const is_mesh::TetrahedronKey& tid);

        double volume_destination(const is_mesh::TetrahedronKey& tid);

        double volume_destination(const is_mesh::SimplexSet<is_mesh::NodeKey>& nids);

        double signed_volume_destination(const is_mesh::SimplexSet<is_mesh::NodeKey>& nids);

        vec3 barycenter(const is_mesh::TetrahedronKey& tid);

        vec3 barycenter_destination(const is_mesh::TetrahedronKey& tid);

        double quality(const is_mesh::TetrahedronKey& tid);

        double min_angle(const is_mesh::FaceKey& fid);

        double max_angle(const is_mesh::FaceKey& fid);

        double quality(const is_mesh::FaceKey& fid);

        double quality(const is_mesh::EdgeKey& eid);

        is_mesh::FaceKey largest_face(const is_mesh::SimplexSet<is_mesh::FaceKey>& fids);

        is_mesh::EdgeKey shortest_edge(const is_mesh::SimplexSet<is_mesh::EdgeKey>& eids);

        is_mesh::EdgeKey longest_edge(const is_mesh::SimplexSet<is_mesh::EdgeKey>& eids);

        double min_quality(const is_mesh::SimplexSet<is_mesh::TetrahedronKey>& tids);

        double min_quality(const is_mesh::SimplexSet<is_mesh::FaceKey>& fids, const vec3& pos);

        double min_quality(const is_mesh::SimplexSet<is_mesh::FaceKey>& fids, const vec3& pos_old, const vec3& pos_new);

        void min_quality(const is_mesh::SimplexSet<is_mesh::FaceKey>& fids, const vec3& pos_old, const vec3& pos_new, double& min_q_old, double& min_q_new);
        
        
    private:

        void check_consistency(const is_mesh::SimplexSet<is_mesh::NodeKey>& nids, is_mesh::SimplexSet<is_mesh::NodeKey>& polygon);
        
        ////////////////////////
        // DOCUMENT FUNCTIONS //
        ////////////////////////
    public:

        double compute_avg_edge_length();

        double cos_dihedral_angle(const is_mesh::FaceKey& f1, const is_mesh::FaceKey& f2);

        double dihedral_angle(const is_mesh::FaceKey& f1, const is_mesh::FaceKey& f2);

        std::vector<double> cos_dihedral_angles(const is_mesh::TetrahedronKey& tid);

        double min_cos_dihedral_angle(const is_mesh::TetrahedronKey& t);

        double min_dihedral_angle(const is_mesh::TetrahedronKey& t);

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
