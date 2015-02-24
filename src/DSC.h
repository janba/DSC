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
#include "geometry.h"
#include "mesh_io.h"
#include "util.h"


struct parameters {
    
    // Thresholds on the quality of edges
    double DEG_EDGE_QUALITY = 0.1;
    double MIN_EDGE_QUALITY = 0.5;
    
    // Thresholds on the quality of faces.
    double DEG_FACE_QUALITY = 0.0005;
    double MIN_FACE_QUALITY = 0.015;
    
    // Thresholds on the quality of tetrahedra.
    double DEG_TET_QUALITY = 0.02;
    double MIN_TET_QUALITY = 0.3;
    
    // Thresholds on the length of edges relative to average edge length.
    double MIN_LENGTH = 0.0;
    double MAX_LENGTH = 2.;
    
    // Thresholds on the area of faces.
    double MIN_AREA = 0.2;
    double MAX_AREA = 5.;
    
    // Thresholds on the volume of tetrahedra.
    double MIN_VOLUME = 0.2;
    double MAX_VOLUME = INFINITY;


    parameters() {
    }

    parameters(double DEG_EDGE_QUALITY, double MIN_EDGE_QUALITY, double DEG_FACE_QUALITY, double MIN_FACE_QUALITY, double DEG_TET_QUALITY, double MIN_TET_QUALITY, double MIN_LENGTH, double MAX_LENGTH, double MIN_AREA, double MAX_AREA, double MIN_VOLUME, double MAX_VOLUME)
            : DEG_EDGE_QUALITY(DEG_EDGE_QUALITY),
              MIN_EDGE_QUALITY(MIN_EDGE_QUALITY),
              DEG_FACE_QUALITY(DEG_FACE_QUALITY),
              MIN_FACE_QUALITY(MIN_FACE_QUALITY),
              DEG_TET_QUALITY(DEG_TET_QUALITY),
              MIN_TET_QUALITY(MIN_TET_QUALITY),
              MIN_LENGTH(MIN_LENGTH),
              MAX_LENGTH(MAX_LENGTH),
              MIN_AREA(MIN_AREA),
              MAX_AREA(MAX_AREA),
              MIN_VOLUME(MIN_VOLUME),
              MAX_VOLUME(MAX_VOLUME) {
    }
};

namespace DSC {
    class DeformableSimplicialComplex
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
        std::shared_ptr<is_mesh::ISMesh> is_mesh_ptr;
        is_mesh::ISMesh& is_mesh;
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

        /// SimplicialComplex constructor.
        DeformableSimplicialComplex(std::shared_ptr<is_mesh::ISMesh> ismesh);

        virtual ~DeformableSimplicialComplex();

    public:

        std::shared_ptr<is_mesh::ISMesh> get_shared_is_mesh();

        is_mesh::ISMesh & get_is_mesh();

        void set_avg_edge_length(double avg_edge_length = 0.);

        void set_parameters(parameters pars_);

        parameters get_parameters();

        DEPRECATED // Use add_design_domain
        void set_design_domain(is_mesh::Geometry *geometry) { add_design_domain(std::shared_ptr<is_mesh::Geometry>(geometry)); }

        void add_design_domain(std::shared_ptr<is_mesh::Geometry> geometry);

        // set sub domain elements which are modified
        void set_subdomain(std::shared_ptr<is_mesh::Geometry> subdomain);

        void clear_subdomain();

        std::shared_ptr<is_mesh::Geometry> get_subdomain();

        virtual void set_labels(const is_mesh::Geometry& geometry, int label);

    private:
        // For debugging!
        void print(const is_mesh::NodeKey& n);
        
        /////////////////////////
        // ATTRIBUTE FUNCTIONS //
        /////////////////////////
        
    protected:
        /**
        * key exists and is non boundary
        */
        virtual bool is_unsafe_editable(const is_mesh::NodeKey& nid);
        /**
        * key exists and is non boundary
        */
        virtual bool is_unsafe_editable(const is_mesh::EdgeKey& eid);
        /**
        * key exists and is non boundary
        */
        virtual bool is_unsafe_editable(const is_mesh::FaceKey& fid);
        /**
        * key exists
        */
        virtual bool is_unsafe_editable(const is_mesh::TetrahedronKey& tid);
        /**
        * key exists and is non interface and non boundary
        */
        virtual bool is_safe_editable(const is_mesh::NodeKey& nid);
        /**
        * key exists and is non interface and non boundary
        */
        virtual bool is_safe_editable(const is_mesh::EdgeKey& eid);
        /**
        * key exists and is non interface and non boundary
        */
        virtual bool is_safe_editable(const is_mesh::FaceKey& fid);
        /**
        * key exists
        */
        virtual bool is_safe_editable(const is_mesh::TetrahedronKey& tid);

    public:
        virtual bool is_movable(const is_mesh::NodeKey& nid);

    protected:
        /**
        * Sets the position of node n.
        */
        void set_pos(const is_mesh::NodeKey& nid, const vec3& p);

    public:
        inline static std::string header_version() { return "0.0.13"; };

        static std::string lib_version();


        /**
        * Sets the destination where the node n is moved to when deform() is called.
        */
        virtual void set_destination(const is_mesh::NodeKey& nid, const vec3& dest);

        /////////////
        // GETTERS //
        /////////////
    public:
        // Returns nodes.
        is_mesh::NodeIterator nodes() const;

        // Returns edges.
        is_mesh::EdgeIterator edges() const;

        // Returns faces.
        is_mesh::FaceIterator faces() const;

        // Returns tetrahedra.
        is_mesh::TetrahedronIterator tetrahedra() const;

        double get_min_tet_quality() const;

        double get_deg_tet_quality() const;

        double get_deg_face_quality() const;

        double get_min_face_quality() const;

        double get_avg_edge_length() const;

        // possible limited by subdomain
        const is_mesh::MultipleGeometry & get_design_domain() const;

        /////////////
        // GETTERS //
        /////////////

    private:
        
        //////////////////////////////
        // TOPOLOGICAL EDGE REMOVAL //
        //////////////////////////////
        /**
        * Build a table K for the dynamic programming method by Klincsek (see Shewchuk "Two Discrete Optimization Algorithms
        * for the Topological Improvement of Tetrahedral Meshes" article for details).
        */
        double build_table(const is_mesh::EdgeKey& e, const is_mesh::SimplexSet<is_mesh::NodeKey>& polygon, std::vector<std::vector<int>>& K);

        is_mesh::NodeKey get_next(const is_mesh::NodeKey& nid, is_mesh::SimplexSet<is_mesh::EdgeKey>& eids);

        is_mesh::SimplexSet<is_mesh::NodeKey> get_polygon(is_mesh::SimplexSet<is_mesh::EdgeKey>& eids);
        
        std::vector<is_mesh::SimplexSet<is_mesh::NodeKey>> get_polygons(const is_mesh::EdgeKey& eid);

        void flip_23_recursively(const is_mesh::SimplexSet<is_mesh::NodeKey>& polygon, const is_mesh::NodeKey& n1, const is_mesh::NodeKey& n2, std::vector<std::vector<int>>& K, int i, int j);

        void topological_edge_removal(const is_mesh::SimplexSet<is_mesh::NodeKey>& polygon, const is_mesh::NodeKey& n1, const is_mesh::NodeKey& n2, std::vector<std::vector<int>>& K);
        /**
        * Attempt to remove edge e by mesh reconnection using the dynamic programming method by Klincsek (see Shewchuk "Two Discrete Optimization Algorithms
        * for the Topological Improvement of Tetrahedral Meshes" article for details).
        */
        bool topological_edge_removal(const is_mesh::EdgeKey& eid);

        void topological_boundary_edge_removal(const is_mesh::SimplexSet<is_mesh::NodeKey>& polygon1, const is_mesh::SimplexSet<is_mesh::NodeKey>& polygon2, const is_mesh::EdgeKey& eid, std::vector<std::vector<int>>& K1, std::vector<std::vector<int>>& K2);

        bool topological_boundary_edge_removal(const is_mesh::EdgeKey& eid);

        void topological_edge_removal();
        
        //////////////////////////////
        // TOPOLOGICAL FACE REMOVAL //
        //////////////////////////////

        is_mesh::SimplexSet<is_mesh::EdgeKey> test_neighbour(const is_mesh::FaceKey& f, const is_mesh::NodeKey& a, const is_mesh::NodeKey& b, const is_mesh::NodeKey& u, const is_mesh::NodeKey& w, double& q_old, double& q_new);
        /**
        * Attempt to remove the faces sandwiched between the apices of f using multi-face removal. The face f is used as a starting point.
        */
        bool topological_face_removal(const is_mesh::FaceKey& f);
        /**
        * Attempt to remove the faces sandwiched between the nodes apex1 and apex2 using multi-face removal.
        * The face which intersects with the line segment |apex1 apex2| is used as a starting point.
        */
        bool topological_face_removal(const is_mesh::NodeKey& apex1, const is_mesh::NodeKey& apex2);
        /**
        * Improve tetrahedra quality by the topological operation (re-connection) edge removal. It do so only for tetrahedra of quality lower than MIN_TET_QUALITY.
        * For details see "Two Discrete Optimization Algorithms for the Topological Improvement of Tetrahedral Meshes" by Shewchuk.
        */
        void topological_face_removal();
        
        ////////////////
        // THICKENING //
        ////////////////
        /**
        * Splits all interface edges with a volume greater than MAX_EDGE_LENGTH by inserting a vertex.
        */
        void thickening_interface();
        /**
        * Splits all tetrahedra with a volume greater than MAX_TET_VOLUME by inserting a vertex.
        */
        void thickening();
        
        //////////////
        // THINNING //
        //////////////

        void thinning_interface();
        /**
        * Collapses all edges which is shorter than MIN_EDGE_LENGTH. However, the collapse is only performed if it does not alter the interface and it improves the minimum quality of the tetrahedra in the star of the edge.
        */
        void thinning();
        
        /////////////////////////
        // REMOVE DEGENERACIES //
        /////////////////////////
        /**
        * Attempt to remove edges with lower quality than DEG_EDGE_QUALITY by collapsing them.
        */
        void remove_degenerate_edges();

        void remove_degenerate_faces();

        void remove_degenerate_tets();
        
        //////////////////////////////////
        // REMOVE LOW QUALITY SIMPLICES //
        //////////////////////////////////
        /**
        * Attempt to remove edges with worse quality than MIN_EDGE_QUALITY by safely collapsing them.
        */
        void remove_low_quality_edges();
        /**
        * Attempt to remove the cap f by splitting the longest edge and collapsing it with cap's apex.
        */
        bool remove_cap(const is_mesh::FaceKey& fid);
        /**
        * Attempt to remove the cap f by splitting the longest edge and collapsing it with cap's apex.
        */
        bool remove_needle(const is_mesh::FaceKey& fid);
        /**
        * Attempt to remove the face f by first determining whether it's a cap or a needle.
        */
        bool remove_face(const is_mesh::FaceKey& f);
        /**
        * Attempts to remove degenerate faces (faces with a minimum angle smaller than MIN_ANGLE).
        */
        void remove_low_quality_faces();

        ///////////////
        // SMOOTHING //
        ///////////////
    private:
        /**
        * Performs Laplacian smoothing if it improves the minimum tetrahedron quality locally.
        */
        bool smart_laplacian(const is_mesh::NodeKey& nid, double alpha = 1.);

        /**
        * Performs smart_laplacian smoothing on all editable nodes
        */
        void smooth();
        
        ///////////////////
        // FIX FUNCTIONS //
        ///////////////////

        // Performs an laplacian smoothing of editable nodes
        // Optionally improve mesh quality using face/edge flipping and removing degenerates
        void fix_complex(bool optimizeMeshStructure = true);

        void resize_complex();
        
        ////////////////////
        // MOVE FUNCTIONS //
        ////////////////////
    public:
        /**
        * Moves all the vertices to their destination which can be set by the set_destination() function.
        */
        void deform(int num_steps = 10, bool optimizeMeshStructure = true);
        
    private:
        /**
        * Tries moving the node n to the new position new_pos. Returns true if it succeeds.
        */
        bool move_vertex(const is_mesh::NodeKey & n);
        
        
    public:
        /**
        * Returns the intersection point (= pos + t*(destination-pos)) with the link of the node n and
        * when moving the node n to the new position destination.
        */
        double intersection_with_link(const is_mesh::NodeKey & n, const vec3& destination);
        
        ///////////
        // FLIPS //
        ///////////
    private:
        /**
        * Returns whether the interface or boundary faces in fids are flat.
        */
        bool is_flat(const is_mesh::SimplexSet<is_mesh::FaceKey>& fids);
        /**
        * Returns whether it is possible to flip the edge e or not, i.e. whether the edge is a feature edge
        * (it is not a feature edge if its neighborhood is sufficiently flat).
        */
        bool is_flippable(const is_mesh::EdgeKey & eid);
        /**
        * Returns whether it is possible to flip the edge e or not.
        */
        bool precond_flip_edge(const is_mesh::EdgeKey& eid, const is_mesh::FaceKey& f1, const is_mesh::FaceKey& f2);
        
        ////////////
        // SPLITS //
        ////////////
    public:
        /**
        * Split a tetrahedron t and returns the new node which is positioned at the barycenter of the vertices of t.
        */
        is_mesh::NodeKey split(const is_mesh::TetrahedronKey& tid);
        /**
        * Split a face f and returns the new node which is positioned at the barycenter of the vertices of f.
        */
        is_mesh::NodeKey split(const is_mesh::FaceKey& fid);
        /**
        * Split an edge e and returns the new node which is placed at the middle of e.
        */
        is_mesh::NodeKey split(const is_mesh::EdgeKey& eid);

        ///////////////
        // COLLAPSES //
        ///////////////
    private:
        bool is_collapsable(const is_mesh::EdgeKey& eid, const is_mesh::NodeKey& nid, bool safe);
        /**
        * Collapses the edge e and places the new node at the most optimal position of the position of either end node or their barycenter.
        * If the parameter safe is true, the method if the nodes of edge e are editable, i.e. not a part of the interface, and will therefore not change the interface.
        * Returns whether the collapse was successful.
        */
        bool collapse(const is_mesh::EdgeKey& eid, bool safe = true);

        bool collapse(is_mesh::SimplexSet<is_mesh::EdgeKey>& eids, bool safe);

        bool collapse(const is_mesh::FaceKey& fid, bool safe = true);

        bool collapse(const is_mesh::TetrahedronKey& tid, bool safe = true);
        
        //////////////////////
        // GETTER FUNCTIONS //
        //////////////////////
    public:

        std::vector<vec3> get_interface_face_positions();
        /**
        Returns the normal to interface face fid.
        */
        vec3 get_normal(const is_mesh::FaceKey& fid);
        /**
        Returns the normal to interface node n.
        */
        vec3 get_normal(const is_mesh::NodeKey& nid);
        /**
        * Calculates the average position of the nodes in the simplex set nids.
        * If interface is true, the average position is only calculated among the nodes which are interface.
        */
        vec3 get_barycenter(const is_mesh::SimplexSet<is_mesh::NodeKey>& nids, bool interface = false);
        /**
        * Calculates the average position of the neighbouring nodes to node n.
        * If interface is true, the average position is only calculated among the neighbouring nodes which are interface.
        */
        vec3 get_barycenter(const is_mesh::NodeKey& nid, bool interface = false);
        
        ///////////////////////
        // UTILITY FUNCTIONS //
        ///////////////////////
        
    public:
        DEPRECATED // use Edge::length()
        double length(const is_mesh::EdgeKey& eid);

        DEPRECATED // use Edge::length_destination()
        double length_destination(const is_mesh::EdgeKey& eid);

        DEPRECATED // use Face::area()
        double area(const is_mesh::FaceKey& fid);

        DEPRECATED // use Face::area_destination()
        double area_destination(const is_mesh::FaceKey& fid);

        DEPRECATED // use Tetrahedron::volume()
        double volume(const is_mesh::TetrahedronKey& tid);

        DEPRECATED // use Tetrahedron::volume_destination()
        double volume_destination(const is_mesh::TetrahedronKey& tid);

        DEPRECATED // use ISMesh::volume_destination(nids)
        double volume_destination(const is_mesh::SimplexSet<is_mesh::NodeKey>& nids);

        DEPRECATED // use ISMesh::signed_volume_destination(nids)
        double signed_volume_destination(const is_mesh::SimplexSet<is_mesh::NodeKey>& nids);

        DEPRECATED // use Tetrahedron::barycenter()
        vec3 barycenter(const is_mesh::TetrahedronKey& tid);

        DEPRECATED // use Tetrahedron::barycenter_destination()
        vec3 barycenter_destination(const is_mesh::TetrahedronKey& tid);

        DEPRECATED // use Tetrahedron::quality()
        double quality(const is_mesh::TetrahedronKey& tid);

        DEPRECATED // use Face::min_angle()
        double min_angle(const is_mesh::FaceKey& fid);

        DEPRECATED // use Face::max_angle()
        double max_angle(const is_mesh::FaceKey& fid);

        DEPRECATED // use Face::quality()
        double quality(const is_mesh::FaceKey& fid);

        double quality(const is_mesh::EdgeKey& eid);

        /**
        * Returns the largest face in the simplex set.
        */
        is_mesh::FaceKey largest_face(const is_mesh::SimplexSet<is_mesh::FaceKey>& fids);
        /**
        * Returns the shortest edge in the simplex set.
        */
        is_mesh::EdgeKey shortest_edge(const is_mesh::SimplexSet<is_mesh::EdgeKey>& eids);
        /**
        * Returns the longest edge in the simplex set.
        */
        is_mesh::EdgeKey longest_edge(const is_mesh::SimplexSet<is_mesh::EdgeKey>& eids);
        /**
        * Returns the minimum quality of the tetrahedra in simplex set s.
        */
        double min_quality(const is_mesh::SimplexSet<is_mesh::TetrahedronKey>& tids);
        /**
        * Returns the minimum tetrahedral quality of a node with position pos. The faces in the link of the node should be passed in fids.
        */
        double min_quality(const is_mesh::SimplexSet<is_mesh::FaceKey>& fids, const vec3& pos);
        /**
        * Returns the new minimum tetrahedral quality when moving a node from old_pos to new_pos. The faces in the link of the node should be passed in fids.
        */
        double min_quality(const is_mesh::SimplexSet<is_mesh::FaceKey>& fids, const vec3& pos_old, const vec3& pos_new);
        /**
        * Returns the old (in min_q_old) and new (in min_q_new) minimum tetrahedral quality of when moving a node from old_pos to new_pos. The faces in the link of the node should be passed in fids.
        */
        void min_quality(const is_mesh::SimplexSet<is_mesh::FaceKey>& fids, const vec3& pos_old, const vec3& pos_new, double& min_q_old, double& min_q_new);
        
        
    private:
        /**
        * Check if the sequence of vertices in polygon is consistent with positive orientation of tetrahedra in the mesh
        * with respect to the ordered pair of vertices in vv. If not, reverse the order of vertices in polygon.
        */
        void check_consistency(const is_mesh::SimplexSet<is_mesh::NodeKey>& nids, is_mesh::SimplexSet<is_mesh::NodeKey>& polygon);
        
        ////////////////////////
        // DOCUMENT FUNCTIONS //
        ////////////////////////
    public:

        double compute_avg_edge_length();
        /**
        * Returns the cosine to the dihedral angle between face f1 and face f2.
        */
        double cos_dihedral_angle(const is_mesh::FaceKey& f1, const is_mesh::FaceKey& f2);
        /**
        * Returns the dihedral angle between face f1 and face f2.
        */
        double dihedral_angle(const is_mesh::FaceKey& f1, const is_mesh::FaceKey& f2);

        std::vector<double> cos_dihedral_angles(const is_mesh::TetrahedronKey& tid);
        /**
        * Returns the cosine of the minimum dihedral angle between the faces of tetrahedron t.
        */
        double min_cos_dihedral_angle(const is_mesh::TetrahedronKey& t);
        /**
        * Returns the minimum dihedral angle between the faces of tetrahedron t.
        */
        double min_dihedral_angle(const is_mesh::TetrahedronKey& t);

        void get_qualities(std::vector<int>& histogram, double& min_quality);

        /**
        * Calculates the dihedral angles in the SimplicialComplex and returns these in a histogram,
        * along with the minimum and maximum dihedral angles.
        */
        void get_dihedral_angles(std::vector<int> & histogram, double & min_angle, double & max_angle);

        double min_quality();
        /// Counts the total number of nodes and the number of nodes on the interface(s).
        void count_nodes(int & total, int & object);
        /// Counts the total number of edges and the number of edges on the interface(s).
        void count_edges(int & total, int & object);
        /// Counts the total number of faces and the number of faces on the interface(s).
        void count_faces(int & total, int & object);
        /// Counts the total number of tetrahedra and the number of tetrahedra in the object(s).
        void count_tetrahedra(int & total, int & object);
        
    public:
        void test_split_collapse();

        void test_flip23_flip32();

        void test_flip44();

        void test_flip22();
    };
    
}
