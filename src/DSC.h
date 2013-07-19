#ifndef DSC_H
#define DSC_H

#include "is_mesh_API.h"

#include "util.h"

#include "printing.h"

template <typename MT>
struct quality_struct
{
    typename MT::real_type quality;
    typename MT::real_type slope;
    typename MT::vector3_type grad_quality;
    
    bool operator<(const quality_struct& tqi) const
    {
        return quality < tqi.quality || (quality == tqi.quality && slope < tqi.slope);
    }
};
        

template <typename MT>
class DeformableSimplicialComplex : public ISMesh<MT>
{
    typedef typename MT::real_type      T;
    typedef typename MT::vector3_type   V;
    typedef typename MT::vector4_type   V4;
    
    typedef ISMesh<MT> Complex;
public:
    
    typedef typename Complex::node_key      node_key;
    typedef typename Complex::edge_key      edge_key;
    typedef typename Complex::face_key      face_key;
    typedef typename Complex::tet_key       tet_key;
    typedef typename Complex::simplex_set   simplex_set;
    
private:
    
    T DEG_ANGLE;
    T MIN_ANGLE;
    T COS_MIN_ANGLE;
    
    T DEG_TET_VOLUME;
    T MIN_TET_VOLUME;
    T MAX_TET_VOLUME;
    
    T DEG_EDGE_LENGTH;
    T MIN_EDGE_LENGTH;
    T MAX_EDGE_LENGTH;
    
    T DEG_TET_QUALITY;
    T MIN_TET_QUALITY;
    
    T AVG_EDGE_LENGTH;
    T MIN_DEFORMATION;
    
    //old:
    T MIN_DIHEDRAL_ANGLE;
    
    T FLIP_EDGE_TET_FLATNESS;
    
    
#define MAX_COS_FACE_ANGLE 0.996
#define MAX_COS 0.998
    
    int step_no;
    
    //////////////////////////
    // INITIALIZE FUNCTIONS //
    //////////////////////////
    
public:
    
    /// SimplicialComplex constructor.
    DeformableSimplicialComplex(std::vector<T> & points, std::vector<int> & tets, std::vector<int> & tet_labels):
        ISMesh<MT>(points, tets, tet_labels)
    {
        step_no = 0;
        
        T ie_min, ie_avg;
        calc_interface_edge_length(ie_min, ie_avg);
        AVG_EDGE_LENGTH = ie_avg;
        
        MIN_EDGE_LENGTH = 0.5 * AVG_EDGE_LENGTH;
        DEG_EDGE_LENGTH = 0.1 * AVG_EDGE_LENGTH;
        
        DEG_ANGLE = 5.*M_PI/180.;
        MIN_ANGLE = 10.*M_PI/180.;
        
        MIN_DEFORMATION = 0.25 * AVG_EDGE_LENGTH;
        
        MIN_DIHEDRAL_ANGLE = 5.*M_PI/180.;
        DEG_TET_QUALITY = 0.01;
        MIN_TET_QUALITY = 0.3;
        FLIP_EDGE_TET_FLATNESS = 0.995;
        
        T vol_avg = AVG_EDGE_LENGTH*AVG_EDGE_LENGTH*AVG_EDGE_LENGTH*sqrt(2.)/12.;
        MIN_TET_VOLUME = 0.5*vol_avg;
        MAX_TET_VOLUME = 10.*vol_avg;
    }
    
    DeformableSimplicialComplex() {}
    
private:
    
    /// Calculates the minimum and average length of edges on the interface.
    void calc_interface_edge_length(T& ie_min, T& ie_avg)
    {
        ie_min = INFINITY;
        ie_avg = 0.;
        int cnt = 0;
        for(auto eit = Complex::edges_begin(); eit != Complex::edges_end(); eit++)
        {
            if (eit->is_interface())
            {
                T l = length(eit.key());
                ie_min = std::min(l, ie_min);
                ie_avg += l;
                ++cnt;
            }
        }
        assert(cnt != 0);
        ie_avg /= cnt;
    }
    
    /////////////////////////
    // ATTRIBUTE FUNCTIONS //
    /////////////////////////
public:
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
    
    int get_label(const tet_key& t)
    {
        return Complex::get_label(t);
    }
    
    /**
     * Returns the position of node n.
     */
    V get_pos(const node_key& n)
    {
        V p = Complex::get(n).get_pos();
        assert(!MT::is_nan(p[0]) && !MT::is_nan(p[1]) && !MT::is_nan(p[2]));
        return p;
    }
    
    /// Returns the positions of the nodes of edge e in verts.
    void get_pos(const edge_key & e, std::vector<V>& verts)
    {
        verts = std::vector<V>(2);
        std::vector<node_key> nodes(2);
        Complex::get_nodes(e, nodes);
        for (int k = 0; k < 2; ++k)
        {
            verts[k] = get_pos(nodes[k]);
        }
    }
    
    /// Returns the positions of the nodes of face f in verts.
    void get_pos(const face_key & f, std::vector<V>& verts)
    {
        verts = std::vector<V>(3);
        Complex::orient_face(f);
        std::vector<node_key> nodes(3);
        Complex::get_nodes(f, nodes);
        for (int k = 0; k < 3; ++k)
        {
            verts[k] = get_pos(nodes[k]);
        }
    }
    
    /// Returns the positions of the nodes of tetrahedron t in verts.
    void get_pos(const tet_key& t, std::vector<V>& verts)
    {
        verts = std::vector<V>(4);
        std::vector<node_key> nodes(4);
        Complex::get_nodes(t, nodes);
        for (int k = 0; k < 4; ++k)
        {
            verts[k] = get_pos(nodes[k]);
        }
    }
    
protected:
    /**
     * Sets the position of node n.
     */
    void set_pos(const node_key& n, V p)
    {
        Complex::get(n).set_pos(p);
    }
    
public:
    V get_destination(const node_key& n)
    {
        return Complex::get(n).get_destination();
    }
    
    /**
     * Sets the destination where the node n is moved to when deform() is called.
     */
    void set_destination(const node_key& n, V p)
    {
        Complex::get(n).set_destination(p);
    }  
    
    /////////////
    // GETTERS //
    /////////////
public:
    
    V get_center() const
    {
        return V(0., 0., 0.);
    }
    
    T get_min_tet_quality() const
    {
        return MIN_TET_QUALITY;
    }
    
    T get_deg_tet_quality() const
    {
        return DEG_TET_QUALITY;
    }
    
    T get_min_dihedral_angle() const
    {
        return MIN_DIHEDRAL_ANGLE;
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
    T build_table(const edge_key& e, const std::vector<node_key>& polygon, std::vector<std::vector<int>>& K)
    {
        std::vector<V> verts;
        get_pos(e, verts);
        
        const int m = (int) polygon.size();
        
        auto Q = std::vector<std::vector<T>>(m-1, std::vector<T>(m, 0.));
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
                    T q2 = Util::quality<MT>(get_pos(polygon[i]), get_pos(polygon[k]), get_pos(polygon[j]), verts[1]);
                    T q1 = Util::quality<MT>(get_pos(polygon[k]), get_pos(polygon[i]), get_pos(polygon[j]), verts[0]);
                    T q = std::min(q1, q2);
                    if (k < j-1)
                    {
                        q = std::min(q, Q[k][j]);
                    }
                    if (k > i+1)
                    {
                        q = std::min(q, Q[i][k]);
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
        check_consistency(e, polygon);
        return polygon;
    }
    
    void flip_23_recursively(std::vector<node_key>& polygon, const node_key& n1, const node_key& n2, std::vector<std::vector<int>>& K, int i, int j)
    {
        if(j >= i+2)
        {
            int k = K[i][j];
            flip_23_recursively(polygon, n1, n2, K, i, k);
            flip_23_recursively(polygon, n1, n2, K, k, j);
            Complex::flip_23(Complex::get_face(n1, n2, polygon[k]));
        }
    }
    
    void edge_removal(std::vector<node_key>& polygon, const node_key& n1, const node_key& n2, std::vector<std::vector<int>>& K)
    {
        const int m = (int) polygon.size();
        int k = K[0][m-1];
        flip_23_recursively(polygon, n1, n2, K, 0, k);
        flip_23_recursively(polygon, n1, n2, K, k, m-1);
        Complex::flip_32(Complex::get_edge(n1, n2));
    }
    
    /**
     * Attempt to remove edge e by mesh reconnection using the dynamic programming method by Klincsek (see Shewchuk "Two Discrete Optimization Algorithms
     * for the Topological Improvement of Tetrahedral Meshes" article for details).
     */
    bool edge_removal(edge_key const & e)
    {
        simplex_set st_e;
        Complex::star(e, st_e);
        
        std::vector<node_key> polygon = get_polygon(e);
        std::vector<node_key> nodes;
        Complex::get_nodes(e, nodes);
        
        std::vector<std::vector<int>> K;
        T q_new = build_table(e, polygon, K);
        
        if (q_new > min_quality(e))
        {
            edge_removal(polygon, nodes[0], nodes[1], K);
            return true;
        }
        return false;
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
            return Complex::NULL_FACE;
        }
        
        simplex_set st_f, cl_st_f;
        Complex::star(f, st_f);
        Complex::closure(st_f, cl_st_f);
        
        st_e.difference(cl_st_f);
        assert(st_e.size_faces() == 1);
        return *st_e.faces_begin();
    }
    
    std::vector<edge_key> test_neighbour(const face_key& f, const node_key& a, const node_key& b, node_key& u, node_key& w, T& q_old, T& q_new)
    {
        edge_key e = Complex::get_edge(u,w);
        face_key g = get_neighbour(f, e);
        T q = Util::quality<MT>(get_pos(a), get_pos(b), get_pos(u), get_pos(w));
        
        if(g != Complex::NULL_FACE && !is_boundary(e) && !is_interface(e))
        {
            node_key v = Complex::get_apex(g, e);
            T V_uv = Util::signed_volume<MT>(get_pos(a), get_pos(b), get_pos(u), get_pos(v));
            T V_vw = Util::signed_volume<MT>(get_pos(a), get_pos(b), get_pos(v), get_pos(w));
            T V_wu = Util::signed_volume<MT>(get_pos(a), get_pos(b), get_pos(w), get_pos(u));
            
            if((V_uv > 0. && V_vw > 0.) || (V_vw > 0. && V_wu > 0.) || (V_wu > 0. && V_uv > 0.))
            {
                q_old = std::min(Util::quality<MT>(get_pos(a), get_pos(u), get_pos(v), get_pos(w)),
                                 Util::quality<MT>(get_pos(u), get_pos(v), get_pos(w), get_pos(b)));
                
                T q_uv_old, q_uv_new, q_vw_old, q_vw_new;
                auto uv_edges = test_neighbour(g, a, b, u, v, q_uv_old, q_uv_new);
                auto vw_edges = test_neighbour(g, a, b, v, w, q_vw_old, q_vw_new);
                
                q_old = std::min(std::min(q_old, q_uv_old), q_vw_old);
                q_new = std::min(q_uv_new, q_vw_new);
                
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
        V x = get_pos(a) - get_pos(nodes[0]);
        V y = get_pos(nodes[1]) - get_pos(nodes[0]);
        V z = get_pos(nodes[2]) - get_pos(nodes[0]);
        T val = MT::dot(x, cross(y,z));
        assert(val != 0.);
        if(val > 0.)
        {
            node_key t = nodes[0];
            nodes[0] = nodes[2];
            nodes[2] = t;
        }
    }
    
    bool remove_multi_face(const face_key& f)
    {
        std::vector<node_key> apices;
        Complex::get_apices(f, apices);
        std::vector<node_key> nodes;
        Complex::get_nodes(f, nodes);
        orient_cc(apices[0], nodes);
        
        T q_01_new, q_01_old, q_12_new, q_12_old, q_20_new, q_20_old;
        auto e01 = test_neighbour(f, apices[0], apices[1], nodes[0], nodes[1], q_01_old, q_01_new);
        auto e12 = test_neighbour(f, apices[0], apices[1], nodes[1], nodes[2], q_12_old, q_12_new);
        auto e20 = test_neighbour(f, apices[0], apices[1], nodes[2], nodes[0], q_20_old, q_20_new);
        
        T q_old = std::min(std::min(std::min(min_quality(f), q_01_old), q_12_old), q_20_old);
        T q_new = std::min(std::min(q_01_new, q_12_new), q_20_new);
        
        if(q_new > q_old)
        {
            node_key n = Complex::flip_23(f);
            assert(n != Complex::NULL_NODE);
            for(auto &e : e01)
            {
                n = Complex::flip_32(e);
                assert(n != Complex::NULL_NODE);
            }
            for(auto &e : e12)
            {
                n = Complex::flip_32(e);
                assert(n != Complex::NULL_NODE);
            }
            for(auto &e : e20)
            {
                n = Complex::flip_32(e);
                assert(n != Complex::NULL_NODE);
            }
            return true;
        }
        return false;
    }
    
    /**
    * Attempt to remove face f using multi-face removal.
    */
    bool multi_face_removal(const face_key& f)
    {
        std::vector<node_key> apices;
        Complex::get_apices(f, apices);
        
        simplex_set lk_n1, lk_n2;
        Complex::link(apices[0], lk_n1);
        Complex::link(apices[1], lk_n2);
        lk_n1.intersection(lk_n2);
        for(auto f = lk_n1.faces_begin(); f != lk_n1.faces_end(); f++)
        {
            std::vector<node_key> nodes;
            Complex::get_nodes(*f, nodes);
            orient_cc(apices[1], nodes);
            
            T t = Util::intersection_ray_triangle<MT>(get_pos(apices[0]), get_pos(apices[1]), get_pos(nodes[0]), get_pos(nodes[1]), get_pos(nodes[2]));
            if(0. < t && t < 1.)
            {
                if(remove_multi_face(*f))
                {
                    return true;
                }
            }
        }
        return false;
    }
    
    /**
     * The helper function for optimal_multi_face_remove().
     * TODO: Sanity check.
     */
    T omfr_quality_helper(int n, std::vector<face_key> & faces, std::map<face_key, T> & face_quality_values, std::vector<bool> & subset, std::vector<V> & vv)
    {
        T quality = 1e99;
        simplex_set multi_face, mf_boundary;
        
        for (int i = 0; i < n; ++i)
            if (subset[i] == false)
            {
                if (face_quality_values[faces[i]] < quality)
                    quality = face_quality_values[faces[i]];
            }
            else
                multi_face.insert(faces[i]);
        
        Complex::boundary_2manifold(multi_face, mf_boundary);
        T q = edge_remove_quality(mf_boundary, vv);
        
        if (q < quality) quality = q;
        return quality;
    }
    
    /**
     * Returns the faces in the set candidate which is connected to the face f. Connected means there is no interface simplices between f and the candidate face and all faces between f and the candidate face is in the candidate set.
     */
    void connected_component(const face_key& f, simplex_set& candidate, simplex_set& cc)
    {
        cc.insert(f);
        simplex_set cl_f, st_cl_f;
        Complex::closure(f, cl_f);
        Complex::star(cl_f, st_cl_f);
        
        for (auto fit = st_cl_f.faces_begin(); fit != st_cl_f.faces_end(); fit++)
        {
            if(Complex::exists(*fit) && candidate.contains(*fit) && !cc.contains(*fit) && !is_interface(*fit))
            {
                edge_key e = Complex::get_edge(*fit, f);
                if(e != Complex::NULL_EDGE && Complex::exists(e) && !is_interface(e))
                {
                    connected_component(*fit, candidate, cc);
                }
            }
        }
    }
    
    /**
     * Finds the biggest multi-face containing f, which can be removed.
     */
    void removable_multi_face(const face_key& f, simplex_set& multi_face)
    {
        std::vector<node_key> apices;
        Complex::get_apices(f, apices);
        
        simplex_set link0, link1;
        Complex::link(apices[0], link0);
        Complex::link(apices[1], link1);
        link0.intersection(link1);
        
        connected_component(f, link0, multi_face);
    }    
    
    /**
     * Attempt multi-face removal (a reconnection method), which finds the optimal subset (via dynamic programming)
     * of the multi-face to be swapped for an edge connecting the apices (in vv).
     * TODO: Sanity-check
     */
    void optimal_multi_face_removal(const face_key& f, simplex_set& new_simplices)
    {
        std::vector<node_key> apices;
        Complex::get_apices(f, apices);
        
        simplex_set multi_face;
        removable_multi_face(f, multi_face);
        
        if(multi_face.size_faces() < 2)
        {
            return;
        }
        
        simplex_set st_mf;
        Complex::star(multi_face, st_mf);
        T q_old = min_quality(st_mf);
        
        std::vector<V> verts(2);
        verts[0] = get_pos(apices[0]);
        verts[1] = get_pos(apices[1]);
        
        unsigned int n = static_cast<unsigned int>(multi_face.size_faces());
        std::map<face_key, double> face_quality_values;
        std::vector<face_key> faces(n);
        
        int i = 0;
        face_key seed_face;
        bool seeded = false;
        int seed_index = -1;
        
        // Find the face intersected by the line connecting vv[0] and vv[1] to seed the dynamic programming algorithm.
        
        for (auto fit = multi_face.faces_begin(); fit != multi_face.faces_end(); fit++)
        {
            face_quality_values[*fit] = min_quality(*fit);
            faces[i] = *fit;
            if (!seeded)
            {
                simplex_set multi_face, mf_boundary;
                multi_face.insert(*fit);
                Complex::boundary_2manifold(multi_face, mf_boundary);
                double q = edge_remove_quality(mf_boundary, verts);
                if (q >= 0)
                {
                    seed_face = *fit;
                    seed_index = i;
                    seeded = true;
                }
            }
            ++i;
        }
        
        if (!seeded) return;
        
        // Find an optimal subset of the multi-face for the removal (using dynamic programming).
        int max_quality_index = -1;
        T max_quality = -1.0;
        std::vector<std::vector<bool> > mf_subsets;
        
        std::vector<bool> seed_set(n);
        seed_set[seed_index] = true;
        mf_subsets.push_back(seed_set);
        max_quality_index = 0;
        max_quality = omfr_quality_helper(n, faces, face_quality_values, seed_set, verts);
        unsigned int prev = 0, next = 1;
        
        for (unsigned int i = 2; i <= n; ++i) // size of the created subset
        {
            for (unsigned int j = prev; j < next; ++j) // index of the set I'm expanding
            {
                for (unsigned int k = 0; k < n; ++k) // potential new face
                {
                    if (!mf_subsets[j][k])
                    {
                        std::vector<bool> subset(n);
                        for (unsigned int l = 0; l < n; ++l)
                            subset[l] = mf_subsets[j][l] || (l == k);
                        if (next == mf_subsets.size() || Util::compare(n, subset, mf_subsets[mf_subsets.size()-1]))
                        {
                            bool connected = false;
                            for (unsigned int l = 0; l < n; ++l)
                            {
                                if (mf_subsets[j][l])
                                {
                                    if (Complex::get_edge(faces[k], faces[l]) != Complex::NULL_EDGE)
                                        connected = true;
                                }
                            }
                            
                            if (connected)
                            {
                                mf_subsets.push_back(subset);
                                T q = omfr_quality_helper(n, faces, face_quality_values, seed_set, verts);
                                if (q > max_quality)
                                {
                                    max_quality = q;
                                    max_quality_index = static_cast<int>(mf_subsets.size())-1;
                                }
                            }
                        }
                    }
                }
            }
            prev = next;
            next = static_cast<int>(mf_subsets.size());
        }
        
        // Attempt the removal if multi-face removal would improve quality.
        if (max_quality > q_old && max_quality > 0)
        {
            std::vector<bool> subset = mf_subsets[max_quality_index];
            simplex_set multi_face;
            for (unsigned int l = 0; l < n; ++l)
            {
                if (subset[l])
                    multi_face.insert(faces[l]);
            }
            
            Complex::multi_face_remove(multi_face, new_simplices);
            cleanup_set(new_simplices);
        }
    }
    
    /**
     * Attempt to remove face f using multi-face removal.
     */
    void multi_face_removal_old(const face_key& f)
    {
        simplex_set st_f;
        Complex::star(f, st_f);
        int label = get_label(*st_f.tetrahedra_begin());
        
        simplex_set new_simplices;
        optimal_multi_face_removal(f, new_simplices);
        
        for (auto tit = new_simplices.tetrahedra_begin(); tit != new_simplices.tetrahedra_end(); tit++)
        {
            set_label(*tit, label);
        }
        
        simplex_set ns_cl;
        Complex::closure(new_simplices, ns_cl);
        update(ns_cl);
    }
    
    /**
     * Improve tetrahedra by topological operations (re-connection): edge removal, multi-face removal,
     * multi-face retriangulation. It do so only for tetrahedra of quality lower than MIN_TET_QUALITY.
     */
    void topological_pass()
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
        int i = 0, j = 0;
        for (auto &t : tets)
        {
            if (Complex::exists(t) && quality(t) < MIN_TET_QUALITY)
            {
                simplex_set cl_t;
                Complex::closure(t, cl_t);
                
                for (auto eit = cl_t.edges_begin(); eit != cl_t.edges_end(); eit++)
                {
                    if (Complex::exists(*eit) && !is_interface(*eit) && !is_boundary(*eit))
                    {
                        if(edge_removal(*eit))
                        {
                            i++;
                        }
                        j++;
                    }
                }
            }
        }
        std::cout << "Topological edge removals: " << i << " / " << j << std::endl;
        
        // Attempt to remove each face of each remaining tetrahedron in tets using multi-face removal.
        // Accept if it increases the minimum quality locally.
        i = 0, j = 0;
        for (auto &t : tets)
        {
            if (Complex::exists(t) && quality(t) < MIN_TET_QUALITY)
            {
                simplex_set cl_t;
                Complex::closure(t, cl_t);
                
                for (auto fit = cl_t.faces_begin(); fit != cl_t.faces_end(); fit++)
                {
                    if (Complex::exists(*fit) && !is_interface(*fit) && !is_boundary(*fit))
                    {
                        if(multi_face_removal(*fit))
                        {
                            i++;
                        }
                        j++;
                    }
                }
            }
        }
        std::cout << "Topological face removals: " << i << " / " << j << std::endl;
        
        Complex::garbage_collect();
    }
    
    ////////////////////
    // EDGE FLIP PASS //
    ////////////////////
    
    void interface_edge_flip_pass()
    {
        std::vector<edge_key> flippable_edges;
        
        for (auto eit = Complex::edges_begin(); eit != Complex::edges_end(); eit++)
        {
            if (eit->is_interface() || eit->is_boundary())
            {
                flippable_edges.push_back(eit.key());
            }
        }
        
        for (auto e : flippable_edges)
        {
            if (Complex::exists(e) && is_flippable(e))
            {
                flip(e);
            }
        }
    }
    
    
    
    /////////////////////
    // THICKENING PASS //
    /////////////////////
    
    /**
     * Splits all tetrahedra with a volume greater than MAX_TET_VOLUME by inserting a vertex.
     */
    void thickening_pass()
    {
        std::vector<tet_key> tetrahedra;
        for (auto tit = Complex::tetrahedra_begin(); tit != Complex::tetrahedra_end(); tit++)
        {
            if (volume(tit.key()) > MAX_TET_VOLUME)
            {
                tetrahedra.push_back(tit.key());
            }
        }
        int i = 0, j = 0;
        for(auto &t : tetrahedra)
        {
            if (Complex::exists(t) && volume(t) > MAX_TET_VOLUME)
            {
                if(split(t) != Complex::NULL_NODE)
                {
                    i++;
                }
                j++;
            }
        }
        std::cout << "Thickening pass splits: " << i << "/" << j << std::endl;
    }
    
    ///////////////////
    // THINNING PASS //
    ///////////////////
    /**
     * Collapses all edges which is shorter than MIN_EDGE_LENGTH. However, the collapse is only performed if it does not alter the interface and it improves the minimum quality of the tetrahedra in the star of the edge.
     */
    void thinning_pass()
    {
        std::vector<tet_key> tetrahedra;
        for (auto tit = Complex::tetrahedra_begin(); tit != Complex::tetrahedra_end(); tit++)
        {
            if (volume(tit.key()) < MIN_TET_VOLUME)
            {
                tetrahedra.push_back(tit.key());
            }
        }
        int i = 0, j = 0;
        for(auto &t : tetrahedra)
        {
            if (Complex::exists(t) && volume(t) < MIN_TET_VOLUME)
            {
                if(collapse(t))
                {
                    i++;
                }
                j++;
            }
        }
        std::cout << "Thinning pass collapses: " << i << "/" << j << std::endl;
    }
    
    ///////////////////////
    // RELABEL TETS PASS //
    ///////////////////////
    
    /**
     * Relabel all tetrahedra in the SimplicialComplex mesh when
     * a) their quality is lower than TET_RELABEL_THRESHOLD and the minimum sine of their dihedral angles is lower than TET_RELABEL_THRESHOLD
     * b) they fulfill criteria for relabelling (described in detail in relabel_condition()).
     * Compliant with multiple phases.
     */
    bool relabel_tets()
    {
        int counter = 0;
        for (auto tit = Complex::tetrahedra_begin(); tit != Complex::tetrahedra_end(); tit++)
        {
            if(relabel_tet(tit.key()))
            {
                ++counter;
            }
        }
        
        if (counter > 0)
        {
            std::cout << "Relabeled " << counter << " tets" << std::endl;
        }
        
        return (counter > 0);
    }
    
    /////////////////
    // SMOOTH PASS //
    /////////////////
    
    void smooth()
    {
        for (auto nit = Complex::nodes_begin(); nit != Complex::nodes_end(); nit++)
        {
            if (Complex::exists(nit.key()) && !is_boundary(nit.key()) && !is_interface(nit.key()))
            {
                if (!smart_laplacian(nit.key()))
                {
//                    if (min_quality(nit.key()) < MIN_TET_QUALITY)
//                    {
//                        freitag_smoothing(nit.key());
//                    }
                }
            }
        }
    }
    
    //////////////////////////////
    // REMOVE DEGENERACIES PASS //
    //////////////////////////////
    /**
     * Attempt to remove edges shorter than DEG_EDGE_LENGTH by collapsing them.
     */
    void remove_degenerate_edges()
    {
        std::list<edge_key> edges;
        for (auto eit = Complex::edges_begin(); eit != Complex::edges_end(); eit++)
        {
            if (length(eit.key()) < DEG_EDGE_LENGTH)
            {
                edges.push_back(eit.key());
            }
        }
        int i = 0, j = 0;
        for(auto e : edges)
        {
            if(Complex::exists(e) && length(e) < DEG_EDGE_LENGTH)
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
            if(min_angle(fit.key()) < DEG_ANGLE)
            {
                faces.push_back(fit.key());
            }
        }
        
        int i = 0, j = 0;
        for (auto &f : faces)
        {
            if (Complex::exists(f) && min_angle(f) < DEG_ANGLE)
            {
                if(collapse(f, false))
                {
                    i++;
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
                j++;
            }
        }
        std::cout << "Removed " << i <<"/"<< j << " degenerate tets" << std::endl;
        Complex::garbage_collect();
    }
    
    //////////////////////////////////
    // REMOVE LOW QUALITY SIMPLICES //
    //////////////////////////////////
    
    bool remove_edge(edge_key& e)
    {
        node_key n = collapse(e);
        if(n == Complex::NULL_NODE)
        {
            n = collapse(e, false);
        }
        return n != Complex::NULL_NODE;
    }
    
    /**
     * Attempt to remove edges shorter than MIN_EDGE_LENGTH by safely collapsing them.
     */
    void remove_edges()
    {
        std::list<edge_key> edges;
        for (auto eit = Complex::edges_begin(); eit != Complex::edges_end(); eit++)
        {
            if (length(eit.key()) < DEG_EDGE_LENGTH)
            {
                edges.push_back(eit.key());
            }
        }
        int i = 0, j = 0;
        for(auto e : edges)
        {
            if(Complex::exists(e) && length(e) < DEG_EDGE_LENGTH)
            {
                if(remove_edge(e))
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
        std::vector<V> verts;
        get_pos(e, verts);
        V p = Util::project<MT>(get_pos(apex), verts[0], verts[1]);
        
        // Split longest edge
        node_key n = split(e);
        set_pos(n, p);
        
        // Collapse new edge
        edge_key e_rem = Complex::get_edge(apex, n);
        return collapse(e_rem, false) != Complex::NULL_NODE;
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
        return collapse(e, false) != Complex::NULL_NODE;
    }
    
    /**
     * Attempt to remove the face f by first determining whether it's a cap or a needle.
     */
    bool remove_face(const face_key& f)
    {
        if(max_angle(f) > M_PI - MIN_ANGLE)
        {
            return remove_cap(f);
        }
        return remove_needle(f);
    }
    
    /**
     * Attempts to remove degenerate faces (faces with a minimum angle smaller than DEG_ANGLE).
     */
    void remove_faces()
    {
        std::list<face_key> faces;
        
        for (auto fit = Complex::faces_begin(); fit != Complex::faces_end(); fit++)
        {
            if(min_angle(fit.key()) < DEG_ANGLE)
            {
               faces.push_back(fit.key());
            }
        }
        
        int i = 0, j = 0;
        for (auto &f : faces)
        {
            if (Complex::exists(f) && min_angle(f) < DEG_ANGLE)
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
        return collapse(e, false) != Complex::NULL_NODE;
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
        std::vector<V> verts;
        get_pos(f, verts);
        V p = Util::project<MT>(get_pos(apex), verts);
        
        // Split the face
        node_key n = split(f);
        set_pos(n, p);
        
        // Collapse edge
        edge_key e = Complex::get_edge(n, apex);
        return collapse(e, false) != Complex::NULL_NODE;
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
            if(collapse(e, false) != Complex::NULL_NODE)
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
            if(collapse(e, false) != Complex::NULL_NODE)
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
        std::vector<V> verts;
        get_pos(f, verts);
        V proj_apex = Util::project<MT>(get_pos(apex), verts);
        
        // Find barycentric coordinates
        std::vector<T> barycentric_coords = Util::barycentric_coords<MT>(proj_apex, verts[0], verts[1], verts[2]);
        
        if(barycentric_coords[0] > 0.2 && barycentric_coords[1] > 0.2 && barycentric_coords[2] > 0.2) // The tetrahedron is a cap
        {
            return remove_cap(t);
        }
        else if(barycentric_coords[0] < -0.2 || barycentric_coords[1] < -0.2 || barycentric_coords[2] < -0.2) // The tetrahedron is a sliver
        {
            return remove_sliver(t);
        }
        
        T mean_dist = 0.;
        for(V &p : verts)
        {
            mean_dist += (p-proj_apex).length()/3.;
        }
        int close = 0;
        for(V &p : verts)
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
     * Attempt to remove tetrahedra with quality lower than DEG_TET_QUALITY.
     */
    void remove_tets()
    {
        std::vector<tet_key> tets;
        
        for (auto tit = Complex::tetrahedra_begin(); tit != Complex::tetrahedra_end(); tit++)
        {
            if (quality(tit.key()) < DEG_TET_QUALITY)
            {
                tets.push_back(tit.key());
            }
        }
        int i = 0, j=0;
        for (auto &tet : tets)
        {
            if (Complex::exists(tet) && quality(tet) < DEG_TET_QUALITY)
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
    
    ///////////////////
    // FIX FUNCTIONS //
    ///////////////////
    
    void fix_complex()
    {
        simplex_set region;
        
        print_out("Smooth.");
        smooth();
        validity_check();
        
        print_out("Topological pass.");
        topological_pass();
        validity_check();
        
        print_out("Edge flip pass.");
        interface_edge_flip_pass();
        validity_check();
        
        print_out("Remove low quality.");
        remove_tets();
        validity_check();
        remove_faces();
        validity_check();
        remove_edges();
        validity_check();
        
//        print_out("Remove degeneracies.");
//        remove_degenerate_tets();
//        validity_check();
//        remove_degenerate_faces();
//        validity_check();
//        remove_degenerate_edges();
//        validity_check();
        
        print_out("Smooth.");
        smooth();
        validity_check();
        
//        print_out("Relabel tets.");
//        bool relabeled = relabel_tets();
//        validity_check();
//        
//        if (relabeled) {
//            smooth();
//        }
    }
    
    void resize_complex()
    {
        print_out("Thickening pass.");
        thickening_pass();
        validity_check();
        
        print_out("Thinning pass.");
        thinning_pass();
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
                if (nit->is_interface() && Complex::exists(nit.key()))
                {
                    if(!move_vertex(nit.key(), get_destination(nit.key())))
                    {
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
        
        if(step_no%5 == 0)
        {
//            resize_complex();
        }
        
        Complex::garbage_collect();
        ++step_no;
    }
    
private:
    
    /**
     * Tries moving the node n to the new position new_pos. Returns true if it succeeds.
     */
    bool move_vertex(const node_key & n, const typename MT::vector3_type & new_pos)
    {
        V pos = get_pos(n);
        T l = MT::length(new_pos - pos);
        
        if (l < EPSILON) // The vertex is not moved
        {
            return true;
        }
        
        T t = intersection_with_link(n, new_pos);
        
        l = std::max(std::min(l*t - l*MIN_DEFORMATION, l), 0.);
        set_pos(n, pos + l*normalize(new_pos - pos));
        
        if (MT::length(new_pos - pos) < EPSILON)
        {
            return true;
        }
        return false;
    } // move_vertex
    
    
    /**
     * Returns the intersection point (= pos + t*(new_pos-pos)) with the link of the node n and
     * when moving the node n to the new position new_pos.
     */
    T intersection_with_link(const node_key & n, const V& new_pos)
    {
        V pos = get_pos(n);
        simplex_set ln;
        Complex::link(n,ln);
        T min_t = INFINITY;
        std::vector<V> face_pos;
        for(auto fit = ln.faces_begin(); fit != ln.faces_end(); fit++)
        {
            get_pos(*fit, face_pos);
            T t = Util::intersection_ray_plane<MT>(pos, new_pos, face_pos[0], face_pos[1], face_pos[2]);
            if (0. <= t)
            {
                min_t = std::min(t, min_t);
            }
        }
        assert(min_t < INFINITY);
        return min_t;
    }
    
    ///////////////
    // SMOOTHING //
    ///////////////
private:
    
    
    /**
     * Performs Laplacian smoothing if it improves the minimum tetrahedron quality locally.
     */
    bool smart_laplacian(const node_key& n, T alpha = 1.)
    {
        T q_old = min_quality(n);
        V old_pos = get_pos(n);
        V avg_pos = get_barycenter(n);
        set_pos(n, old_pos + alpha * (avg_pos - old_pos));
        
        T q_new = min_quality(n);
        if (q_new < q_old)
        {
            set_pos(n, old_pos);
            return false;
        }
        return true;
    }
    
    /// Find the active set for Freitag smoothing.
    T get_quality_info(const node_key & n, std::vector<quality_struct<MT>> & full_quality_info, std::list<quality_struct<MT>> & active)
    {
        simplex_set st_n;
        Complex::star(n, st_n);
        
        T min_quality = INFINITY;
        
        for (auto tit = st_n.tetrahedra_begin(); tit != st_n.tetrahedra_end(); tit++)
        {
            quality_struct<MT> quality_info;
            
            std::vector<node_key> nodes;
            Complex::get_nodes(*tit, nodes);
            
            int alpha = 0;
            for (int k = 0; k < 4; ++k)
            {
                if (nodes[k] == n)
                {
                    alpha = k;
                    break;
                }
            }
            assert(alpha != -1 || !"n is not a vertex of t");
            
            std::vector<V> verts;
            get_pos(*tit, verts);
            quality_info.quality = quality(*tit);
            quality_info.grad_quality = Util::grad_rms_length<MT>(verts, alpha);
            
            full_quality_info.push_back(quality_info);
            
            if ((min_quality - quality_info.quality) >= EPSILON)
            {
                active.clear();
                min_quality = quality_info.quality;
            }
            
            if ((quality_info.quality - min_quality) < EPSILON)
            {
                active.push_back(quality_info);
            }
        }
        
        return min_quality;
    }
    
    /**
     * Helper function for the vertex position optimizers. Finds the approximate step-length
     * based on the Taylor expansion of the objective functions.
     *
     * @param full_quality_info     (Multi-valued) objective function values with gradients.
     * @param predicted_quality     Predicted value of the objective function.
     *
     * @return  Estimated step-length.
     */
    T initial_step_length(std::vector<quality_struct<MT>> & full_quality_info, T & predicted_quality)
    {
        std::sort(&full_quality_info[0], &full_quality_info[0] + full_quality_info.size());
        std::list<quality_struct<MT>> active_elements;
        std::list<T> x_values;
        active_elements.push_front(full_quality_info[0]);
        
        T x_max = 0, y_max = full_quality_info[0].quality;
        
        x_values.push_front(x_max);
        
        for (unsigned int i = 1; i < full_quality_info.size(); ++i)
        {
            quality_struct<MT> active_element = *(active_elements.begin());
            
            if (full_quality_info[i].quality == active_element.quality ||
                full_quality_info[i].slope   >= active_element.slope)
                continue;
            
            T dy, x, y;
            
            do
            {
                active_element = *(active_elements.begin());
                
                dy = full_quality_info[i].quality - active_element.quality;
                x = dy / (active_element.slope - full_quality_info[i].slope);
                y = active_element.quality + x * active_element.slope;
                assert(!MT::is_nan(x) && !MT::is_nan(y));
                
                if (x < *(x_values.begin()))
                {
                    active_elements.pop_front();
                    x_values.pop_front();
                }
                else
                    break;
                
            } while (!x_values.empty());
            
            if (active_element.slope > 0)
            {
                active_elements.push_front(full_quality_info[i]);
                x_values.push_front(x);
                y_max = y;
            }
        }
        
        predicted_quality = y_max;
        return *(x_values.begin());
    }
    
    T calc_step_length(const node_key &n, T init_alpha, V s, T min_q_old, T min_q, T min_q_pre, T min_step)
    {
        T alpha = init_alpha;
        
        T quality_prev = min_q_old;
        V pos = get_pos(n);
        int i = 0;
        while (abs(alpha) > min_step)
        {
            set_pos(n, pos + alpha * s);
            min_q = min_quality(n);
            if (min_q < quality_prev)
            {
                if (i > 0)
                {
                    alpha *= 2.0;
                    break;
                }
                else
                {
                    alpha /= 2.0;
                    continue;
                }
            }
            if (min_q - min_q_old > 0.6 * (min_q_pre - min_q_old))
            {
                break;
            }
            alpha /= 2.0;
            quality_prev = min_q;
            ++i;
        }
        if (i == 0)
        {
            alpha = -init_alpha;
            
            T quality_prev = min_q_old;
            
            while (abs(alpha) > min_step)
            {
                set_pos(n, pos + alpha * s);
                min_q = min_quality(n);
                if (min_q < quality_prev)
                {
                    alpha *= 2.0;
                    break;
                }
                if (min_q - min_q_old > 0.6 * (min_q_pre - min_q_old))
                    break;
                alpha /= 2.0;
                quality_prev = min_q;
                ++i;
            }
        }
        return alpha;
    }
    
public:
    void freitag_smoothing(node_key const & n, T min_step = 1e-4, T min_improvement = 1e-5, int max_iter = 3)
    {
        std::list<quality_struct<MT>> active_set = std::list<quality_struct<MT>>();
        std::vector<quality_struct<MT>> full_quality_info = std::vector<quality_struct<MT>>();
        std::vector<V> g;
        
        T min_q = min_quality(n);
        T alpha = INFINITY;
        T min_q_old;
        V s;
        int iter = 0;
        do
        {
            min_q_old = min_q;
            
            full_quality_info.clear();
            active_set.clear();
            min_q = get_quality_info(n, full_quality_info, active_set);
            if (min_q < 0 || active_set.size() == 0)
            {
                return;
            }
            
            g.clear();
            for (auto &it : active_set)
            {
                g.push_back(it.grad_quality);
            }
            
            s = Util::min_convex_hull_point<MT>(g);
            if (MT::length(s) < EPSILON)
            {
                return;
            }
            
            T r = INFINITY;
            for (auto &it : active_set)
            {
                r = std::min(r, MT::dot(s, it.grad_quality));
            }
            
            for(auto &q : full_quality_info)
            {
                q.slope = MT::dot(s, q.grad_quality);
            }
            
            T predicted_quality = 0.;
            T init_alpha = initial_step_length(full_quality_info, predicted_quality);
            alpha = calc_step_length(n, init_alpha, s, min_q_old, min_q, predicted_quality, min_step);
            if (alpha < min_step)
            {
                return;
            }
            V pos = get_pos(n);
            set_pos(n, pos + alpha * s);
            min_q = min_quality(n);
            
            if (min_q < min_q_old)
            {
                set_pos(n, pos);
                break;
            }
            ++iter;
            
        } while (iter <= max_iter && alpha > min_step && std::abs(min_q - min_q_old) > min_improvement);
    }
    
    ////////////////
    // RELABELING //
    ////////////////
private:
    /**
     * Precondition for relable_tet function.
     * Returns true when
     * a) the quality of tetrahedron t is lower than DEG_TET_QUALITY.
     * b) the biggest face f of tetrahedron t lies on the interface.
     * c) the area of f is at least 0.48 of the total area of all the faces of t.
     * d) the vertex opposite to f lies on the interface.
     */
    bool precond_relabel(const tet_key& t, face_key & f)
    {
        if (quality(t) > DEG_TET_QUALITY)
        {
            return false;
        }
        
        // Find the face (f) of t that lies on the interface and has the biggest area (max_area).
        // Compute the total area of all faces of t (total_area).
        T max_area = 0.;
        T total_area = 0.;
        simplex_set cl_t;
        Complex::closure(t, cl_t);
        for (auto fit = cl_t.faces_begin(); fit != cl_t.faces_end(); fit++)
        {
            T a = area(*fit);
            total_area += a;
            
            if (is_interface(*fit) && !is_boundary(*fit) && a > max_area)
            {
                max_area = a;
                f = *fit;
            }
        }
        
        // Relabel the tetrahedron if it's biggest face lies on the interface, it's at least 0.48 of the sum of all it's faces' areas
        // and if the vertex opposite to this face also lies on the interface.
        if (max_area > 0.48*total_area)
        {
            std::vector<V> verts;
            get_pos(f, verts);
            
            simplex_set lk_f;
            Complex::link(f, lk_f);
            for(auto nit = lk_f.nodes_begin(); nit != lk_f.nodes_end(); nit++)
            {
                if (cl_t.contains(*nit) && is_interface(*nit) &&
                    Util::distance<MT>(get_pos(*nit), verts[0], verts[1], verts[2]) < MIN_EDGE_LENGTH)
                {
                    return true;
                }
            }
        }
        return false;
    }
    
    /**
     * Relabels tetrahedron t if precond_relabel() returns true.
     */
    bool relabel_tet(const tet_key& t)
    {
        face_key f;
        
        if (precond_relabel(t, f))
        {
            int label = -1;
            simplex_set st_f;
            Complex::star(f, st_f);
            for (auto tit = st_f.tetrahedra_begin(); tit != st_f.tetrahedra_end(); tit++)
            {
                if(*tit != t)
                {
                    label = get_label(*tit);
                }
            }
            
            if (label != -1 && label != get_label(t))
            {
                set_label(t, label);
                simplex_set cl_t;
                Complex::closure(t, cl_t);
                update(cl_t);
                return true;
            }
        }
        return false;
    }
    
    
    ///////////
    // FLIPS //
    ///////////
private:
    /**
     * Returns whether it is possible to flip the edge e or not, i.e.: whether the edge does not fulfill Delaunay criterion, and whether
     * it is not a feature edge (it is not a feature edge if its neighborhood is sufficiently flat).
     */
    bool is_flippable(const edge_key & e)
    {
        std::vector<V> verts;
        get_pos(e, verts);
        
        simplex_set ln_e;
        Complex::link(e, ln_e);
        for (auto nit = ln_e.nodes_begin(); nit != ln_e.nodes_end(); nit++)
        {
            if (is_interface(*nit) || is_boundary(*nit))
            {
                verts.push_back(get_pos(*nit));
            }
        }
        
        if (verts.size() != 4)
        {
            return false;
        }
        
        T alpha0 = Util::angle<MT>(verts[2], verts[0], verts[1]);
        T alpha1 = Util::angle<MT>(verts[3], verts[0], verts[1]);
        
        if (alpha0 + alpha1 <= M_PI)
        {
            return false;
        }
        
        T flatness = Util::calc_flatness<MT>(verts[0], verts[1], verts[2], verts[3]);
        
        if (flatness >= FLIP_EDGE_TET_FLATNESS)
        {
            return true;
        }
        
        return false;
    }
    
    /**
     * Flips the edge e (which is a special case of the edge remove operation in the embedding mesh).
     * Relabels the tetrahedra accordingly so that the interface mesh geometry does not change.
     */
    void flip(const edge_key & e)
    {
        // Find the pair of nodes to be connected by a new edge as a result of edge flip.
        simplex_set lk_e;
        Complex::link(e, lk_e);
        std::vector<node_key> nodes;
        for (auto nit = lk_e.nodes_begin(); nit != lk_e.nodes_end(); nit++)
        {
            if (is_interface(*nit) || is_boundary(*nit)) {
                nodes.push_back(*nit);
            }
        }
        
        assert(nodes.size() == 2);
        if(Complex::get_edge(nodes[0], nodes[1]) != Complex::NULL_EDGE) // Check that there does not already exist an edge between the pair of nodes.
        {
            return;
        }
        
        node_key n_new = split(e);
        
        edge_key e_c = Complex::get_edge(nodes[0], n_new);
        assert(e_c != Complex::NULL_EDGE);
        V p = get_pos(nodes[0]);
        V p_new = get_destination(nodes[0]);
        
        if(precond_collapse(e_c, p))
        {
            collapse(e_c, p, p_new);
        }
    }
    
    ////////////
    // SPLITS //
    ////////////
public:
    /**
     * Split a tetrahedron t and returns the new node which is positioned at the barycenter of the vertices of t.
     */
    node_key split(const tet_key& t)
    {
        std::vector<V> verts;
        get_pos(t, verts);
        V p = Util::barycenter<MT>(verts[0], verts[1], verts[2], verts[3]);
        
        node_key n = Complex::split(t);
        set_pos(n, p);
        set_destination(n, p);
        
        simplex_set st_n;
        Complex::star(n, st_n);
        st_n.insert(n);
        update(st_n);
        return n;
    }
    
    /**
     * Split a face f and returns the new node which is positioned at the barycenter of the vertices of f.
     */
    node_key split(const face_key & f)
    {   
        std::vector<node_key> nodes;
        Complex::get_nodes(f, nodes);
        V p = Util::barycenter<MT>(get_pos(nodes[0]), get_pos(nodes[1]), get_pos(nodes[2]));
        V p_new = p;
        if(is_interface(f))
        {
            p_new = Util::barycenter<MT>(get_destination(nodes[0]), get_destination(nodes[1]), get_destination(nodes[2]));
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
        std::vector<node_key> nodes;
        Complex::get_nodes(e, nodes);
        V p = Util::barycenter<MT>(get_pos(nodes[0]), get_pos(nodes[1]));
        V p_new = p;
        if(is_interface(e))
        {
            p_new = Util::barycenter<MT>(get_destination(nodes[0]), get_destination(nodes[1]));
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
    node_key collapse(edge_key& e, const V& p, const V& p_new)
    {
        node_key n_new = Complex::collapse(e);
        
        if (Complex::exists(n_new) && n_new != Complex::NULL_NODE)
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
    bool precond_collapse(const edge_key& e, const V& p)
    {
        std::vector<node_key> nodes;
        Complex::get_nodes(e, nodes);
        
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
     * Returns the minimum quality of neighbouring tetrahedra if the edge e is collapsed and the resulting node is moved to v_new.
     */
    T min_quality(const edge_key& e, const V& v_new)
    {
        std::vector<node_key> nodes;
        Complex::get_nodes(e, nodes);
        
        simplex_set lk_n0, lk_n1;
        Complex::link(nodes[0], lk_n0);
        Complex::link(nodes[1], lk_n1);
        
        simplex_set st_e, cl_st_e;
        Complex::star(e, st_e);
        Complex::closure(st_e, cl_st_e);
        lk_n0.add(lk_n1);
        lk_n0.difference(cl_st_e);
        return min_quality(lk_n0, v_new);
    }
    
    /**
     * Collapses the edge e and places the new node at the most optimal position of the position of either end node or their barycenter.
     * If the parameter safe is true, the method if the nodes of edge e are editable, i.e. not a part of the interface, and will therefore not change the interface.
     * If non of the nodes are editable or precond_collapse returns false, the method returns NULL_NODE.
     */
    node_key collapse(edge_key& e, bool safe = true)
    {
        if (!Complex::exists(e) || e == Complex::NULL_EDGE)
        {
            return Complex::NULL_NODE;
        }
        std::vector<node_key> nodes;
        Complex::get_nodes(e, nodes);
        
        V p_opt, p_new_opt;
        T q_max = 0.;
        
        if (!is_boundary(nodes[0]) && !is_boundary(nodes[1]) && (!safe || (!is_interface(nodes[0]) && !is_interface(nodes[1]))))
        {
            V p = Util::barycenter<MT>(get_pos(nodes[0]), get_pos(nodes[1]));
            T q = min_quality(e, p);
            if (precond_collapse(e, p) && q > q_max)
            {
                p_new_opt = Util::barycenter<MT>(get_destination(nodes[0]), get_destination(nodes[1]));
                p_opt = p;
                q_max = q;
            }
        }
        
        if (!is_boundary(nodes[0]) && (!safe || !is_interface(nodes[0])))
        {
            V p = get_pos(nodes[1]);
            T q = min_quality(e, p);
            
            if (precond_collapse(e, p) && q > q_max)
            {
                p_new_opt = get_destination(nodes[1]);
                p_opt = p;
                q_max = q;
            }
        }
        
        if (!is_boundary(nodes[1]) && (!safe || !is_interface(nodes[1])))
        {
            V p = get_pos(nodes[0]);
            T q = min_quality(e, p);
            
            if (precond_collapse(e, p) && q > q_max)
            {
                p_new_opt = get_destination(nodes[0]);
                p_opt = p;
                q_max = q;
            }
        }
        T q = std::min(min_quality(nodes[0]), min_quality(nodes[1]));
        if((!safe && q_max > 0.) || (safe && q_max > q))
        {
            return collapse(e, p_opt, p_new_opt);
        }
        return Complex::NULL_NODE;
    }
    
    bool collapse(const face_key& f, bool safe = true)
    {
        simplex_set cl_f;
        Complex::closure(f, cl_f);
        while(cl_f.size_edges() > 0)
        {
            edge_key e = shortest_edge(cl_f);
            if(collapse(e, safe) != Complex::NULL_NODE)
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
            if(collapse(e, safe) != Complex::NULL_NODE)
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
     Returns the normal to interface face f.
     */
    V get_normal(const face_key & f)
    {
        std::vector<V> verts;
        get_pos(f, verts);
        return Util::normal_direction<MT>(verts[0], verts[1], verts[2]);
    }
    
    /**
     Returns the normal to interface node n.
     */
    V get_normal(const node_key & n)
    {
        simplex_set st;
        Complex::star(n, st);
        
        V result(0.);
        for (auto fit = st.faces_begin(); fit != st.faces_end(); fit++)
        {
            if (is_interface(*fit))
            {
                result += get_normal(*fit);
            }
        }
        if (MT::length(result) < EPSILON) {
            return V(0.);
        }
#ifdef DEBUG
        assert(!MT::is_nan(result[0]) && !MT::is_nan(result[1]) && !MT::is_nan(result[2]));
#endif
        return MT::normalize(result);
    }
    
    /**
     * Calculates the average position of the neighbouring nodes to node n.
     * If interface is true, the average position is only calculated among the neighbouring nodes which are interface.
     */
    V get_barycenter(const node_key& n, bool interface = false)
    {
        simplex_set lk_n;
        Complex::link(n, lk_n);
        
        V avg_pos(0.);
        int i = 0;
        for (auto nit = lk_n.nodes_begin(); nit != lk_n.nodes_end(); nit++)
        {
            if (!interface || is_interface(*nit))
            {
                avg_pos += get_pos(*nit);
                i++;
            }
        }
        assert(i != 0);
        return avg_pos / static_cast<T>(i);
    }
    
    ///////////////////////
    // UTILITY FUNCTIONS //
    ///////////////////////
    
public:
    
    T length(const edge_key& e)
    {
        std::vector<V> verts;
        get_pos(e,verts);
        return MT::length(verts[0] - verts[1]);
    }
    
    T area(const face_key& f)
    {
        std::vector<V> verts;
        get_pos(f,verts);
        return Util::area<MT>(verts);
    }
    
    T volume(const tet_key& t)
    {
        std::vector<V> verts;
        get_pos(t,verts);
        return Util::volume<MT>(verts[0], verts[1], verts[2], verts[3]);
    }
    
    T signed_volume(const tet_key& t)
    {
        std::vector<V> verts;
        get_pos(t,verts);
        return Util::signed_volume<MT>(verts[0], verts[1], verts[2], verts[3]);
    }
    
    T quality(const tet_key& t)
    {
        std::vector<V> verts;
        get_pos(t, verts);
        return Util::quality<MT>(verts[0], verts[1], verts[2], verts[3]);
    }
    
    T min_angle(const face_key& f)
    {
        std::vector<V> verts;
        get_pos(f, verts);
        return Util::min_angle<MT>(verts[0], verts[1], verts[2]);
    }
    
    T max_angle(const face_key& f)
    {
        std::vector<V> verts;
        get_pos(f, verts);
        return Util::max_angle<MT>(verts[0], verts[1], verts[2]);
    }
    
    /**
     * Returns the largest face in the simplex set.
     */
    face_key largest_face(simplex_set& set)
    {
        T max_a = -INFINITY;
        face_key max_f;
        for(auto f = set.faces_begin(); f != set.faces_end(); f++)
        {
            T a = area(*f);
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
        T min_l = INFINITY;
        edge_key min_e;
        for(auto e = set.edges_begin(); e != set.edges_end(); e++)
        {
            T l = length(*e);
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
        T max_l = -INFINITY;
        edge_key max_e;
        for(auto e = set.edges_begin(); e != set.edges_end(); e++)
        {
            T l = length(*e);
            if(l > max_l)
            {
                max_l = l;
                max_e = *e;
            }
        }
        return max_e;
    }
    
    /// Checks whether the interface patch around n is 2-manifold
    //    bool is_2manifold(const node_key & n)
    //    {
    //        typename SimplicialComplex<MT>::simplex_set_type st_n, interface_simplices;
    //        SimplicialComplex<MT>::star(n, st_n);
    //
    //        std::vector<int> labels;
    //        int l;
    //        for (typename SimplicialComplex<MT>::simplex_set_type::tetrahedron_set_iterator tit = st_n.tetrahedra_begin(); tit != st_n.tetrahedra_end(); tit++)
    //        {
    //            l = SimplicialComplex<MT>::find_tetrahedron(*tit).label;
    //            if(find(labels.begin(), labels.end(), l) == labels.end())
    //                labels.push_back(l);
    //        }
    //
    //        if (labels.size() != 2) return false;
    //
    //        typename SimplicialComplex<MT>::simplex_set_type::edge_set_iterator eit = st_n.edges_begin();
    //        while (eit != st_n.edges_end())
    //        {
    //            if (SimplicialComplex<MT>::find_edge(*eit).is_interface() ||
    //                (!(SimplicialComplex<MT>::find_node(n).is_interface()) && SimplicialComplex<MT>::find_edge(*eit).is_boundary()))
    //                interface_simplices.insert(*eit);
    //            ++eit;
    //        }
    //
    //        typename SimplicialComplex<MT>::simplex_set_type::face_set_iterator fit = st_n.faces_begin();
    //        while (fit != st_n.faces_end())
    //        {
    //            if (SimplicialComplex<MT>::find_face(*fit).is_interface() ||
    //                (!(SimplicialComplex<MT>::find_node(n).is_interface()) && SimplicialComplex<MT>::find_face(*fit).is_boundary()))
    //                interface_simplices.insert(*fit);
    //            ++fit;
    //        }
    //
    //        typename SimplicialComplex<MT>::face_key f = *(interface_simplices.faces_begin());
    //        do
    //        {
    //            interface_simplices.erase(f);
    //            typename SimplicialComplex<MT>::edge_key e;
    //            typename SimplicialComplex<MT>::simplex_set_type::face_set_iterator fit = interface_simplices.faces_begin();
    //            while (fit != interface_simplices.faces_end())
    //            {
    //                if (SimplicialComplex<MT>::get_intersection(f, *fit, e))
    //                    if (interface_simplices.contains(e))
    //                        break;
    //                ++fit;
    //            }
    //            if (fit != interface_simplices.faces_end())
    //            {
    //                f = *fit;
    //                interface_simplices.erase(e);
    //            }
    //            else
    //                break;
    //        } while(1);
    //
    //        if (interface_simplices.size_faces() > 0) return false;
    //
    //        return true;
    //    } // is_2manifold
    
    /**
     * Returns the minimum quality of the tetrahedra in simplex set s.
     */
    T min_quality(simplex_set& set)
    {
        T q_min = INFINITY;
        for (auto tit = set.tetrahedra_begin(); tit != set.tetrahedra_end(); tit++)
        {
            q_min = std::min(quality(*tit), q_min);
        }
        return q_min;
    }
    
    /**
     * Returns the minimum quality of the tetrahedra spanned by the vertices of the faces in s and pos.
     */
    T min_quality(simplex_set& set, const V& pos)
    {
        T min_q = INFINITY;
        std::vector<V> verts;
        for(auto fit = set.faces_begin(); fit != set.faces_end(); fit++)
        {
            get_pos(*fit, verts);
            min_q = std::min(min_q, std::abs(Util::quality<MT>(verts[0], verts[1], verts[2], pos)));
        }
        return min_q;
    }
    
    /**
     * Returns the minimum quality among tetrahedra in the star of n.
     */
    T min_quality(const node_key& n)
    {
        simplex_set st_n;
        Complex::star(n, st_n);
        return min_quality(st_n);
    }
    
    /**
     * Returns the minimum quality among tetrahedra in the star of each end node of the edge e.
     */
    T min_quality(const edge_key& e)
    {
        simplex_set st_e;
        Complex::star(e, st_e);
        return min_quality(st_e);
    }
    
    /**
     * Returns minimum quality among the tetrahedra adjacent to face f.
     */
    T min_quality(const face_key& f)
    {
        simplex_set st_f;
        Complex::star(f, st_f);
        return min_quality(st_f);
    }
    
    // The minimum quality among tetrahedra that would be produced by connecting apices with the edges on mf_boundary.
    T edge_remove_quality(simplex_set& mf_boundary, std::vector<V> & apices)
    {
        T q_min = 1e99;
        unsigned int n = static_cast<unsigned int>(mf_boundary.size_nodes());
        
        assert ( n > 2 || !"Too few vertices in the multi-face boundary!" );
        
        std::vector<node_key> polygon(n);
        
        sort_vertices(mf_boundary, polygon);
        check_consistency(apices, polygon);
        
        for (unsigned int i = 0; i < n; ++i)
        {
            T q = Util::quality<MT>(apices[1], apices[0], get_pos(polygon[i]), get_pos(polygon[(i+1)%n]));
            if (q < q_min) q_min = q;
        }
        
        return q_min;
    }
    
    /// Computes the difference in the volume enclosed by the interface after collapsing edge e.
    T volume_difference(edge_key & e, node_key & n0, node_key & n1)
    {
        T result = 0.0;
        simplex_set st0, st1, lk0, st_u, patch, cl_patch;
        
        Complex::star(n0, st0);
        Complex::link(n0, lk0);
        Complex::star(n1, st1);
        st0.insert(n0);
        st1.insert(n1);
        st_u.add(st0);
        st_u.add(st1);
        
        
        for (auto fit = st_u.faces_begin(); fit != st_u.faces_end(); fit++)
        {
            if (is_interface(e))
            {
                if (is_interface(*fit))
                    patch.insert(*fit);
            }
            else if (is_boundary(e))
            {
                if (is_boundary(*fit))
                    patch.insert(*fit);
            }
        }
        
        Complex::closure(patch, cl_patch);
        cl_patch.difference(st_u);
        cl_patch.difference(lk0);
        
        std::vector<V> pos;
        for (auto eit = cl_patch.edges_begin(); eit != cl_patch.edges_end(); eit++)
        {
            get_pos(*eit, pos);
            result += Util::volume<MT>(get_pos(n0), get_pos(n1), pos[0], pos[1]);
        }
        
        return result;
    }
    
private:
    
    /**
     * Check if the sequence of vertices in polygon is consistent with positive orientation of tetrahedra in the mesh
     * with respect to the ordered pair of vertices in vv. If not, reverse the order of vertices in polygon.
     */
    void check_consistency(std::vector<V> & vv, std::vector<node_key> & polygon)
    {
        unsigned int n = static_cast<unsigned int>(polygon.size());
        std::vector<V> vp(n);
        
        for (unsigned int i = 0; i < n; ++i)
        {
            vp[i] = get_pos(polygon[i]);
        }
        
        T sum = 0;
        for (unsigned int i = 0; i < n; ++i)
        {
            sum += Util::signed_volume<MT>(vv[0], vv[1], vp[(i+1)%n], vp[i]);
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
    
    /**
     * Check if the sequence of vertices in polygon is consistent with positive orientation of tetrahedra in the mesh
     * with respect to edge. If not, reverse the order of vertices in polygon.
     */
    void check_consistency(edge_key const & e, std::vector<node_key> & polygon)
    {
        std::vector<V> verts;
        get_pos(e, verts);
        std::vector<V> vv(2);
        vv[1] = verts[0]; // VERTEX FLIP
        vv[0] = verts[1];
        
        check_consistency(vv, polygon);
    }
    
    
    /// Check whether they are no inverted tetrahedra in the Complex::
    bool simplicial_complex_criterion_check()
    {
        for (auto tit = Complex::tetrahedra_begin(); tit != Complex::tetrahedra_end(); tit++)
        {
            if (0. > signed_volume(tit.key()))
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
            if (0. > signed_volume(*tit))
            {
                return true;
            }
        }
        return false;
    }
    
    /**
     * Returns whether any of the tetrahedra in the star of n is inverted.
     */
    bool inverted(const node_key& n)
    {
        simplex_set st_n;
        Complex::star(n, st_n);
        return inverted(st_n);
    }
    
    /**
     * Returns whether any of the tetrahedra in the simplex set will invert if node n is moved to p_new.
     */
    bool will_invert(const node_key& n, const V p_new, simplex_set& set)
    {
        V p = get_pos(n);
        std::vector<V> verts;
        for(auto fit = set.faces_begin(); fit != set.faces_end(); fit++)
        {
            get_pos(*fit, verts);
            T vol1 = Util::signed_volume<MT>(verts[0], verts[1], verts[2], p);
            T vol2 = Util::signed_volume<MT>(verts[0], verts[1], verts[2], p_new);
            if(Util::sgn<MT>(vol1) !=  Util::sgn<MT>(vol2))
            {
                return true;
            }
        }
        return false;
    }
    
    /**
     * Returns whether any of the tetrahedra in the star of n will invert if node n is moved to p_new.
     */
//    bool will_invert(const node_key& n, const V p_new)
//    {
//        simplex_set st_n, cl_st_n;
//        Complex::star(n, st_n);
//        Complex::closure(st_n, cl_st_n);
//        return will_invert(n, p_new, cl_st_n);
//    }
    
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
        
        for (auto eit = Complex::edges_begin(); eit != Complex::edges_end(); eit++)
        {
            valid = valid & Complex::exists(eit.key());
        }
        assert(valid);
        
        for (auto nit = Complex::nodes_begin(); nit != Complex::nodes_end(); nit++)
        {
            valid = valid & Complex::exists(nit.key());
        }
        assert(valid);
    }
    
    
    /// Remove the simplices that are no longer represented in the mesh from set.
    void cleanup_set(simplex_set& set)
    {
        simplex_set unused_simplices;
        
        for (auto nit = set.nodes_begin(); nit != set.nodes_end(); nit++)
        {
            if (!Complex::exists(*nit))
                unused_simplices.insert(*nit);
        }
        
        for (auto eit = set.edges_begin(); eit != set.edges_end(); eit++)
        {
            if (!Complex::exists(*eit))
                unused_simplices.insert(*eit);
        }
        
        for (auto fit = set.faces_begin(); fit != set.faces_end(); fit++)
        {
            if (!Complex::exists(*fit))
                unused_simplices.insert(*fit);
        }
        
        for (auto tit = set.tetrahedra_begin(); tit != set.tetrahedra_end(); tit++)
        {
            if (!Complex::exists(*tit))
            {
                unused_simplices.insert(*tit);
            }
        }
        
        set.difference(unused_simplices);
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
                    std::vector<node_key> nodes;
                    Complex::get_nodes(*eit, nodes);
                    
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
    void extract_interface(std::vector<typename MT::vector3_type> & verts, std::vector<int> & indices)
    {
        std::map<node_key, int> vert_index;
        
        // Extract vertices
        for (auto nit = Complex::nodes_begin(); nit != Complex::nodes_end(); nit++)
        {
            if (nit->is_interface())
            {
                verts.push_back(nit->v);
                vert_index[nit.key()] = verts.size();
            }
        }
        
        // Extract faces
        for (auto fit = Complex::faces_begin(); fit != Complex::faces_end(); fit++)
        {
            if (fit->is_interface())
            {
                std::vector<node_key> nodes(3);
                Complex::get(fit.key(), nodes);
                
                indices.push_back(vert_index[nodes[0]]);
                indices.push_back(vert_index[nodes[1]]);
                indices.push_back(vert_index[nodes[2]]);
            }
        }
    } // extract_interface
    
    void extract_tet_mesh(std::vector<typename MT::vector3_type> & points, std::vector< std::vector<int> > & tets)
    {
        std::map<node_key, int> indices;
        Complex::garbage_collect();
        
        int counter = 0;
        for (auto nit = Complex::nodes_begin(); nit != Complex::nodes_end(); nit++)
        {
            indices[nit.key()] = counter;
            points.push_back(nit->v);
            ++counter;
        }
        
        for (auto tit = Complex::tetrahedra_begin(); tit != Complex::tetrahedra_end(); tit++)
        {
            std::vector<node_key> nodes(4);
            std::vector<int> tet;
            Complex::get(tit.key(), nodes);
            
            for (int i = 0; i < nodes.size(); i++)
                tet.push_back(indices[nodes[i]]);
            
            tet.push_back(tit->label);
            tets.push_back(tet);
        }
    } // extract_tet_mesh
    
    /**
     * Returns the cosine to the dihedral angle between face f1 and face f2.
     */
    T cos_dihedral_angle(const face_key& f1, const face_key& f2)
    {
        V n0 = get_normal(f1);
        V n1 = get_normal(f2);
        T angle = MT::dot(n0, n1);
        assert(angle <= 1.);
        assert(angle >= -1.);
        return angle;
    }
    
    /**
     * Returns the dihedral angle between face f1 and face f2.
     */
    T dihedral_angle(const face_key& f1, const face_key& f2)
    {
        return acos(cos_dihedral_angle(f1, f2));
    }
    
    /**
     * Returns the minimum dihedral angle between the faces of tetrahedron t.
     */
    T min_cos_dihedral_angle(const tet_key& t)
    {
        T min_angle = -1.;
        simplex_set cl_t;
        Complex::closure(t, cl_t);
        for(auto fit1 = cl_t.faces_begin(); fit1 != cl_t.faces_end(); fit1++)
        {
            for(auto fit2 = fit1; fit2 != cl_t.faces_end(); fit2++)
            {
                if(*fit1 != *fit2)
                {
                    min_angle = std::max(min_angle, cos_dihedral_angle(*fit1, *fit2));
                }
            }
        }
        return min_angle;
    }
    
    T min_dihedral_angle(const tet_key& t)
    {
        return acos(min_cos_dihedral_angle(t));
    }
    
    void get_qualities(std::vector<int>& histogram, T& min_quality)
    {
        min_quality = INFINITY;
        
        histogram = std::vector<int>(180);
        for (int i = 0; i < 100; ++i)
        {
            histogram[i] = 0;
        }
        
        for (auto tit = Complex::tetrahedra_begin(); tit != Complex::tetrahedra_end(); tit++)
        {
            T q = quality(tit.key());
            min_quality = std::min(min_quality, q);
            histogram[(int)floor(q*100.)] += 1;
        }
    }
    
    /**
     * Calculates the dihedral angles in the SimplicialComplex and returns these in a histogram,
     * along with the minimum and maximum dihedral angles.
     */
    void get_dihedral_angles(std::vector<int> & histogram, T & min_angle, T & max_angle)
    {
        max_angle = -INFINITY, min_angle = INFINITY;
        
        histogram = std::vector<int>(180);
        for (int i = 0; i < 180; ++i)
        {
            histogram[i] = 0;
        }
        
        for (auto tit = Complex::tetrahedra_begin(); tit != Complex::tetrahedra_end(); tit++)
        {
            simplex_set cl_t;
            Complex::closure(tit.key(), cl_t);
            
            for(auto fit1 = cl_t.faces_begin(); fit1 != cl_t.faces_end(); fit1++)
            {
                for(auto fit2 = fit1; fit2 != cl_t.faces_end(); fit2++)
                {
                    if(*fit1 != *fit2)
                    {
                        T angle = dihedral_angle(*fit1, *fit2)*180./M_PI;
                        min_angle = std::min(min_angle, angle);
                        max_angle = std::max(max_angle, angle);
                        histogram[(int)floor(angle)] += 1;
                    }
                }
            }
        }
    }
    
    T min_quality()
    {
        T min_q;
        for (auto tit = Complex::tetrahedra_begin(); tit != Complex::tetrahedra_end(); tit++)
        {
            min_q = std::min(min_q, quality(tit.key()));
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
            if (tit->label != 0)
                object++;
        }
    }
    
};
        
        // SIMPLICIAL_COMPLEX_H
#endif
        
