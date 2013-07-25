#ifndef DSC_H
#define DSC_H

#include "is_mesh_API.h"
#include "util.h"
#include "printing.h"

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
    
    // Thresholds on the dihedral angle between two neighbouring faces
    T DEG_ANGLE;
    T MIN_ANGLE;
    
    // Thresholds on the quality of tetrahedra
    T DEG_TET_QUALITY;
    T MIN_TET_QUALITY;
    
    // Thresholds on the volume of tetrahedra
    T DEG_TET_VOLUME;
    T MIN_TET_VOLUME;
    T MAX_TET_VOLUME;
    
    // Thresholds on the length of edges
    T DEG_EDGE_LENGTH;
    T MIN_EDGE_LENGTH;
    T MAX_EDGE_LENGTH;
    
    // User defined parameters
    T AVG_EDGE_LENGTH;
    T MIN_DEFORMATION;

    T FLIP_EDGE_INTERFACE_FLATNESS;
        
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
        
        DEG_EDGE_LENGTH = 0.1 * AVG_EDGE_LENGTH;
        MIN_EDGE_LENGTH = 0.5 * AVG_EDGE_LENGTH;
        MAX_EDGE_LENGTH = 2. * AVG_EDGE_LENGTH;
        
        DEG_ANGLE = 5.*M_PI/180.;
        MIN_ANGLE = 10.*M_PI/180.;
        
        MIN_DEFORMATION = 0.25 * AVG_EDGE_LENGTH;
        
        DEG_TET_QUALITY = 0.01;
        MIN_TET_QUALITY = 0.3;
        FLIP_EDGE_INTERFACE_FLATNESS = 0.995;
        
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
    
    template<typename key>
    bool is_crossing(const key& k)
    {
        return Complex::is_crossing(k);
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
    
    T get_min_angle() const
    {
        return MIN_ANGLE;
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
    
    void get_half_polygon(const edge_key& e, std::vector<node_key>& polygon1, std::vector<node_key>& polygon2)
    {
        std::vector<node_key> nodes;
        Complex::get_nodes(e, nodes);
        std::vector<node_key> polygon = get_polygon(e);
        int i1 = -1, i2 = -1;
        for(int i = 0; i < polygon.size(); i++)
        {
            face_key f = Complex::get_face(nodes[0], nodes[1], polygon[i]);
            assert(f != Complex::NULL_FACE);
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
        assert(i2 != -1);
        
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
        assert(polygon1.size() > 2);
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
        T q_new = build_table(e, polygon, K);
        
        if (q_new > min_quality(e))
        {
            std::vector<node_key> nodes;
            Complex::get_nodes(e, nodes);
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
        
        if(m2 > 2)
        {
            k = K2[0][m2-1];
            flip_23_recursively(polygon2, n1, n2, K2, 0, k);
            flip_23_recursively(polygon2, n1, n2, K2, k, m2-1);
        }
        
        // Find the faces to flip about.
        face_key f1 = Complex::get_face(n1, n2, polygon1.front());
        face_key f2 = Complex::get_face(n1, n2, polygon1.back());
        
        if(m2 <= 2) {
            Complex::flip_22(f1, f2);
        }
        else {
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
        T q_new = build_table(e, polygon1, K1);
        
        if(polygon2.size() > 2)
        {
            q_new = std::min(q_new, build_table(e, polygon2, K2));
        }
        
        if (q_new > min_quality(e))
        {
            std::vector<node_key> nodes;
            Complex::get_nodes(e, nodes);
            topological_boundary_edge_removal(polygon1, polygon2, nodes[0], nodes[1], K1, K2);
            return true;
        }
        return false;
    }
    
    /**
     * Returns whether it is possible to flip the edge e or not, i.e. whether the edge is not a feature edge 
     * (it is not a feature edge if its neighborhood is sufficiently flat).
     */
    bool is_topological_removable(const edge_key & e)
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
        assert(faces.size() == 2);
        
        T angle = cos_dihedral_angle(e, faces[0], faces[1]);
        if(angle > FLIP_EDGE_INTERFACE_FLATNESS)
        {
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
        int i = 0, j = 0;
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
                        if(!is_interface(*eit) && !is_boundary(*eit))
                        {
                            if(topological_edge_removal(*eit))
                            {
                                i++;
                            }
                        }
                        else if(is_topological_removable(*eit))
                        {
                            if(topological_boundary_edge_removal(*eit))
                            {
                                i++;
                            }
                        }
                        j++;
                    }
                }
            }
        }
        std::cout << "Topological edge removals: " << i << "/" << j << std::endl;
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
    
    /**
     * Attempt to remove the faces sandwiched between the apices of f using multi-face removal. The face f is used as a starting point.
     */
    bool topological_face_removal(const face_key& f)
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
    * Attempt to remove the faces sandwiched between the nodes apex1 and apex2 using multi-face removal. 
    * The face which intersects with the line segment |apex1 apex2| is used as a starting point.
    */
    bool topological_face_removal(const node_key& apex1, const node_key& apex2)
    {        
        simplex_set lk_n1, lk_n2;
        Complex::link(apex1, lk_n1);
        Complex::link(apex2, lk_n2);
        lk_n1.intersection(lk_n2);
        for(auto f = lk_n1.faces_begin(); f != lk_n1.faces_end(); f++)
        {
            if(!is_boundary(*f) && !is_interface(*f))
            {
                std::vector<node_key> nodes;
                Complex::get_nodes(*f, nodes);
                orient_cc(apex2, nodes);
                
                T t = Util::intersection_ray_triangle<MT>(get_pos(apex1), get_pos(apex2), get_pos(nodes[0]), get_pos(nodes[1]), get_pos(nodes[2]));
                if(0. < t && t < 1.)
                {
                    if(topological_face_removal(*f))
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
                    if (Complex::exists(*fit) && !is_interface(*fit) && !is_boundary(*fit))
                    {
                        std::vector<node_key> apices;
                        Complex::get_apices(*fit, apices);
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
            if (is_interface(eit.key()) && length(eit.key()) > MAX_EDGE_LENGTH)
            {
                edges.push_back(eit.key());
            }
        }
        int i = 0, j = 0;
        for(auto &e : edges)
        {
            if (Complex::exists(e) && is_interface(e) && length(e) > MAX_EDGE_LENGTH)
            {
                if(split(e) != Complex::NULL_NODE)
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
        std::cout << "Thinning collapses: " << i << "/" << j << std::endl;
    }
    
    /////////////////////////
    // REMOVE DEGENERACIES //
    /////////////////////////
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
     * Attempt to remove edges shorter than MIN_EDGE_LENGTH by safely collapsing them.
     */
    void remove_edges()
    {
        std::list<edge_key> edges;
        for (auto eit = Complex::edges_begin(); eit != Complex::edges_end(); eit++)
        {
            if (length(eit.key()) < MIN_EDGE_LENGTH)
            {
                edges.push_back(eit.key());
            }
        }
        int i = 0, j = 0;
        for(auto e : edges)
        {
            if(Complex::exists(e) && length(e) < MIN_EDGE_LENGTH)
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
        std::vector<V> verts;
        get_pos(e, verts);
        V p = Util::project<MT>(get_pos(apex), verts[0], verts[1]);
        
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
        if(max_angle(f) > M_PI - MIN_ANGLE)
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
            if(min_angle(fit.key()) < MIN_ANGLE)
            {
               faces.push_back(fit.key());
            }
        }
        
        int i = 0, j = 0;
        for (auto &f : faces)
        {
            if (Complex::exists(f) && min_angle(f) < MIN_ANGLE)
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
        std::vector<V> verts;
        get_pos(f, verts);
        V p = Util::project<MT>(get_pos(apex), verts);
        
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
    
    void smooth()
    {
        int i = 0, j = 0;
        for (auto nit = Complex::nodes_begin(); nit != Complex::nodes_end(); nit++)
        {
            if (Complex::exists(nit.key()) && !is_boundary(nit.key()) && !is_interface(nit.key()))
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
                Complex::set_label(t, label);
                return true;
            }
        }
        return false;
    }
    
    /**
     * Relabel all tetrahedra in the SimplicialComplex mesh when
     * a) their quality is lower than TET_RELABEL_THRESHOLD and the minimum sine of their dihedral angles is lower than TET_RELABEL_THRESHOLD
     * b) they fulfill criteria for relabelling (described in detail in relabel_condition()).
     * Compliant with multiple phases.
     */
    bool relabel_tets()
    {
        int i = 0, j = 0;
        for (auto tit = Complex::tetrahedra_begin(); tit != Complex::tetrahedra_end(); tit++)
        {
            if(relabel_tet(tit.key()))
            {
                i++;
            }
            j++;
        }
        std::cout << "Relabeled tetrahedra: " << i << "/" << j << std::endl;
        
        return (i > 0);
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
        
        print_out("Topological removals.");
        topological_edge_removal();
        validity_check();
        topological_face_removal();
        validity_check();
        
//        print_out("Low quality removal.");
//        remove_tets();
//        validity_check();
//        remove_faces();
//        validity_check();
//        remove_edges();
//        validity_check();
        
        print_out("Degeneracy removal.");
        remove_degenerate_tets();
        validity_check();
        remove_degenerate_faces();
        validity_check();
        remove_degenerate_edges();
        validity_check();
        
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
        print_out("Thickening interface pass.");
        thickening_interface();
        validity_check();
        
//        print_out("Thinning interface pass.");
//        thinning_interface();
//        validity_check();
        
//        print_out("Thickening pass.");
//        thickening();
//        validity_check();
        
        print_out("Thinning pass.");
        thinning();
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
                if (Complex::exists(nit.key()) && !nit->is_crossing() && nit->is_interface())
                {
                    if(!move_vertex(nit.key()))
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
        
        resize_complex();
        
        Complex::garbage_collect();
        ++step_no;
    }
    
private:
    
    /**
     * Tries moving the node n to the new position new_pos. Returns true if it succeeds.
     */
    bool move_vertex(const node_key & n)
    {
        V pos = get_pos(n);
        V new_pos = get_destination(n);
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
    }
    
    
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
        std::vector<node_key> nodes;
        for (auto nit = ln_e.nodes_begin(); nit != ln_e.nodes_end(); nit++)
        {
            if (is_interface(*nit) || is_boundary(*nit))
            {
                verts.push_back(get_pos(*nit));
                nodes.push_back(*nit);
            }
        }
        
        if (verts.size() != 4)
        {
            return false;
        }
        if(Complex::get_edge(nodes[0], nodes[1]) != Complex::NULL_EDGE) // Check that there does not already exist an edge between the pair of nodes.
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
        
        if (flatness >= FLIP_EDGE_INTERFACE_FLATNESS)
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
        
        node_key n_new = split(e);
        
        edge_key e_c = Complex::get_edge(nodes[0], n_new);
        assert(e_c != Complex::NULL_EDGE);
        V p = get_pos(nodes[0]);
        V p_new = get_destination(nodes[0]);
        
        if(precond_collapse(e_c, p))
        {
            node_key n = collapse(e_c, p, p_new);
            return n != Complex::NULL_NODE;
        }
        return false;
    }
    
    /**
     * Flips the edge e (which is a special case of the edge remove operation in the embedding mesh).
     * Relabels the tetrahedra accordingly so that the interface mesh geometry does not change.
     */
    bool border_flip(const edge_key & e)
    {
        if(!is_interface(e) && !is_boundary(e))
        {
            return false;
        }
        // Find the pair of nodes to be connected by a new edge as a result of edge flip.
        simplex_set st_e;
        Complex::star(e, st_e);
        std::vector<face_key> faces;
        std::cout << std::endl;
        for (auto fit = st_e.faces_begin(); fit != st_e.faces_end(); fit++)
        {
            if (is_interface(*fit) || is_boundary(*fit))
            {
                std::cout << "Face: " << *fit << std::endl;
                std::vector<V> verts;
                get_pos(*fit, verts);
                for(auto &v : verts)
                {
                    std::cout << v << std::endl;
                }
                faces.push_back(*fit);
            }
        }
        
        if(faces.size() != 2)
        {
            return false;
        }
        
        node_key n = Complex::flip_44(faces[0], faces[1]);
        validity_check();
        return n != Complex::NULL_NODE;
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
                std::vector<node_key> nodes;
                Complex::get_nodes(fit.key(), nodes);
                
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
            std::vector<node_key> nodes;
            std::vector<int> tet;
            Complex::get_nodes(tit.key(), nodes);
            
            for (int i = 0; i < nodes.size(); i++)
                tet.push_back(indices[nodes[i]]);
            
            tet.push_back(tit->label);
            tets.push_back(tet);
        }
    } // extract_tet_mesh
    
    /**
     * Returns the cosine to the dihedral angle between face f1 and face f2 at edge e.
     */
    T cos_dihedral_angle(const edge_key& e, const face_key& f1, const face_key& f2)
    {
        std::vector<V> verts;
        get_pos(e, verts);
        V c = get_pos(Complex::get_apex(f1, e));
        V d = get_pos(Complex::get_apex(f2, e));
        return Util::cos_dihedral_angle<MT>(verts[0], verts[1], c, d);
    }
    
    /**
     * Returns the dihedral angle between face f1 and face f2.
     */
    T dihedral_angle(const edge_key& e, const face_key& f1, const face_key& f2)
    {
        return acos(cos_dihedral_angle(e, f1, f2));
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
                    min_angle = std::max(min_angle, cos_dihedral_angle(Complex::get_edge(*fit1, *fit2),*fit1, *fit2));
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
                        T angle = dihedral_angle(Complex::get_edge(*fit1, *fit2),*fit1, *fit2)*180./M_PI;
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
        
