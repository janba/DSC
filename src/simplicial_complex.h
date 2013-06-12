#ifndef SIMPLICIAL_COMPLEX_H
#define SIMPLICIAL_COMPLEX_H

/** 2013 Marek K. Misztal and Mark Viinblad Jensen **/

#include <is_mesh/is_mesh.h>
#include <is_mesh/io/is_mesh_lists_read.h>

#include "attributes.h"
#include "util.h"

#include "printing.h"

/**
 * Data structure used by klincsek_triangulation() (indices for triangle's vertices and triangle's quality).
 */
template <typename MT>
struct kt_chunk
{
    typename MT::real_type quality;
    int i;
    int j;
    int k;
};

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
class SimplicialComplex
{
    typedef typename MT::real_type      T;
    typedef typename MT::vector3_type   V;
    typedef typename MT::vector4_type   V4;
    
    typedef typename OpenTissue::is_mesh::t4mesh< default_node_traits<MT>, default_tetrahedron_traits<MT>, default_edge_traits<MT>, default_face_traits<MT> > Mesh;
    
public:
    
    typedef default_node_traits<MT>        node_traits;
    typedef default_edge_traits<MT>        edge_traits;
    typedef default_face_traits<MT>        face_traits;
    typedef default_tetrahedron_traits<MT> tetrahedron_traits;
    
    typedef typename Mesh::node_key_type    node_key_type;
    typedef typename Mesh::edge_key_type    edge_key_type;
    typedef typename Mesh::face_key_type    face_key_type;
    typedef typename Mesh::tetrahedron_key_type tetrahedron_key_type;
    
    const node_key_type NULL_NODE;
    const edge_key_type NULL_EDGE;
    const face_key_type NULL_FACE;
    const tetrahedron_key_type NULL_TETRAHEDRON;
    
    typedef typename Mesh::simplex_set_type simplex_set;
    
    typedef std::map<node_key_type, V> movement_map;
    
private:
    Mesh mesh;
    
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
    SimplicialComplex(std::vector<T> & points, std::vector<int> & tets, std::vector<int> & tet_labels):
    NULL_NODE(-1), NULL_EDGE(-1), NULL_FACE(-1), NULL_TETRAHEDRON(-1)
    {
        vectors_read(points, tets, mesh);
        step_no = 0;
        init(tet_labels);
        
        T ie_min, ie_avg;
        calc_interface_edge_length(ie_min, ie_avg);
        AVG_EDGE_LENGTH = ie_avg;
        
        MIN_EDGE_LENGTH = 0.5 * AVG_EDGE_LENGTH;
        
        MIN_DEFORMATION = 0.25 * AVG_EDGE_LENGTH;
        
        MIN_DIHEDRAL_ANGLE = 5.*M_PI/180.;
        DEG_TET_QUALITY = 0.01;
        MIN_TET_QUALITY = 0.3;
        FLIP_EDGE_TET_FLATNESS = 0.995;
        
        T vol_avg = AVG_EDGE_LENGTH*sqrt(2.)/12;
        MIN_TET_VOLUME = 0.25*vol_avg;
        MAX_TET_VOLUME = 2.*vol_avg;
        
        check_validity();
        check_validity();
    }
    
    SimplicialComplex() {}
    
private:    
    
    /// Calculates the minimum and average length of edges on the interface.
    void calc_interface_edge_length(T& ie_min, T& ie_avg)
    {
        ie_min = INFINITY;
        ie_avg = 0.;
        int cnt = 0;
        for(auto eit = edges_begin(); eit != edges_end(); eit++)
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
    
    /**
     * Label all tetrahedra according to tet_labels and perform an initial update
     * of flags and attributes of all simplices
     */
    void init(std::vector<int> & tet_labels)
    {
        // Update and label all tetrahedra
        for (auto tit = tetrahedra_begin(); tit != tetrahedra_end(); tit++)
        {
            tit->label = tet_labels[tit.key()];
        }
        
        // Update all faces
        for (auto fit = faces_begin(); fit != faces_end(); fit++)
        {
            update_flag(fit.key());
        }
        
        // Update all edges
        for (auto eit = edges_begin(); eit != SimplicialComplex<MT>::edges_end(); eit++)
        {
            update_flag(eit.key());
        }
        
        // Update all nodes
        for (auto nit = nodes_begin(); nit != nodes_end(); nit++)
        {
            update_flag(nit.key());
        }
    }
    
    //////////////////////
    // UPDATE FUNCTIONS //
    //////////////////////
    /**
     * Updates the flags (is interface, is boundary, is locked) of simplices in set.
     */
    void update(simplex_set & set)
    {
        // Update faces
        for (auto fit = set.faces_begin(); fit != set.faces_end(); fit++)
        {
            if (mesh.exists(*fit))
            {
                update_flag(*fit);
            }
        }
        
        // Update edges
        for (auto eit = set.edges_begin(); eit != set.edges_end(); eit++)
        {
            if (mesh.exists(*eit))
            {
                update_flag(*eit);
            }
        }
        
        // Update nodes
        for (auto nit = set.nodes_begin(); nit != set.nodes_end(); nit++)
        {
            if (mesh.exists(*nit))
            {
                update_flag(*nit);
            }
        }
    } // update
    
    void update_flag(const face_key_type & f)
    {
        get(f).set_interface(false);
        get(f).set_boundary(false);
        get(f).set_locked(false);
        
        simplex_set st_f;
        SimplicialComplex<MT>::star(f, st_f);
        if (st_f.size_tetrahedra() == 1)
        {
            // On the boundary
            get(f).set_boundary(true);
            get(f).set_locked(true);
            if (get(*(st_f.tetrahedra_begin())).label != 0)
            {
                get(f).set_interface(true);
            }
        }
        else if(st_f.size_tetrahedra() == 2)
        {
            auto tit = st_f.tetrahedra_begin();
            int label0 = get(*tit).label;   ++tit;
            int label1 = get(*tit).label;
            if (label0 != label1)
            {
                // On the interface
                get(f).set_interface(true);
                get(f).set_locked(true);
            }
        }
    } // update_flag(face_key_type)
    
    void update_flag(const edge_key_type & e)
    {
        get(e).set_boundary(false);
        get(e).set_interface(false);
        get(e).set_locked(false);
        
        simplex_set ste;
        SimplicialComplex<MT>::star(e, ste);
        
        for (auto efit = ste.faces_begin(); efit != ste.faces_end(); efit++)
        {
            if (mesh.exists(*efit))
            {
                if (get(*efit).is_boundary())
                {
                    get(e).set_boundary(true);
                    get(e).set_locked(true);
                }
                else if (get(*efit).is_interface())
                {
                    get(e).set_interface(true);
                    get(e).set_locked(true);
                }
            }
        }
    } // update_flage(edge_key_type)
    
    void update_flag(const node_key_type & n)
    {
        get(n).set_interface(false);
        get(n).set_boundary(false);
        get(n).set_locked(false);
        
        simplex_set st_n;
        SimplicialComplex<MT>::star(n, st_n);
        for (auto neit = st_n.edges_begin(); neit != st_n.edges_end(); neit++)
        {
            if (mesh.exists(*neit))
            {
                if (get(*neit).is_interface())
                {
                    get(n).set_interface(true);
                    get(n).set_locked(true);
                }
                if (get(*neit).is_boundary())
                {
                    get(n).set_boundary(true);
                    get(n).set_locked(true);
                }
            }
        }
    } // update_flag(node_key_type)
    
    ////////////////////
    // MESH FUNCTIONS //
    ////////////////////
public:
    typename Mesh::node_iterator nodes_begin()
    {
        return mesh.nodes_begin();
    }
    
    typename Mesh::node_iterator nodes_end()
    {
        return mesh.nodes_end();
    }
    
    typename Mesh::edge_iterator edges_begin()
    {
        return mesh.edges_begin();
    }
    
    typename Mesh::edge_iterator edges_end()
    {
        return mesh.edges_end();
    }
    
    typename Mesh::face_iterator faces_begin()
    {
        return mesh.faces_begin();
    }
    
    typename Mesh::face_iterator faces_end()
    {
        return mesh.faces_end();
    }
    
    typename Mesh::tetrahedron_iterator tetrahedra_begin()
    {
        return mesh.tetrahedra_begin();
    }
    
    typename Mesh::tetrahedron_iterator tetrahedra_end()
    {
        return mesh.tetrahedra_end();
    }
    
    
    typename Mesh::node_type & get(const node_key_type& k)
    {
        return mesh.find_node(k);
    }
    
    typename Mesh::edge_type & get(const edge_key_type& k)
    {
        return mesh.find_edge(k);
    }
    
    typename Mesh::face_type & get(const face_key_type& k)
    {
        return mesh.find_face(k);
    }
    
    typename Mesh::tetrahedron_type & get(const tetrahedron_key_type& k)
    {
        return mesh.find_tetrahedron(k);
    }
    
    void get_nodes(const edge_key_type& e, std::vector<node_key_type>& nodes)
    {
        nodes = std::vector<node_key_type>(2);
        mesh.vertices(e, nodes);
    }
    
    void get_nodes(const face_key_type& f, std::vector<node_key_type>& nodes)
    {
        nodes = std::vector<node_key_type>(3);
        mesh.vertices(f, nodes);
    }
    
    void get_nodes(const tetrahedron_key_type& t, std::vector<node_key_type>& nodes)
    {
        nodes = std::vector<node_key_type>(4);
        mesh.vertices(t, nodes);
    }
    
    edge_key_type get_edge(const node_key_type& n1, const node_key_type& n2)
    {
        simplex_set st1, st2;
        star(n1, st1);
        star(n2, st2);
        st1.intersection(st2);
        
        if (st1.size_edges() != 1)
        {
            return NULL_EDGE;
        }
        return *(st1.edges_begin());
    }
    
    edge_key_type get_edge(const face_key_type& f1, const face_key_type& f2)
    {
        simplex_set cl1, cl2;
        closure(f1, cl1);
        closure(f2, cl2);
        cl1.intersection(cl2);
        
        if (cl1.size_edges() != 1)
        {
            return NULL_EDGE;
        }
        return *(cl1.edges_begin());
    }
    
    void get_apices(const face_key_type& f, std::vector<node_key_type>& apices)
    {
        apices = std::vector<node_key_type>(0);
        simplex_set lk_f;
        link(f, lk_f);
        assert(lk_f.size_nodes() == 2);
        for(auto nit = lk_f.nodes_begin(); nit != lk_f.nodes_end(); nit++)
        {
            apices.push_back(*nit);
        }
    }
    
    template<typename Key>
    void star(const Key& k, simplex_set& set)
    {
        mesh.star(k, set);
    }
    
    template<typename Key>
    void link(const Key& k, simplex_set& set)
    {
        mesh.link(k, set);
    }
    
    template<typename Key>
    void closure(const Key& k, simplex_set& set)
    {
        mesh.closure(k, set);
    }
    
    //    void closure_set(const simplex_set& set_, simplex_set& set)
    //    {
    //        mesh.closure(set_, set);
    //    }
    
    //    bool is_editable(const node_key_type& n)
    //    {
    //        return mesh.exists(n);
    //    }
    
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
    
    //////////////////////
    // TOPOLOGICAL PASS //
    //////////////////////
    
    /**
     * Finds the triangulation of the input polygon "sandwiched" between two vertices in vv,
     * which maximizes the worst tetrahedron's quality (using dynamic programming -- Klincsek's algorithm).
     */
    T klincsek_triangulation(const std::vector<node_key_type>& polygon, std::vector<node_key_type>& new_edges, std::vector<V>& vv)
    {
        // Build an array for the dynamic programming method
        T max_double = INFINITY;
        T q;
        int n = (int) polygon.size();
        
        std::vector<std::vector<kt_chunk<MT> > > kt_array;
        kt_array.resize(n-1);
        
        for (int i=0; i<n-1; ++i)
        {
            kt_array[i].resize(n);
            kt_array[i][i+1].quality = max_double;
            kt_array[i][i+1].i = i;
            kt_array[i][i+1].j = i+1;
        }
        
        for (int i=n-3; i>=0; --i)
        {
            for (int j=i+2; j<n; ++j)
            {
                for (int k=i+1; k<j; ++k)
                {
                    q = Util::quality<MT>(get_pos(polygon[i]), get_pos(polygon[k]), get_pos(polygon[j]), vv[0], vv[1]);
                    if (k<j-1)
                    {
                        q = std::min(q, kt_array[k][j].quality);
                    }
                    if (k>i+1)
                    {
                        q = std::min(q, kt_array[i][k].quality);
                    }

                    if (k==i+1 || q>kt_array[i][j].quality)
                    {
                        kt_array[i][j].quality = q;
                        kt_array[i][j].k = k;
                        kt_array[i][j].i = i;
                        kt_array[i][j].j = j;
                    }
                }
            }
        }
        
        // Find optimal triangulation of the polygon
        std::list<kt_chunk<MT> > kt_stack;
        kt_stack.push_front(kt_array[0][kt_array[0][n-1].k]);
        kt_stack.push_front(kt_array[kt_array[0][n-1].k][n-1]);
        
        q = kt_array[0][n-1].quality;
        
        while (!kt_stack.empty())
        {
            kt_chunk<MT> kt = kt_stack.front();
            kt_stack.pop_front();
            
            if (kt.j - kt.i != 1)
            {
                new_edges.push_back(polygon[kt.i]);
                new_edges.push_back(polygon[kt.j]);
            
                kt_stack.push_front(kt_array[kt.i][kt.k]);
                kt_stack.push_front(kt_array[kt.k][kt.j]);
            }
        }
        
        return q;
    }
    
    /**
     * Attempt to remove edge e (the topological operation -- mesh reconnection).
     */
    void edge_removal(edge_key_type const & e)
    {
        simplex_set st_e;
        star(e, st_e);

        T min_q = min_quality(st_e);
        int label = mesh.find_tetrahedron(*(st_e.tetrahedra_begin())).label;
                
        std::vector<V> vv;
        get_pos(e, vv);
        std::vector<V> verts(2);
        verts[0] = vv[1];// VERTEX FLIP (due to strange choice in vertices())
        verts[1] = vv[0];
        
        simplex_set lk_e;
        link(e, lk_e);
        std::vector<node_key_type> polygon, new_edges;
        
        sort_vertices(lk_e, polygon);
        check_consistency(e, polygon);
        
        T q = klincsek_triangulation(polygon, new_edges, verts);
        
        if (q > min_q && q > 0.)
        {
            simplex_set new_simplices;
            mesh.remove_edge(e, new_edges, new_simplices);
            
            for (auto tit = new_simplices.tetrahedra_begin(); tit != new_simplices.tetrahedra_end(); tit++)
            {
                mesh.find_tetrahedron(*tit).label = label;
            }
            
            simplex_set cl_ns;
            mesh.closure(new_simplices, cl_ns);
            update(cl_ns);
        }
    }
    
    /**
     * The helper function for optimal_multi_face_remove().
     * TODO: Sanity check.
     */
    T omfr_quality_helper(int n, std::vector<face_key_type> & faces, std::map<face_key_type, T> & face_quality_values, std::vector<bool> & subset, std::vector<V> & vv)
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
        
        mesh.boundary_2manifold(multi_face, mf_boundary);
        T q = edge_remove_quality(mf_boundary, vv);
        
        if (q < quality) quality = q;
        return quality;
    }
    
    
    /**
     * Returns the boundary of the multi-face in mf_boundary.
     */
    void multi_face_boundary(simplex_set& multi_face, simplex_set& mf_boundary)
    {
        simplex_set cl_mf, e_mf;
        mesh.closure(multi_face, cl_mf);
        for(auto eit = cl_mf.edges_begin(); eit != cl_mf.edges_end(); eit++)
        {
            simplex_set st_e;
            star(*eit, st_e);
            int i = 0;
            for(auto fit = st_e.faces_begin(); fit != st_e.faces_end(); fit++)
            {
                if(cl_mf.contains(*fit))
                {
                    i++;
                }
            }
            assert(i <= 2);
            if(i == 1)
            {
                e_mf.insert(*eit);
            }
        }
        mesh.closure(e_mf, mf_boundary);
    }
    
    /**
     * Returns the faces in the set candidate which is connected to the face f. Connected means there is no interface simplices between f and the candidate face and all faces between f and the candidate face is in the candidate set.
     */
    void connected_component(const face_key_type& f, simplex_set& candidate, simplex_set& cc)
    {
        cc.insert(f);
        simplex_set cl_f, st_cl_f;
        closure(f, cl_f);
        mesh.star(cl_f, st_cl_f);
        
        for (auto fit = st_cl_f.faces_begin(); fit != st_cl_f.faces_end(); fit++)
        {
            if(mesh.exists(*fit) && candidate.contains(*fit) && !cc.contains(*fit) && !get(*fit).is_interface())
            {
                edge_key_type e = get_edge(*fit, f);
                if(e != NULL_EDGE && mesh.exists(e) && !get(e).is_interface())
                {
                    connected_component(*fit, candidate, cc);
                }
            }
        }
    }
    
    /**
     * Finds the biggest multi-face containing f, which can be removed.
     */
    void removable_multi_face(const face_key_type& f, simplex_set& multi_face)
    {
        std::vector<node_key_type> apices;
        get_apices(f, apices);
        
        simplex_set link0, link1;
        link(apices[0], link0);
        link(apices[1], link1);
        link0.intersection(link1);
        
        connected_component(f, link0, multi_face);
    }
    
    /**
     * Attempt multi-face retriangulation on the initial set given in removed-faces.
     */
    void optimal_multi_face_retriangulation(const face_key_type& f, simplex_set & new_simplices)
    {
        simplex_set multi_face;
        removable_multi_face(f, multi_face);
        
        simplex_set mf_boundary;
        multi_face_boundary(multi_face, mf_boundary);
        
        if (mf_boundary.size_nodes() != mf_boundary.size_edges()) // if mf_boundary is not homomorphic to a circle/segment
        {
            simplex_set new_mf;
            mesh.find_min_multi_face(f, multi_face, new_mf);
            multi_face.intersection(new_mf);
            
            mf_boundary.clear();
            multi_face_boundary(multi_face, mf_boundary);
        }
        
        if (multi_face.size_faces() > 1)
        {
            std::vector<node_key_type> apices;
            get_apices(f, apices);
            
            std::vector<V> verts(2);
            verts[0] = get_pos(apices[0]);
            verts[1] = get_pos(apices[1]);
            
            simplex_set st_mf;
            mesh.star(multi_face, st_mf);
            T q_old = min_quality(st_mf);
            T q = -1.;
            
            std::vector<node_key_type> polygon, new_edges;
            sort_vertices(mf_boundary, polygon);
            check_consistency(verts, polygon);
            
            q = klincsek_triangulation(polygon, new_edges, verts);
            
            // If Klincsek's algorithm found a better triangulation, proceed with multi-face retriangulation
            if (q > q_old + EPSILON && q > 0.)
            {
                simplex_set new_multi_face;
                mesh.multi_face_retriangulation(multi_face, new_edges, new_multi_face, new_simplices);
                multi_face = new_multi_face;
                
                assert (new_simplices.size_tetrahedra() == 2 * multi_face.size_faces());
            }
        }
    }
    
    /**
     * Attempt to remove face f using multi-face retriangulation.
     */
    void multi_face_retriangulation(const face_key_type& f)
    {
        simplex_set st_f;
        star(f, st_f);
        int label = get(*st_f.tetrahedra_begin()).label;
        
        // Attempt to remove the multi-face by multi-face retriangulation
        simplex_set new_simplices;
        optimal_multi_face_retriangulation(f, new_simplices);
        
        for (auto tit = new_simplices.tetrahedra_begin(); tit != new_simplices.tetrahedra_end(); tit++)
        {
            get(*tit).label = label;
        }
        
        simplex_set ns_cl;
        mesh.closure(new_simplices, ns_cl);
        update(ns_cl);
    }
    
    
    /**
     * Attempt multi-face removal (a reconnection method), which finds the optimal subset (via dynamic programming)
     * of the multi-face to be swapped for an edge connecting the apices (in vv).
     * TODO: Sanity-check
     */
    void optimal_multi_face_removal(const face_key_type& f, simplex_set& new_simplices)
    {
        std::vector<node_key_type> apices;
        get_apices(f, apices);
        
        simplex_set multi_face;
        removable_multi_face(f, multi_face);
        
        if(multi_face.size_faces() < 2)
        {
            return;
        }
        
        simplex_set st_mf;
        mesh.star(multi_face, st_mf);
        T q_old = min_quality(st_mf);
        
        std::vector<V> verts(2);
        verts[0] = get_pos(apices[0]);
        verts[1] = get_pos(apices[1]);
        
        unsigned int n = static_cast<unsigned int>(multi_face.size_faces());
        std::map<face_key_type, double> face_quality_values;
        std::vector<face_key_type> faces(n);
        
        int i = 0;
        face_key_type seed_face;
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
                mesh.boundary_2manifold(multi_face, mf_boundary);
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
                                    typename SimplicialComplex<MT>::edge_key_type ek;
                                    if (mesh.get_intersection(faces[k], faces[l], ek))
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
            
            mesh.multi_face_remove(multi_face, new_simplices);
            cleanup_set(new_simplices);
        }
    }
    
    /**
     * Attempt to remove face f using multi-face removal.
     */
    void multi_face_removal(const face_key_type& f)
    {
        simplex_set st_f;
        star(f, st_f);
        int label = get(*st_f.tetrahedra_begin()).label;
        
        simplex_set new_simplices;
        optimal_multi_face_removal(f, new_simplices);
        
        for (auto tit = new_simplices.tetrahedra_begin(); tit != new_simplices.tetrahedra_end(); tit++)
        {
            get(*tit).label = label;
        }
        
        simplex_set ns_cl;
        mesh.closure(new_simplices, ns_cl);
        update(ns_cl);
    }
    
    /**
     * Improve tetrahedra by topological operations (re-connection): edge removal, multi-face removal,
     * multi-face retriangulation. It do so only for tetrahedra of quality lower than MIN_TET_QUALITY.
     */
    void topological_pass()
    {
        std::vector<tetrahedron_key_type> tets;
        for (auto tit = tetrahedra_begin(); tit != tetrahedra_end(); tit++)
        {
            if (quality(tit.key()) < MIN_TET_QUALITY)
            {
                tets.push_back(tit.key());
            }
        }
        
        // Attempt to remove each edge of each tetrahedron in tets. Accept if it increases the minimum quality locally.
        for (auto &t : tets)
        {
            if (mesh.exists(t) && quality(t) < MIN_TET_QUALITY)
            {
                simplex_set cl_t;
                closure(t, cl_t);
                
                for (auto eit = cl_t.edges_begin(); eit != cl_t.edges_end(); eit++)
                {
                    if (mesh.exists(*eit) && !get(*eit).is_interface() && !get(*eit).is_boundary())
                    {
                        edge_removal(*eit);
                    }
                }
            }
        }
        
        // Attempt to remove each face of each remaining tetrahedron in tets using multi-face retriangulation.
        // Accept if it increases the minimum quality locally.
        for (auto &t : tets)
        {
            if (mesh.exists(t) && quality(t) < MIN_TET_QUALITY)
            {
                simplex_set cl_t;
                closure(t, cl_t);
                
                for (auto fit = cl_t.faces_begin(); fit != cl_t.faces_end(); fit++)
                {
                    if (mesh.exists(*fit) && !get(*fit).is_interface() && !get(*fit).is_boundary())
                    {
                        multi_face_retriangulation(*fit);
                    }
                }
            }
        }
        
        // Attempt to remove each face of each remaining tetrahedron in tets using multi-face removal.
        // Accept if it increases the minimum quality locally.
        for (auto &t : tets)
        {
            if (mesh.exists(t) && quality(t) < MIN_TET_QUALITY)
            {
                simplex_set cl_t;
                closure(t, cl_t);
                
                for (auto fit = cl_t.faces_begin(); fit != cl_t.faces_end(); fit++)
                {
                    if (mesh.exists(*fit) && !get(*fit).is_interface() && !get(*fit).is_boundary())
                    {
                        multi_face_removal(*fit);
                    }
                }
            }
        }
        
        mesh.garbage_collect();
    }

    ////////////////////
    // EDGE FLIP PASS //
    ////////////////////
    
    void interface_edge_flip_pass()
    {
        std::vector<edge_key_type> flippable_edges;
        
        for (auto eit = edges_begin(); eit != edges_end(); eit++)
        {
            if (eit->is_interface() || eit->is_boundary())
            {
                flippable_edges.push_back(eit.key());
            }
        }
        
        for (auto e : flippable_edges)
        {
            if (is_flippable(e))
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
        std::vector<tetrahedron_key_type> tetrahedra;
        for (auto tit = mesh.tetrahedra_begin(); tit != mesh.tetrahedra_end(); tit++)
        {
            if (volume(tit.key()) > MAX_TET_VOLUME)
            {
                tetrahedra.push_back(tit.key());
            }
        }
        
        for(auto &tit : tetrahedra)
        {
            if (mesh.exists(tit))
            {
                split_tetrahedron(tit);
            }
        }
    }
    
    ///////////////////
    // THINNING PASS //
    ///////////////////
    /**
     * Collapses all edges which is shorter than MIN_EDGE_LENGTH. However, the collapse is only performed if it does not alter the interface and it improves the minimum quality of the tetrahedra in the star of the edge.
     */
    void thinning_pass()
    {
        std::vector<edge_key_type> edges;
        for (auto eit = edges_begin(); eit != edges_end(); eit++)
        {
            if (!eit->is_interface() && !eit->is_boundary())
            {
                edges.push_back(eit.key());
            }
        }
        
        std::vector<node_key_type> nodes;
        for (auto &e : edges)
        {
            if (mesh.exists(e) && length(e) < MIN_EDGE_LENGTH)
            {
                node_key_type n = safe_collapse(e);
                if (n != NULL_NODE)
                {
                    freitag_smoothing(n);
                }
            }
        }
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
        for (auto tit = SimplicialComplex<MT>::tetrahedra_begin(); tit != SimplicialComplex<MT>::tetrahedra_end(); tit++)
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
        for (auto nit = nodes_begin(); nit != nodes_end(); nit++)
        {
            if (!get(nit.key()).is_boundary() && !get(nit.key()).is_interface())
            {
                if (!smart_laplacian(nit.key()))
                {
                    if (min_quality(nit.key()) < MIN_TET_QUALITY)
                    {
                        freitag_smoothing(nit.key());
                    }
                }
            }
        }
    }
    
    //////////////////////////////
    // REMOVE DEGENERACIES PASS //
    //////////////////////////////
    
    /**
     * Attempt to remove edges shorter than threshold, provided that the volume enclosed by the interface changes within 0.0025 * threshold^3 tolerance
     */
    void remove_degenerate_edges()
    {
        std::list<edge_key_type> edge_list;
        T bnd_threshold = 0.25;
        T volume_int_threshold = 0.0025 * MIN_EDGE_LENGTH * MIN_EDGE_LENGTH * MIN_EDGE_LENGTH;
        T volume_bnd_threshold = 0.0025 * bnd_threshold * bnd_threshold * bnd_threshold;
        
        for (auto eit = mesh.edges_begin(); eit != mesh.edges_end(); eit++)
        {
            if (eit->is_interface() && length(eit.key()) < MIN_EDGE_LENGTH)
            {
                edge_list.push_back(eit.key());
            }
            else if (eit->is_boundary() && !eit->is_interface() && length(eit.key()) < bnd_threshold)
            {
                edge_list.push_back(eit.key());
            }
        }
        
        for (auto eit = edge_list.begin(); eit != edge_list.end(); eit++)
        {
            if (mesh.exists(*eit))
            {
                T l = length(*eit);
                if (mesh.find_edge(*eit).is_interface() && l < MIN_EDGE_LENGTH)
                {
                    simplex_set st_e;
                    mesh.star(*eit, st_e);
                    
                    
                    bool removable = false;
                    for (auto tit = st_e.tetrahedra_begin(); tit != st_e.tetrahedra_end(); tit++)
                    {
                        if (mesh.find_tetrahedron(*tit).label > 0)
                            removable = true;
                    }
                    
                    if (!removable) continue;
                    
                    node_key_type n0, n1;
                    typename Mesh::edge_type::boundary_list bnd_e = get(*eit).get_boundary();
                    auto beit = bnd_e->begin();
                    n0 = *beit; ++beit;
                    n1 = *beit;
                    
                    T v0 = volume_difference(*eit, n0, n1);
                    T v1 = volume_difference(*eit, n1, n0);
                    
                    if (v0 < v1)
                    {
                        if (v0 < volume_int_threshold)
                            unsafe_collapse(n0, n1);
                    }
                    else
                    {
                        if (v1 < volume_int_threshold)
                            unsafe_collapse(n1, n0);
                    }
                }
                else if (mesh.find_edge(*eit).is_boundary() && l < bnd_threshold)
                {
                    node_key_type n0, n1;
                    typename Mesh::edge_type::boundary_list bnd_e = mesh.find_edge(*eit).get_boundary();
                    auto beit = bnd_e->begin();
                    n0 = *beit; ++beit;
                    n1 = *beit;
                    
                    T v0 = volume_difference(*eit, n0, n1);
                    T v1 = volume_difference(*eit, n1, n0);
                    
                    if (v0 < v1)
                    {
                        if (v0 < volume_bnd_threshold)
                            unsafe_collapse(n0, n1);
                    }
                    else
                    {
                        if (v1 < volume_bnd_threshold)
                            unsafe_collapse(n1, n0);
                    }
                }
            }
        }
        
        mesh.garbage_collect();
    }
    
    /**
     * Attempt to remove the face f by first determining whether it's a cap and, if that's the case
     * splitting the longest edge at the projection of the opposite vertex and collapsing it with cap's apex.
     */
    void remove_degenerate_face(face_key_type const & f, std::vector<node_key_type> & nodes, std::vector<V> & verts, int i)
    {        
        T ratio = 1.03; // if the ratio between the lengths of two longest edges of f is less than this, it's treated as a needle.
        
        V d1 = verts[(i+1)%3] - verts[i];
        V d2 = verts[(i+2)%3] - verts[i];
        node_key_type n1 = nodes[(i+1)%3];
        node_key_type n2 = nodes[(i+2)%3];
        double l1 = d1.length(),
        l2 = d2.length();
        
        //    int cap_tip;
        V cap_projection;
        
        if (l1 > l2)
        {
            //        cap_tip = (i+2)%3;
            d1.normalize();
            cap_projection = verts[i] + d1 * MT::dot(d1, d2);
        }
        else
        {
            //        cap_tip = (i+1)%3;
            d2.normalize();
            cap_projection = verts[i] + d2 * MT::dot(d1, d2);
        }
        
        if ((l1 < ratio * l2) && (l2 < ratio * l1))
        {
            if (!(mesh.find_node(n2).is_interface()))
            {
                unsafe_collapse(n2,n1);
            }
            else
            {
                unsafe_collapse(n1,n2);
            }
        }
    }
    
    /**
     * Attempt to remove faces with at least one angle whose cosine is greater than MAX_COS_FACE_ANGLE.
     */
    void remove_degenerate_faces()
    {        
        std::list<face_key_type> face_list;
        
        for (auto fit = faces_begin(); fit != faces_end(); fit++)
            face_list.push_back(fit.key());
        
        for (auto it = face_list.begin(); it != face_list.end(); it++)
        {
            if (mesh.exists(*it))
            {
                std::vector<node_key_type> nodes(3);
                get_nodes(*it, nodes);
                std::vector<V> verts(3);
                get_pos(*it, verts);
                
                bool added = false;
                for (int i = 0; i < 3; ++i)
                {
                    V d1 = verts[(i+1)%3] - verts[i];
                    V d2 = verts[(i+2)%3] - verts[i];
                    d1.normalize();
                    d2.normalize();
                    if (MT::dot(d1, d2) > MAX_COS_FACE_ANGLE && (!added))
                        remove_degenerate_face(*it, nodes, verts, i);
                }
            }
        }
    }
    
    /**
     * Remove a degenerate tetrahedron of a type "sliver" by splitting the longest edges at their projected intersection point
     * and collapsing the newly created vertices together. Return true if successful.
     */
    bool remove_sliver(tetrahedron_key_type const & t, std::vector<node_key_type> & nodes, std::vector<V> & verts, int i, int j, int k, int l)
    {
        
        edge_key_type e1 = get_edge(nodes[i], nodes[j]);
        edge_key_type e2 = get_edge(nodes[k], nodes[l]);
        assert(e1 != NULL_EDGE);
        assert(e2 != NULL_EDGE);
        
        simplex_set st_e1, st_e2, interface_set;
        star(e1, st_e1);
        star(e2, st_e2);
        st_e2.add(st_e1);
        mesh.closure(st_e2, interface_set);
        
        mesh.set_undo_mark(interface_set);
        node_key_type n1 = mesh.split_edge(e1);
        node_key_type n2 = mesh.split_edge(e2);
        assert(n1 != NULL_NODE);
        assert(n2 != NULL_NODE);
        node_key_type result = unsafe_collapse(n1, n2);
        if (result == NULL_NODE)
        {
            mesh.undo();
            return false;
        }
        else
        {
            mesh.commit();
            cleanup_set(interface_set);
            update(interface_set);
            return true;
        }
    }
    
    /**
     * Remove a degenerate tetrahedron of a type "cap" by splitting the face opposite cap's apex and collapsing cap's apex with the newly created vertex.
     * Return true if successful.
     */
    bool remove_cap(const tetrahedron_key_type & t, std::vector<node_key_type> & nodes, std::vector<V> & verts, std::vector<T> & b_coords, V & projection, int j)
    {
        face_key_type f_split = -1;
        simplex_set bnd_t;
        mesh.boundary(t, bnd_t);
        for (auto fit = bnd_t.faces_begin(); fit != bnd_t.faces_end(); fit++)
        {
            simplex_set bnd_f;
            mesh.boundary(*fit, bnd_f);
            if (!(bnd_f.contains(nodes[j])))
                f_split = *fit;
        }
        
        //TODO: Possibly replace with exact collision point computation
        //    bool collision = false;
        //    if (mesh.find_face(f_split).is_interface() && mesh.find_node(nodes[j]).is_interface())
        //        collision = true;
        
        assert (mesh.exists(f_split));
        simplex_set st_f, interface_set;
        mesh.star(f_split, st_f);
        st_f.insert(f_split);
        mesh.closure(st_f, interface_set);
        mesh.set_undo_mark(interface_set);
        node_key_type nn = split_face(f_split);
        node_key_type result = unsafe_collapse(nodes[j], nn);
        if (result == NULL_NODE)
        {
            mesh.undo();
            return false;
        }
        else
        {
            mesh.commit();
            cleanup_set(interface_set);
            update(interface_set);
            return true;
        }
    }
    
    /**
     * Remove a degenerate tetrahedron of a type "wedge" by collapsing the shortest edge.
     * Return true if successful.
     */
    bool remove_wedge(const tetrahedron_key_type & t, std::vector<node_key_type> & nodes, std::vector<V> & verts, int i, int j, int k)
    {
        simplex_set si, sj;
        star(nodes[i], si);
        star(nodes[j], sj);
        si.intersection(sj);
        assert (si.size_edges() == 1);
        simplex_set st_e, interface_set;
        star(*(si.edges_begin()), st_e);
        st_e.insert(*(si.edges_begin()));
        mesh.closure(st_e, interface_set);
        mesh.set_undo_mark(interface_set);
        
        //    bool collision = false;
        //    if (mesh.find_edge(*(si.edges_begin())).is_interface() && mesh.find_node(nodes[k]).is_interface())
        //        collision = true;
        
        V vij = verts[j]-verts[i];
        double p1 = vij.length();
        vij.normalize();
        double p2 = dot(vij, verts[k]-verts[i]);
        node_key_type result;
        if (p2 < 0)
        {
            result = unsafe_collapse(nodes[k], nodes[i]);
        }
        else if (p2 > p1)
        {
            result = unsafe_collapse(nodes[k], nodes[j]);
        }
        else
        {
            node_key_type nn = split_edge(*(si.edges_begin()));
            result = unsafe_collapse(nodes[k], nn);
        }
        
        if (result == NULL_NODE)
        {
            mesh.undo();
            return false;
        }
        else
        {
            mesh.commit();
            cleanup_set(interface_set);
            update(interface_set);
            return true;
        }
    }
    
    /**
     * Destroy degenerate (nearly flat) tetrahedron t by splits and collapses.
     * This function detects what type of degeneracy tetrahedron t is (sliver, cap or wedge)
     * and selects appropriate degeneracy removal routine.
     */
    inline void remove_degenerate_tet(const tetrahedron_key_type & t)
    {
        // Find the largest face
        simplex_set cl_t;
        closure(t, cl_t);
        face_key_type f = largest_face(cl_t);
        std::vector<V> verts;
        get_pos(f, verts);
        
        // Find the apex
        simplex_set cl_f;
        closure(f, cl_f);
        cl_t.difference(cl_f);
        node_key_type apex = *cl_t.nodes_begin();
        
        // Project the apex
        V proj_apex = Util::project<MT>(get_pos(apex), verts);
        
        std::vector<T> barycentric_coords(3);
        Util::get_barycentric_coords<MT>(proj_apex, verts[0], verts[1], verts[2], barycentric_coords);

        std::vector<T> cosines(6);
        Util::get_cosines<MT>(proj_apex, verts, cosines);
        
        // Determine the position of the apex with regard to the edges of the biggest face.
        std::vector<int> inside(3);
        for (int i = 0; i < 3; i++)
        {
            if (barycentric_coords[i] > EPSILON)
            {
                if (cosines[2*i] <= MAX_COS && cosines[2*i+1] <= MAX_COS)
                {
                    inside[i] = 1;
                }
            }
            else if (barycentric_coords[i] < -EPSILON)
            {
                if (cosines[2*i] <= MAX_COS && cosines[2*i+1] <= MAX_COS)
                {
                    inside[i] = -1;
                }
            }
            else
                inside[i] = 0;
        }
        
        // Select appropriate degeneracy removal routine based on the location of the apex with regard to the largest face.
        
        std::vector<node_key_type> nodes;
        get_nodes(f, nodes);
        nodes.push_back(apex);
        verts.push_back(get_pos(apex));
        
        if (inside[0] == 1 && inside[1] == 1 && inside[2] == 1)
            remove_cap(t, nodes, verts, barycentric_coords, proj_apex, 3);
        else if (inside[0] == -1 && inside[1] ==  1 && inside[2] ==  1)
            remove_sliver(t, nodes, verts, 0, 1, 2, 3);
        else if (inside[0] ==  1 && inside[1] == -1 && inside[2] ==  1)
            remove_sliver(t, nodes, verts, 1, 2, 0, 3);
        else if (inside[0] ==  1 && inside[1] ==  1 && inside[2] == -1)
            remove_sliver(t, nodes, verts, 2, 0, 1, 3);
        else if (inside[0] ==  0 && inside[1] ==  1 && inside[2] ==  1)
            remove_wedge(t, nodes, verts, 0, 1, 3);
        else if (inside[0] ==  1 && inside[1] ==  0 && inside[2] ==  1)
            remove_wedge(t, nodes, verts, 1, 2, 3);
        else if (inside[0] ==  1 && inside[1] ==  1 && inside[2] ==  0)
            remove_wedge(t, nodes, verts, 2, 0, 3);
    }
    
    /**
     * Attempt to remove tetrahedra with quality lower than MIN_QUALITY.
     */
    void remove_degenerate_tets()
    {
        std::vector<tetrahedron_key_type> degenerated_tets;
        
        for (auto tit = tetrahedra_begin(); tit != tetrahedra_end(); tit++)
        {
            if (quality(tit.key()) <= DEG_TET_QUALITY)
            {
                degenerated_tets.push_back(tit.key());
            }
        }
        
        for (auto &tet : degenerated_tets)
        {
            if (mesh.exists(tet) && quality(tet) <= DEG_TET_QUALITY)
            {
                remove_degenerate_tet(tet);
            }
        }
        
        mesh.garbage_collect();
    }
    
    void fix_complex()
    {
        simplex_set region;
        
        print_out("Smooth.");
        smooth();
        
        print_out("Topological pass.");
        topological_pass();
        validity_check();
        
        print_out("Edge flip pass.");
        interface_edge_flip_pass();
        validity_check();
        
        print_out("Remove degeneracies.");
        remove_degenerate_tets();
        remove_degenerate_faces();
        remove_degenerate_edges();
        validity_check();
        
        print_out("Smooth.");
        smooth();
        validity_check();
        
        print_out("Relabel tets.");
        bool relabeled = relabel_tets();
        validity_check();
        
        if (relabeled) {
            smooth();
        }
    }
    
    void resize_complex()
    {
        print_out("Thinning pass.");
        thinning_pass();
        validity_check();
        
        print_out("Thickening pass.");
        thickening_pass();
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
        bool finished = false;
        int step = 0;
        while (!finished && step < num_steps)
        {
            std::cout << "Move vertices step " << step << std::endl;
            finished = true;
            for (auto nit = nodes_begin(); nit != nodes_end(); nit++)
            {
                if (nit->is_interface() && mesh.exists(nit.key()))
                {
                    finished = finished & move_vertex(nit.key(), nit->get_destination());
                }
            }
            
            fix_complex();
            
            ++step;
        }
        
        if(step_no%5 == 0)
        {
            resize_complex();
        }
        
        mesh.garbage_collect();
        ++step_no;
    }
    
private:
    
    /**
     * Tries moving the node n to the new position new_pos. Returns true if it succeeds.
     */
    bool move_vertex(const node_key_type & n, const typename MT::vector3_type & new_pos)
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
    T intersection_with_link(const node_key_type & n, const V& new_pos)
    {
        V pos = get_pos(n);
        simplex_set ln;
        link(n,ln);
        T min_t = INFINITY;
        std::vector<V> face_pos;
        for(auto fit = ln.faces_begin(); fit != ln.faces_end(); fit++)
        {
            get_pos(*fit, face_pos);
            T t = Util::intersection<MT>(pos, new_pos, face_pos);
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
    bool smart_laplacian(const node_key_type& n, T alpha = 1.)
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
    T get_quality_info(const node_key_type & n, std::vector<quality_struct<MT>> & full_quality_info, std::list<quality_struct<MT>> & active)
    {
        simplex_set st_n;
        star(n, st_n);
        
        T min_quality = INFINITY;
        
        for (auto tit = st_n.tetrahedra_begin(); tit != st_n.tetrahedra_end(); tit++)
        {
            quality_struct<MT> quality_info;
            
            std::vector<node_key_type> nodes;
            get_nodes(*tit, nodes);
            
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
    
    T calc_step_length(const node_key_type &n, T init_alpha, V s, T min_q_old, T min_q, T min_q_pre, T min_step)
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
    void freitag_smoothing(node_key_type const & n, T min_step = 1e-4, T min_improvement = 1e-5, int max_iter = 3)
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
    bool precond_relabel(const tetrahedron_key_type& t, face_key_type & f)
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
        closure(t, cl_t);
        for (auto fit = cl_t.faces_begin(); fit != cl_t.faces_end(); fit++)
        {
            T a = area(*fit);
            total_area += a;
            
            if (get(*fit).is_interface() && !get(*fit).is_boundary() && a > max_area)
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
            link(f, lk_f);
            for(auto nit = lk_f.nodes_begin(); nit != lk_f.nodes_end(); nit++)
            {
                if (cl_t.contains(*nit) && get(*nit).is_interface() &&
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
    bool relabel_tet(const tetrahedron_key_type& t)
    {
        face_key_type f;
        
        if (precond_relabel(t, f))
        {
            int label = -1;
            simplex_set st_f;
            star(f, st_f);
            for (auto tit = st_f.tetrahedra_begin(); tit != st_f.tetrahedra_end(); tit++)
            {
                if(*tit != t)
                {
                    label = get(*tit).label;
                }
            }
            
            if (label != -1 && label != get(t).label)
            {
                get(t).label = label;
                simplex_set cl_t;
                closure(t, cl_t);
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
    bool is_flippable(const edge_key_type & e)
    {
        std::vector<V> verts;
        get_pos(e, verts);
        
        simplex_set ln_e;
        link(e, ln_e);
        for (auto nit = ln_e.nodes_begin(); nit != ln_e.nodes_end(); nit++)
        {
            if (get(*nit).is_interface())
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
    void flip(const edge_key_type & e)
    {
        // Find the pair of nodes to be connected by a new edge as a result of edge flip.
        simplex_set lk_e;
        link(e, lk_e);
        std::vector<node_key_type> nodes;
        for (auto nit = lk_e.nodes_begin(); nit != lk_e.nodes_end(); nit++)
        {
            if (get(*nit).is_interface() || get(*nit).is_boundary()) {
                nodes.push_back(*nit);
            }
        }
        
        assert(nodes.size() == 2);
        if(get_edge(nodes[0], nodes[1]) == NULL_EDGE) // Check that there does not already exist an edge between the pair of nodes.
        {
            return;
        }
        
        node_key_type n_new = split_edge(e);
        if (n_new == NULL_NODE) {
            return;
        }
        
        edge_key_type e_c = get_edge(nodes[0], n_new);
        assert(e_c != NULL_EDGE);
        V p = get_pos(nodes[0]);
        V p_new = get(nodes[0]).get_destination();
        
        node_key_type n = collapse_edge(e_c , nodes[0], n_new, p);
        if(n != NULL_NODE)
        {
            set_destination(n, p_new);
        }
    }
    
    ////////////
    // SPLITS //
    ////////////
public:
    /**
     * Split a tetrahedron t and returns the new node which is positioned at the barycenter of the vertices of t.
     */
    node_key_type split_tetrahedron(tetrahedron_key_type& t)
    {
        std::vector<V> verts;
        get_pos(t, verts);
        int label = get(t).label;
        V p = Util::barycenter<MT>(verts[0], verts[1], verts[2], verts[3]);
        
        node_key_type n = mesh.split_tetrahedron(t);
        set_pos(n, p);
        set_destination(n, p);
        
        simplex_set st_n;
        star(n, st_n);
        st_n.insert(n);
        
        for (auto tit = st_n.tetrahedra_begin(); tit != st_n.tetrahedra_end(); tit++)
        {
            get(*tit).label = label;
        }
        update(st_n);
        return n;
    }
    
    /**
     * Split a face f and returns the new node which is positioned at the barycenter of the vertices of f.
     */
    node_key_type split_face(const face_key_type & f)
    {
        std::map<tetrahedron_key_type, int> tt;
        simplex_set st_f;
        star(f, st_f);
        
        for (auto tit = st_f.tetrahedra_begin(); tit != st_f.tetrahedra_end(); tit++)
        {
            tt[*tit] = get(*tit).label;
        }
        
        std::vector<V> verts;
        get_pos(f, verts);
        V p = Util::barycenter<MT>(verts[0], verts[1], verts[2]);
        
        std::map<tetrahedron_key_type, tetrahedron_key_type> new_tets;
        node_key_type n = mesh.split_face_helper(f, new_tets);
        set_pos(n, p);
        set_destination(n, p);
        
        for(auto it = new_tets.begin(); it != new_tets.end(); it++)
        {
            get(it->first).label = tt[it->second];
        }
        
        simplex_set st_n;
        star(n, st_n);
        st_n.insert(n);
        
        update(st_n);
        return n;
    }
    
    /**
     * Split an edge e and returns the new node which is placed at the middle of e.
     */
    node_key_type split_edge(const edge_key_type & e)
    {
        std::vector<V> verts;
        get_pos(e, verts);
        V p = Util::barycenter<MT>(verts[0], verts[1]);
        
        std::map<tetrahedron_key_type, int> tt;
        simplex_set st_e;
        star(e, st_e);
        
        for (auto tit = st_e.tetrahedra_begin(); tit != st_e.tetrahedra_end(); tit++)
        {
            tt[*tit] = get(*tit).label;
        }
        
        std::map<tetrahedron_key_type, tetrahedron_key_type> new_tets;
        node_key_type n = mesh.split_edge_helper(e, new_tets);
        set_pos(n, p);
        set_destination(n, p);
        
        for (auto it = new_tets.begin(); it != new_tets.end(); it++)
        {
            get(it->first).label = tt[it->second];
        }
        
        simplex_set st_n;
        star(n, st_n);
        st_n.insert(n);
        
        update(st_n);
        return n;
    }
    
    ///////////////
    // COLLAPSES //
    ///////////////
private:
    /**
     * Collapses the edge e by moving the node n2 to n1. Returns n1 if successful, otherwise NULL_NODE.
     */
    node_key_type collapse_edge(edge_key_type& e, node_key_type const & n1, node_key_type const & n2, const V& p)
    {
        node_key_type n_new = mesh.edge_collapse_helper(e, n1, n2);
        
        if (n_new != NULL_NODE)
        {
            set_pos(n_new, p);
            set_destination(n_new, p);
            
            simplex_set st, st_cl;
            star(n_new, st);
            mesh.closure(st, st_cl);
            update(st_cl);
        }
        return n_new;
    }
    
    /**
     * Returns true if the collapse of the edge e results in a minimum tetrahedron quality above MIN_TET_QUALITY or in an improvement of the minimum tetrahedron quality.
     * The merged nodes are assumed moved to v_new after the collapse.
     */
    bool precond_collapse(const edge_key_type& e, const V& v_new)
    {
        if (e == NULL_EDGE)
        {
            return false;
        }
        std::vector<node_key_type> nodes;
        get_nodes(e, nodes);
        
        if(inverted(nodes[0]) || inverted(nodes[1])) // Should not be necessary!
        {
            return false;
        }
            
        simplex_set lk_n0, lk_n1, lk_e;
        link(nodes[0], lk_n0);
        link(nodes[1], lk_n1);
        
        simplex_set st_e, cl_st_e;
        star(e, st_e);
        mesh.closure(st_e, cl_st_e);
        lk_n0.add(lk_n1);
        lk_n0.difference(cl_st_e);
        
        if(will_invert(nodes[0], v_new, lk_n0) || will_invert(nodes[1], v_new, lk_n0))
        {
            return false;
        }
        return true;
    }
    
    /**
     * Returns the increase (positive) or decrease in quality if the edge e is collapsed and the resulting node is moved to v_new.
     */
    T quality_improvement(const edge_key_type& e, const V& v_new)
    {
        std::vector<node_key_type> nodes;
        get_nodes(e, nodes);
        
        simplex_set lk_n0, lk_n1, lk_e;
        link(nodes[0], lk_n0);
        link(nodes[1], lk_n1);
        
        simplex_set st_e, cl_st_e;
        star(e, st_e);
        mesh.closure(st_e, cl_st_e);
        lk_n0.add(lk_n1);
        lk_n0.difference(cl_st_e);
        T q_new = min_quality(lk_n0, v_new);
        
        simplex_set st_n0, st_n1;
        star(nodes[0], st_n0);
        star(nodes[1], st_n1);
        st_n0.add(st_n1);
        
        T q = min_quality(st_n0);
        return q_new - q;
    }
    
    /**
     * Checks if the nodes of edge e are editable, i.e. neither a part of the interface nor the boundary. If both are editable, unsafe_collapse is called.
     * If one is editable then this node is moved to the other if precond_collapse returns true.
     * If non of the nodes are editable or precond_collapse returns false, the method returns NULL_NODE.
     */
    node_key_type safe_collapse(edge_key_type& e)
    {
        std::vector<node_key_type> nodes;
        get_nodes(e, nodes);
        node_key_type n1 = nodes[0];
        node_key_type n2 = nodes[1];
        
        bool editable1 = !get(n1).is_interface() && !get(n1).is_boundary();
        bool editable2 = !get(n2).is_interface() && !get(n2).is_boundary();
        
        if(editable1 || editable2)
        {
            V p;
            if (editable1 && editable2)
            {
                p = Util::barycenter<MT>(get_pos(n1), get_pos(n2));
                
            }
            else if (editable1)
            {
                p = get_pos(n2);
            }
            else if (editable2)
            {
                p = get_pos(n1);
            }
            
            if (precond_collapse(e, p) && quality_improvement(e, p) > 0.)
            {
                return collapse_edge(e, n1, n2, p);
            }
        }
        return NULL_NODE;
    }
    
    /**
     * n1 and n2 are moved to the barycenter of their positions and the edge between them are collapsed if precond_collapse returns true.
     * NOTE: This method allows to move interface vertices. Use safe_collapse to make sure the interface is not changed.
     */
    node_key_type unsafe_collapse(const node_key_type & n1, const node_key_type & n2)
    {
        edge_key_type e = get_edge(n1, n2);
        
        V p = Util::barycenter<MT>(get_pos(n1), get_pos(n2));
        V p_new = Util::barycenter<MT>(get(n1).get_destination(), get(n2).get_destination());
        
        if (precond_collapse(e, p))
        {
            node_key_type n_new = collapse_edge(e, n1, n2, p);
            if (n_new != NULL_NODE)
            {
                set_destination(n_new, p_new);
            }
            return n_new;
        }
        return NULL_NODE;
    }
    
    
    //////////////////////
    // GETTER FUNCTIONS //
    //////////////////////
public:
    
    /**
     Returns the normal to interface face f.
     */
    V get_normal(const face_key_type & f)
    {
        std::vector<V> verts;
        get_pos(f, verts);
        return Util::normal_direction<MT>(verts[0], verts[1], verts[2]);
    }
    
    /**
     Returns the normal to interface node n.
     */
    V get_normal(const node_key_type & n)
    {
        simplex_set st;
        star(n, st);
        
        V result(0.);
        for (auto fit = st.faces_begin(); fit != st.faces_end(); fit++)
        {
            if (get(*fit).is_interface())
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
    V get_barycenter(const node_key_type& n, bool interface = false)
    {        
        simplex_set lk_n;
        link(n, lk_n);
        
        V avg_pos(0.);
        int i = 0;
        for (auto nit = lk_n.nodes_begin(); nit != lk_n.nodes_end(); nit++)
        {
            if (!interface || get(*nit).is_interface())
            {
                avg_pos += get_pos(*nit);
                i++;
            }
        }
        assert(i != 0);
        return avg_pos / static_cast<T>(i);
    }
    
    /**
     * Returns the position of node n.
     */
    V get_pos(const node_key_type& n)
    {
        V p = get(n).get_pos();
        assert(!MT::is_nan(p[0]) && !MT::is_nan(p[1]) && !MT::is_nan(p[2]));
        return p;
    }
    
    /// Returns the positions of the nodes of edge e in verts.
    void get_pos(const edge_key_type & e, std::vector<V>& verts)
    {
        verts = std::vector<V>(2);
        std::vector<node_key_type> nodes(2);
        get_nodes(e, nodes);
        for (int k = 0; k < 2; ++k)
        {
            verts[k] = get_pos(nodes[k]);
        }
    }
    
    /// Returns the positions of the nodes of face f in verts.
    void get_pos(const face_key_type & f, std::vector<V>& verts)
    {
        verts = std::vector<V>(3);
        orient_face(f);
        std::vector<node_key_type> nodes(3);
        get_nodes(f, nodes);
        for (int k = 0; k < 3; ++k)
        {
            verts[k] = get_pos(nodes[k]);
        }
    }
    
    /// Returns the positions of the nodes of tetrahedron t in verts.
    void get_pos(const tetrahedron_key_type& t, std::vector<V>& verts)
    {
        verts = std::vector<V>(4);
        std::vector<node_key_type> nodes(4);
        get_nodes(t, nodes);
        for (int k = 0; k < 4; ++k)
        {
            verts[k] = get_pos(nodes[k]);
        }
    }
    
    //////////////////////
    // SETTER FUNCTIONS //
    //////////////////////
private:
    /**
     * Sets the position of node n.
     */
    void set_pos(const node_key_type& n, V p)
    {
        get(n).set_pos(p);
    }
    
public:
    /**
     * Sets the destination where the node n is moved to when deform() is called.
     */
    void set_destination(const node_key_type& n, V p)
    {
        get(n).set_destination(p);
    }
    
    ///////////////////////
    // UTILITY FUNCTIONS //
    ///////////////////////
    
public:
    
    T length(const edge_key_type& e)
    {
        std::vector<V> verts;
        get_pos(e,verts);
        return MT::length(verts[0] - verts[1]);
    }
    
    T area(const face_key_type& f)
    {
        std::vector<V> verts;
        get_pos(f,verts);
        return Util::area<MT>(verts);
    }
    
    T volume(const tetrahedron_key_type& t)
    {
        std::vector<V> verts;
        get_pos(t,verts);
        return Util::volume<MT>(verts[0], verts[1], verts[2], verts[3]);
    }
    
    T signed_volume(const tetrahedron_key_type& t)
    {
        std::vector<V> verts;
        get_pos(t,verts);
        return Util::signed_volume<MT>(verts[0], verts[1], verts[2], verts[3]);
    }
    
    T quality(const tetrahedron_key_type& t)
    {
        std::vector<V> verts;
        get_pos(t, verts);
        return Util::quality<MT>(verts[0], verts[1], verts[2], verts[3]);
    }
    
    /**
     * Returns the largest face in the simplex set.
     */
    face_key_type largest_face(simplex_set& set)
    {
        T max_a = -INFINITY;
        face_key_type max_f;
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
    
    /// Checks whether the interface patch around n is 2-manifold
    //    bool is_2manifold(const node_key_type & n)
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
    //        typename SimplicialComplex<MT>::face_key_type f = *(interface_simplices.faces_begin());
    //        do
    //        {
    //            interface_simplices.erase(f);
    //            typename SimplicialComplex<MT>::edge_key_type e;
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
    T min_quality(const node_key_type& n)
    {
        simplex_set st_n;
        star(n, st_n);
        return min_quality(st_n);
    }
    
    /**
     * Returns minimum quality among the tetrahedra adjacent to face f.
     */
    T min_quality(const face_key_type& f)
    {
        simplex_set st_f;
        star(f, st_f);
        return min_quality(st_f);
    }
    
    // The minimum quality among tetrahedra that would be produced by connecting apices with the edges on mf_boundary.
    T edge_remove_quality(simplex_set& mf_boundary, std::vector<V> & apices)
    {
        T q_min = 1e99;
        unsigned int n = static_cast<unsigned int>(mf_boundary.size_nodes());
        
        assert ( n > 2 || !"Too few vertices in the multi-face boundary!" );
        
        std::vector<node_key_type> polygon(n);
        
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
    T volume_difference(edge_key_type & e, node_key_type & n0, node_key_type & n1)
    {        
        T result = 0.0;
        simplex_set st0, st1, lk0, st_u, patch, cl_patch;
        
        star(n0, st0);
        link(n0, lk0);
        star(n1, st1);
        st0.insert(n0);
        st1.insert(n1);
        st_u.add(st0);
        st_u.add(st1);
        
        
        for (auto fit = st_u.faces_begin(); fit != st_u.faces_end(); fit++)
        {
            if (get(e).is_interface())
            {
                if (get(*fit).is_interface())
                    patch.insert(*fit);
            }
            else if (get(e).is_boundary())
            {
                if (get(*fit).is_boundary())
                    patch.insert(*fit);
            }
        }
        
        mesh.closure(patch, cl_patch);
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
     * Ensures consistent orientation of all faces to the two tetrahedra which are in the star of f.
     */
    void orient_face(const face_key_type& f)
    {
        if (get(f).is_interface())
        {
            simplex_set st_f;
            star(f, st_f);
            int label = -100;
            for (auto tit = st_f.tetrahedra_begin(); tit != st_f.tetrahedra_end(); tit++)
            {
                int tl = get(*tit).label;
                
                if (tl > label)
                {
                    mesh.orient_faces_consistently(*tit);
                }
                label = tl;
            }
        }
    }
    
    /**
     * Check if the sequence of vertices in polygon is consistent with positive orientation of tetrahedra in the mesh
     * with respect to the ordered pair of vertices in vv. If not, reverse the order of vertices in polygon.
     */
    void check_consistency(std::vector<V> & vv, std::vector<node_key_type> & polygon)
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
                node_key_type temp = polygon[i];
                polygon[i] = polygon[n-1-i];
                polygon[n-1-i] = temp;
            }
        }
    }
    
    /**
     * Check if the sequence of vertices in polygon is consistent with positive orientation of tetrahedra in the mesh
     * with respect to edge. If not, reverse the order of vertices in polygon.
     */
    void check_consistency(edge_key_type const & e, std::vector<node_key_type> & polygon)
    {
        std::vector<V> verts;
        get_pos(e, verts);
        std::vector<V> vv(2);
        vv[1] = verts[0]; // VERTEX FLIP
        vv[0] = verts[1];
        
        check_consistency(vv, polygon);
    }
    
    
    /// Check whether they are no inverted tetrahedra in the mesh.
    bool simplicial_complex_criterion_check()
    {
        for (auto tit = tetrahedra_begin(); tit != tetrahedra_end(); tit++)
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
    bool inverted(const node_key_type& n)
    {
        simplex_set st_n;
        star(n, st_n);
        return inverted(st_n);
    }
    
    /**
     * Returns whether any of the tetrahedra in the simplex set will invert if node n is moved to p_new.
     */
    bool will_invert(const node_key_type& n, const V p_new, simplex_set& set)
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
    bool will_invert(const node_key_type& n, const V p_new)
    {
        simplex_set st_n;
        star(n);
        return will_invert(n, p_new, st_n);
    }
    
    void validity_check()
    {
        bool valid = simplicial_complex_criterion_check();
        assert(valid);
        for (auto tit = tetrahedra_begin(); tit != tetrahedra_end(); tit++)
        {
            valid = valid & mesh.exists(tit.key());
        }
        assert(valid);
        
        for (auto fit = faces_begin(); fit != faces_end(); fit++)
        {
            valid = valid & mesh.exists(fit.key());
        }
        assert(valid);
        
        for (auto eit = edges_begin(); eit != edges_end(); eit++)
        {
            valid = valid & mesh.exists(eit.key());
        }
        assert(valid);
        
        for (auto nit = nodes_begin(); nit != nodes_end(); nit++)
        {
            valid = valid & mesh.exists(nit.key());
        }
        assert(valid);
    }
    
    void check_validity()
    {
        std::cout << "Validity check" << std::endl;
        
        bool valid = true;
        
        for (auto fit = faces_begin(); fit != faces_end(); fit++)
        {
            simplex_set st_f;
            mesh.star(fit.key(), st_f);
            if (st_f.size_tetrahedra() > 2)
            {
                std::map<tetrahedron_key_type, simplex_set> t_bnds;
                
                std::cout << fit.key() << ": " << std::endl;
                auto tit = st_f.tetrahedra_begin();
                while (tit != st_f.tetrahedra_end())
                {
                    simplex_set bnd_t;
                    mesh.boundary(*tit, bnd_t);
                    t_bnds[*tit] = bnd_t;
                    
                    auto ftit = bnd_t.nodes_begin();
                    while (ftit != bnd_t.nodes_end())
                    {
                        std::cout << *ftit << " ";
                        ++ftit;
                    }
                    std::cout << std::endl;
                    ++tit;
                }
                
                while (!t_bnds.empty())
                {
                    typename std::map<tetrahedron_key_type, simplex_set>::iterator it;
                    it = t_bnds.begin();
                    tetrahedron_key_type t = it->first;
                    simplex_set bnd0;
                    bnd0.add(it->second);
                    ++it;
                    while (it != t_bnds.end())
                    {
                        simplex_set bnd1;
                        bnd1.add(it->second);
                        bnd1.difference(bnd0);
                        if (bnd1.size_faces() == 0)
                        {
                            mesh.unsafe_remove(t); // FUCKED
                            break;
                        }
                        ++it;
                    }
                    t_bnds.erase(t_bnds.begin());
                }
                
                valid = false;
            }
        }
        
        if (!valid)
            std::cout << "Input mesh invalid" << std::endl;
        else
            std::cout << "Input mesh valid" << std::endl;
    } // check_validity
    
    
    /// Remove the simplices that are no longer represented in the mesh from set.
    void cleanup_set(simplex_set& set)
    {
        simplex_set unused_simplices;
        
        for (auto nit = set.nodes_begin(); nit != set.nodes_end(); nit++)
        {
            if (!mesh.exists(*nit))
                unused_simplices.insert(*nit);
        }
        
        for (auto eit = set.edges_begin(); eit != set.edges_end(); eit++)
        {
            if (!mesh.exists(*eit))
                unused_simplices.insert(*eit);
        }
        
        for (auto fit = set.faces_begin(); fit != set.faces_end(); fit++)
        {
            if (!mesh.exists(*fit))
                unused_simplices.insert(*fit);
        }
        
        for (auto tit = set.tetrahedra_begin(); tit != set.tetrahedra_end(); tit++)
        {
            if (!mesh.exists(*tit))
            {
                unused_simplices.insert(*tit);
            }
        }
        
        set.difference(unused_simplices);
    }
    
    /**
     * Sort the vertices from set according to their connectivity and returned the sorted vertices in sorted_vertices.
     */
    void sort_vertices(simplex_set& set, std::vector<node_key_type>& sorted_vertices)
    {
        sorted_vertices = std::vector<node_key_type>();
        sorted_vertices.push_back(*(set.nodes_begin()));
        
        std::map<edge_key_type, bool> edge_used;
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
                    std::vector<node_key_type> nodes;
                    get_nodes(*eit, nodes);
                    
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
        std::map<node_key_type, int> vert_index;
        
        // Extract vertices
        for (auto nit = nodes_begin(); nit != nodes_end(); nit++)
        {
            if (nit->is_interface())
            {
                verts.push_back(nit->v);
                vert_index[nit.key()] = verts.size();
            }
        }
        
        // Extract faces
        for (auto fit = SimplicialComplex<MT>::faces_begin(); fit != SimplicialComplex<MT>::faces_end(); fit++)
        {
            if (fit->is_interface())
            {
                std::vector<node_key_type> nodes(3);
                get(fit.key(), nodes);
                
                indices.push_back(vert_index[nodes[0]]);
                indices.push_back(vert_index[nodes[1]]);
                indices.push_back(vert_index[nodes[2]]);
            }
        }
    } // extract_interface
    
    void extract_tet_mesh(std::vector<typename MT::vector3_type> & points, std::vector< std::vector<int> > & tets)
    {
        std::map<node_key_type, int> indices;
        mesh.garbage_collect();
        
        int counter = 0;
        for (auto nit = nodes_begin(); nit != nodes_end(); nit++)
        {
            indices[nit.key()] = counter;
            points.push_back(nit->v);
            ++counter;
        }
        
        for (auto tit = tetrahedra_begin(); tit != tetrahedra_end(); tit++)
        {
            std::vector<node_key_type> nodes(4);
            std::vector<int> tet;
            get(tit.key(), nodes);
            
            for (int i = 0; i < nodes.size(); i++)
                tet.push_back(indices[nodes[i]]);
            
            tet.push_back(tit->label);
            tets.push_back(tet);
        }
    } // extract_tet_mesh
    
    /**
     * Returns the cosine to the dihedral angle between face f1 and face f2.
     */
    T cos_dihedral_angle(const face_key_type& f1, const face_key_type& f2)
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
    T dihedral_angle(const face_key_type& f1, const face_key_type& f2)
    {
        return acos(cos_dihedral_angle(f1, f2));
    }
    
    /**
     * Returns the minimum dihedral angle between the faces of tetrahedron t.
     */
    T min_cos_dihedral_angle(const tetrahedron_key_type& t)
    {
        T min_angle = -1.;
        simplex_set cl_t;
        closure(t, cl_t);
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
    
    T min_dihedral_angle(const tetrahedron_key_type& t)
    {
        return acos(min_cos_dihedral_angle(t));
    }
    
    /**
     * Calculates the dihedral angles in the SimplicialComplex and returns these in a histogram,
     * along with the minimum and maximum dihedral angles.
     */
    void calc_dihedral_angles(std::vector<int> & histogram, T & min_angle, T & max_angle)
    {
        max_angle = -INFINITY, min_angle = INFINITY;
        
        histogram = std::vector<int>(180);
        for (int i = 0; i < 180; ++i)
        {
            histogram[i] = 0;
        }
        
        for (auto tit = tetrahedra_begin(); tit != tetrahedra_end(); tit++)
        {
            simplex_set cl_t;
            closure(tit.key(), cl_t);
            
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
        for (auto tit = tetrahedra_begin(); tit != tetrahedra_end(); tit++)
        {
            min_q = std::min(min_q, quality(tit.key()));
        }
        return min_q;
    }
    
    /// Counts the total number of nodes and the number of nodes on the interface(s).
    void count_nodes(int & total, int & object)
    {
        total = 0, object = 0;
        for (auto nit = nodes_begin(); nit != nodes_end(); nit++)
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
        for (auto eit = edges_begin(); eit != edges_end(); eit++)
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
        for (auto fit = faces_begin(); fit != faces_end(); fit++)
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
        for (auto tit = tetrahedra_begin(); tit != tetrahedra_end(); tit++)
        {
            total++;
            if (tit->label != 0)
                object++;
        }
    }
    
};

// SIMPLICIAL_COMPLEX_H
#endif
