#include "DSC.h"

using namespace is_mesh;


namespace DSC {

    DeformableSimplicialComplex::DeformableSimplicialComplex(std::vector<vec3> & points, std::vector<int> & tets, const std::vector<int>& tet_labels)
            : is_mesh::ISMesh(points, tets, tet_labels) {
                pars = {0.1, 0.5, 0.0005, 0.015, 0.02, 0.3, 0., 2., 0.2, 5., 0.2, INFINITY};
                set_avg_edge_length();
            }


    DeformableSimplicialComplex::~DeformableSimplicialComplex() {

    }

    void DeformableSimplicialComplex::set_avg_edge_length(double avg_edge_length) {
        if(avg_edge_length <= 0.)
        {
            AVG_LENGTH = compute_avg_edge_length();
        }
        else {
            AVG_LENGTH = avg_edge_length;
        }

        AVG_AREA = 0.5*std::sqrt(3./4.)*AVG_LENGTH*AVG_LENGTH;
        AVG_VOLUME = (sqrt(2.)/12.)*AVG_LENGTH*AVG_LENGTH*AVG_LENGTH;
    }

    void DeformableSimplicialComplex::set_parameters(parameters pars_) {
        pars = pars_;
    }

    void DeformableSimplicialComplex::set_design_domain(is_mesh::Geometry *geometry) {
        design_domain.add_geometry(geometry);
    }

    void DeformableSimplicialComplex::set_labels(const is_mesh::Geometry& geometry, int label) {
        for (auto tit = tetrahedra_begin(); tit != tetrahedra_end(); tit++) {
            is_mesh::SimplexSet<is_mesh::NodeKey> nids = get_nodes(tit.key());
            if(geometry.is_inside(Util::barycenter(get(nids[0]).get_pos(), get(nids[1]).get_pos(), get(nids[2]).get_pos(), get(nids[3]).get_pos())))
            {
                set_label(tit.key(), label);
            }
        }
    }

    void DeformableSimplicialComplex::print(const node_key& n) {
        std::cout << "Node: " << n << std::endl;
        vec3 p = get_pos(n);
        vec3 d = get(n).get_destination();
        std::cout << "P = " << p[0] << ", " << p[1] << ", " << p[2] << std::endl;
        std::cout << "D = " << d[0] << ", " << d[1] << ", " << d[2]  << std::endl;

        std::cout << "\nStar_edges = [";
        for(auto e : get_edges(n))
        {
            auto verts = get_pos(get_nodes(e));
            vec3 p1 = verts[0];
            vec3 p2 = verts[1];
            std::cout << p1[0] << ", " << p1[1] << ", " << p1[2] << "; " << std::endl;
            std::cout << p2[0] << ", " << p2[1] << ", " << p2[2] << "; " << std::endl;
        }
        std::cout << "];" << std::endl;

        std::cout << "\nStar_Iedges = [";
        for(auto e : get_edges(n))
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
                auto verts = get_pos(get_nodes(e));
                vec3 p1 = verts[0];
                vec3 p2 = verts[1];
                std::cout << p1[0] << ", " << p1[1] << ", " << p1[2] << "; " << std::endl;
                std::cout << p2[0] << ", " << p2[1] << ", " << p2[2] << "; " << std::endl;
            }
        }
        std::cout << "];" << std::endl;
    }

    bool DeformableSimplicialComplex::is_unsafe_editable(const node_key& nid) {
        return exists(nid) && !get(nid).is_boundary();
    }

    bool DeformableSimplicialComplex::is_unsafe_editable(const edge_key& eid) {
        return exists(eid) && !get(eid).is_boundary();
    }

    bool DeformableSimplicialComplex::is_unsafe_editable(const face_key& fid) {
        return exists(fid) && !get(fid).is_boundary();
    }

    bool DeformableSimplicialComplex::is_unsafe_editable(const tet_key& tid) {
        return exists(tid);
    }

    bool DeformableSimplicialComplex::is_safe_editable(const node_key& nid) {
        return is_unsafe_editable(nid) && !get(nid).is_interface();
    }

    bool DeformableSimplicialComplex::is_safe_editable(const edge_key& eid) {
        return is_unsafe_editable(eid) && !get(eid).is_interface();
    }

    bool DeformableSimplicialComplex::is_safe_editable(const face_key& fid) {
        return is_unsafe_editable(fid) && !get(fid).is_interface();
    }

    bool DeformableSimplicialComplex::is_safe_editable(const tet_key& tid) {
        return is_unsafe_editable(tid);
    }

    bool DeformableSimplicialComplex::is_movable(const node_key& nid) {
        return is_unsafe_editable(nid) && get(nid).is_interface() && !get(nid).is_crossing();
    }

    void DeformableSimplicialComplex::set_pos(const node_key& nid, const vec3& p) {
        get(nid).set_pos(p);
        if(!is_movable(nid))
        {
            get(nid).set_destination(p);
        }
    }

    void DeformableSimplicialComplex::set_destination(const node_key& nid, const vec3& dest) {
        if(is_movable(nid))
        {
            vec3 p = get_pos(nid);
            vec3 vec = dest - p;
            design_domain.clamp_vector(p, vec);
            get(nid).set_destination(p + vec);
        }
        else {
            get(nid).set_destination(get(nid).get_pos());
        }
    }

    vec3 DeformableSimplicialComplex::get_center() const {
        return vec3(0.);
    }

    double DeformableSimplicialComplex::get_min_tet_quality() const {
        return pars.MIN_TET_QUALITY;
    }

    double DeformableSimplicialComplex::get_deg_tet_quality() const {
        return pars.DEG_TET_QUALITY;
    }

    double DeformableSimplicialComplex::get_deg_face_quality() const {
        return pars.DEG_FACE_QUALITY;
    }

    double DeformableSimplicialComplex::get_min_face_quality() const {
        return pars.MIN_FACE_QUALITY;
    }

    double DeformableSimplicialComplex::get_avg_edge_length() const {
        return AVG_LENGTH;
    }

    const is_mesh::MultipleGeometry &DeformableSimplicialComplex::get_design_domain() const {
        return design_domain;
    }

    double DeformableSimplicialComplex::build_table(const edge_key& e, const is_mesh::SimplexSet<node_key>& polygon, std::vector<std::vector<int>>& K) {
        is_mesh::SimplexSet<node_key> nids = get_nodes(e);

        const int m = (int) polygon.size();

        auto Q = std::vector<std::vector<double>>(m-1, std::vector<double>(m, 0.));
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
                    double q2 = Util::quality(get_pos(polygon[i]), get_pos(polygon[k]), get_pos(nids[0]), get_pos(polygon[j]));
                    double q1 = Util::quality(get_pos(polygon[k]), get_pos(polygon[i]), get_pos(nids[1]), get_pos(polygon[j]));
                    double q = Util::min(q1, q2);
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

    is_mesh::NodeKey DeformableSimplicialComplex::get_next(const is_mesh::NodeKey& nid, is_mesh::SimplexSet<is_mesh::EdgeKey>& eids) {
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

    is_mesh::SimplexSet<is_mesh::NodeKey> DeformableSimplicialComplex::get_polygon(is_mesh::SimplexSet<is_mesh::EdgeKey>& eids) {
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

    std::vector<is_mesh::SimplexSet<is_mesh::NodeKey>> DeformableSimplicialComplex::get_polygons(const is_mesh::EdgeKey& eid)
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

        std::vector<is_mesh::SimplexSet<is_mesh::NodeKey>> polygons;
        is_mesh::SimplexSet<is_mesh::EdgeKey> m_eids = get_edges(get_faces(eid));
        for(auto& tids : tid_groups)
        {
            is_mesh::SimplexSet<is_mesh::EdgeKey> eids = get_edges(tids) - m_eids;
            is_mesh::SimplexSet<is_mesh::NodeKey> polygon = get_polygon(eids);
            check_consistency(get_nodes(eid), polygon);
            polygons.push_back(polygon);
        }

        struct {
            bool operator()(const is_mesh::SimplexSet<is_mesh::NodeKey>& a, const is_mesh::SimplexSet<is_mesh::NodeKey>& b)
            {
                return a.size() > b.size();
            }
        } compare;
        std::sort(polygons.begin(), polygons.end(), compare);

        return polygons;
    }

    void DeformableSimplicialComplex::flip_23_recursively(const is_mesh::SimplexSet<node_key>& polygon, const node_key& n1, const node_key& n2, std::vector<std::vector<int>>& K, int i, int j) {
        if(j >= i+2)
        {
            int k = K[i][j];
            flip_23_recursively(polygon, n1, n2, K, i, k);
            flip_23_recursively(polygon, n1, n2, K, k, j);
            flip_23(get_face(n1, n2, polygon[k]));
        }
    }

    void DeformableSimplicialComplex::topological_edge_removal(const is_mesh::SimplexSet<node_key>& polygon, const node_key& n1, const node_key& n2, std::vector<std::vector<int>>& K) {
        const int m = static_cast<int>(polygon.size());
        int k = K[0][m-1];
        flip_23_recursively(polygon, n1, n2, K, 0, k);
        flip_23_recursively(polygon, n1, n2, K, k, m-1);
        flip_32(get_edge(n1, n2));
    }

    bool DeformableSimplicialComplex::topological_edge_removal(const edge_key& eid) {
        std::vector<is_mesh::SimplexSet<node_key>> polygon = get_polygons(eid);
#ifdef DEBUG
        assert(polygon.size() == 1 && polygon.front().size() > 2);
#endif

        std::vector<std::vector<int>> K;
        double q_new = build_table(eid, polygon.front(), K);

        if (q_new > min_quality(get_tets(eid)))
        {
            const is_mesh::SimplexSet<node_key>& nodes = get_nodes(eid);
            topological_edge_removal(polygon.front(), nodes[0], nodes[1], K);
            return true;
        }
        return false;
    }

    void DeformableSimplicialComplex::topological_boundary_edge_removal(const is_mesh::SimplexSet<node_key>& polygon1, const is_mesh::SimplexSet<node_key>& polygon2, const edge_key& eid, std::vector<std::vector<int>>& K1, std::vector<std::vector<int>>& K2) {
        auto nids = get_nodes(eid);
        const int m1 = static_cast<int>(polygon1.size());
        const int m2 = static_cast<int>(polygon2.size());
        int k = K1[0][m1-1];
        flip_23_recursively(polygon1, nids[0], nids[1], K1, 0, k);
        flip_23_recursively(polygon1, nids[0], nids[1], K1, k, m1-1);

        if(m2 <= 2) {
            // Find the faces to flip about.
            face_key f1 = get_face(nids[0], nids[1], polygon1.front());
            face_key f2 = get_face(nids[0], nids[1], polygon1.back());
#ifdef DEBUG
            assert(get(f1).is_boundary() && get(f2).is_boundary());
#endif

            if(precond_flip_edge(get_edge(f1, f2), f1, f2))
            {
                flip_22(f1, f2);
            }
        }
        else {
            k = K2[0][m2-1];
            flip_23_recursively(polygon2, nids[0], nids[1], K2, 0, k);
            flip_23_recursively(polygon2, nids[0], nids[1], K2, k, m2-1);

            // Find the faces to flip about.
            face_key f1 = get_face(nids[0], nids[1], polygon1.front());
            face_key f2 = get_face(nids[0], nids[1], polygon1.back());

            if(precond_flip_edge(get_edge(f1, f2), f1, f2))
            {
                flip_44(f1, f2);
            }
        }
    }

    bool DeformableSimplicialComplex::topological_boundary_edge_removal(const edge_key& eid) {
        std::vector<is_mesh::SimplexSet<node_key>> polygons = get_polygons(eid);

        if(polygons.size() > 2 || polygons[0].size() <= 2)
        {
            return false;
        }

        std::vector<std::vector<int>> K1, K2;
        double q_new = build_table(eid, polygons[0], K1);

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

    void DeformableSimplicialComplex::topological_edge_removal() {
        std::vector<tet_key> tets;
        for (auto tit = tetrahedra_begin(); tit != tetrahedra_end(); tit++)
        {
            if (quality(tit.key()) < pars.MIN_TET_QUALITY)
            {
                tets.push_back(tit.key());
            }
        }

        // Attempt to remove each edge of each tetrahedron in tets. Accept if it increases the minimum quality locally.
        int i = 0, j = 0, k = 0;
        for (auto &t : tets)
        {
            if (is_unsafe_editable(t) && quality(t) < pars.MIN_TET_QUALITY)
            {
                for (auto e : get_edges(t))
                {
                    if(is_safe_editable(e))
                    {
                        if(topological_edge_removal(e))
                        {
                            i++;
                            break;
                        }
                    }
                    else if(exists(e) && (get(e).is_interface() || get(e).is_boundary()) && is_flippable(e))
                    {
                        if(topological_boundary_edge_removal(e))
                        {
                            k++;
                            break;
                        }
                    }
                }
                j++;
            }
        }
#ifdef DEBUG
        std::cout << "Topological edge removals: " << i + k << "/" << j << " (" << k << " at interface)" << std::endl;
#endif
        garbage_collect();
    }

    is_mesh::SimplexSet<is_mesh::EdgeKey> DeformableSimplicialComplex::test_neighbour(const face_key& f, const node_key& a, const node_key& b, const node_key& u, const node_key& w, double& q_old, double& q_new) {
        is_mesh::EdgeKey e = get_edge(u,w);
        is_mesh::SimplexSet<face_key> g_set = get_faces(e) - get_faces(get_tets(f));
        double q = Util::quality(get_pos(a), get_pos(b), get_pos(w), get_pos(u));

        if(g_set.size() == 1 && is_safe_editable(e))
        {
            face_key g = g_set.front();
            node_key v = (get_nodes(g) - get_nodes(e)).front();
            double V_uv = Util::signed_volume(get_pos(a), get_pos(b), get_pos(v), get_pos(u));
            double V_vw = Util::signed_volume(get_pos(a), get_pos(b), get_pos(w), get_pos(v));
            double V_wu = Util::signed_volume(get_pos(a), get_pos(b), get_pos(u), get_pos(w));

            if((V_uv >= EPSILON && V_vw >= EPSILON) || (V_vw >= EPSILON && V_wu >= EPSILON) || (V_wu >= EPSILON && V_uv >= EPSILON))
            {
                q_old = Util::min(Util::quality(get_pos(a), get_pos(u), get_pos(w), get_pos(v)),
                        Util::quality(get_pos(u), get_pos(v), get_pos(b), get_pos(w)));

                double q_uv_old, q_uv_new, q_vw_old, q_vw_new;
                is_mesh::SimplexSet<edge_key> uv_edges = test_neighbour(g, a, b, u, v, q_uv_old, q_uv_new);
                is_mesh::SimplexSet<edge_key> vw_edges = test_neighbour(g, a, b, v, w, q_vw_old, q_vw_new);

                q_old = Util::min(Util::min(q_old, q_uv_old), q_vw_old);
                q_new = Util::min(q_uv_new, q_vw_new);

                if(q_new > q_old || q_new > q)
                {
                    is_mesh::SimplexSet<edge_key> edges = {get_edge(f, g)};
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

    bool DeformableSimplicialComplex::topological_face_removal(const face_key& f) {
        is_mesh::SimplexSet<node_key> nids = get_nodes(f);
        is_mesh::SimplexSet<node_key> apices = get_nodes(get_tets(f)) - nids;
        this->orient_cc(apices[0], nids);

        double q_01_new, q_01_old, q_12_new, q_12_old, q_20_new, q_20_old;
        is_mesh::SimplexSet<edge_key> e01 = test_neighbour(f, apices[0], apices[1], nids[0], nids[1], q_01_old, q_01_new);
        is_mesh::SimplexSet<edge_key> e12 = test_neighbour(f, apices[0], apices[1], nids[1], nids[2], q_12_old, q_12_new);
        is_mesh::SimplexSet<edge_key> e20 = test_neighbour(f, apices[0], apices[1], nids[2], nids[0], q_20_old, q_20_new);

        double q_old = Util::min(Util::min(Util::min(min_quality(get_tets(f)), q_01_old), q_12_old), q_20_old);
        double q_new = Util::min(Util::min(q_01_new, q_12_new), q_20_new);

        if(q_new > q_old)
        {
            flip_23(f);
            for(auto &e : e01)
            {
                flip_32(e);
            }
            for(auto &e : e12)
            {
                flip_32(e);
            }
            for(auto &e : e20)
            {
                flip_32(e);
            }
            return true;
        }
        return false;
    }

    bool DeformableSimplicialComplex::topological_face_removal(const node_key& apex1, const node_key& apex2) {
        is_mesh::SimplexSet<face_key> fids = get_faces(get_tets(apex1)) & get_faces(get_tets(apex2));
        vec3 p = get_pos(apex1);
        vec3 ray = get_pos(apex2) - p;
        for(auto f : fids)
        {
            if(is_safe_editable(f))
            {
                auto nids = get_nodes(f);
                this->orient_cc(apex2, nids);

                double t = Util::intersection_ray_triangle(p, ray, get_pos(nids[0]), get_pos(nids[1]), get_pos(nids[2]));
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

    void DeformableSimplicialComplex::topological_face_removal() {
        std::vector<tet_key> tets;
        for (auto tit = tetrahedra_begin(); tit != tetrahedra_end(); tit++)
        {
            if (quality(tit.key()) < pars.MIN_TET_QUALITY)
            {
                tets.push_back(tit.key());
            }
        }

        // Attempt to remove each face of each remaining tetrahedron in tets using multi-face removal.
        // Accept if it increases the minimum quality locally.
        int i = 0, j = 0;
        for (auto &t : tets)
        {
            if (is_unsafe_editable(t) && quality(t) < pars.MIN_TET_QUALITY)
            {
                for (auto f : get_faces(t))
                {
                    if (is_safe_editable(f))
                    {
                        auto apices = get_nodes(get_tets(f)) - get_nodes(f);
                        if(topological_face_removal(apices[0], apices[1]))
                        {
                            i++;
                            break;
                        }
                    }
                }
                j++;
            }
        }
#ifdef DEBUG
        std::cout << "Topological face removals: " << i << "/" << j << std::endl;
#endif

        garbage_collect();
    }

    void DeformableSimplicialComplex::thickening_interface() {
        if(pars.MAX_LENGTH == INFINITY)
        {
            return;
        }

        std::vector<edge_key> edges;
        for (auto eit = edges_begin(); eit != edges_end(); eit++)
        {
            if ((eit->is_interface() || eit->is_boundary()) && length(eit.key()) > pars.MAX_LENGTH*AVG_LENGTH)
            {
                edges.push_back(eit.key());
            }
        }
        int i = 0;
        for(auto &e : edges)
        {
            if (exists(e) && length(e) > pars.MAX_LENGTH*AVG_LENGTH && !is_flat(get_faces(e)))
            {
                split(e);
                i++;
            }
        }
#ifdef DEBUG
        std::cout << "Thickening interface splits: " << i << std::endl;
#endif
    }

    void DeformableSimplicialComplex::thickening() {
        if(pars.MAX_VOLUME == INFINITY)
        {
            return;
        }

        std::vector<tet_key> tetrahedra;
        for (auto tit = tetrahedra_begin(); tit != tetrahedra_end(); tit++)
        {
            if (volume(tit.key()) > pars.MAX_VOLUME*AVG_VOLUME)
            {
                tetrahedra.push_back(tit.key());
            }
        }
        int i = 0;
        for(auto &t : tetrahedra)
        {
            if (is_unsafe_editable(t) && volume(t) > pars.MAX_VOLUME*AVG_VOLUME)
            {
                split(t);
                i++;
            }
        }
#ifdef DEBUG
        std::cout << "Thickening splits: " << i << std::endl;
#endif
    }

    void DeformableSimplicialComplex::thinning_interface() {
        if(pars.MIN_LENGTH <= 0.)
        {
            return;
        }

        std::vector<edge_key> edges;
        for (auto eit = edges_begin(); eit != edges_end(); eit++)
        {
            if ((eit->is_interface() || eit->is_boundary()) && length(eit.key()) < pars.MIN_LENGTH*AVG_LENGTH)
            {
                edges.push_back(eit.key());
            }
        }
        int i = 0, j = 0;
        for(auto &e : edges)
        {
            if (exists(e) && length(e) < pars.MIN_LENGTH*AVG_LENGTH)
            {
                if(collapse(e))
                {
                    i++;
                }
                j++;
            }
        }
#ifdef DEBUG
        std::cout << "Thinning interface collapses: " << i << "/" << j << std::endl;
#endif
    }

    void DeformableSimplicialComplex::thinning() {
        if(pars.MIN_VOLUME <= 0.)
        {
            return;
        }

        std::vector<tet_key> tetrahedra;
        for (auto tit = tetrahedra_begin(); tit != tetrahedra_end(); tit++)
        {
            if (volume(tit.key()) < pars.MIN_VOLUME*AVG_VOLUME)
            {
                tetrahedra.push_back(tit.key());
            }
        }
        int i = 0, j = 0;
        for(auto &t : tetrahedra)
        {
            if (is_unsafe_editable(t) && volume(t) < pars.MIN_VOLUME*AVG_VOLUME)
            {
                if(collapse(t))
                {
                    i++;
                }
                j++;
            }
        }
#ifdef DEBUG
        std::cout << "Thinning collapses: " << i << "/" << j << std::endl;
#endif
    }

    void DeformableSimplicialComplex::remove_degenerate_edges() {
        std::list<edge_key> edges;
        for (auto eit = edges_begin(); eit != edges_end(); eit++)
        {
            if (quality(eit.key()) < pars.DEG_EDGE_QUALITY)
            {
                edges.push_back(eit.key());
            }
        }
        int i = 0, j = 0;
        for(auto e : edges)
        {
            if(exists(e) && quality(e) < pars.DEG_EDGE_QUALITY && !collapse(e))
            {
                if(collapse(e, false))
                {
                    i++;
                }
                j++;
            }
        }
#ifdef DEBUG
        std::cout << "Removed " << i <<"/"<< j << " degenerate edges" << std::endl;
#endif
        garbage_collect();
    }

    void DeformableSimplicialComplex::remove_degenerate_faces() {
        std::list<face_key> faces;

        for (auto fit = faces_begin(); fit != faces_end(); fit++)
        {
            if(quality(fit.key()) < pars.DEG_FACE_QUALITY)
            {
                faces.push_back(fit.key());
            }
        }

        int i = 0, j = 0;
        for (auto &f : faces)
        {
            if (exists(f) && quality(f) < pars.DEG_FACE_QUALITY && !collapse(f))
            {
                if(collapse(f, false))
                {
                    i++;
                }
                else {
                    edge_key e = longest_edge(get_edges(f));
                    if(length(e) > AVG_LENGTH)
                    {
                        split(e);
                    }
                }
                j++;
            }
        }
#ifdef DEBUG
        std::cout << "Removed " << i <<"/"<< j << " degenerate faces" << std::endl;
#endif
        garbage_collect();
    }

    void DeformableSimplicialComplex::remove_degenerate_tets() {
        std::vector<tet_key> tets;

        for (auto tit = tetrahedra_begin(); tit != tetrahedra_end(); tit++)
        {
            if (quality(tit.key()) < pars.DEG_TET_QUALITY)
            {
                tets.push_back(tit.key());
            }
        }
        int i = 0, j = 0;
        for (auto &t : tets)
        {
            if (exists(t) && quality(t) < pars.DEG_TET_QUALITY && !collapse(t))
            {
                if(collapse(t, false))
                {
                    i++;
                }
                else {
                    edge_key e = longest_edge(get_edges(t));
                    if(length(e) > AVG_LENGTH)
                    {
                        split(e);
                    }
                }
                j++;
            }
        }
#ifdef DEBUG
        std::cout << "Removed " << i <<"/"<< j << " degenerate tets" << std::endl;
#endif
        garbage_collect();
    }

    void DeformableSimplicialComplex::remove_edges() {
        std::list<edge_key> edges;
        for (auto eit = edges_begin(); eit != edges_end(); eit++)
        {
            if (quality(eit.key()) < pars.MIN_EDGE_QUALITY)
            {
                edges.push_back(eit.key());
            }
        }
        int i = 0, j = 0;
        for(auto e : edges)
        {
            if(is_unsafe_editable(e) && quality(e) < pars.MIN_EDGE_QUALITY)
            {
                if(collapse(e))
                {
                    i++;
                }
                j++;
            }
        }
#ifdef DEBUG
        std::cout << "Removed " << i <<"/"<< j << " low quality edges" << std::endl;
#endif
        garbage_collect();
    }

    bool DeformableSimplicialComplex::remove_cap(const face_key& fid) {
        // Find longest edge
        edge_key eid = longest_edge(get_edges(fid));

        // Find apex
        node_key apex = (get_nodes(fid) - get_nodes(eid)).front();
        // Find the projected position of the apex
        auto verts = get_pos(get_nodes(eid));
        vec3 p = Util::project_point_line(get_pos(apex), verts[1], verts[0]);

        // Split longest edge
        // node_key n = split(eid, p, p);
        throw std::string{"Unsupported"};

        // Collapse new edge
        //edge_key e_rem = get_edge(apex, n);
        //return collapse(e_rem);

    }

    bool DeformableSimplicialComplex::remove_needle(const face_key& fid) {
        // Find shortest edge
        edge_key e = shortest_edge(get_edges(fid));

        // Remove edge
        return collapse(e);
    }

    bool DeformableSimplicialComplex::remove_face(const face_key& f) {
        if(max_angle(f) > 0.9*M_PI)
        {
            return remove_cap(f);
        }
        return remove_needle(f);
    }

    void DeformableSimplicialComplex::remove_faces() {
        std::list<face_key> faces;

        for (auto fit = faces_begin(); fit != faces_end(); fit++)
        {
            if(quality(fit.key()) < pars.MIN_FACE_QUALITY)
            {
                faces.push_back(fit.key());
            }
        }

        int i = 0, j = 0;
        for (auto &f : faces)
        {
            if (is_unsafe_editable(f) && quality(f) < pars.MIN_FACE_QUALITY)
            {
                if(remove_face(f))
                {
                    i++;
                }
                j++;
            }
        }
#ifdef DEBUG
        std::cout << "Removed " << i <<"/"<< j << " low quality faces" << std::endl;
#endif
        garbage_collect();
    }

    bool DeformableSimplicialComplex::remove_sliver(const tet_key& tid) {
        is_mesh::SimplexSet<edge_key> eids = get_edges(tid);
        edge_key e1 = longest_edge(eids);
        eids -= e1;
        edge_key e2 = longest_edge(eids);

        //node_key n1 = split(e1);
        //node_key n2 = split(e2);

        //edge_key e = get_edge(n1, n2);
        //return collapse(e);
        throw std::string("Not implemented");
    }

    bool DeformableSimplicialComplex::remove_cap(const tet_key& tid) {
        // Find the largest face
        face_key fid = largest_face(get_faces(tid));

        // Find the apex
        node_key apex = (get_nodes(tid) - get_nodes(fid)).front();

        // Project the apex
        auto verts = get_pos(get_nodes(fid));
        vec3 p = Util::project_point_plane(get_pos(apex), verts[0], verts[1], verts[2]);

        // Split the face
//        node_key n = split(fid, p, p);

        // Collapse edge
//        edge_key e = get_edge(n, apex);
//        return collapse(e);
        throw std::string{"Not implemented"};
    }

    bool DeformableSimplicialComplex::remove_wedge(const tet_key& tid) {
        is_mesh::SimplexSet<edge_key> eids = get_edges(tid);
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
        //        closure(t, cl_t);
        //        edge_key e1 = longest_edge(cl_t);
        //        cl_t.erase(e1);
        //        edge_key e2 = longest_edge(cl_t);
        //
        //        node_key n1 = split(e1);
        //        node_key n2 = split(e2);
        //
        //        edge_key e = get_edge(n1, n2);
        //        return collapse(e);
    }

    bool DeformableSimplicialComplex::remove_needle(const tet_key& tid) {
        is_mesh::SimplexSet<edge_key> eids = get_edges(tid);
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

    bool DeformableSimplicialComplex::remove_tet(const tet_key& tid) {
        // Find the largest face
        is_mesh::SimplexSet<face_key> fids = get_faces(tid);
        face_key fid = largest_face(fids);
        is_mesh::SimplexSet<node_key> nids = get_nodes(fid);

        // Find the apex
        node_key apex = (get_nodes(tid) - nids).front();

        // Project the apex
        auto verts = get_pos(nids);
        vec3 proj_apex = Util::project_point_plane(get_pos(apex), verts[0], verts[1], verts[2]);

        // Find barycentric coordinates
        std::vector<double> barycentric_coords = Util::barycentric_coords(proj_apex, verts[0], verts[1], verts[2]);

        if(barycentric_coords[0] > 0.2 && barycentric_coords[1] > 0.2 && barycentric_coords[2] > 0.2) // The tetrahedron is a cap
        {
            return remove_cap(tid);
        }
        else if(barycentric_coords[0] < -0.2 || barycentric_coords[1] < -0.2 || barycentric_coords[2] < -0.2) // The tetrahedron is a sliver
        {
            return remove_sliver(tid);
        }

        double mean_dist = 0.;
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

    void DeformableSimplicialComplex::remove_tets() {
        std::vector<tet_key> tets;

        for (auto tit = tetrahedra_begin(); tit != tetrahedra_end(); tit++)
        {
            if (quality(tit.key()) < pars.MIN_TET_QUALITY)
            {
                tets.push_back(tit.key());
            }
        }
        int i = 0, j=0;
        for (auto &tet : tets)
        {
            if (is_unsafe_editable(tet) && quality(tet) < pars.MIN_TET_QUALITY)
            {
                if(remove_tet(tet))
                {
                    i++;
                }
                j++;
            }
        }
#ifdef DEBUG
        std::cout << "Removed " << i <<"/"<< j << " low quality tets" << std::endl;
#endif
        garbage_collect();
    }

    bool DeformableSimplicialComplex::smart_laplacian(const node_key& nid, double alpha) {
        is_mesh::SimplexSet<tet_key> tids = get_tets(nid);
        is_mesh::SimplexSet<face_key> fids = get_faces(tids) - get_faces(nid);

        vec3 old_pos = get_pos(nid);
        vec3 avg_pos = get_barycenter(get_nodes(fids));
        vec3 new_pos = old_pos + alpha * (avg_pos - old_pos);

        double q_old, q_new;
        min_quality(fids, old_pos, new_pos, q_old, q_new);
        if(q_new > pars.MIN_TET_QUALITY || q_new > q_old)
        {
            set_pos(nid, new_pos);
            return true;
        }
        return false;
    }

    void DeformableSimplicialComplex::smooth() {
        int i = 0, j = 0;
        for (auto nit = nodes_begin(); nit != nodes_end(); nit++)
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
#ifdef DEBUG
        std::cout << "Smoothed: " << i << "/" << j << std::endl;
#endif
    }

    void DeformableSimplicialComplex::fix_complex() {
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

    void DeformableSimplicialComplex::resize_complex() {
        thickening_interface();

        thinning_interface();

        thickening();

        thinning();

        fix_complex();
    }

    void DeformableSimplicialComplex::deform(int num_steps) {
#ifdef DEBUG
        validity_check();
        std::cout << std::endl << "********************************" << std::endl;
#endif
        int missing;
        int step = 0;
        do {
            std::cout << "\n\tMove vertices step " << step << std::endl;
            missing = 0;
            int movable = 0;
            for (auto nit = nodes_begin(); nit != nodes_end(); nit++)
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
            std::cout << "\tVertices missing to be moved: " << missing <<"/" << movable << std::endl;
            fix_complex();
#ifdef DEBUG
            validity_check();
#endif
            ++step;
        } while (missing > 0 && step < num_steps);

        resize_complex();

        garbage_collect();
        for (auto nit = nodes_begin(); nit != nodes_end(); nit++)
        {
            nit->set_destination(nit->get_pos());
        }
#ifdef DEBUG
        validity_check();
#endif
    }

    bool DeformableSimplicialComplex::move_vertex(const node_key & n) {
        vec3 pos = get_pos(n);
        vec3 destination = get(n).get_destination();
        double l = Util::length(destination - pos);

        if (l < 1e-4*AVG_LENGTH) // The vertex is not moved
        {
            return true;
        }

        double max_l = l*intersection_with_link(n, destination) - 1e-4 * AVG_LENGTH;
        l = Util::max(Util::min(0.5*max_l, l), 0.);
        set_pos(n, pos + l*Util::normalize(destination - pos));

        if (Util::length(destination - get_pos(n)) < 1e-4*AVG_LENGTH)
        {
            return true;
        }
        return false;
    }

    double DeformableSimplicialComplex::intersection_with_link(const node_key & n, const vec3& destination) {
        vec3 pos = get_pos(n);
        vec3 ray = destination - pos;

        double min_t = INFINITY;
        auto fids = get_faces(get_tets(n)) - get_faces(n);
        for(auto f : fids)
        {
            auto face_pos = get_pos(get_nodes(f));
            double t = Util::intersection_ray_plane(pos, ray, face_pos[0], face_pos[1], face_pos[2]);
            if (0. <= t)
            {
                min_t = Util::min(t, min_t);
            }
        }
#ifdef DEBUG
        assert(min_t < INFINITY);
#endif
        return min_t;
    }

    bool DeformableSimplicialComplex::is_flat(const is_mesh::SimplexSet<face_key>& fids) {
        for (const face_key& f1 : fids) {
            if (get(f1).is_interface() || get(f1).is_boundary())
            {
                vec3 normal1 = get_normal(f1);
                for (const face_key& f2 : fids) {
                    if (f1 != f2 && (get(f2).is_interface() || get(f2).is_boundary()))
                    {
                        vec3 normal2 = get_normal(f2);
                        if(std::abs(dot(normal1, normal2)) < FLIP_EDGE_INTERFACE_FLATNESS)
                        {
                            return false;
                        }
                    }
                }
            }
        }
        return true;
    }

    bool DeformableSimplicialComplex::is_flippable(const edge_key & eid) {
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
#ifdef DEBUG
        assert(new_e_nids.size() == 2);
#endif

        // Check that there does not already exist an edge.
        if(get_edge(new_e_nids[0], new_e_nids[1]).is_valid())
        {
            return false;
        }

        // Check that the edge is not a feature edge if it is a part of the interface or boundary.
        if(get(eid).is_interface() || get(eid).is_boundary())
        {
            return is_flat(fids);
        }

        return true;
    }

    bool DeformableSimplicialComplex::precond_flip_edge(const edge_key& eid, const face_key& f1, const face_key& f2) {
        is_mesh::SimplexSet<node_key> e_nids = get_nodes(eid);
        is_mesh::SimplexSet<node_key> new_e_nids = (get_nodes(f1) + get_nodes(f2)) - e_nids;
        is_mesh::SimplexSet<node_key> apices = (get_nodes(get_faces(eid)) - e_nids) - new_e_nids;
#ifdef DEBUG
        assert(e_nids.size() == 2);
        assert(new_e_nids.size() == 2);
#endif

        // Check that there does not already exist an edge.
        if(get_edge(new_e_nids[0], new_e_nids[1]).is_valid())
        {
            return false;
        }

        vec3 p = get_pos(new_e_nids[0]);
        vec3 r = get_pos(new_e_nids[1]) - p;
        vec3 a = get_pos(e_nids[0]);
        vec3 b = get_pos(e_nids[1]);

        for (node_key n : apices) {
            vec3 c = get_pos(n);
            double t = Util::intersection_ray_plane(p, r, a, b, c);
            if(t > 0. && t < 1.)
            {
                std::vector<double> coords = Util::barycentric_coords(p + t*r, c, a, b);
                if(coords[0] > EPSILON && coords[1] > EPSILON && coords[2] >= 0.)
                {
                    return true;
                }
            }
        }

        return false;
    }

    void DeformableSimplicialComplex::split(const tet_key& tid) {
        is_mesh::SimplexSet<edge_key> eids = get_edges(tid);
        edge_key eid = longest_edge(eids);
        split(eid);
    }

    void DeformableSimplicialComplex::split(const face_key& fid) {
        is_mesh::SimplexSet<edge_key> eids = get_edges(fid);
        edge_key eid = longest_edge(eids);
        split(eid);
    }

    void DeformableSimplicialComplex::split(const EdgeKey& eid) {
        auto verts = get_pos(get_nodes(eid));
        vec3 pos = Util::barycenter(verts[0], verts[1]);
        vec3 destination = pos;
        if(get(eid).is_interface())
        {
            auto nids = get_nodes(eid);
            destination = Util::barycenter(get(nids[0]).get_destination(), get(nids[1]).get_destination());
        }

        //split(eid, pos, destination);
        throw std::string{"Not implemented"};
    }

    bool DeformableSimplicialComplex::is_collapsable(const edge_key& eid, const node_key& nid, bool safe) {
        if(safe)
        {
            if(is_safe_editable(nid))
            {
                return true;
            }
        }
        else {
            if(is_unsafe_editable(nid))
            {
                return true;
            }
        }
        if(get(eid).is_boundary() || get(eid).is_interface())
        {
            return is_flat(get_faces(nid));
        }
        return false;
    }

    bool DeformableSimplicialComplex::collapse(const EdgeKey& eid, bool safe) {
        is_mesh::SimplexSet<node_key> nids = get_nodes(eid);
        bool n0_is_editable = is_collapsable(eid, nids[0], safe);
        bool n1_is_editable = is_collapsable(eid, nids[1], safe);

        if (!n0_is_editable && !n1_is_editable)
        {
            return false;
        }
        std::vector<double> test_weights;
        if (!n0_is_editable || !n1_is_editable)
        {
            test_weights = {0.};
            if(!n0_is_editable)
            {
                nids.swap();
            }
        }
        else {
            test_weights = {0., 0.5, 1.};
        }

        is_mesh::SimplexSet<tet_key> e_tids = get_tets(eid);
        is_mesh::SimplexSet<face_key> fids0 = get_faces(get_tets(nids[0]) - e_tids) - get_faces(nids[0]);
        is_mesh::SimplexSet<face_key> fids1 = get_faces(get_tets(nids[1]) - e_tids) - get_faces(nids[1]);

        double q_max = -INFINITY;
        double weight;
        for (double w : test_weights)
        {
            vec3 p = (1.-w) * get(nids[1]).get_pos() + w * get(nids[0]).get_pos();
            double q = Util::min(min_quality(fids0, get(nids[0]).get_pos(), p), min_quality(fids1, get(nids[1]).get_pos(), p));

            if (q > q_max && ((!get(nids[0]).is_interface() && !get(nids[1]).is_interface()) || design_domain.is_inside(p)))
            {
                q_max = q;
                weight = w;
            }
        }

        if(q_max > EPSILON)
        {
            if(!safe || q_max > Util::min(min_quality(get_tets(nids[0]) + get_tets(nids[1])), pars.MIN_TET_QUALITY) + EPSILON)
            {
                throw std::string{"Not implemented"};
                //collapse(eid, nids[1], weight);
                return true;
            }
        }
        return false;
    }

    bool DeformableSimplicialComplex::collapse(is_mesh::SimplexSet<edge_key>& eids, bool safe) {
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

    bool DeformableSimplicialComplex::collapse(const face_key& fid, bool safe) {
        is_mesh::SimplexSet<edge_key> eids = get_edges(fid);
        return collapse(eids, safe);
    }

    bool DeformableSimplicialComplex::collapse(const tet_key& tid, bool safe) {
        is_mesh::SimplexSet<edge_key> eids = get_edges(tid);
        return collapse(eids, safe);
    }

    std::vector<vec3> DeformableSimplicialComplex::get_interface_face_positions() {
        std::vector<vec3> verts;
        for (auto fit = faces_begin(); fit != faces_end(); fit++) {
            if(fit->is_interface())
            {
                for(auto n : get_nodes(fit.key()))
                {
                    verts.push_back(get_pos(n));
                }
            }
        }
        return verts;
    }

    vec3 DeformableSimplicialComplex::get_normal(const face_key& fid) {
        auto pos = get_pos(this->get_sorted_nodes(fid));
        return Util::normal_direction(pos[0], pos[1], pos[2]);
    }

    vec3 DeformableSimplicialComplex::get_normal(const node_key& nid) {
        vec3 result(0.);
        for (auto f : get_faces(nid))
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

    vec3 DeformableSimplicialComplex::get_barycenter(const is_mesh::SimplexSet<node_key>& nids, bool interface) {
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
        return avg_pos / static_cast<double>(i);
    }

    vec3 DeformableSimplicialComplex::get_barycenter(const node_key& nid, bool interface) {
        if(interface && !get(nid).is_interface())
        {
            return get_pos(nid);
        }

        is_mesh::SimplexSet<node_key> nids = get_nodes(get_tets(nid)) - nid;
        return get_barycenter(nids, interface);
    }

    double DeformableSimplicialComplex::length(const edge_key& eid) {
        is_mesh::SimplexSet<node_key> nids = get_nodes(eid);
        return Util::length(get_pos(nids[0]) - get_pos(nids[1]));
    }

    double DeformableSimplicialComplex::length_destination(const edge_key& eid) {
        is_mesh::SimplexSet<node_key> nids = get_nodes(eid);
        return Util::length(get(nids[0]).get_destination() - get(nids[1]).get_destination());
    }

    double DeformableSimplicialComplex::area(const face_key& fid) {
        is_mesh::SimplexSet<node_key> nids = get_nodes(fid);
        return Util::area(get_pos(nids[0]), get_pos(nids[1]), get_pos(nids[2]));
    }

    double DeformableSimplicialComplex::area_destination(const face_key& fid) {
        is_mesh::SimplexSet<node_key> nids = get_nodes(fid);
        return Util::area(get(nids[0]).get_destination(), get(nids[1]).get_destination(), get(nids[2]).get_destination());
    }

    double DeformableSimplicialComplex::volume(const tet_key& tid) {
        is_mesh::SimplexSet<node_key> nids = get_nodes(tid);
        return Util::volume(get_pos(nids[0]), get_pos(nids[1]), get_pos(nids[2]), get_pos(nids[3]));
    }

    double DeformableSimplicialComplex::volume_destination(const tet_key& tid) {
        is_mesh::SimplexSet<node_key> nids = get_nodes(tid);
        return Util::volume(get(nids[0]).get_destination(), get(nids[1]).get_destination(), get(nids[2]).get_destination(), get(nids[3]).get_destination());
    }

    double DeformableSimplicialComplex::volume_destination(const is_mesh::SimplexSet<node_key>& nids) {
        return Util::volume(get(nids[0]).get_destination(), get(nids[1]).get_destination(), get(nids[2]).get_destination(), get(nids[3]).get_destination());
    }

    double DeformableSimplicialComplex::signed_volume_destination(const is_mesh::SimplexSet<node_key>& nids) {
        return Util::signed_volume(get(nids[0]).get_destination(), get(nids[1]).get_destination(), get(nids[2]).get_destination(), get(nids[3]).get_destination());
    }

    vec3 DeformableSimplicialComplex::barycenter(const tet_key& tid) {
        is_mesh::SimplexSet<node_key> nids = get_nodes(tid);
        return Util::barycenter(get(nids[0]).get_pos(), get(nids[1]).get_pos(), get(nids[2]).get_pos(), get(nids[3]).get_pos());
    }

    vec3 DeformableSimplicialComplex::barycenter_destination(const tet_key& tid) {
        is_mesh::SimplexSet<node_key> nids = get_nodes(tid);
        return Util::barycenter(get(nids[0]).get_destination(), get(nids[1]).get_destination(), get(nids[2]).get_destination(), get(nids[3]).get_destination());
    }

    double DeformableSimplicialComplex::quality(const tet_key& tid) {
        is_mesh::SimplexSet<node_key> nids = get_nodes(tid);
        return std::abs(Util::quality(get_pos(nids[0]), get_pos(nids[1]), get_pos(nids[2]), get_pos(nids[3])));
    }

    double DeformableSimplicialComplex::min_angle(const face_key& fid) {
        is_mesh::SimplexSet<node_key> nids = get_nodes(fid);
        return Util::min_angle(get_pos(nids[0]), get_pos(nids[1]), get_pos(nids[2]));
    }

    double DeformableSimplicialComplex::max_angle(const face_key& fid) {
        is_mesh::SimplexSet<node_key> nids = get_nodes(fid);
        return Util::max_angle(get_pos(nids[0]), get_pos(nids[1]), get_pos(nids[2]));
    }

    double DeformableSimplicialComplex::quality(const face_key& fid) {
        is_mesh::SimplexSet<node_key> nids = get_nodes(fid);
        auto angles = Util::cos_angles(get_pos(nids[0]), get_pos(nids[1]), get_pos(nids[2]));
        double worst_a = -INFINITY;
        for(auto a : angles)
        {
            worst_a = std::max(worst_a, std::abs(a));
        }
        return 1. - worst_a;
    }

    double DeformableSimplicialComplex::quality(const edge_key& eid) {
        return length(eid)/AVG_LENGTH;
    }

    FaceKey DeformableSimplicialComplex::largest_face(const is_mesh::SimplexSet<FaceKey>& fids) {
        double max_a = -INFINITY;
        face_key max_f;
        for(auto f : fids)
        {
            double a = area(f);
            if(a > max_a)
            {
                max_a = a;
                max_f = f;
            }
        }
        return max_f;
    }

    EdgeKey DeformableSimplicialComplex::shortest_edge(const is_mesh::SimplexSet<EdgeKey>& eids) {
        double min_l = INFINITY;
        EdgeKey min_e;
        for(auto e : eids)
        {
            double l = length(e);
            if(l < min_l)
            {
                min_l = l;
                min_e = e;
            }
        }
        return min_e;
    }

    EdgeKey DeformableSimplicialComplex::longest_edge(const is_mesh::SimplexSet<EdgeKey>& eids) {
        double max_l = -INFINITY;
        EdgeKey max_e;
        for(auto e : eids)
        {
            double l = length(e);
            if(l > max_l)
            {
                max_l = l;
                max_e = e;
            }
        }
        return max_e;
    }

    double DeformableSimplicialComplex::min_quality(const is_mesh::SimplexSet<TetrahedronKey>& tids) {
        double q_min = INFINITY;
        for (auto t : tids)
        {
            q_min = Util::min(quality(t), q_min);
        }
        return q_min;
    }

    double DeformableSimplicialComplex::min_quality(const is_mesh::SimplexSet<FaceKey>& fids, const vec3& pos) {
        double min_q = INFINITY;
        for (auto f : fids)
        {
            is_mesh::SimplexSet<NodeKey> nids = get_nodes(f);
            min_q = Util::min(min_q, std::abs(Util::quality(get_pos(nids[0]), get_pos(nids[1]), get_pos(nids[2]), pos)));
        }
        return min_q;
    }

    double DeformableSimplicialComplex::min_quality(const is_mesh::SimplexSet<FaceKey>& fids, const vec3& pos_old, const vec3& pos_new) {
        double min_q = INFINITY;
        for (auto f : fids)
        {
            is_mesh::SimplexSet<NodeKey> nids = get_nodes(f);
            if(Util::sign(Util::signed_volume(get_pos(nids[0]), get_pos(nids[1]), get_pos(nids[2]), pos_old)) !=
                    Util::sign(Util::signed_volume(get_pos(nids[0]), get_pos(nids[1]), get_pos(nids[2]), pos_new)))
            {
                return -INFINITY;
            }
            min_q = Util::min(min_q, std::abs(Util::quality(get_pos(nids[0]), get_pos(nids[1]), get_pos(nids[2]), pos_new)));
        }
        return min_q;
    }

    void DeformableSimplicialComplex::min_quality(const is_mesh::SimplexSet<FaceKey>& fids, const vec3& pos_old, const vec3& pos_new, double& min_q_old, double& min_q_new) {
        min_q_old = INFINITY;
        min_q_new = INFINITY;
        for (auto f : fids)
        {
            SimplexSet<NodeKey> nids = get_nodes(f);
            if(Util::sign(Util::signed_volume(get_pos(nids[0]), get_pos(nids[1]), get_pos(nids[2]), pos_old)) !=
                    Util::sign(Util::signed_volume(get_pos(nids[0]), get_pos(nids[1]), get_pos(nids[2]), pos_new)))
            {
                min_q_old = INFINITY;
                min_q_new = -INFINITY;
                break;
            }
            min_q_old = Util::min(min_q_old, std::abs(Util::quality(get_pos(nids[0]), get_pos(nids[1]), get_pos(nids[2]), pos_old)));
            min_q_new = Util::min(min_q_new, std::abs(Util::quality(get_pos(nids[0]), get_pos(nids[1]), get_pos(nids[2]), pos_new)));
        }
    }

    void DeformableSimplicialComplex::check_consistency(const is_mesh::SimplexSet<NodeKey>& nids, is_mesh::SimplexSet<NodeKey>& polygon) {
        unsigned int n = static_cast<unsigned int>(polygon.size());

        double sum = 0;
        for (unsigned int i = 0; i < n; ++i)
        {
            sum += Util::signed_volume(get_pos(nids[0]), get_pos(nids[1]), get_pos(polygon[i]), get_pos(polygon[(i+1)%n]));
        }

        if (sum < 0.)
        {
            for (unsigned int i = 0; i < n/2; ++i)
            {
                polygon.swap(i, n-1-i);
            }
        }
    }

    double DeformableSimplicialComplex::compute_avg_edge_length() {
        double avg_edge_length = 0.;
        int N = 0;
        for (auto eit = edges_begin(); eit != edges_end(); eit++) {
            if(eit->is_interface())
            {
                avg_edge_length += length(eit.key());
                N++;
            }
        }
        if (N > 0) {
            avg_edge_length /= static_cast<double>(N);
        }
        return avg_edge_length;
    }

    double DeformableSimplicialComplex::cos_dihedral_angle(const face_key& f1, const face_key& f2) {
        auto nids1 = get_nodes(f1);
        auto nids2 = get_nodes(f2);
        is_mesh::SimplexSet<node_key> nids = nids1 & nids2;
        is_mesh::SimplexSet<node_key> apices = (nids1 + nids2) - nids;

        return Util::cos_dihedral_angle(get_pos(nids[0]), get_pos(nids[1]), get_pos(apices[0]), get_pos(apices[1]));
    }

    double DeformableSimplicialComplex::dihedral_angle(const face_key& f1, const face_key& f2) {
        return acos(cos_dihedral_angle(f1, f2));
    }

    std::vector<double> DeformableSimplicialComplex::cos_dihedral_angles(const tet_key& tid) {
        auto verts = get_pos(get_nodes(tid));
        std::vector<double> angles;
        std::vector<int> apices;
        for (unsigned int i = 0; i < verts.size(); i++) {
            for (unsigned int j = 0; j < verts.size(); j++) {
                if(i < j)
                {
                    apices.clear();
                    for (unsigned int k = 0; k < verts.size(); k++) {
                        if(k != i && k != j)
                        {
                            apices.push_back(k);
                        }
                    }
                    angles.push_back(Util::cos_dihedral_angle(verts[i], verts[j], verts[apices[0]], verts[apices[1]]));
                }
            }
        }
        return angles;
    }

    double DeformableSimplicialComplex::min_cos_dihedral_angle(const tet_key& t) {
        double min_angle = -1.;
        std::vector<double> angles = cos_dihedral_angles(t);
        for(auto a : angles)
        {
            min_angle = Util::max(min_angle, a);
        }
        return min_angle;
    }

    double DeformableSimplicialComplex::min_dihedral_angle(const tet_key& t) {
        return acos(min_cos_dihedral_angle(t));
    }

    void DeformableSimplicialComplex::get_qualities(std::vector<int>& histogram, double& min_quality) {
        min_quality = INFINITY;

        histogram = std::vector<int>(100);
        for (int i = 0; i < 100; ++i)
        {
            histogram[i] = 0;
        }

        for (auto tit = tetrahedra_begin(); tit != tetrahedra_end(); tit++)
        {
            double q = quality(tit.key());
            min_quality = Util::min(min_quality, q);
            int index = static_cast<int>(floor(q*100.));
#ifdef DEBUG
            assert(index < 100 && index >= 0);
#endif
            histogram[index] += 1;
        }
    }

    void DeformableSimplicialComplex::get_dihedral_angles(std::vector<int> & histogram, double & min_angle, double & max_angle) {
        max_angle = -INFINITY, min_angle = INFINITY;

        histogram = std::vector<int>(180);
        for (int i = 0; i < 180; ++i)
        {
            histogram[i] = 0;
        }

        for (auto tit = tetrahedra_begin(); tit != tetrahedra_end(); tit++)
        {
            std::vector<double> angles = cos_dihedral_angles(tit.key());
            for(auto cos_a : angles)
            {
                double a = acos(cos_a)*180./M_PI;
                min_angle = Util::min(min_angle, a);
                max_angle = Util::max(max_angle, a);
                histogram[(int)floor(a)] += 1;
            }
        }
    }

    double DeformableSimplicialComplex::min_quality() {
        double min_q = INFINITY;
        for (auto tit = tetrahedra_begin(); tit != tetrahedra_end(); tit++)
        {
            min_q = Util::min(min_q, quality(tit.key()));
        }
        return min_q;
    }

    void DeformableSimplicialComplex::count_nodes(int & total, int & object) {
        total = 0, object = 0;
        for (auto nit = nodes_begin(); nit != nodes_end(); nit++)
        {
            total++;
            if (nit->is_interface())
            {
                object++;
            }
        }
    }

    void DeformableSimplicialComplex::count_edges(int & total, int & object) {
        total = 0, object = 0;
        for (auto eit = edges_begin(); eit != edges_end(); eit++)
        {
            total++;
            if (eit->is_interface())
            {
                object++;
            }
        }
    }

    void DeformableSimplicialComplex::count_faces(int & total, int & object) {
        total = 0, object = 0;
        for (auto fit = faces_begin(); fit != faces_end(); fit++)
        {
            total++;
            if (fit->is_interface())
            {
                object++;
            }
        }
    }

    void DeformableSimplicialComplex::count_tetrahedra(int & total, int & object) {
        total = 0, object = 0;
        for (auto tit = tetrahedra_begin(); tit != tetrahedra_end(); tit++)
        {
            total++;
            if (tit->label() != 0)
            {
                object++;
            }
        }
    }

    void DeformableSimplicialComplex::test_split_collapse() {
        is_mesh::SimplexSet<edge_key> eids;
        for (auto eit = edges_begin(); eit != edges_end(); eit++)
        {
            auto neighbours = get_edges(get_faces(eit.key()));
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

        //            int j = 0;
        //            std::cout << "Split test # = " << eids.size();
        //            is_mesh::SimplexSet<edge_key> new_eids;
        //            std::vector<node_key> old_nids;
        //            for (auto e : eids) {
        //                auto nids = get_nodes(e);
        //                auto new_nid = split(e);
        //                auto new_eid = (get_edges(nids) & get_edges(new_nid)) - e;
        //                assert(new_eid.size() == 1);
        //                new_eids += new_eid[0];
        //                auto old_nid = get_nodes(new_eid) - new_nid;
        //                assert(old_nid.size() == 1);
        //                old_nids.push_back(old_nid.front());
        //                j++;
        //                if(j%1000 == 0)
        //                {
        //                    std::cout << ".";
        //                }
        //            }
        //            std::cout << " DONE" << std::endl;
        //            garbage_collect();
        //            validity_check();
        //
        //            std::cout << "Collapse test # = " << new_eids.size();
        //            j = 0;
        //            for (unsigned int i = 0; i < new_eids.size(); i++) {
        //                assert(exists(new_eids[i]));
        //                auto nids = get_nodes(new_eids[i]);
        //                collapse(new_eids[i], old_nids[i], 0.);
        //                assert(nids[0].is_valid());
        //
        //                j++;
        //                if(j%1000 == 0)
        //                {
        //                    std::cout << ".";
        //                }
        //            }
        std::cout << " DONE" << std::endl;
        garbage_collect();
        validity_check();
    }

    void DeformableSimplicialComplex::test_flip23_flip32() {
        is_mesh::SimplexSet<face_key> fids;
        for (auto fit = faces_begin(); fit != faces_end(); fit++)
        {
            if(is_safe_editable(fit.key()))
            {
                auto nids = get_nodes(fit.key());
                nids += get_nodes(get_tets(fit.key()));
                double t = Util::intersection_ray_triangle(get_pos(nids[3]), get_pos(nids[4]) - get_pos(nids[3]), get_pos(nids[0]), get_pos(nids[1]), get_pos(nids[2]));

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
            assert(exists(f));
            auto new_eid = flip_23(f);
            assert(new_eid.is_valid());
            new_eids += new_eid;
            i++;
            if(i%1000 == 0)
            {
                std::cout << ".";
            }
        }
        std::cout << " DONE" << std::endl;
        garbage_collect();
        validity_check();

        i=0;
        std::cout << "Flip 3-2 test # = " << new_eids.size();
        for (auto e : new_eids) {
            auto new_fid = flip_32(e);
            assert(new_fid.is_valid());
            i++;
            if(i%1000 == 0)
            {
                std::cout << ".";
            }
        }
        std::cout << " DONE" << std::endl;
        garbage_collect();
        validity_check();
    }

    void DeformableSimplicialComplex::test_flip44() {
        is_mesh::SimplexSet<edge_key> eids;
        for (auto eit = edges_begin(); eit != edges_end(); eit++)
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
                flip_44(flip_fids[0], flip_fids[1]);
                i++;
                if(i%100 == 0)
                {
                    std::cout << ".";
                }
            }
            std::cout << " DONE" << std::endl;
            garbage_collect();
            validity_check();
        }
    }

    void DeformableSimplicialComplex::test_flip22() {
        is_mesh::SimplexSet<edge_key> eids;
        for (auto eit = edges_begin(); eit != edges_end(); eit++)
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
                assert(exists(e));
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
                flip_22(flip_fids[0], flip_fids[1]);
                i++;
                if(i%10 == 0)
                {
                    std::cout << ".";
                }
            }
            std::cout << " DONE" << std::endl;
            garbage_collect();
            validity_check();
        }
    }
}
