#include "DSC.h"

using namespace is_mesh;
using namespace std;


namespace DSC {

    DeformableSimplicialComplex::DeformableSimplicialComplex(vector<vec3> & points, vector<int> & tets, const vector<int>& tet_labels)
            : is_mesh_ptr(std::make_shared<ISMesh>(points, tets, tet_labels)), is_mesh(*is_mesh_ptr)
    {
        pars = {0.1, 0.5, 0.0005, 0.015, 0.02, 0.3, 0., 2., 0.2, 5., 0.2, INFINITY};
        set_avg_edge_length();
        gcListenerId = is_mesh.add_gc_listener([&](const GarbageCollectDeletions& gc){
            on_gc(gc);
        });
    }

    DeformableSimplicialComplex::~DeformableSimplicialComplex() {
        is_mesh.remove_gc_listener(gcListenerId);
    }

    void DeformableSimplicialComplex::set_avg_edge_length(double avg_edge_length) {
        if(avg_edge_length <= 0.)
        {
            AVG_LENGTH = compute_avg_edge_length();
        }
        else {
            AVG_LENGTH = avg_edge_length;
        }

        AVG_AREA = 0.5*sqrt(3./4.)*AVG_LENGTH*AVG_LENGTH;
        AVG_VOLUME = (sqrt(2.)/12.)*AVG_LENGTH*AVG_LENGTH*AVG_LENGTH;
    }

    void DeformableSimplicialComplex::set_parameters(parameters pars_) {
        pars = pars_;
    }

    void DeformableSimplicialComplex::add_design_domain(std::shared_ptr<Geometry> geometry) {
        design_domain.add_geometry(geometry);
    }

    void DeformableSimplicialComplex::set_labels(const Geometry& geometry, int label) {
        for (auto & tit : tetrahedra()) {
            SimplexSet<NodeKey> nids = is_mesh.get_nodes(tit.key());
            if(geometry.is_inside(Util::barycenter(is_mesh.get(nids[0]).get_pos(), is_mesh.get(nids[1]).get_pos(), is_mesh.get(nids[2]).get_pos(), is_mesh.get(nids[3]).get_pos())))
            {
                is_mesh.set_label(tit.key(),label);
            }
        }
    }

    void DeformableSimplicialComplex::print(const NodeKey& n) {
        cout << "Node: " << n << endl;
        vec3 p = is_mesh.get(n).get_pos();
        vec3 d = is_mesh.get(n).get_destination();
        cout << "P = " << p[0] << ", " << p[1] << ", " << p[2] << endl;
        cout << "D = " << d[0] << ", " << d[1] << ", " << d[2]  << endl;

        cout << "\nStar_edges = [";
        for(auto e : is_mesh.get(n).edge_keys())
        {
            auto verts = is_mesh.get_pos(is_mesh.get(e).node_keys());
            vec3 p1 = verts[0];
            vec3 p2 = verts[1];
            cout << p1[0] << ", " << p1[1] << ", " << p1[2] << "; " << endl;
            cout << p2[0] << ", " << p2[1] << ", " << p2[2] << "; " << endl;
        }
        cout << "];" << endl;

        cout << "\nStar_Iedges = [";
        for(auto e : is_mesh.get(n).edge_keys())
        {
            if(is_mesh.get(e).is_interface())
            {
                auto verts = is_mesh.get_pos(is_mesh.get(e).node_keys());
                vec3 p1 = verts[0];
                vec3 p2 = verts[1];
                cout << p1[0] << ", " << p1[1] << ", " << p1[2] << "; " << endl;
                cout << p2[0] << ", " << p2[1] << ", " << p2[2] << "; " << endl;
            }
        }
        cout << "];" << endl;

        cout << "\nedges = [";
        auto eids = is_mesh.get_edges(is_mesh.get_tets(n)) - is_mesh.get(n).edge_keys();
        for(auto e : eids)
        {
            auto verts = is_mesh.get_pos(is_mesh.get(e).node_keys());
            vec3 p1 = verts[0];
            vec3 p2 = verts[1];
            cout << p1[0] << ", " << p1[1] << ", " << p1[2] << "; " << endl;
            cout << p2[0] << ", " << p2[1] << ", " << p2[2] << "; " << endl;
        }
        cout << "];" << endl;

        cout << "\nIedges = [";
        for(auto e : eids)
        {
            if(is_mesh.get(e).is_interface())
            {
                auto verts = is_mesh.get_pos(is_mesh.get(e).node_keys());
                vec3 p1 = verts[0];
                vec3 p2 = verts[1];
                cout << p1[0] << ", " << p1[1] << ", " << p1[2] << "; " << endl;
                cout << p2[0] << ", " << p2[1] << ", " << p2[2] << "; " << endl;
            }
        }
        cout << "];" << endl;
    }

    bool DeformableSimplicialComplex::is_unsafe_editable(const NodeKey& nid) {
        return is_mesh.exists(nid) && !is_mesh.get(nid).is_boundary();
    }

    bool DeformableSimplicialComplex::is_unsafe_editable(const EdgeKey& eid) {
        return is_mesh.exists(eid) && !is_mesh.get(eid).is_boundary();
    }

    bool DeformableSimplicialComplex::is_unsafe_editable(const FaceKey& fid) {
        return is_mesh.exists(fid) && !is_mesh.get(fid).is_boundary();
    }

    bool DeformableSimplicialComplex::is_unsafe_editable(const TetrahedronKey& tid) {
        return is_mesh.exists(tid);
    }

    bool DeformableSimplicialComplex::is_safe_editable(const NodeKey& nid) {
        return is_unsafe_editable(nid) && !is_mesh.get(nid).is_interface();
    }

    bool DeformableSimplicialComplex::is_safe_editable(const EdgeKey& eid) {
        return is_unsafe_editable(eid) && !is_mesh.get(eid).is_interface();
    }

    bool DeformableSimplicialComplex::is_safe_editable(const FaceKey& fid) {
        return is_unsafe_editable(fid) && !is_mesh.get(fid).is_interface();
    }

    bool DeformableSimplicialComplex::is_safe_editable(const TetrahedronKey& tid) {
        return is_unsafe_editable(tid);
    }

    bool DeformableSimplicialComplex::is_movable(const NodeKey& nid) {
        return is_unsafe_editable(nid) && is_mesh.get(nid).is_interface() && !is_mesh.get(nid).is_crossing();
    }

    void DeformableSimplicialComplex::set_pos(const NodeKey& nid, const vec3& p) {
        auto & n = is_mesh.get(nid);
        n.set_pos(p);
        if(!is_movable(nid))
        {
            n.set_destination(p);
        }
    }

    void DeformableSimplicialComplex::set_destination(const NodeKey& nid, const vec3& dest) {
        auto &n = is_mesh.get(nid);
        if(is_movable(nid))
        {
            vec3 p = n.get_pos();
            vec3 vec = dest - p;
            design_domain.clamp_vector(p, vec);
            n.set_destination(p + vec);
        }
        else {
            n.set_destination(n.get_pos());
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

    const MultipleGeometry &DeformableSimplicialComplex::get_design_domain() const {
        return design_domain;
    }

    double DeformableSimplicialComplex::build_table(const EdgeKey& e, const SimplexSet<NodeKey>& polygon, vector<vector<int>>& K) {
        SimplexSet<NodeKey> nids = is_mesh.get(e).node_keys();

        vector<vec3> poly_pos = is_mesh.get_pos(polygon);

        const int m = (int) polygon.size();

        auto Q = vector<vector<double>>(m-1, vector<double>(m, 0.));
        K = vector<vector<int>>(m-1, vector<int>(m, 0) );

        for(int i = 0; i < m-1; i++)
        {
            Q[i][i+1] = INFINITY;
        }

        auto polyn0 = is_mesh.get(nids[0]).get_pos();
        auto polyn1 = is_mesh.get(nids[1]).get_pos();
        for (int i = m-3; i >= 0; i--)
        {
            auto polyi = poly_pos[i];
            for (int j = i+2; j < m; j++)
            {
                auto polyj = poly_pos[j];
                for (int k = i+1; k < j; k++)
                {
                    auto polyk = poly_pos[k];
                    double q2 = Util::quality(polyi, polyk, polyn0, polyj);
                    double q1 = Util::quality(polyk, polyi, polyn1, polyj);
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

    NodeKey DeformableSimplicialComplex::get_next(const NodeKey& nid, SimplexSet<EdgeKey>& eids) {
        for (auto e : eids)
        {
            auto n = is_mesh.get(e).node_keys() - nid;
            if(n.size() == 1)
            {
                eids -= e;
                return n.front();
            }
        }
        return NodeKey();
    }

    SimplexSet<NodeKey> DeformableSimplicialComplex::get_polygon(SimplexSet<EdgeKey>& eids) {
        SimplexSet<NodeKey> polygon = {is_mesh.get(eids[0]).node_keys().front()};
        NodeKey next_nid;
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

    vector<SimplexSet<NodeKey>> DeformableSimplicialComplex::get_polygons(const EdgeKey& eid)
    {
        vector<SimplexSet<TetrahedronKey>> tid_groups;
        for (auto t : is_mesh.get_tets(eid))
        {
            bool found = false;
            for(auto& tids : tid_groups)
            {
                if(is_mesh.get(t).label() == is_mesh.get(tids.front()).label())
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

        vector<SimplexSet<NodeKey>> polygons;
        SimplexSet<EdgeKey> m_eids = is_mesh.get_edges(is_mesh.get(eid).face_keys());
        for(auto& tids : tid_groups)
        {
            SimplexSet<EdgeKey> eids = is_mesh.get_edges(tids) - m_eids;
            SimplexSet<NodeKey> polygon = get_polygon(eids);
            check_consistency(is_mesh.get(eid).node_keys(), polygon);
            polygons.push_back(polygon);
        }

        struct {
            bool operator()(const SimplexSet<NodeKey>& a, const SimplexSet<NodeKey>& b)
            {
                return a.size() > b.size();
            }
        } compare;
        sort(polygons.begin(), polygons.end(), compare);

        return polygons;
    }

    void DeformableSimplicialComplex::flip_23_recursively(const SimplexSet<NodeKey>& polygon, const NodeKey& n1, const NodeKey& n2, vector<vector<int>>& K, int i, int j) {
        if(j >= i+2)
        {
            int k = K[i][j];
            flip_23_recursively(polygon, n1, n2, K, i, k);
            flip_23_recursively(polygon, n1, n2, K, k, j);
            is_mesh.flip_23(is_mesh.get_face(n1, n2, polygon[k]));
        }
    }

    void DeformableSimplicialComplex::topological_edge_removal(const SimplexSet<NodeKey>& polygon, const NodeKey& n1, const NodeKey& n2, vector<vector<int>>& K) {
        const int m = static_cast<int>(polygon.size());
        int k = K[0][m-1];
        flip_23_recursively(polygon, n1, n2, K, 0, k);
        flip_23_recursively(polygon, n1, n2, K, k, m-1);
        is_mesh.flip_32(is_mesh.get_edge(n1, n2));
    }

    bool DeformableSimplicialComplex::topological_edge_removal(const EdgeKey& eid) {
        vector<SimplexSet<NodeKey>> polygon = get_polygons(eid);
#ifdef DEBUG
        assert(polygon.size() == 1 && polygon.front().size() > 2);
#endif

        vector<vector<int>> K;
        double q_new = build_table(eid, polygon.front(), K);

        if (q_new > min_quality(is_mesh.get_tets(eid)))
        {
            const SimplexSet<NodeKey>& nodes = is_mesh.get(eid).node_keys();
            topological_edge_removal(polygon.front(), nodes[0], nodes[1], K);
            return true;
        }
        return false;
    }

    void DeformableSimplicialComplex::topological_boundary_edge_removal(const SimplexSet<NodeKey>& polygon1, const SimplexSet<NodeKey>& polygon2, const EdgeKey& eid, vector<vector<int>>& K1, vector<vector<int>>& K2) {
        auto nids = is_mesh.get(eid).node_keys();
        const int m1 = static_cast<int>(polygon1.size());
        const int m2 = static_cast<int>(polygon2.size());
        int k = K1[0][m1-1];
        flip_23_recursively(polygon1, nids[0], nids[1], K1, 0, k);
        flip_23_recursively(polygon1, nids[0], nids[1], K1, k, m1-1);

        if(m2 <= 2) {
            // Find the faces to flip about.
            FaceKey f1 = is_mesh.get_face(nids[0], nids[1], polygon1.front());
            FaceKey f2 = is_mesh.get_face(nids[0], nids[1], polygon1.back());
#ifdef DEBUG
            assert(is_mesh.get(f1).is_boundary() && is_mesh.get(f2).is_boundary());
#endif

            if(precond_flip_edge(is_mesh.get_edge(f1, f2), f1, f2))
            {
                is_mesh.flip_22(f1, f2);
            }
        }
        else {
            k = K2[0][m2-1];
            flip_23_recursively(polygon2, nids[0], nids[1], K2, 0, k);
            flip_23_recursively(polygon2, nids[0], nids[1], K2, k, m2-1);

            // Find the faces to flip about.
            FaceKey f1 = is_mesh.get_face(nids[0], nids[1], polygon1.front());
            FaceKey f2 = is_mesh.get_face(nids[0], nids[1], polygon1.back());

            if(precond_flip_edge(is_mesh.get_edge(f1, f2), f1, f2))
            {
                is_mesh.flip_44(f1, f2);
            }
        }
    }

    bool DeformableSimplicialComplex::topological_boundary_edge_removal(const EdgeKey& eid) {
        vector<SimplexSet<NodeKey>> polygons = get_polygons(eid);

        if(polygons.size() > 2 || polygons[0].size() <= 2)
        {
            return false;
        }

        vector<vector<int>> K1, K2;
        double q_new = build_table(eid, polygons[0], K1);

        if(polygons.size() == 2 && polygons[1].size() > 2)
        {
            q_new = Util::min(q_new, build_table(eid, polygons[1], K2));
        }
        else
        {
            polygons.push_back({polygons[0].front(), polygons[0].back()});
        }

        if (q_new > min_quality(is_mesh.get_tets(eid)))
        {
            topological_boundary_edge_removal(polygons[0], polygons[1], eid, K1, K2);
            return true;
        }
        return false;
    }

    void DeformableSimplicialComplex::topological_edge_removal() {
        vector<TetrahedronKey> tets;
        for (auto & tit : tetrahedra())
        {
            if (tit.quality() < pars.MIN_TET_QUALITY)
            {
                tets.push_back(tit.key());
            }
        }

        // Attempt to remove each edge of each tetrahedron in tets. Accept if it increases the minimum quality locally.
        int i = 0, j = 0, k = 0;
        for (auto &t : tets)
        {
            if (is_unsafe_editable(t) && is_mesh.get(t).quality() < pars.MIN_TET_QUALITY)
            {
                for (auto e : is_mesh.get_edges(t))
                {
                    if(is_safe_editable(e))
                    {
                        if(topological_edge_removal(e))
                        {
                            i++;
                            break;
                        }
                    }
                    else if(is_mesh.exists(e) && (is_mesh.get(e).is_interface() || is_mesh.get(e).is_boundary()) && is_flippable(e))
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
        cout << "Topological edge removals: " << i + k << "/" << j << " (" << k << " at interface)" << endl;
#endif
        is_mesh.garbage_collect();
    }

    SimplexSet<EdgeKey> DeformableSimplicialComplex::test_neighbour(const FaceKey& f, const NodeKey& a, const NodeKey& b, const NodeKey& u, const NodeKey& w, double& q_old, double& q_new) {
        EdgeKey e = is_mesh.get_edge(u,w);
        SimplexSet<FaceKey> g_set = is_mesh.get(e).face_keys() - is_mesh.get_faces(is_mesh.get(f).tet_keys());
        double q = Util::quality(is_mesh.get(a).get_pos(), is_mesh.get(b).get_pos(), is_mesh.get(w).get_pos(), is_mesh.get(u).get_pos());

        if(g_set.size() == 1 && is_safe_editable(e))
        {
            FaceKey g = g_set.front();
            NodeKey v = (is_mesh.get_nodes(g) - is_mesh.get(e).node_keys()).front();
            double V_uv = Util::signed_volume(is_mesh.get(a).get_pos(), is_mesh.get(b).get_pos(), is_mesh.get(v).get_pos(), is_mesh.get(u).get_pos());
            double V_vw = Util::signed_volume(is_mesh.get(a).get_pos(), is_mesh.get(b).get_pos(), is_mesh.get(w).get_pos(), is_mesh.get(v).get_pos());
            double V_wu = Util::signed_volume(is_mesh.get(a).get_pos(), is_mesh.get(b).get_pos(), is_mesh.get(u).get_pos(), is_mesh.get(w).get_pos());

            if((V_uv >= EPSILON && V_vw >= EPSILON) || (V_vw >= EPSILON && V_wu >= EPSILON) || (V_wu >= EPSILON && V_uv >= EPSILON))
            {
                q_old = Util::min(Util::quality(is_mesh.get(a).get_pos(), is_mesh.get(u).get_pos(), is_mesh.get(w).get_pos(), is_mesh.get(v).get_pos()),
                        Util::quality(is_mesh.get(u).get_pos(), is_mesh.get(v).get_pos(), is_mesh.get(b).get_pos(), is_mesh.get(w).get_pos()));

                double q_uv_old, q_uv_new, q_vw_old, q_vw_new;
                SimplexSet<EdgeKey> uv_edges = test_neighbour(g, a, b, u, v, q_uv_old, q_uv_new);
                SimplexSet<EdgeKey> vw_edges = test_neighbour(g, a, b, v, w, q_vw_old, q_vw_new);

                q_old = Util::min(Util::min(q_old, q_uv_old), q_vw_old);
                q_new = Util::min(q_uv_new, q_vw_new);

                if(q_new > q_old || q_new > q)
                {
                    SimplexSet<EdgeKey> edges = {is_mesh.get_edge(f, g)};
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

    bool DeformableSimplicialComplex::topological_face_removal(const FaceKey& f) {
        SimplexSet<NodeKey> nids = is_mesh.get_nodes(f);
        SimplexSet<NodeKey> apices = is_mesh.get_nodes(is_mesh.get(f).tet_keys()) - nids;
        is_mesh.orient_cc(apices[0], nids);

        double q_01_new, q_01_old, q_12_new, q_12_old, q_20_new, q_20_old;
        SimplexSet<EdgeKey> e01 = test_neighbour(f, apices[0], apices[1], nids[0], nids[1], q_01_old, q_01_new);
        SimplexSet<EdgeKey> e12 = test_neighbour(f, apices[0], apices[1], nids[1], nids[2], q_12_old, q_12_new);
        SimplexSet<EdgeKey> e20 = test_neighbour(f, apices[0], apices[1], nids[2], nids[0], q_20_old, q_20_new);

        double q_old = Util::min(Util::min(Util::min(min_quality(is_mesh.get(f).tet_keys()), q_01_old), q_12_old), q_20_old);
        double q_new = Util::min(Util::min(q_01_new, q_12_new), q_20_new);

        if(q_new > q_old)
        {
            is_mesh.flip_23(f);
            for(auto &e : e01)
            {
                is_mesh.flip_32(e);
            }
            for(auto &e : e12)
            {
                is_mesh.flip_32(e);
            }
            for(auto &e : e20)
            {
                is_mesh.flip_32(e);
            }
            return true;
        }
        return false;
    }

    bool DeformableSimplicialComplex::topological_face_removal(const NodeKey& apex1, const NodeKey& apex2) {
        SimplexSet<FaceKey> fids = is_mesh.get_faces(is_mesh.get_tets(apex1)) & is_mesh.get_faces(is_mesh.get_tets(apex2));
        vec3 p = is_mesh.get(apex1).get_pos();
        vec3 ray = is_mesh.get(apex2).get_pos() - p;
        for(auto f : fids)
        {
            if(is_safe_editable(f))
            {
                auto nids = is_mesh.get_nodes(f);
                is_mesh.orient_cc(apex2, nids);

                double t = Util::intersection_ray_triangle(p, ray, is_mesh.get(nids[0]).get_pos(), is_mesh.get(nids[1]).get_pos(), is_mesh.get(nids[2]).get_pos());
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
        vector<TetrahedronKey> tets;
        for (auto & tit : tetrahedra())
        {
            if (tit.quality() < pars.MIN_TET_QUALITY)
            {
                tets.push_back(tit.key());
            }
        }

        // Attempt to remove each face of each remaining tetrahedron in tets using multi-face removal.
        // Accept if it increases the minimum quality locally.
        int i = 0, j = 0;
        for (auto &t : tets)
        {
            if (is_unsafe_editable(t) && is_mesh.get(t).quality() < pars.MIN_TET_QUALITY)
            {
                for (auto f : is_mesh.get(t).face_keys())
                {
                    if (is_safe_editable(f))
                    {
                        auto apices = is_mesh.get_nodes(is_mesh.get(f).tet_keys()) - is_mesh.get_nodes(f);
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
        cout << "Topological face removals: " << i << "/" << j << endl;
#endif

        is_mesh.garbage_collect();
    }

    void DeformableSimplicialComplex::thickening_interface() {
        if(pars.MAX_LENGTH == INFINITY)
        {
            return;
        }

        vector<EdgeKey> edges;
        for (auto & eit : DeformableSimplicialComplex::edges())
        {
            if ((eit.is_interface() || eit.is_boundary()) && eit.length() > pars.MAX_LENGTH*AVG_LENGTH)
            {
                edges.push_back(eit.key());
            }
        }
        int i = 0;
        for(auto &e : edges)
        {
            if (is_mesh.exists(e) && is_mesh.get(e).length() > pars.MAX_LENGTH*AVG_LENGTH && !is_flat(is_mesh.get(e).face_keys()))
            {
                split(e);
                i++;
            }
        }
#ifdef DEBUG
        cout << "Thickening interface splits: " << i << endl;
#endif
    }

    void DeformableSimplicialComplex::thickening() {
        if(pars.MAX_VOLUME == INFINITY)
        {
            return;
        }

        vector<TetrahedronKey> tetrahedra;
        for (auto & tit : DeformableSimplicialComplex::tetrahedra())
        {
            if (tit.volume() > pars.MAX_VOLUME*AVG_VOLUME)
            {
                tetrahedra.push_back(tit.key());
            }
        }
        int i = 0;
        for(auto &t : tetrahedra)
        {
            if (is_unsafe_editable(t) && is_mesh.get(t).volume() > pars.MAX_VOLUME*AVG_VOLUME)
            {
                split(t);
                i++;
            }
        }
#ifdef DEBUG
        cout << "Thickening splits: " << i << endl;
#endif
    }

    void DeformableSimplicialComplex::thinning_interface() {
        if(pars.MIN_LENGTH <= 0.)
        {
            return;
        }

        vector<EdgeKey> edges;
        for (auto & eit : DeformableSimplicialComplex::edges())
        {
            if ((eit.is_interface() || eit.is_boundary()) && eit.length() < pars.MIN_LENGTH*AVG_LENGTH)
            {
                edges.push_back(eit.key());
            }
        }
        int i = 0, j = 0;
        for(auto &e : edges)
        {
            if (is_mesh.exists(e) && is_mesh.get(e).length() < pars.MIN_LENGTH*AVG_LENGTH)
            {
                if(collapse(e))
                {
                    i++;
                }
                j++;
            }
        }
#ifdef DEBUG
        cout << "Thinning interface collapses: " << i << "/" << j << endl;
#endif
    }

    void DeformableSimplicialComplex::thinning() {
        if(pars.MIN_VOLUME <= 0.)
        {
            return;
        }

        vector<TetrahedronKey> tetrahedra;
        for (auto & tit : DeformableSimplicialComplex::tetrahedra())
        {
            if (tit.volume() < pars.MIN_VOLUME*AVG_VOLUME)
            {
                tetrahedra.push_back(tit.key());
            }
        }
        int i = 0, j = 0;
        for(auto &t : tetrahedra)
        {
            if (is_unsafe_editable(t) && is_mesh.get(t).volume() < pars.MIN_VOLUME*AVG_VOLUME)
            {
                if(collapse(t))
                {
                    i++;
                }
                j++;
            }
        }
#ifdef DEBUG
        cout << "Thinning collapses: " << i << "/" << j << endl;
#endif
    }

    void DeformableSimplicialComplex::remove_degenerate_edges() {
        list<EdgeKey> edges;
        for (auto & eit : DeformableSimplicialComplex::edges())
        {
            if (quality(eit.key()) < pars.DEG_EDGE_QUALITY)
            {
                edges.push_back(eit.key());
            }
        }
        int i = 0, j = 0;
        for(auto e : edges)
        {
            if(is_mesh.exists(e) && quality(e) < pars.DEG_EDGE_QUALITY && !collapse(e))
            {
                if(collapse(e, false))
                {
                    i++;
                }
                j++;
            }
        }
#ifdef DEBUG
        cout << "Removed " << i <<"/"<< j << " degenerate edges" << endl;
#endif
        is_mesh.garbage_collect();
    }

    void DeformableSimplicialComplex::remove_degenerate_faces() {
        list<FaceKey> faces;

        for (auto & fit : DeformableSimplicialComplex::faces())
        {
            if(fit.quality() < pars.DEG_FACE_QUALITY)
            {
                faces.push_back(fit.key());
            }
        }

        int i = 0, j = 0;
        for (auto &f : faces)
        {
            if (is_mesh.exists(f) && is_mesh.get(f).quality() < pars.DEG_FACE_QUALITY && !collapse(f))
            {
                if(collapse(f, false))
                {
                    i++;
                }
                else {
                    EdgeKey e = longest_edge(is_mesh.get(f).edge_keys());
                    if(is_mesh.get(e).length() > AVG_LENGTH)
                    {
                        split(e);
                    }
                }
                j++;
            }
        }
#ifdef DEBUG
        cout << "Removed " << i <<"/"<< j << " degenerate faces" << endl;
#endif
        is_mesh.garbage_collect();
    }

    void DeformableSimplicialComplex::remove_degenerate_tets() {
        vector<TetrahedronKey> tets;

        for (auto & tit : DeformableSimplicialComplex::tetrahedra())
        {
            if (tit.quality() < pars.DEG_TET_QUALITY)
            {
                tets.push_back(tit.key());
            }
        }
        int i = 0, j = 0;
        for (auto &t : tets)
        {
            if (is_mesh.exists(t) && is_mesh.get(t).quality() < pars.DEG_TET_QUALITY && !collapse(t))
            {
                if(collapse(t, false))
                {
                    i++;
                }
                else {
                    EdgeKey e = longest_edge(is_mesh.get_edges(t));
                    if(is_mesh.get(e).length() > AVG_LENGTH)
                    {
                        split(e);
                    }
                }
                j++;
            }
        }
#ifdef DEBUG
        cout << "Removed " << i <<"/"<< j << " degenerate tets" << endl;
#endif
        is_mesh.garbage_collect();
    }

    void DeformableSimplicialComplex::remove_low_quality_edges() {
        list<EdgeKey> edges;
        for (auto & eit : DeformableSimplicialComplex::edges())
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
        cout << "Removed " << i <<"/"<< j << " low quality edges" << endl;
#endif
        is_mesh.garbage_collect();
    }

    bool DeformableSimplicialComplex::remove_cap(const FaceKey& fid) {
        // Find longest edge
        EdgeKey eid = is_mesh.get(fid).longest_edge();

        // Find apex
        NodeKey apex = (is_mesh.get_nodes(fid) - is_mesh.get(eid).node_keys()).front();
        // Find the projected position of the apex
        auto verts = is_mesh.get_pos(is_mesh.get(eid).node_keys());
        vec3 p = Util::project_point_line(is_mesh.get(apex).get_pos(), verts[1], verts[0]);

        // Split longest edge
        NodeKey n = is_mesh.split(eid, p, p);

        // Collapse new edge
        EdgeKey e_rem = is_mesh.get_edge(apex, n);
        return collapse(e_rem);
    }

    bool DeformableSimplicialComplex::remove_needle(const FaceKey& fid) {
        // Find shortest edge
        EdgeKey e = shortest_edge(is_mesh.get(fid).edge_keys());

        // Remove edge
        return collapse(e);
    }

    bool DeformableSimplicialComplex::remove_face(const FaceKey& f) {
        if(is_mesh.get(f).max_angle() > 0.9*M_PI)
        {
            return remove_cap(f);
        }
        return remove_needle(f);
    }

    void DeformableSimplicialComplex::remove_low_quality_faces() {
        list<FaceKey> faces;

        for (auto & fit : DeformableSimplicialComplex::faces())
        {
            if(fit.quality() < pars.MIN_FACE_QUALITY)
            {
                faces.push_back(fit.key());
            }
        }

        int i = 0, j = 0;
        for (auto &f : faces)
        {
            if (is_unsafe_editable(f) && is_mesh.get(f).quality() < pars.MIN_FACE_QUALITY)
            {
                if(remove_face(f))
                {
                    i++;
                }
                j++;
            }
        }
#ifdef DEBUG
        cout << "Removed " << i <<"/"<< j << " low quality faces" << endl;
#endif
        is_mesh.garbage_collect();
    }

    bool DeformableSimplicialComplex::smart_laplacian(const NodeKey& nid, double alpha) {
        SimplexSet<TetrahedronKey> tids = is_mesh.get_tets(nid);
        SimplexSet<FaceKey> fids = is_mesh.get_faces(tids) - is_mesh.get_faces(nid);

        vec3 old_pos = is_mesh.get(nid).get_pos();
        vec3 avg_pos = get_barycenter(is_mesh.get_nodes(fids));
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

        for (auto & nit : nodes())
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
        cout << "Smoothed: " << i << "/" << j << endl;
#endif
    }

    void DeformableSimplicialComplex::fix_complex() {
        smooth();

        topological_edge_removal();
        topological_face_removal();

        //            remove_tets();
        //            remove_low_quality_faces();
        //            remove_low_quality_edges();

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
        is_mesh.validity_check();
        cout << endl << "********************************" << endl;
#endif
        int missing;
        int step = 0;
        do {
#ifdef DEBUG
            cout << "\n\tMove vertices step " << step << endl;
#endif
            missing = 0;
            int movable = 0;
            for (auto & nit : nodes())
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
#ifdef DEBUG
            cout << "\tVertices missing to be moved: " << missing <<"/" << movable << endl;
#endif
            fix_complex();
#ifdef DEBUG
            is_mesh.validity_check();
#endif
            ++step;
        } while (missing > 0 && step < num_steps);

        resize_complex();

        is_mesh.garbage_collect();
        for (auto & nit : is_mesh.nodes()) // reset all nodes
        {
            nit.set_destination(nit.get_pos());
        }
#ifdef DEBUG
        is_mesh.validity_check();
#endif
    }

    bool DeformableSimplicialComplex::move_vertex(const NodeKey & n) {
        Node & node = is_mesh.get(n);
        vec3 pos = node.get_pos();
        vec3 destination = node.get_destination();
        vec3 move_direction = destination - pos;
        double l = Util::length(move_direction);
        move_direction = move_direction * (1.0/l);

        if (l < 1e-4*AVG_LENGTH) // The vertex is not moved
        {
            return true;
        }

        double max_l = l*intersection_with_link(n, destination) - 1e-4 * AVG_LENGTH;
        l = Util::max(Util::min(0.5*max_l, l), 0.);
        set_pos(n, pos + l*move_direction);

        if (Util::length(destination - node.get_pos()) < 1e-4*AVG_LENGTH)
        {
            return true;
        }
        return false;
    }

    double DeformableSimplicialComplex::intersection_with_link(const NodeKey & n, const vec3& destination) {
        vec3 pos = is_mesh.get(n).get_pos();
        vec3 ray = destination - pos;

        double min_t = INFINITY;
        auto fids = is_mesh.get_faces(is_mesh.get_tets(n)) - is_mesh.get_faces(n);
        for(auto f : fids)
        {
            auto face_pos = is_mesh.get_pos(is_mesh.get_nodes(f));
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

    bool DeformableSimplicialComplex::is_flat(const SimplexSet<FaceKey>& fids) {
        for (const FaceKey& fk1 : fids) {
            Face& f1 = is_mesh.get(fk1);
            if (f1.is_interface() || f1.is_boundary())
            {
                vec3 normal1 = get_normal(fk1);
                for (const FaceKey& f2 : fids) {
                    if (fk1 != f2 && (is_mesh.get(f2).is_interface() || is_mesh.get(f2).is_boundary()))
                    {
                        vec3 normal2 = get_normal(f2);
                        if(abs(dot(normal1, normal2)) < FLIP_EDGE_INTERFACE_FLATNESS)
                        {
                            return false;
                        }
                    }
                }
            }
        }
        return true;
    }

    bool DeformableSimplicialComplex::is_flippable(const EdgeKey & eid) {
        SimplexSet<FaceKey> fids;
        for(auto f : is_mesh.get(eid).face_keys())
        {
            auto & ff = is_mesh.get(f);
            if (ff.is_interface() || ff.is_boundary())
            {
                fids += f;
            }
        }
        if (fids.size() != 2)
        {
            return false;
        }

        SimplexSet<NodeKey> e_nids = is_mesh.get(eid).node_keys();
        SimplexSet<NodeKey> new_e_nids = (is_mesh.get_nodes(fids[0]) + is_mesh.get_nodes(fids[1])) - e_nids;
#ifdef DEBUG
        assert(new_e_nids.size() == 2);
#endif

        // Check that there does not already exist an edge.
        if(is_mesh.get_edge(new_e_nids[0], new_e_nids[1]).is_valid())
        {
            return false;
        }

        // Check that the edge is not a feature edge if it is a part of the interface or boundary.
        if(is_mesh.get(eid).is_interface() || is_mesh.get(eid).is_boundary())
        {
            return is_flat(fids);
        }

        return true;
    }

    bool DeformableSimplicialComplex::precond_flip_edge(const EdgeKey& eid, const FaceKey& f1, const FaceKey& f2) {
        SimplexSet<NodeKey> e_nids = is_mesh.get(eid).node_keys();
        SimplexSet<NodeKey> new_e_nids = (is_mesh.get_nodes(f1) + is_mesh.get_nodes(f2)) - e_nids;
        SimplexSet<NodeKey> apices = (is_mesh.get_nodes(is_mesh.get(eid).face_keys()) - e_nids) - new_e_nids;
#ifdef DEBUG
        assert(e_nids.size() == 2);
        assert(new_e_nids.size() == 2);
#endif

        // Check that there does not already exist an edge.
        if(is_mesh.get_edge(new_e_nids[0], new_e_nids[1]).is_valid())
        {
            return false;
        }

        vec3 p = is_mesh.get(new_e_nids[0]).get_pos();
        vec3 r = is_mesh.get(new_e_nids[1]).get_pos() - p;
        vec3 a = is_mesh.get(e_nids[0]).get_pos();
        vec3 b = is_mesh.get(e_nids[1]).get_pos();

        for (NodeKey n : apices) {
            vec3 c = is_mesh.get(n).get_pos();
            double t = Util::intersection_ray_plane(p, r, a, b, c);
            if(t > 0. && t < 1.)
            {
                vector<double> coords = Util::barycentric_coords(p + t*r, c, a, b);
                if(coords[0] > EPSILON && coords[1] > EPSILON && coords[2] >= 0.)
                {
                    return true;
                }
            }
        }

        return false;
    }

    NodeKey DeformableSimplicialComplex::split(const TetrahedronKey& tid) {
        SimplexSet<EdgeKey> eids = is_mesh.get_edges(tid);
        EdgeKey eid = longest_edge(eids);
        return split(eid);
    }

    NodeKey DeformableSimplicialComplex::split(const FaceKey& fid) {
        SimplexSet<EdgeKey> eids = is_mesh.get(fid).edge_keys();
        EdgeKey eid = longest_edge(eids);
        return split(eid);
    }

    NodeKey DeformableSimplicialComplex::split(const EdgeKey& eid) {
        auto nids = is_mesh.get(eid).node_keys();
        auto verts = is_mesh.get_pos(nids);
        vec3 pos = Util::barycenter(verts[0], verts[1]);
        vec3 destination = pos;
        if(is_mesh.get(eid).is_interface())
        {
            destination = Util::barycenter(is_mesh.get(nids[0]).get_destination(), is_mesh.get(nids[1]).get_destination());
        }

        return is_mesh.split(eid, pos, destination);
    }

    bool DeformableSimplicialComplex::is_collapsable(const EdgeKey& eid, const NodeKey& nid, bool safe) {
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
        if(is_mesh.get(eid).is_boundary() || is_mesh.get(eid).is_interface())
        {
            return is_flat(is_mesh.get_faces(nid));
        }
        return false;
    }

    bool DeformableSimplicialComplex::collapse(const EdgeKey& eid, bool safe) {
        SimplexSet<NodeKey> nids = is_mesh.get(eid).node_keys();
        bool n0_is_editable = is_collapsable(eid, nids[0], safe);
        bool n1_is_editable = is_collapsable(eid, nids[1], safe);

        if (!n0_is_editable && !n1_is_editable)
        {
            return false;
        }
        vector<double> test_weights;
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

        SimplexSet<TetrahedronKey> e_tids = is_mesh.get_tets(eid);
        SimplexSet<FaceKey> fids0 = is_mesh.get_faces(is_mesh.get_tets(nids[0]) - e_tids) - is_mesh.get_faces(nids[0]);
        SimplexSet<FaceKey> fids1 = is_mesh.get_faces(is_mesh.get_tets(nids[1]) - e_tids) - is_mesh.get_faces(nids[1]);

        double q_max = -INFINITY;
        double weight;
        for (double w : test_weights)
        {
            vec3 p = (1.-w) * is_mesh.get(nids[1]).get_pos() + w * is_mesh.get(nids[0]).get_pos();
            double q = Util::min(min_quality(fids0, is_mesh.get(nids[0]).get_pos(), p), min_quality(fids1, is_mesh.get(nids[1]).get_pos(), p));

            if (q > q_max && ((!is_mesh.get(nids[0]).is_interface() && !is_mesh.get(nids[1]).is_interface()) || design_domain.is_inside(p)))
            {
                q_max = q;
                weight = w;
            }
        }

        if(q_max > EPSILON)
        {
            if(!safe || q_max > Util::min(min_quality(is_mesh.get_tets(nids[0]) + is_mesh.get_tets(nids[1])), pars.MIN_TET_QUALITY) + EPSILON)
            {
                is_mesh.collapse(eid, nids[1], weight);
                return true;
            }
        }
        return false;
    }

    bool DeformableSimplicialComplex::collapse(SimplexSet<EdgeKey>& eids, bool safe) {
        while(eids.size() > 0)
        {
            EdgeKey e = shortest_edge(eids);
            if(collapse(e, safe))
            {
                return true;
            }
            eids -= e;
        }
        return false;
    }

    bool DeformableSimplicialComplex::collapse(const FaceKey& fid, bool safe) {
        SimplexSet<EdgeKey> eids = is_mesh.get(fid).edge_keys();
        return collapse(eids, safe);
    }

    bool DeformableSimplicialComplex::collapse(const TetrahedronKey& tid, bool safe) {
        SimplexSet<EdgeKey> eids = is_mesh.get_edges(tid);
        return collapse(eids, safe);
    }

    vector<vec3> DeformableSimplicialComplex::get_interface_face_positions() {
        vector<vec3> verts;
        for (auto & fit : faces()) {
            if(fit.is_interface())
            {
                for(auto n : is_mesh.get_nodes(fit.key()))
                {
                    verts.push_back(is_mesh.get(n).get_pos());
                }
            }
        }
        return verts;
    }

    vec3 DeformableSimplicialComplex::get_normal(const FaceKey& fid) {
        auto pos = is_mesh.get_pos(is_mesh.get_sorted_nodes(fid));
        return Util::normal_direction(pos[0], pos[1], pos[2]);
    }

    vec3 DeformableSimplicialComplex::get_normal(const NodeKey& nid) {
        vec3 result(0.);
        for (auto f : is_mesh.get_faces(nid))
        {
            if (is_mesh.get(f).is_interface())
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

    vec3 DeformableSimplicialComplex::get_barycenter(const SimplexSet<NodeKey>& nids, bool interface) {
        vec3 avg_pos(0.);
        int i = 0;
        for (auto n : nids)
        {
            if (!interface || is_mesh.get(n).is_interface())
            {
                avg_pos += is_mesh.get(n).get_pos();
                i++;
            }
        }
#ifdef DEBUG
        assert(i != 0);
#endif
        return avg_pos / static_cast<double>(i);
    }

    vec3 DeformableSimplicialComplex::get_barycenter(const NodeKey& nid, bool interface) {
        if(interface && !is_mesh.get(nid).is_interface())
        {
            return is_mesh.get(nid).get_pos();
        }

        SimplexSet<NodeKey> nids = is_mesh.get_nodes(is_mesh.get_tets(nid)) - nid;
        return get_barycenter(nids, interface);
    }

    double DeformableSimplicialComplex::length(const EdgeKey& eid) {
        return is_mesh.get(eid).length();
    }

    double DeformableSimplicialComplex::length_destination(const EdgeKey& eid) {
        return is_mesh.get(eid).length_destination();
    }

    double DeformableSimplicialComplex::area(const FaceKey& fid) {
        return is_mesh.get(fid).area();
    }

    double DeformableSimplicialComplex::area_destination(const FaceKey& fid) {
        return is_mesh.get(fid).area_destination();
    }

    double DeformableSimplicialComplex::volume(const TetrahedronKey& tid) {
        return is_mesh.get(tid).volume();
    }

    double DeformableSimplicialComplex::volume_destination(const TetrahedronKey& tid) {
        return is_mesh.get(tid).volume_destination();
    }

    double DeformableSimplicialComplex::volume_destination(const SimplexSet<NodeKey>& nids) {
        return is_mesh.volume_destination(nids);
    }

    double DeformableSimplicialComplex::signed_volume_destination(const SimplexSet<NodeKey>& nids) {
        return is_mesh.signed_volume_destination(nids);
    }

    vec3 DeformableSimplicialComplex::barycenter(const TetrahedronKey& tid) {
        return is_mesh.get(tid).barycenter();
    }

    vec3 DeformableSimplicialComplex::barycenter_destination(const TetrahedronKey& tid) {
        return is_mesh.get(tid).barycenter_destination();
    }

    double DeformableSimplicialComplex::quality(const TetrahedronKey& tid) {
        return is_mesh.get(tid).quality();
    }

    double DeformableSimplicialComplex::min_angle(const FaceKey& fid) {
        return is_mesh.get(fid).min_angle();
    }

    double DeformableSimplicialComplex::max_angle(const FaceKey& fid) {
        return is_mesh.get(fid).max_angle();
    }

    double DeformableSimplicialComplex::quality(const FaceKey& fid) {
        return is_mesh.get(fid).quality();
    }

    double DeformableSimplicialComplex::quality(const EdgeKey& eid) {
        return is_mesh.get(eid).length()/AVG_LENGTH;
    }

    FaceKey DeformableSimplicialComplex::largest_face(const SimplexSet<FaceKey>& fids) {
        double max_a = -INFINITY;
        FaceKey max_f;
        for(auto f : fids)
        {
            double a = is_mesh.get(f).area();
            if(a > max_a)
            {
                max_a = a;
                max_f = f;
            }
        }
        return max_f;
    }

    EdgeKey DeformableSimplicialComplex::shortest_edge(const SimplexSet<EdgeKey>& eids) {
        double min_l = INFINITY;
        EdgeKey min_e;
        for(auto e : eids)
        {
            double l = is_mesh.get(e).length();
            if(l < min_l)
            {
                min_l = l;
                min_e = e;
            }
        }
        return min_e;
    }

    EdgeKey DeformableSimplicialComplex::longest_edge(const SimplexSet<EdgeKey>& eids) {
        double max_l = -INFINITY;
        EdgeKey max_e;
        for(auto e : eids)
        {
            double l = is_mesh.get(e).length();
            if(l > max_l)
            {
                max_l = l;
                max_e = e;
            }
        }
        return max_e;
    }

    double DeformableSimplicialComplex::min_quality(const SimplexSet<TetrahedronKey>& tids) {
        double q_min = INFINITY;
        for (auto t : tids)
        {
            q_min = Util::min(is_mesh.get(t).quality(), q_min);
        }
        return q_min;
    }

    double DeformableSimplicialComplex::min_quality(const SimplexSet<FaceKey>& fids, const vec3& pos) {
        double min_q = INFINITY;
        for (auto f : fids)
        {
            SimplexSet<NodeKey> nids = is_mesh.get_nodes(f);
            min_q = Util::min(min_q, abs(Util::quality(is_mesh.get(nids[0]).get_pos(), is_mesh.get(nids[1]).get_pos(), is_mesh.get(nids[2]).get_pos(), pos)));
        }
        return min_q;
    }

    double DeformableSimplicialComplex::min_quality(const SimplexSet<FaceKey>& fids, const vec3& pos_old, const vec3& pos_new) {
        double min_q = INFINITY;
        for (auto f : fids)
        {
            SimplexSet<NodeKey> nids = is_mesh.get_nodes(f);
            vector<vec3> pos = is_mesh.get_pos(nids);
            if(Util::sign(Util::signed_volume(pos[0], pos[1], pos[2], pos_old)) !=
                    Util::sign(Util::signed_volume(pos[0], pos[1], pos[2], pos_new)))
            {
                return -INFINITY;
            }
            min_q = Util::min(min_q, abs(Util::quality(pos[0], pos[1], pos[2], pos_new)));
        }
        return min_q;
    }

    void DeformableSimplicialComplex::min_quality(const SimplexSet<FaceKey>& fids, const vec3& pos_old, const vec3& pos_new, double& min_q_old, double& min_q_new) {
        min_q_old = INFINITY;
        min_q_new = INFINITY;
        for (auto f : fids)
        {
            SimplexSet<NodeKey> nids = is_mesh.get_nodes(f);
            vector<vec3> pos = is_mesh.get_pos(nids);
            if(Util::sign(Util::signed_volume(pos[0], pos[1], pos[2], pos_old)) !=
                    Util::sign(Util::signed_volume(pos[0], pos[1], pos[2], pos_new)))
            {
                min_q_old = INFINITY;
                min_q_new = -INFINITY;
                break;
            }
            min_q_old = Util::min(min_q_old, abs(Util::quality(pos[0], pos[1], pos[2], pos_old)));
            min_q_new = Util::min(min_q_new, abs(Util::quality(pos[0], pos[1], pos[2], pos_new)));
        }
    }

    void DeformableSimplicialComplex::check_consistency(const SimplexSet<NodeKey>& nids, SimplexSet<NodeKey>& polygon) {
        unsigned int n = static_cast<unsigned int>(polygon.size());
        vector<vec3> pos = is_mesh.get_pos(nids);
        vector<vec3> poly_pos = is_mesh.get_pos(polygon);
        double sum = 0;
        for (unsigned int i = 0; i < n; ++i)
        {
            sum += Util::signed_volume(pos[0], pos[1], poly_pos[i], poly_pos[(i+1)%n]);
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
        for (auto & eit : edges()) {
            if(eit.is_interface())
            {
                avg_edge_length += eit.length();
                N++;
            }
        }
        if (N > 0) {
            avg_edge_length /= static_cast<double>(N);
        }
        return avg_edge_length;
    }

    double DeformableSimplicialComplex::cos_dihedral_angle(const FaceKey& f1, const FaceKey& f2) {
        auto nids1 = is_mesh.get_nodes(f1);
        auto nids2 = is_mesh.get_nodes(f2);
        SimplexSet<NodeKey> nids = nids1 & nids2;
        SimplexSet<NodeKey> apices = (nids1 + nids2) - nids;

        return Util::cos_dihedral_angle(is_mesh.get(nids[0]).get_pos(), is_mesh.get(nids[1]).get_pos(), is_mesh.get(apices[0]).get_pos(), is_mesh.get(apices[1]).get_pos());
    }

    double DeformableSimplicialComplex::dihedral_angle(const FaceKey& f1, const FaceKey& f2) {
        return acos(cos_dihedral_angle(f1, f2));
    }

    vector<double> DeformableSimplicialComplex::cos_dihedral_angles(const TetrahedronKey& tid) {
        auto verts = is_mesh.get_pos(is_mesh.get_nodes(tid));
        vector<double> angles;
        vector<int> apices;
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

    double DeformableSimplicialComplex::min_cos_dihedral_angle(const TetrahedronKey& t) {
        double min_angle = -1.;
        vector<double> angles = cos_dihedral_angles(t);
        for(auto a : angles)
        {
            min_angle = Util::max(min_angle, a);
        }
        return min_angle;
    }

    double DeformableSimplicialComplex::min_dihedral_angle(const TetrahedronKey& t) {
        return acos(min_cos_dihedral_angle(t));
    }

    void DeformableSimplicialComplex::get_qualities(vector<int>& histogram, double& min_quality) {
        min_quality = INFINITY;

        histogram = vector<int>(100);
        for (int i = 0; i < 100; ++i)
        {
            histogram[i] = 0;
        }

        for (auto & tit : tetrahedra())
        {
            double q = tit.quality();
            min_quality = Util::min(min_quality, q);
            int index = static_cast<int>(floor(q*100.));
#ifdef DEBUG
            assert(index < 100 && index >= 0);
#endif
            histogram[index] += 1;
        }
    }

    void DeformableSimplicialComplex::get_dihedral_angles(vector<int> & histogram, double & min_angle, double & max_angle) {
        max_angle = -INFINITY, min_angle = INFINITY;

        histogram = vector<int>(180);
        for (int i = 0; i < 180; ++i)
        {
            histogram[i] = 0;
        }

        for (auto & tit : tetrahedra())
        {
            vector<double> angles = cos_dihedral_angles(tit.key());
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
        for (auto & tit : tetrahedra())
        {
            min_q = Util::min(min_q, tit.quality());
        }
        return min_q;
    }

    void DeformableSimplicialComplex::count_nodes(int & total, int & object) {
        total = 0, object = 0;
        for (auto & nit : nodes())
        {
            total++;
            if (nit.is_interface())
            {
                object++;
            }
        }
    }

    void DeformableSimplicialComplex::count_edges(int & total, int & object) {
        total = 0, object = 0;
        for (auto & e : edges()){
            total++;
            if (e.is_interface()){
                object++;
            }
        }
    }

    void DeformableSimplicialComplex::count_faces(int & total, int & object) {
        total = 0, object = 0;
        for (auto & f : faces()){
            total++;
            if (f.is_interface())
            {
                object++;
            }
        }
    }

    void DeformableSimplicialComplex::count_tetrahedra(int & total, int & object) {
        total = 0, object = 0;
        for (auto & t : tetrahedra()){
            total++;
            if (t.label() != 0)
            {
                object++;
            }
        }
    }

    void DeformableSimplicialComplex::test_split_collapse() {
        SimplexSet<EdgeKey> eids;
        for (auto & eit : edges())
        {
            auto neighbours = is_mesh.get_edges(eit.face_keys());
            bool ok = true;
            for(auto e : neighbours)
            {
                if (eids.contains(e))
                {
                    ok = false;
                }
            }
            if (ok)
            {
                eids += eit.key();
            }
        }

        cout << " DONE" << endl;
        is_mesh.garbage_collect();
        is_mesh.validity_check();
    }

    void DeformableSimplicialComplex::test_flip23_flip32() {
        SimplexSet<FaceKey> fids;
        for (auto & fit : faces())
        {
            if (is_safe_editable(fit.key()))
            {
                auto nids = is_mesh.get_nodes(fit.key());
                nids += is_mesh.get_nodes(fit.tet_keys());
                double t = Util::intersection_ray_triangle(is_mesh.get(nids[3]).get_pos(), is_mesh.get(nids[4]).get_pos() - is_mesh.get(nids[3]).get_pos(), is_mesh.get(nids[0]).get_pos(), is_mesh.get(nids[1]).get_pos(), is_mesh.get(nids[2]).get_pos());

                auto neighbours = is_mesh.get_faces(fit.tet_keys());
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

        cout << "Flip 2-3 test # = " << fids.size();
        SimplexSet<EdgeKey> new_eids;
        int i = 0;
        for (auto f : fids) {
            assert(is_mesh.exists(f));
            auto new_eid = is_mesh.flip_23(f);
            assert(new_eid.is_valid());
            new_eids += new_eid;
            i++;
            if(i%1000 == 0)
            {
                cout << ".";
            }
        }
        cout << " DONE" << endl;
        is_mesh.garbage_collect();
        is_mesh.validity_check();

        i=0;
        cout << "Flip 3-2 test # = " << new_eids.size();
        for (auto e : new_eids) {
            auto new_fid = is_mesh.flip_32(e);
            assert(new_fid.is_valid());
            i++;
            if(i%1000 == 0)
            {
                cout << ".";
            }
        }
        cout << " DONE" << endl;
        is_mesh.garbage_collect();
        is_mesh.validity_check();
    }

    void DeformableSimplicialComplex::test_flip44() {
        SimplexSet<EdgeKey> eids;
        for (auto & eit : edges())
        {

            if(is_unsafe_editable(eit.key()) && eit.is_interface() && eit.face_keys().size() == 4)
            {
                auto neighbours = is_mesh.get_edges(is_mesh.get_tets(eit.key()));
                bool ok = true;
                for(auto e : neighbours)
                {
                    if(eids.contains(e))
                    {
                        ok = false;
                    }
                }
                SimplexSet<FaceKey> flip_fids;
                for(auto f : eit.face_keys())
                {
                    if(is_mesh.get(f).is_interface())
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
            cout << "Flip 4-4 test # = " << eids.size();
            int i = 0;
            for (auto e : eids) {
                SimplexSet<FaceKey> flip_fids;
                for(auto f : is_mesh.get(e).face_keys())
                {
                    if(is_mesh.get(f).is_interface())
                    {
                        flip_fids += f;
                    }
                }
                assert(flip_fids.size() == 2);
                assert(is_mesh.get(e).face_keys().size() == 4);
                is_mesh.flip_44(flip_fids[0], flip_fids[1]);
                i++;
                if(i%100 == 0)
                {
                    cout << ".";
                }
            }
            cout << " DONE" << endl;
            is_mesh.garbage_collect();
            is_mesh.validity_check();
        }
    }

    void DeformableSimplicialComplex::test_flip22() {
        SimplexSet<EdgeKey> eids;
        for (auto & eit : edges())
        {
            if(eit.is_boundary() && eit.face_keys().size() == 3)
            {
                auto neighbours = is_mesh.get_edges(is_mesh.get_tets(eit.key()));
                bool ok = true;
                for(auto e : neighbours)
                {
                    if(eids.contains(e))
                    {
                        ok = false;
                    }
                }
                SimplexSet<FaceKey> flip_fids;
                for(auto f : eit.face_keys())
                {
                    if(is_mesh.get(f).is_boundary())
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
            cout << "Flip 2-2 test # = " << eids.size();
            int i = 0;
            for (auto e : eids) {
                assert(is_mesh.exists(e));
                auto fids = is_mesh.get(e).face_keys();
                assert(fids.size() == 3);
                SimplexSet<FaceKey> flip_fids;
                for(auto f : fids)
                {
                    if(is_mesh.get(f).is_boundary())
                    {
                        flip_fids += f;
                    }
                }

                assert(flip_fids.size() == 2);
                is_mesh.flip_22(flip_fids[0], flip_fids[1]);
                i++;
                if(i%10 == 0)
                {
                    cout << ".";
                }
            }
            cout << " DONE" << endl;
            is_mesh.garbage_collect();
            is_mesh.validity_check();
        }
    }

    is_mesh::ISMesh &DeformableSimplicialComplex::get_is_mesh() {
        return is_mesh;
    }

    void DeformableSimplicialComplex::set_subdomain(std::shared_ptr<Subdomain> subdomain) {
        this->subdomain = subdomain;
    }

    void DeformableSimplicialComplex::clear_subdomain() {
        set_subdomain(nullptr);
    }

    void DeformableSimplicialComplex::on_gc(const is_mesh::GarbageCollectDeletions& gc){
        if (subdomain){
            subdomain->rebuild();
        }
    }

    std::shared_ptr<Subdomain> DeformableSimplicialComplex::get_subdomain() {
        return subdomain;
    }

    is_mesh::NodeIterator DeformableSimplicialComplex::nodes() const {
        return is_mesh.nodes(subdomain?&(subdomain->nodes()): nullptr);
    }

    is_mesh::EdgeIterator DeformableSimplicialComplex::edges() const {
        return is_mesh.edges(subdomain?&subdomain->edges(): nullptr);
    }

    is_mesh::FaceIterator DeformableSimplicialComplex::faces() const {
        return is_mesh.faces(subdomain?&subdomain->faces(): nullptr);
    }

    is_mesh::TetrahedronIterator DeformableSimplicialComplex::tetrahedra() const {
        return is_mesh.tetrahedra(subdomain?&subdomain->tetrahedra(): nullptr);
    }
}
