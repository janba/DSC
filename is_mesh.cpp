#include "is_mesh.h"

using namespace std;

namespace is_mesh{
    ISMesh::ISMesh(vector<vec3> & points, vector<int> & tets, const vector<int>& tet_labels) {
        m_node_kernel = new kernel<node_type, NodeKey>();
        m_edge_kernel = new kernel<edge_type, EdgeKey>();
        m_face_kernel = new kernel<face_type, FaceKey>();
        m_tetrahedron_kernel = new kernel<tetrahedron_type, TetrahedronKey>();

        create(points, tets);
        init_flags(tet_labels);
        validity_check();
    }


    ISMesh::~ISMesh() {
        delete m_tetrahedron_kernel;
        delete m_face_kernel;
        delete m_edge_kernel;
        delete m_node_kernel;
    }

    unsigned int ISMesh::get_no_nodes() const {
        return static_cast<unsigned int>(m_node_kernel->size());
    }

    unsigned int ISMesh::get_no_edges() const {
        return static_cast<unsigned int>(m_edge_kernel->size());
    }

    unsigned int ISMesh::get_no_faces() const {
        return static_cast<unsigned int>(m_face_kernel->size());
    }

    unsigned int ISMesh::get_no_tets() const {
        return static_cast<unsigned int>(m_tetrahedron_kernel->size());
    }

    NodeIterator<NodeAttributes> ISMesh::nodes() const {
        return NodeIterator<NodeAttributes>{m_node_kernel};
    }

    typename kernel<Node, NodeKey>::iterator ISMesh::nodes_begin() {
        return m_node_kernel->begin();
    }

    typename kernel<Node, NodeKey>::iterator ISMesh::nodes_end() {
        return m_node_kernel->end();
    }

    EdgeIterator<EdgeAttributes> ISMesh::edges() const {
        return EdgeIterator<EdgeAttributes>{m_edge_kernel};
    }

    typename kernel<Edge, EdgeKey>::iterator ISMesh::edges_begin() {
        return m_edge_kernel->begin();
    }

    typename kernel<Edge, EdgeKey>::iterator ISMesh::edges_end() {
        return m_edge_kernel->end();
    }

    FaceIterator<FaceAttributes> ISMesh::faces() const {
        return FaceIterator<FaceAttributes>{m_face_kernel};
    }

    typename kernel<Face, FaceKey>::iterator ISMesh::faces_begin() {
        return m_face_kernel->begin();
    }

    typename kernel<Face, FaceKey>::iterator ISMesh::faces_end() {
        return m_face_kernel->end();
    }

    TetrahedronIterator<TetAttributes> ISMesh::tetrahedra() const {
        return TetrahedronIterator<TetAttributes>{m_tetrahedron_kernel};
    }

    typename kernel<Tetrahedron, TetrahedronKey>::iterator ISMesh::tetrahedra_begin() {
        return m_tetrahedron_kernel->begin();
    }

    typename kernel<Tetrahedron, TetrahedronKey>::iterator ISMesh::tetrahedra_end() {
        return m_tetrahedron_kernel->end();
    }

    int ISMesh::get_label(const TetrahedronKey& t) {
        return get(t).label();
    }

    void ISMesh::set_label(const TetrahedronKey& tid, int label) {
        get(tid).label(label);
        SimplexSet<TetrahedronKey> tids = {tid};
        update(tids);
    }

    bool ISMesh::create(const vector<vec3>& points, const vector<int>& tets) {
        map<edge_key, int> edge_map;
        map<face_key, int> face_map;

        for (vec3 p : points)
        {
            insert_node(p);
        }

        for (unsigned int j = 0; j < tets.size()/4; ++j)
        {
            int idx[4];
            idx[0] = tets[4*j];
            idx[1] = tets[4*j+1];
            idx[2] = tets[4*j+2];
            idx[3] = tets[4*j+3];

            int edges[6];
            edges[0] = create_edge(idx[0], idx[1], edge_map, *this); //edge 01
            edges[1] = create_edge(idx[0], idx[2], edge_map, *this); //edge 02
            edges[2] = create_edge(idx[0], idx[3], edge_map, *this); //edge 03
            edges[3] = create_edge(idx[1], idx[2], edge_map, *this); //edge 12
            edges[4] = create_edge(idx[1], idx[3], edge_map, *this); //edge 13
            edges[5] = create_edge(idx[2], idx[3], edge_map, *this); //edge 23

            int faces[4];
            faces[0] = create_face(edges[3], edges[5], edges[4], face_map, *this); //12-23-31
            faces[1] = create_face(edges[1], edges[5], edges[2], face_map, *this); //02-23-30
            faces[2] = create_face(edges[0], edges[4], edges[2], face_map, *this); //01-13-30
            faces[3] = create_face(edges[0], edges[3], edges[1], face_map, *this); //01-12-20

            insert_tetrahedron( faces[0], faces[1], faces[2], faces[3] );
        }

        return true;
    }

    void ISMesh::init_flags(const vector<int>& tet_labels) {
        if(tet_labels.size() > 0)
        {
            for (auto tit = tetrahedra_begin(); tit != tetrahedra_end(); tit++)
            {
                tit->label(tet_labels[tit.key()]);
            }
        }

        for (auto fit = faces_begin(); fit != faces_end(); fit++)
        {
            update_flag(fit.key());
        }

        for (auto eit = edges_begin(); eit != edges_end(); eit++)
        {
            update_flag(eit.key());
        }

        for (auto nit = nodes_begin(); nit != nodes_end(); nit++)
        {
            update_flag(nit.key());
        }
    }

    void ISMesh::update(const SimplexSet<TetrahedronKey>& tids) {
        // Update faces
        SimplexSet<FaceKey> fids = get_faces(tids);
        for (auto f : fids)
        {
            if (exists(f))
            {
                update_flag(f);
            }
        }

        // Update edges
        SimplexSet<EdgeKey> eids = get_edges(fids);
        for (auto e : eids)
        {
            if (exists(e))
            {
                update_flag(e);
            }
        }

        // Update nodes
        SimplexSet<NodeKey> nids = get_nodes(eids);
        for (auto n : nids)
        {
            if (exists(n))
            {
                update_flag(n);
            }
        }
    }

    void ISMesh::update_flag(const FaceKey & f) {
        set_interface(f, false);
        set_boundary(f, false);

        SimplexSet<TetrahedronKey> tids = get_tets(f);
        if (tids.size() == 1)
        {
            // On the boundary
            set_boundary(f, true);
            if (get_label(tids.front()) != 0)
            {
                set_interface(f, true);
            }
        }
        else if(tids.size() == 2)
        {
            if (get_label(tids.front()) != get_label(tids.back()))
            {
                // On the interface
                set_interface(f, true);
            }
        }
    }

    void ISMesh::update_flag(const EdgeKey & e) {
        set_boundary(e, false);
        set_interface(e, false);
        set_crossing(e, false);

        int i = 0;
        for (auto f : get_faces(e))
        {
            if (exists(f))
            {
                if (get(f).is_boundary())
                {
                    set_boundary(e, true);
                }
                if (get(f).is_interface())
                {
                    set_interface(e, true);
                    i++;
                }
            }
        }
        if(i > 2)
        {
            set_crossing(e, true);
        }
    }

    void ISMesh::connected_component(SimplexSet<TetrahedronKey>& tids, const TetrahedronKey& tid) {
                int label = get_label(tid);
                tids -= tid;

                for(auto f : get_faces(tid))
                {
                    if(get_tets(f).size() == 2)
                    {
                        TetrahedronKey tid2 = get_tet(tid, f);
                        if(tids.contains(tid2) && label == get_label(tid2))
                        {
                            connected_component(tids, tid2);
                        }
                    }
                }
            }

    bool ISMesh::crossing(const NodeKey& n) {
        SimplexSet<TetrahedronKey> tids = get_tets(n);

        int c = 0;
        while (tids.size() > 0)
        {
            if(c == 2)
            {
                return true;
            }
            TetrahedronKey tid = tids.front();
            connected_component(tids, tid);
            c++;
        }
        return false;
    }

    void ISMesh::update_flag(const NodeKey & n) {
        set_interface(n, false);
        set_boundary(n, false);
        set_crossing(n, false);

        for (auto e : get_edges(n))
        {
            if (exists(e))
            {
                if (get(e).is_interface())
                {
                    set_interface(n, true);
                }
                if (get(e).is_boundary())
                {
                    set_boundary(n, true);
                }
                if (get(e).is_crossing())
                {
                    set_crossing(n, true);
                }
            }
        }
        if(!get(n).is_crossing() && get(n).is_interface() && crossing(n))
        {
            set_crossing(n, true);
        }
    }

    is_mesh::Node &ISMesh::get(const NodeKey& nid) {
        return m_node_kernel->find(nid);
    }

    is_mesh::Edge &ISMesh::get(const EdgeKey& eid) {
        return m_edge_kernel->find(eid);
    }

    is_mesh::Face &ISMesh::get(const FaceKey& fid) {
        return m_face_kernel->find(fid);
    }

    is_mesh::Tetrahedron &ISMesh::get(const TetrahedronKey& tid) {
        return m_tetrahedron_kernel->find(tid);
    }

    const SimplexSet<NodeKey> &ISMesh::get_nodes(const EdgeKey& eid) {
        return get(eid).get_boundary();
    }

    const SimplexSet<EdgeKey> &ISMesh::get_edges(const NodeKey& nid) {
        return get(nid).get_co_boundary();
    }

    const SimplexSet<EdgeKey> &ISMesh::get_edges(const FaceKey& fid) {
        return get(fid).get_boundary();
    }

    const SimplexSet<FaceKey> &ISMesh::get_faces(const EdgeKey& eid) {
        return get(eid).get_co_boundary();
    }

    const SimplexSet<FaceKey> &ISMesh::get_faces(const TetrahedronKey& tid) {
        return get(tid).get_boundary();
    }

    const SimplexSet<TetrahedronKey> &ISMesh::get_tets(const FaceKey& fid) {
        return get(fid).get_co_boundary();
    }

    SimplexSet<NodeKey> ISMesh::get_sorted_nodes(const FaceKey& fid, const TetrahedronKey& tid) {
        SimplexSet<NodeKey> nids = get_nodes(fid);
        NodeKey apex = (get_nodes(tid) - nids).front();
        orient_cc(apex, nids);
        return nids;
    }

    SimplexSet<NodeKey> ISMesh::get_sorted_nodes(const FaceKey& fid) {
        if (get(fid).is_interface())
        {
            int label = -100;
            TetrahedronKey tid;
            for (auto t : get_tets(fid))
            {
                int tl = get_label(t);
                if (tl > label)
                {
                    label = tl;
                    tid = t;
                }
            }
            return get_sorted_nodes(fid, tid);
        }
        else if (get(fid).is_boundary())
        {
            TetrahedronKey tid = get_tets(fid).front();
            return get_sorted_nodes(fid, tid);
        }
        return get_nodes(fid);
    }

    SimplexSet<NodeKey> ISMesh::get_nodes(const FaceKey& fid) {
        const SimplexSet<EdgeKey>& eids = get_edges(fid);
        SimplexSet<NodeKey> nids = get_nodes(eids[0]);
        nids += get_nodes(eids[1]);
        return nids;
    }

    SimplexSet<NodeKey> ISMesh::get_nodes(const TetrahedronKey& tid) {
        const SimplexSet<FaceKey>& fids = get_faces(tid);
        SimplexSet<NodeKey> nids = get_nodes(fids[0]);
        nids += get_nodes(fids[1]);
        return nids;
    }

    SimplexSet<EdgeKey> ISMesh::get_edges(const TetrahedronKey& tid) {
        SimplexSet<EdgeKey> eids;
        for(const FaceKey& f : get_faces(tid))
        {
            eids += get_edges(f);
        }
        return eids;
    }

    SimplexSet<FaceKey> ISMesh::get_faces(const NodeKey& nid) {
        SimplexSet<FaceKey> fids;
        for(const EdgeKey& e : get_edges(nid))
        {
            fids += get_faces(e);
        }
        return fids;
    }

    SimplexSet<TetrahedronKey> ISMesh::get_tets(const NodeKey& nid) {
        SimplexSet<TetrahedronKey> tids;
        for(const EdgeKey& e : get_edges(nid))
        {
            for(const FaceKey& f : get_faces(e))
            {
                tids += get_tets(f);
            }
        }
        return tids;
    }

    SimplexSet<TetrahedronKey> ISMesh::get_tets(const EdgeKey& eid) {
        SimplexSet<TetrahedronKey> tids;
        for(const FaceKey& f : get_faces(eid))
        {
            tids += get_tets(f);
        }
        return tids;
    }

    NodeKey ISMesh::get_node(const EdgeKey& eid1, const EdgeKey& eid2) {
        const SimplexSet<NodeKey>& nids = get_nodes(eid2);
        for (const NodeKey& n : get_nodes(eid1)) {
            if(nids.contains(n))
            {
                return n;
            }
        }
        return NodeKey();
    }

    NodeKey ISMesh::get_node(const EdgeKey& eid, const NodeKey& nid) {
        const SimplexSet<NodeKey>& nids = get_nodes(eid);
        if(nids[0] == nid)
        {
            return nids[1];
        }
        return nids[0];
    }

    EdgeKey ISMesh::get_edge(const NodeKey& nid1, const NodeKey& nid2) {
        const SimplexSet<EdgeKey>& eids = get_edges(nid2);
        for (const EdgeKey& e : get_edges(nid1)) {
            if(eids.contains(e))
            {
                return e;
            }
        }
        return EdgeKey();
    }

    EdgeKey ISMesh::get_edge(const FaceKey& fid1, const FaceKey& fid2) {
        const SimplexSet<EdgeKey>& eids = get_edges(fid2);
        for (const EdgeKey& e : get_edges(fid1)) {
            if(eids.contains(e))
            {
                return e;
            }
        }
        return EdgeKey();
    }

    FaceKey ISMesh::get_face(const NodeKey& nid1, const NodeKey& nid2, const NodeKey& nid3) {
        SimplexSet<FaceKey> fids1 = get_faces(nid1);
        SimplexSet<FaceKey> fids2 = get_faces(nid2);
        for (const FaceKey& f : get_faces(nid3)) {
            if(fids1.contains(f) && fids2.contains(f))
            {
                return f;
            }
        }
        return FaceKey();
    }

    FaceKey ISMesh::get_face(const TetrahedronKey& tid1, const TetrahedronKey& tid2) {
        const SimplexSet<FaceKey>& fids = get_faces(tid2);
        for (const FaceKey& f : get_faces(tid1)) {
            if(fids.contains(f))
            {
                return f;
            }
        }
        return FaceKey();
    }

    TetrahedronKey ISMesh::get_tet(const TetrahedronKey& tid, const FaceKey& fid) {
        for(TetrahedronKey t : get_tets(fid))
        {
            if(t != tid)
            {
                return t;
            }
        }
        return TetrahedronKey();
    }

    vec3 ISMesh::get_pos(const NodeKey& nid) {
        return get(nid).get_pos();
    }

    vector<vec3> ISMesh::get_pos(const SimplexSet<NodeKey>& nids) {
        vector<vec3> verts(nids.size());
        for (unsigned int i = 0; i < nids.size(); i++)
        {
            verts[i] = get_pos(nids[i]);
        }
        return verts;
    }

    bool ISMesh::exists(const TetrahedronKey& t) {
        return m_tetrahedron_kernel->is_valid(t);
    }

    bool ISMesh::exists(const FaceKey& f) {
        return m_face_kernel->is_valid(f);
    }

    bool ISMesh::exists(const EdgeKey& e) {
        return m_edge_kernel->is_valid(e);
    }

    bool ISMesh::exists(const NodeKey& n) {
        return m_node_kernel->is_valid(n);
    }

    bool ISMesh::is_clockwise_order(const NodeKey& nid, SimplexSet<NodeKey>& nids) {
        auto x = get(nid).get_pos() - get(nids[0]).get_pos();
        auto y = get(nids[1]).get_pos() - get(nids[0]).get_pos();
        auto z = get(nids[2]).get_pos() - get(nids[0]).get_pos();
        auto val = dot(x, cross(y,z));

        return val > 0.;
    }

    void ISMesh::orient_cc(const NodeKey& nid, SimplexSet<NodeKey>& nids) {
        if(is_clockwise_order(nid, nids))
        {
            nids.swap();
        }
    }

    bool ISMesh::is_inverted(const TetrahedronKey& tid) {
        for(const FaceKey& f : get_faces(tid))
        {
            const SimplexSet<TetrahedronKey>& tids = get_tets(f);
            if(tids.size() == 2) // Check that f is not a boundary face.
            {
                SimplexSet<NodeKey> nids = get_nodes(f);
                SimplexSet<NodeKey> apices = get_nodes(tids) - nids;
                auto normal = cross(get_pos(nids[0]) - get_pos(nids[2]), get_pos(nids[1]) - get_pos(nids[2]));
                auto d1 = dot(get_pos(apices[0]) - get_pos(nids[2]), normal);
                auto d2 = dot(get_pos(apices[1]) - get_pos(nids[2]), normal);
                if((d1 < 0. && d2 < 0) || (d1 > 0. && d2 > 0.))
                {
                    return true;
                }
            }
        }
        return false;
    }

    bool ISMesh::is_inverted_destination(const TetrahedronKey& tid) {
        for(const FaceKey& f : get_faces(tid))
        {
            const SimplexSet<TetrahedronKey>& tids = get_tets(f);
            if(tids.size() == 2) // Check that f is not a boundary face.
            {
                SimplexSet<NodeKey> nids = get_nodes(f);
                SimplexSet<NodeKey> apices = get_nodes(tids) - nids;
                auto normal = cross(get(nids[0]).get_destination() - get(nids[2]).get_destination(), get(nids[1]).get_destination() - get(nids[2]).get_destination());
                auto d1 = dot(get(apices[0]).get_destination() - get(nids[2]).get_destination(), normal);
                auto d2 = dot(get(apices[1]).get_destination() - get(nids[2]).get_destination(), normal);
                if((d1 < 0. && d2 < 0) || (d1 > 0. && d2 > 0.))
                {
                    return true;
                }
            }
        }
        return false;
    }

    NodeKey ISMesh::insert_node(const vec3& p) {
        auto node = m_node_kernel->create(NodeAttributes(p), this);
        return node.key();
    }

    EdgeKey ISMesh::insert_edge(NodeKey node1, NodeKey node2) {
        auto edge = m_edge_kernel->create(EdgeAttributes(), this);
        //add the new simplex to the co-boundary relation of the boundary simplices
        get(node1).add_co_face(edge.key());
        get(node2).add_co_face(edge.key());
        //set the boundary relation
        edge->add_face(node1);
        edge->add_face(node2);
        return edge.key();
    }

    FaceKey ISMesh::insert_face(EdgeKey edge1, EdgeKey edge2, EdgeKey edge3) {
        auto face = m_face_kernel->create(FaceAttributes(), this);
        //update relations
        get(edge1).add_co_face(face.key());
        get(edge2).add_co_face(face.key());
        get(edge3).add_co_face(face.key());
        face->add_face(edge1);
        face->add_face(edge2);
        face->add_face(edge3);
        return face.key();
    }

    TetrahedronKey ISMesh::insert_tetrahedron(FaceKey face1, FaceKey face2, FaceKey face3, FaceKey face4) {
        auto tetrahedron = m_tetrahedron_kernel->create(TetAttributes(), this);
        //update relations
        get(face1).add_co_face(tetrahedron.key());
        get(face2).add_co_face(tetrahedron.key());
        get(face3).add_co_face(tetrahedron.key());
        get(face4).add_co_face(tetrahedron.key());
        tetrahedron->add_face(face1);
        tetrahedron->add_face(face2);
        tetrahedron->add_face(face3);
        tetrahedron->add_face(face4);

        return tetrahedron.key();
    }

    void ISMesh::remove(const NodeKey& nid) {
        for(auto e : get_edges(nid))
        {
            get(e).remove_face(nid);
        }
        m_node_kernel->erase(nid);
    }

    void ISMesh::remove(const EdgeKey& eid) {
        for(auto f : get_faces(eid))
        {
            get(f).remove_face(eid);
        }
        for(auto n : get_nodes(eid))
        {
            get(n).remove_co_face(eid);
        }
        m_edge_kernel->erase(eid);
    }

    void ISMesh::remove(const FaceKey& fid) {
        for(auto t : get_tets(fid))
        {
            get(t).remove_face(fid);
        }
        for(auto e : get_edges(fid))
        {
            get(e).remove_co_face(fid);
        }
        m_face_kernel->erase(fid);
    }

    void ISMesh::remove(const TetrahedronKey& tid) {
        for(auto f : get_faces(tid))
        {
            get(f).remove_co_face(tid);
        }
        m_tetrahedron_kernel->erase(tid);
    }

    NodeKey ISMesh::merge(const NodeKey& key1, const NodeKey& key2) {
        for(auto e : get_edges(key2))
        {
            connect(key1, e);
        }
        remove(key2);
        return key1;
    }

    void ISMesh::update_split(const NodeKey& nid_new, const NodeKey& nid1, const NodeKey& nid2) {
         throw "Not implemented";
    }

    NodeKey ISMesh::split(const EdgeKey& eid, const vec3& pos, const vec3& destination) {
        auto nids = get_nodes(eid);
        auto fids = get_faces(eid);
        auto tids = get_tets(eid);

        // Split edge
        auto new_nid = insert_node(pos);
        get(new_nid).set_destination(destination);

        disconnect(nids[1], eid);
        connect(new_nid, eid);

        auto new_eid = insert_edge(new_nid, nids[1]);

        // Update faces, create faces
        for (auto f : fids)
        {
            EdgeKey f_eid = (get_edges(f) & get_edges(nids[1])).front();
            disconnect(f_eid, f);

            SimplexSet<NodeKey> new_e_nids = get_nodes(get_edges(f)) - nids[0];
            assert(new_e_nids.size() == 2);
            EdgeKey new_f_eid = insert_edge(new_e_nids[0], new_e_nids[1]);
            connect(new_f_eid, f);

            insert_face(f_eid, new_f_eid, new_eid);
        }

        // Update tetrahedra, create tetrahedra
        SimplexSet<TetrahedronKey> new_tids;
        for (auto t : tids)
        {
            FaceKey t_fid = (get_faces(t) - get_faces(nids[0])).front();
            disconnect(t_fid, t);

            SimplexSet<EdgeKey> new_f_eids = get_edges(get_faces(t)) - get_edges(nids[0]);
            assert(new_f_eids.size() == 3);
            FaceKey new_t_fid = insert_face(new_f_eids[0], new_f_eids[1], new_f_eids[2]);
            connect(new_t_fid, t);

            SimplexSet<FaceKey> t_fids = get_faces(new_eid) & get_faces(get_edges(t_fid));
            assert(t_fids.size() == 2);
            new_tids += insert_tetrahedron(t_fids[0], t_fids[1], new_t_fid, t_fid);
        }

        // Update flags
        for (unsigned int i = 0; i < tids.size(); i++)
        {
            set_label(new_tids[i], get_label(tids[i]));
        }

        update_split(new_nid, nids[0], nids[1]);

        return new_nid;
    }

    void ISMesh::update_collapse(const NodeKey& nid, const NodeKey& nid_removed, double weight) {
        Node& node = get(nid);
        Node& node_removed = get(nid_removed);
        node.set_pos((1.-weight) * node.get_pos() + weight * node_removed.get_pos());
        node.set_destination((1.-weight) * node.get_destination() + weight * node_removed.get_destination());
    }

    void ISMesh::collapse(const EdgeKey& eid, const NodeKey& nid, double weight) {
        NodeKey nid_remove = (get_nodes(eid) - nid).front();
        update_collapse(nid, nid_remove, weight);

        auto fids = get_faces(eid);
        auto tids = get_tets(eid);

        // Remove edge
        remove(eid);

        // Remove faces
        for(auto f : fids)
        {
            SimplexSet<EdgeKey> eids = get_edges(f);
            remove(f);
            if(get_nodes(eids[1]).contains(nid))
            {
                eids.swap();
            }
            merge(eids[0], eids[1]);
        }

        // Remove tetrahedra
        for(auto t : tids)
        {
            SimplexSet<FaceKey> fids = get_faces(t);
            remove(t);
            if(get_nodes(fids[1]).contains(nid))
            {
                fids.swap();
            }
            merge(fids[0], fids[1]);
        }

        // Merge nodes
        merge(nid, nid_remove);

        // Update flags.
        update(get_tets(nid));
    }

    FaceKey ISMesh::flip_32(const EdgeKey& eid) {
        SimplexSet<NodeKey> e_nids = get_nodes(eid);
        SimplexSet<FaceKey> e_fids = get_faces(eid);
#ifdef DEBUG
        assert(e_fids.size() == 3);
#endif
        SimplexSet<TetrahedronKey> e_tids = get_tets(e_fids);
        int label = get_label(e_tids[0]);
#ifdef DEBUG
        assert(e_tids.size() == 3);
        assert(label == get_label(e_tids[1]));
        assert(label == get_label(e_tids[2]));
#endif

        // Remove edge
        remove(eid);

        // Create face
        SimplexSet<EdgeKey> f_eids = get_edges(e_tids) - get_edges(e_fids);
#ifdef DEBUG
        assert(f_eids.size() == 3);
#endif
        FaceKey new_fid = insert_face(f_eids[0], f_eids[1], f_eids[2]);

        // Remove faces
        for(const FaceKey& f : e_fids)
        {
            remove(f);
        }

        // Create tetrahedra
        SimplexSet<FaceKey> exterior_fids = get_faces(e_tids);
        for (const NodeKey& n : e_nids)
        {
            SimplexSet<FaceKey> t_fids = exterior_fids & get_faces(n);
#ifdef DEBUG
            assert(t_fids.size() == 3);
#endif
            insert_tetrahedron(t_fids[0], t_fids[1], t_fids[2], new_fid);
        }

        // Remove tetrahedra
        for(const TetrahedronKey& t : e_tids)
        {
            remove(t);
        }

        // Update flags
#ifdef DEBUG
        assert(get_tets(new_fid).size() == 2);
#endif
        for (auto t : get_tets(new_fid)) {
            set_label(t, label);
        }
        return new_fid;
    }

    EdgeKey ISMesh::flip_23(const FaceKey& fid) {
        SimplexSet<TetrahedronKey> f_tids = get_tets(fid);
#ifdef DEBUG
        assert(f_tids.size() == 2);
#endif
        SimplexSet<EdgeKey> f_eids = get_edges(fid);
#ifdef DEBUG
        assert(f_eids.size() == 3);
#endif
        SimplexSet<NodeKey> f_nids = get_nodes(f_eids);
        int label = get_label(f_tids[0]);
#ifdef DEBUG
        assert(label == get_label(f_tids[1]));
#endif

        // Create edge
        SimplexSet<NodeKey> e_nids = get_nodes(f_tids) - f_nids;
#ifdef DEBUG
        assert(e_nids.size() == 2);
#endif
        EdgeKey new_eid = insert_edge(e_nids[0], e_nids[1]);

        // Create faces
        SimplexSet<EdgeKey> new_fs_eids = get_edges(f_tids) - f_eids;
#ifdef DEBUG
        assert(new_fs_eids.size() == 6);
#endif
        for (const NodeKey& n : f_nids)
        {
            auto new_f_eids = new_fs_eids & get_edges(n);
#ifdef DEBUG
            assert(new_f_eids.size() == 2);
#endif
            insert_face(new_f_eids[0], new_f_eids[1], new_eid);
        }

        // Remove face
        remove(fid);

        // Create tetrahedra
        SimplexSet<FaceKey> new_ts_fids1 = get_faces(f_tids);
        SimplexSet<FaceKey> new_ts_fids2 = get_faces(new_eid);
#ifdef DEBUG
        assert(new_ts_fids1.size() == 6);
        assert(new_ts_fids2.size() == 3);
#endif
        for (const EdgeKey& e : f_eids)
        {
            SimplexSet<FaceKey> new_t_fids = (new_ts_fids1 & get_faces(e)) + (new_ts_fids2 & get_faces(get_nodes(e)));
#ifdef DEBUG
            assert(new_t_fids.size() == 4);
#endif
            insert_tetrahedron(new_t_fids[0], new_t_fids[1], new_t_fids[2], new_t_fids[3]);
        }

        // Remove tetrahedra
        for (auto t : f_tids) {
            remove(t);
        }

        // Update flags
#ifdef DEBUG
        assert(get_tets(new_eid).size() == 3);
#endif
        for (auto t : get_tets(new_eid)) {
            set_label(t, label);
        }
        return new_eid;
    }

    void ISMesh::flip(const EdgeKey& eid, const FaceKey& fid1, const FaceKey& fid2) {
        SimplexSet<FaceKey> fids = {fid1, fid2};
        SimplexSet<NodeKey> e_nids = get_nodes(eid);
        SimplexSet<FaceKey> e_fids = get_faces(eid);
        SimplexSet<TetrahedronKey> e_tids = get_tets(e_fids);

        // Reconnect edge
        SimplexSet<NodeKey> new_e_nids = get_nodes(fids) - e_nids;
#ifdef DEBUG
        assert(new_e_nids.size() == 2);
#endif

        disconnect(e_nids[0], eid);
        disconnect(e_nids[1], eid);
        connect(new_e_nids[0], eid);
        connect(new_e_nids[1], eid);

        // Reconnect faces
        SimplexSet<EdgeKey> swap_eids = {get_edge(e_nids[0], new_e_nids[0]), get_edge(e_nids[1], new_e_nids[1])};
        swap(swap_eids[0], fid1, swap_eids[1], fid2);

#ifdef DEBUG
        assert((e_fids - fids).size() <= 2);
#endif
        for(FaceKey f : (e_fids - fids))
        {
            SimplexSet<EdgeKey> rm_eids = get_edges(f) - eid;
#ifdef DEBUG
            assert(rm_eids.size() == 2);
#endif
            SimplexSet<NodeKey> apex = get_nodes(f) - (new_e_nids + e_nids);
#ifdef DEBUG
            assert(apex.size() == 1);
#endif

            SimplexSet<EdgeKey> add_eids = {get_edge(apex[0], new_e_nids[0]), get_edge(apex[0], new_e_nids[1])};

            disconnect(rm_eids[0], f);
            disconnect(rm_eids[1], f);
            connect(add_eids[0], f);
            connect(add_eids[1], f);

            // Reconnect tetrahedra
            SimplexSet<TetrahedronKey> tids = get_tets(f);
#ifdef DEBUG
            assert(tids.size() == 2);
#endif
            SimplexSet<FaceKey> swap_fids = (get_faces(tids) - e_fids) & get_faces(swap_eids);
#ifdef DEBUG
            assert(swap_fids.size() == 2);
#endif
            swap(swap_fids[0], tids[0], swap_fids[1], tids[1]);
        }

        // Update flags
        update(e_tids);
    }

    void ISMesh::flip_22(const FaceKey& fid1, const FaceKey& fid2) {
        SimplexSet<EdgeKey> eid = get_edges(fid1) & get_edges(fid2);
#ifdef DEBUG
        assert(eid.size() == 1);
        assert(get_faces(eid).size() == 3);
        assert(get_tets(eid).size() == 2);
#endif

        flip(eid[0], fid1, fid2);
    }

    void ISMesh::flip_44(const FaceKey& fid1, const FaceKey& fid2) {
        SimplexSet<EdgeKey> eid = get_edges(fid1) & get_edges(fid2);
#ifdef DEBUG
        assert(eid.size() == 1);
        assert(get_faces(eid).size() == 4);
        assert(get_tets(eid).size() == 4);
#endif

        flip(eid[0], fid1, fid2);
    }

    void ISMesh::garbage_collect() {
        m_node_kernel->garbage_collect();
        m_edge_kernel->garbage_collect();
        m_face_kernel->garbage_collect();
        m_tetrahedron_kernel->garbage_collect();
    }

    void ISMesh::scale(const vec3& s) {
        for (auto nit = nodes_begin(); nit != nodes_end(); nit++) {
            nit->set_pos(s*nit->get_pos());
            nit->set_destination(s*nit->get_destination());
        }
    }

    void ISMesh::extract_surface_mesh(vector<vec3>& points, vector<int>& faces) {
        garbage_collect();

        map<NodeKey, int> indices;
        // Extract vertices
        for (auto nit = nodes_begin(); nit != nodes_end(); nit++)
        {
            if (nit->is_interface())
            {
                points.push_back(nit->get_pos());
                indices[nit.key()] = static_cast<int>(points.size());
            }
        }

        // Extract faces
        for (auto fit = faces_begin(); fit != faces_end(); fit++)
        {
            if (fit->is_interface())
            {
                for (auto &n : get_sorted_nodes(fit.key()))
                {
                    faces.push_back(indices[n]);
                }
            }
        }
    }

    void ISMesh::extract_tet_mesh(vector<vec3>& points, vector<int>& tets, vector<int>& tet_labels) {
        garbage_collect();

        map<NodeKey, int> indices;
        // Extract vertices
        for (auto nit = nodes_begin(); nit != nodes_end(); nit++)
        {
            indices[nit.key()] = static_cast<int>(points.size());
            points.push_back(nit->get_pos());
        }

        for (auto tit = tetrahedra_begin(); tit != tetrahedra_end(); tit++)
        {
            for (auto n : get_nodes(tit.key()))
            {
                tets.push_back(indices[n]);
            }
            tet_labels.push_back(get_label(tit.key()));
        }
    }

    void ISMesh::validity_check() {
        cout << "Testing connectivity of simplicial complex: ";
        for(auto tit = tetrahedra_begin(); tit != tetrahedra_end(); tit++)
        {
            assert(exists(tit.key()));
            // Check faces:
            auto faces = get_faces(tit.key());
            assert(faces.size() == 4);
            for (auto f : faces) {
                assert(exists(f));
                auto cotets = get_tets(f);
                assert((get(f).is_boundary() && cotets.size() == 1) || (!get(f).is_boundary() && cotets.size() == 2));
                assert(find(cotets.begin(), cotets.end(), tit.key()) != cotets.end());
                for (auto f2 : faces) {
                    assert(f == f2 || get_edge(f, f2).is_valid());
                }

                // Check edges:
                auto edges = get_edges(f);
                assert(edges.size() == 3);
                for (auto e : edges)
                {
                    assert(exists(e));
                    auto cofaces = get_faces(e);
                    assert(find(cofaces.begin(), cofaces.end(), f) != cofaces.end());
                    for (auto e2 : edges) {
                        assert(e == e2 || get_node(e, e2).is_valid());
                    }

                    // Check nodes:
                    auto nodes = get_nodes(e);
                    assert(nodes.size() == 2);
                    for (auto n : nodes)
                    {
                        assert(exists(n));
                        auto coedges = get_edges(n);
                        assert(find(coedges.begin(), coedges.end(), e) != coedges.end());
                    }
                }

            }

            assert(get_edges(tit.key()).size() == 6);
            assert(get_nodes(tit.key()).size() == 4);
        }
        cout << "PASSED" << endl;

        cout << "Testing for inverted tetrahedra: ";
        for (auto tit = tetrahedra_begin(); tit != tetrahedra_end(); tit++)
        {
            assert(!is_inverted(tit.key()));
        }
        cout << "PASSED" << endl;

        cout << "Testing for corrupted interface or boundary: ";
        for (auto eit = edges_begin(); eit != edges_end(); eit++)
        {
            int boundary = 0;
            int interface = 0;
            for (auto f : get_faces(eit.key())) {
                if(get(f).is_boundary())
                {
                    boundary++;
                }
                if(get(f).is_interface())
                {
                    interface++;
                }
            }
            assert((eit->is_interface() && interface >= 2) || (!eit->is_interface() && interface == 0)); // Check that the interface is not corrupted
            assert((eit->is_boundary() && boundary == 2) || (!eit->is_boundary() && boundary == 0)); // Check that the boundary is not corrupted
        }
        cout << "PASSED" << endl;
    }
}
