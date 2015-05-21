#include "is_mesh.h"
#include "query.h"

#include <algorithm>

using namespace std;

namespace is_mesh{
    ISMesh::ISMesh(vector<vec3> & points, vector<int> & tets, const vector<int>& tet_labels) {
        create(points, tets);
        validity_check(true);
        init_flags(tet_labels);
        validity_check();
    }

    ISMesh::ISMesh(const  ISMesh &mesh) {
        vector<NodeKey> nodeKeys;
        vector<EdgeKey> edgeKeys;
        vector<FaceKey> faceKeys;
        vector<TetrahedronKey> tetKeys;

        map<NodeKey,NodeKey> nodeToNew;
        map<EdgeKey,EdgeKey> edgeToNew;
        map<FaceKey,FaceKey> faceToNew;
        map<TetrahedronKey,TetrahedronKey> tetToNew;

        // create all primitives
        for (auto & node:mesh.nodes()){
            auto k = m_node_kernel.create(this, node.get_pos()).key();
            nodeKeys.push_back(node.key());
            nodeToNew[node.key()] = k;
        }

        for (auto & edge : mesh.edges()){
            auto k = m_edge_kernel.create(this).key();
            edgeKeys.push_back(edge.key());
            edgeToNew[edge.key()] = k;
        }

        for (auto & face : mesh.faces()){
            auto k = m_face_kernel.create(this).key();
            faceKeys.push_back(face.key());
            faceToNew[face.key()] = k;
        }

        for (auto & tet : mesh.tetrahedra()){
            auto k = m_tetrahedron_kernel.create(this, tet.label()).key();
            tetKeys.push_back(tet.key());
            tetToNew[tet.key()] = k;
        }

        // copy connectivity
        for (auto & node : nodes()){
            const auto & orignal = mesh.get(nodeKeys[(int)node.key()]);
            for (auto & cb : orignal.m_co_boundary){
                if (edgeToNew.find(cb) != edgeToNew.end())
                    node.m_co_boundary += edgeToNew[cb];
            }
        }
        for (auto & edge : edges()){
            const auto & orignal = mesh.get(edgeKeys[(int)edge.key()]);
            for (auto & b : orignal.m_boundary){
                if (nodeToNew.find(b) != nodeToNew.end())
                    edge.m_boundary += nodeToNew[b];
            }
            for (auto & cb : orignal.m_co_boundary){
                if (faceToNew.find(cb) != faceToNew.end())
                    edge.m_co_boundary += faceToNew[cb];
            }
        }
        for (auto & face : faces()){
            const auto & orignal = mesh.get(faceKeys[(int)face.key()]);
            for (auto & b : orignal.m_boundary){
                if (edgeToNew.find(b) != edgeToNew.end())
                    face.m_boundary += edgeToNew[b];
            }
            for (auto & cb : orignal.m_co_boundary){
                if (tetToNew.find(cb) != tetToNew.end())
                    face.m_co_boundary += tetToNew[cb];
            }
        }
        for (auto & tet : tetrahedra()){
            const auto & orignal = mesh.get(tetKeys[(int)tet.key()]);
            for (auto & b : orignal.m_boundary){
                if (faceToNew.find(b) != faceToNew.end())
                    tet.m_boundary += faceToNew[b];
            }
        }
        update_flag();
    }

    ISMesh::~ISMesh() {
    }

    unsigned int ISMesh::get_no_nodes() const {
        return static_cast<unsigned int>(m_node_kernel.size());
    }

    unsigned int ISMesh::get_no_edges() const {
        return static_cast<unsigned int>(m_edge_kernel.size());
    }

    unsigned int ISMesh::get_no_faces() const {
        return static_cast<unsigned int>(m_face_kernel.size());
    }

    unsigned int ISMesh::get_no_tets() const {
        return static_cast<unsigned int>(m_tetrahedron_kernel.size());
    }

    NodeIterator ISMesh::nodes() const {
        return NodeIterator{&m_node_kernel};
    }

    EdgeIterator ISMesh::edges() const {
        return EdgeIterator{&m_edge_kernel};
    }

    FaceIterator ISMesh::faces() const {
        return FaceIterator{&m_face_kernel};
    }

    TetrahedronIterator ISMesh::tetrahedra() const {
        return TetrahedronIterator{&m_tetrahedron_kernel};
    }

    typename kernel<NodeKey,Node>::iterator ISMesh::nodes_begin() {
        return m_node_kernel.begin();
    }

    typename kernel<NodeKey,Node>::iterator ISMesh::nodes_end() {
        return m_node_kernel.end();
    }

    typename kernel<EdgeKey,Edge>::iterator ISMesh::edges_begin() {
        return m_edge_kernel.begin();
    }

    typename kernel<EdgeKey,Edge>::iterator ISMesh::edges_end() {
        return m_edge_kernel.end();
    }

    typename kernel<FaceKey,Face>::iterator ISMesh::faces_begin() {
        return m_face_kernel.begin();
    }

    typename kernel<FaceKey,Face>::iterator ISMesh::faces_end() {
        return m_face_kernel.end();
    }

    typename kernel<TetrahedronKey,Tetrahedron>::iterator ISMesh::tetrahedra_begin() {
        return m_tetrahedron_kernel.begin();
    }

    typename kernel<TetrahedronKey,Tetrahedron>::iterator ISMesh::tetrahedra_end() {
        return m_tetrahedron_kernel.end();
    }

    void ISMesh::set_label(const TetrahedronKey& tid, int newLabel) {
        int oldLabel = get(tid).label();
        get(tid).label(newLabel);
        SimplexSet<TetrahedronKey> tids = {tid};
        update_flag(tids);

        for (auto & l : m_set_label_listeners){
            l.second(tid, oldLabel);
        }
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
            for (auto & tit : tetrahedra())
            {
                tit.label(tet_labels[tit.key()]);
            }
        }

        for (auto & fit : faces())
        {
            update_flag(fit.key());
        }

        for (auto & eit : edges())
        {
            update_flag(eit.key());
        }

        for (auto & nit : nodes())
        {
            update_flag(nit.key());
        }
    }

    void ISMesh::update_flag() {
        // Update faces

        for (auto & f : faces())
        {
            update_flag(f.key());
        }

        // Update edges
        for (auto & e : edges())
        {
            update_flag(e.key());
        }

        // Update nodes
        for (auto & n : nodes())
        {
            update_flag(n.key());
        }
    }

    void ISMesh::update_flag(const SimplexSet<TetrahedronKey>& tids) {
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
        Face & face = get(f);
        face.set_interface(false);

        const SimplexSet<TetrahedronKey> &tids = get(f).tet_keys();
        int size = tids.size();
        if (size == 1)
        {
            if (get(tids.front()).label() != 0)
            {
                face.set_interface(true);
            }
        }
        else
        {

            assert(size == 2);
            if (get(tids.front()).label() != get(tids.back()).label())
            {
                // On the interface
                face.set_interface(true);
            }
        }
    }

    void ISMesh::update_flag(const EdgeKey & e) {
        Edge & edge = get(e);
        edge.set_boundary(false);
        edge.set_interface(false);
        edge.set_crossing(false);

        int i = 0;
        for (auto f : get(e).face_keys())
        {
            if (exists(f))
            {
                if (get(f).is_boundary())
                {
                    edge.set_boundary(true);
                }
                if (get(f).is_interface())
                {
                    edge.set_interface(true);
                    i++;
                }
            }
        }
        if(i > 2)
        {
            edge.set_crossing(true);
        }
    }

    void ISMesh::connected_component(SimplexSet<TetrahedronKey>& tids, const TetrahedronKey& tid) {
                int label = get(tid).label();
                tids -= tid;

                for(auto f : get(tid).face_keys())
                {
                    if(get(f).tet_keys().size() == 2)
                    {
                        TetrahedronKey tid2 = get_tet(tid, f);
                        if(tids.contains(tid2) && label == get(tid2).label())
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
        Node & node = get(n);
        node.set_interface(false);
        node.set_boundary(false);
        node.set_crossing(false);

        for (auto e : get(n).edge_keys())
        {
            if (exists(e))
            {
                if (get(e).is_interface())
                {
                    node.set_interface(true);
                }
                if (get(e).is_boundary())
                {
                    node.set_boundary(true);
                }
                if (get(e).is_crossing())
                {
                    node.set_crossing(true);
                }
            }
        }
        if(!get(n).is_crossing() && get(n).is_interface() && crossing(n))
        {
            node.set_crossing(true);
        }
    }

    is_mesh::Node &ISMesh::get(const NodeKey& nid) {
        return m_node_kernel.find(nid);
    }

    is_mesh::Edge &ISMesh::get(const EdgeKey& eid) {
        return m_edge_kernel.find(eid);
    }

    is_mesh::Face &ISMesh::get(const FaceKey& fid) {
        return m_face_kernel.find(fid);
    }

    is_mesh::Tetrahedron &ISMesh::get(const TetrahedronKey& tid) {
        return m_tetrahedron_kernel.find(tid);
    }

    const is_mesh::Node &ISMesh::get(const NodeKey& nid) const {
        return m_node_kernel.find(nid);
    }

    const is_mesh::Edge &ISMesh::get(const EdgeKey& eid) const {
        return m_edge_kernel.find(eid);
    }

    const is_mesh::Face &ISMesh::get(const FaceKey& fid) const {
        return m_face_kernel.find(fid);
    }

    const is_mesh::Tetrahedron &ISMesh::get(const TetrahedronKey& tid) const  {
        return m_tetrahedron_kernel.find(tid);
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
        Face& f = get(fid);
        if (f.is_interface())
        {
            int label = -100;
            TetrahedronKey tid;
            for (auto t : f.tet_keys())
            {
                int tl = get(t).label();
                if (tl > label)
                {
                    label = tl;
                    tid = t;
                }
            }
            return get_sorted_nodes(fid, tid);
        }
        else if (f.is_boundary())
        {
            TetrahedronKey tid = f.tet_keys().front();
            return get_sorted_nodes(fid, tid);
        }
        return get_nodes(fid);
    }

    SimplexSet<NodeKey> ISMesh::get_nodes(const FaceKey& fid) {
        const SimplexSet<EdgeKey>& eids = get(fid).edge_keys();
        SimplexSet<NodeKey> nids = get(eids[0]).node_keys();
        nids += get(eids[1]).node_keys();
        return nids;
    }

    SimplexSet<NodeKey> ISMesh::get_nodes(const TetrahedronKey& tid) {
        const SimplexSet<FaceKey>& fids = get(tid).face_keys();
        SimplexSet<NodeKey> nids = get_nodes(fids[0]);
        nids += get_nodes(fids[1]);
        return nids;
    }

    SimplexSet<EdgeKey> ISMesh::get_edges(const TetrahedronKey& tid) {
        SimplexSet<EdgeKey> eids;
        for(const FaceKey& f : get(tid).face_keys())
        {
            eids += get(f).edge_keys();
        }
        return eids;
    }

    SimplexSet<FaceKey> ISMesh::get_faces(const NodeKey& nid) {
        SimplexSet<FaceKey> fids;
        for(const EdgeKey& e : get(nid).edge_keys())
        {
            fids += get(e).face_keys();
        }
        return fids;
    }

    SimplexSet<TetrahedronKey> ISMesh::get_tets(const NodeKey& nid) {
        SimplexSet<TetrahedronKey> tids;
        for(const EdgeKey& e : get(nid).edge_keys())
        {
            for(const FaceKey& f : get(e).face_keys())
            {
                tids += get(f).tet_keys();
            }
        }
        return tids;
    }

    SimplexSet<TetrahedronKey> ISMesh::get_tets(const EdgeKey& eid) {
        SimplexSet<TetrahedronKey> tids;
        for(const FaceKey& f : get(eid).face_keys())
        {
            tids += get(f).tet_keys();
        }
        return tids;
    }

    NodeKey ISMesh::get_node(const EdgeKey& eid1, const EdgeKey& eid2) {
        const SimplexSet<NodeKey>& nids = get(eid2).node_keys();
        for (const NodeKey& n : get(eid1).node_keys()) {
            if(nids.contains(n))
            {
                return n;
            }
        }
        return NodeKey();
    }

    NodeKey ISMesh::get_node(const EdgeKey& eid, const NodeKey& nid) {
        const SimplexSet<NodeKey>& nids = get(eid).node_keys();
        if(nids[0] == nid)
        {
            return nids[1];
        }
        return nids[0];
    }

    EdgeKey ISMesh::get_edge(const NodeKey& nid1, const NodeKey& nid2) {
        const SimplexSet<EdgeKey>& eids = get(nid2).edge_keys();
        for (const EdgeKey& e : get(nid1).edge_keys()) {
            if(eids.contains(e))
            {
                return e;
            }
        }
        return EdgeKey();
    }

    EdgeKey ISMesh::get_edge(const FaceKey& fid1, const FaceKey& fid2) {
        const SimplexSet<EdgeKey>& eids = get(fid2).edge_keys();
        for (const EdgeKey& e : get(fid1).edge_keys()) {
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
        const SimplexSet<FaceKey>& fids = get(tid2).face_keys();
        for (const FaceKey& f : get(tid1).face_keys()) {
            if(fids.contains(f))
            {
                return f;
            }
        }
        return FaceKey();
    }

    TetrahedronKey ISMesh::get_tet(const TetrahedronKey& tid, const FaceKey& fid) {
        for(TetrahedronKey t : get(fid).tet_keys())
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
            verts[i] = get(nids[i]).get_pos();
        }
        return verts;
    }

    bool ISMesh::exists(const TetrahedronKey& t) {
        return m_tetrahedron_kernel.is_valid(t);
    }

    bool ISMesh::exists(const FaceKey& f) {
        return m_face_kernel.is_valid(f);
    }

    bool ISMesh::exists(const EdgeKey& e) {
        return m_edge_kernel.is_valid(e);
    }

    bool ISMesh::exists(const NodeKey& n) {
        return m_node_kernel.is_valid(n);
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
        for(const FaceKey& f : get(tid).face_keys())
        {
            const SimplexSet<TetrahedronKey>& tids = get(f).tet_keys();
            if(tids.size() == 2) // Check that f is not a boundary face.
            {
                SimplexSet<NodeKey> nids = get(f).node_keys();
                SimplexSet<NodeKey> apices = get_nodes(tids) - nids;
                auto normal = cross(get(nids[0]).get_pos() - get(nids[2]).get_pos(), get(nids[1]).get_pos() - get(nids[2]).get_pos());
                auto d1 = dot(get(apices[0]).get_pos() - get(nids[2]).get_pos(), normal);
                auto d2 = dot(get(apices[1]).get_pos() - get(nids[2]).get_pos(), normal);
                if((d1 < 0. && d2 < 0) || (d1 > 0. && d2 > 0.))
                {
                    return true;
                }
            }
        }
        return false;
    }

    bool ISMesh::is_inverted_destination(const TetrahedronKey& tid) {
        for(const FaceKey& f : get(tid).face_keys())
        {
            const SimplexSet<TetrahedronKey>& tids = get(f).tet_keys();
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
        auto node = m_node_kernel.create(this,p);
        return node.key();
    }

    EdgeKey ISMesh::insert_edge(NodeKey node1, NodeKey node2) {
        auto edge = m_edge_kernel.create(this);
        //add the new simplex to the co-boundary relation of the boundary simplices
        get(node1).add_co_face(edge.key());
        get(node2).add_co_face(edge.key());
        //set the boundary relation
        edge->add_face(node1);
        edge->add_face(node2);
        return edge.key();
    }

    FaceKey ISMesh::insert_face(EdgeKey edge1, EdgeKey edge2, EdgeKey edge3) {
        auto face = m_face_kernel.create(this);
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
        auto tetrahedron = m_tetrahedron_kernel.create(this);
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
        for(auto e : get(nid).edge_keys())
        {
            get(e).remove_face(nid);
        }
        m_node_kernel.erase(nid);
    }

    void ISMesh::remove(const EdgeKey& eid) {
        for(auto f : get(eid).face_keys())
        {
            get(f).remove_face(eid);
        }
        for(auto n : get(eid).node_keys())
        {
            get(n).remove_co_face(eid);
        }
        m_edge_kernel.erase(eid);
    }

    void ISMesh::remove(const FaceKey& fid) {
        for(auto t : get(fid).tet_keys())
        {
            get(t).remove_face(fid);
        }
        for(auto e : get(fid).edge_keys())
        {
            get(e).remove_co_face(fid);
        }
        m_face_kernel.erase(fid);
    }

    void ISMesh::remove(const TetrahedronKey& tid) {
        for(auto f : get(tid).face_keys())
        {
            get(f).remove_co_face(tid);
        }
        m_tetrahedron_kernel.erase(tid);
    }

    NodeKey ISMesh::merge(const NodeKey& key1, const NodeKey& key2) {
        for(auto e : get(key2).edge_keys())
        {
            connect(key1, e);
        }
        remove(key2);
        return key1;
    }

    void ISMesh::update_split(const NodeKey& nid_new, const NodeKey& nid1, const NodeKey& nid2) {
        for (auto & sl : m_split_listeners){
            sl.second(nid_new, nid1, nid2);
        }
    }

    NodeKey ISMesh::split(const EdgeKey& eid, const vec3& pos, const vec3& destination) {
        Edge &e = get(eid);
        auto nids = e.node_keys();
        auto fids = e.face_keys();
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
            EdgeKey f_eid = (get(f).edge_keys() & get(nids[1]).edge_keys()).front();
            disconnect(f_eid, f);

            SimplexSet<NodeKey> new_e_nids = get_nodes(get(f).edge_keys()) - nids[0];
            assert(new_e_nids.size() == 2);
            EdgeKey new_f_eid = insert_edge(new_e_nids[0], new_e_nids[1]);
            connect(new_f_eid, f);

            insert_face(f_eid, new_f_eid, new_eid);
        }

        // Update tetrahedra, create tetrahedra
        SimplexSet<TetrahedronKey> new_tids;
        for (auto t : tids)
        {
            FaceKey t_fid = (get(t).face_keys() - get_faces(nids[0])).front();
            disconnect(t_fid, t);

            SimplexSet<EdgeKey> new_f_eids = get_edges(get(t).face_keys()) - get(nids[0]).edge_keys();
            assert(new_f_eids.size() == 3);
            FaceKey new_t_fid = insert_face(new_f_eids[0], new_f_eids[1], new_f_eids[2]);
            connect(new_t_fid, t);

            SimplexSet<FaceKey> t_fids = get(new_eid).face_keys() & get_faces(get(t_fid).edge_keys());
            assert(t_fids.size() == 2);
            new_tids += insert_tetrahedron(t_fids[0], t_fids[1], new_t_fid, t_fid);
        }

        // Update flags
        for (unsigned int i = 0; i < tids.size(); i++)
        {
            set_label(new_tids[i], get(tids[i]).label());
        }

        update_split(new_nid, nids[0], nids[1]);

#ifdef DEBUG
        validity_check();
#endif

        return new_nid;
    }

    void ISMesh::update_collapse(const NodeKey& nid, const NodeKey& nid_removed, double weight) {
        Node& node = get(nid);
        Node& node_removed = get(nid_removed);
        node.set_pos((1.-weight) * node.get_pos() + weight * node_removed.get_pos());
        node.set_destination((1.-weight) * node.get_destination() + weight * node_removed.get_destination());

        for (auto & c : m_collapse_listeners){
            c.second(nid, nid_removed, weight);
        }
    }

    void ISMesh::collapse(const EdgeKey& eid, const NodeKey& nid, double weight) {
        NodeKey nid_remove = (get(eid).node_keys() - nid).front();
        update_collapse(nid, nid_remove, weight);

        auto fids = get(eid).face_keys();
        auto tids = get_tets(eid);

        // Remove edge
        remove(eid);

        // Remove faces
        for(auto f : fids)
        {
            SimplexSet<EdgeKey> eids = get(f).edge_keys();
            remove(f);
            if(get(eids[1]).node_keys().contains(nid))
            {
                eids.swap();
            }
            merge(eids[0], eids[1]);
        }

        // Remove tetrahedra
        for(auto t : tids)
        {
            SimplexSet<FaceKey> fids = get(t).face_keys();
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
        update_flag(get_tets(nid));

#ifdef DEBUG
        validity_check();
#endif
    }

    FaceKey ISMesh::flip_32(const EdgeKey& eid) {
        Edge & e = get(eid);
        SimplexSet<NodeKey> e_nids = e.node_keys();
        SimplexSet<FaceKey> e_fids = e.face_keys();

        assert(e_fids.size() == 3);

        SimplexSet<TetrahedronKey> e_tids = get_tets(e_fids);
        int label = get(e_tids[0]).label();
        assert(e_tids.size() == 3);
#ifdef DEBUG
        assert(label == get(e_tids[1]).label());
        assert(label == get(e_tids[2]).label());
#endif

        // Remove edge
        remove(eid);

        // Create face
        SimplexSet<EdgeKey> f_eids = get_edges(e_tids) - get_edges(e_fids);

        assert(f_eids.size() == 3);

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

            assert(t_fids.size() == 3);

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
        for (auto t : get(new_fid).tet_keys()) {
            set_label(t,label);
        }
#ifdef DEBUG
        validity_check();
#endif

        return new_fid;
    }

    EdgeKey ISMesh::flip_23(const FaceKey& fid) {
        SimplexSet<TetrahedronKey> f_tids = get(fid).tet_keys();

        assert(f_tids.size() == 2);

        SimplexSet<EdgeKey> f_eids = get(fid).edge_keys();

        assert(f_eids.size() == 3);

        SimplexSet<NodeKey> f_nids = get_nodes(f_eids);
        int label = get(f_tids[0]).label();
#ifdef DEBUG
        assert(label == get(f_tids[1]).label());
#endif

        // Create edge
        SimplexSet<NodeKey> e_nids = get_nodes(f_tids) - f_nids;

        assert(e_nids.size() == 2);

        EdgeKey new_eid = insert_edge(e_nids[0], e_nids[1]);

        // Create faces
        SimplexSet<EdgeKey> new_fs_eids = get_edges(f_tids) - f_eids;

        assert(new_fs_eids.size() == 6);

        for (const NodeKey& n : f_nids)
        {
            auto new_f_eids = new_fs_eids & get(n).edge_keys();

            assert(new_f_eids.size() == 2);

            insert_face(new_f_eids[0], new_f_eids[1], new_eid);
        }

        // Remove face
        remove(fid);

        // Create tetrahedra
        SimplexSet<FaceKey> new_ts_fids1 = get_faces(f_tids);
        SimplexSet<FaceKey> new_ts_fids2 = get(new_eid).face_keys();

        assert(new_ts_fids1.size() == 6);
        assert(new_ts_fids2.size() == 3);

        for (const EdgeKey& e : f_eids)
        {
            SimplexSet<FaceKey> new_t_fids = (new_ts_fids1 & get(e).face_keys()) + (new_ts_fids2 & get_faces(get(e).node_keys()));

            assert(new_t_fids.size() == 4);

            insert_tetrahedron(new_t_fids[0], new_t_fids[1], new_t_fids[2], new_t_fids[3]);
        }

        // Remove tetrahedra
        for (auto t : f_tids) {
            remove(t);
        }

        // Update flags

        assert(get_tets(new_eid).size() == 3);

        for (auto t : get_tets(new_eid)) {
            set_label(t,label);
        }

        return new_eid;
    }

    void ISMesh::flip(const EdgeKey& eid, const FaceKey& fid1, const FaceKey& fid2) {
        SimplexSet<FaceKey> fids = {fid1, fid2};
        Edge &e = get(eid);
        SimplexSet<NodeKey> e_nids = e.node_keys();
        SimplexSet<FaceKey> e_fids = e.face_keys();
        SimplexSet<TetrahedronKey> e_tids = get_tets(e_fids);

        // Reconnect edge
        SimplexSet<NodeKey> new_e_nids = get_nodes(fids) - e_nids;

        assert(new_e_nids.size() == 2);

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
            SimplexSet<EdgeKey> rm_eids = get(f).edge_keys() - eid;

            assert(rm_eids.size() == 2);

            SimplexSet<NodeKey> apex = get_nodes(f) - (new_e_nids + e_nids);

            assert(apex.size() == 1);

            SimplexSet<EdgeKey> add_eids = {get_edge(apex[0], new_e_nids[0]), get_edge(apex[0], new_e_nids[1])};

            disconnect(rm_eids[0], f);
            disconnect(rm_eids[1], f);
            connect(add_eids[0], f);
            connect(add_eids[1], f);

            // Reconnect tetrahedra
            SimplexSet<TetrahedronKey> tids = get(f).tet_keys();
            assert(tids.size() == 2);

            SimplexSet<FaceKey> swap_fids = (get_faces(tids) - e_fids) & get_faces(swap_eids);
            assert(swap_fids.size() == 2);

            swap(swap_fids[0], tids[0], swap_fids[1], tids[1]);
        }

        // Update flags
        update_flag(e_tids);

    }

    void ISMesh::flip_22(const FaceKey& fid1, const FaceKey& fid2) {
        SimplexSet<EdgeKey> eid = get(fid1).edge_keys() & get(fid2).edge_keys();
        assert(eid.size() == 1);
#ifdef DEBUG
        assert(get_faces(eid).size() == 3);
        assert(get_tets(eid).size() == 2);
#endif

        flip(eid[0], fid1, fid2);

    }

    void ISMesh::flip_44(const FaceKey& fid1, const FaceKey& fid2) {
        SimplexSet<EdgeKey> eid = get(fid1).edge_keys() & get(fid2).edge_keys();
        assert(eid.size() == 1);
#ifdef DEBUG
        assert(get_faces(eid).size() == 4);
        assert(get_tets(eid).size() == 4);
#endif
        flip(eid[0], fid1, fid2);

    }

    double ISMesh::volume_destination(const SimplexSet<NodeKey>& nids) {
        return Util::volume(get(nids[0]).get_destination(), get(nids[1]).get_destination(), get(nids[2]).get_destination(), get(nids[3]).get_destination());
    }

    double ISMesh::signed_volume_destination(const SimplexSet<NodeKey>& nids) {
        return Util::signed_volume(get(nids[0]).get_destination(), get(nids[1]).get_destination(), get(nids[2]).get_destination(), get(nids[3]).get_destination());
    }

    void ISMesh::garbage_collect() {
        GarbageCollectDeletions res {
            m_node_kernel.garbage_collect(),
            m_edge_kernel.garbage_collect(),
            m_face_kernel.garbage_collect(),
            m_tetrahedron_kernel.garbage_collect()};
        for (auto & fn : m_gc_listeners){
            fn.second(res);
        }
    }

    void ISMesh::scale(const vec3& s) {
        for (auto & nit : nodes()) {
            nit.set_pos(s * nit.get_pos());
            nit.set_destination(s * nit.get_destination());
        }
    }

    void ISMesh::extract_surface_mesh(vector<vec3>& points, vector<int>& faces) {
        garbage_collect();

        map<NodeKey, int> indices;
        // Extract vertices
        for (auto & nit : nodes())
        {
            if (nit.is_interface())
            {
                points.push_back(nit.get_pos());
                indices[nit.key()] = static_cast<int>(points.size());
            }
        }

        // Extract faces
        for (auto & fit : ISMesh::faces())
        {
            if (fit.is_interface())
            {
                for (auto &n : get_sorted_nodes(fit.key()))
                {
                    faces.push_back(indices[n]);
                }
            }
        }
    }

    void ISMesh::extract_surface_mesh_debug(vector<vec3>& points, vector<int>& faces) {
        garbage_collect();

        map<NodeKey, int> indices;
        // Extract vertices
        for (auto & nit : nodes())
        {
            points.push_back(nit.get_pos());
            indices[nit.key()] = static_cast<int>(points.size());
        }

        // Extract faces
        for (auto & fit : ISMesh::faces())
        {
            for (auto &n : get_sorted_nodes(fit.key()))
            {
                faces.push_back(indices[n]);
            }
        }
    }

    void ISMesh::extract_tet_mesh(vector<vec3>& points, vector<int>& tets, vector<int>& tet_labels) {
        garbage_collect();

        map<NodeKey, int> indices;
        // Extract vertices
        for (auto & nit : nodes())
        {
            indices[nit.key()] = static_cast<int>(points.size());
            points.push_back(nit.get_pos());
        }

        for (auto & tit : tetrahedra())
        {
            for (auto n : get_nodes(tit.key()))
            {
                tets.push_back(indices[n]);
            }
            tet_labels.push_back(get(tit.key()).label());
        }
    }

    void ISMesh::validity_check(bool skip_boundary_check) {
        cout << "Testing structure:  ";
        for (auto & f : faces()){
            auto tetsPerFace = f.tet_keys().size();
            assert(tetsPerFace >0 && tetsPerFace <=2);
            auto edgesPerFace = f.edge_keys().size();
            assert(edgesPerFace == 3);
        }
        for (auto & e : edges()){
            auto facesPerEdge = e.face_keys().size();
            assert(facesPerEdge >= 2);
            auto nodesPerEdge = e.node_keys().size();
            assert(nodesPerEdge == 2);
        }

        for (auto & n : nodes()){
            auto edgesPerNode = n.edge_keys().size();
            assert(edgesPerNode >= 3);
        }
        cout << "PASSED" << endl;

        cout << "Testing existence: ";
        for (auto & f : faces()){
            assert(exists(f.key()));
            for (auto t : f.tet_keys()){
                assert(exists(t));
            }
            for (auto e : f.edge_keys()){
                assert(exists(e));
            }
        }

        for (auto & f : edges()){
            assert(exists(f.key()));
            for (auto f : f.face_keys()){
                assert(exists(f));
            }
            for (auto n : f.node_keys()){
                assert(exists(n));
            }
        }

        for (auto & n : nodes()){
            assert(exists(n.key()));
            for (auto e : n.edge_keys()){
                assert(exists(e));
            }
        }
        cout << "PASSED" << endl;

        cout << "Testing connectivity of simplicial complex: ";
        for(auto & tit : tetrahedra())
        {
            assert(exists(tit.key()));
            // Check faces:
            auto faces = tit.face_keys();
            assert(faces.size() == 4);
            for (auto f : faces) {
                assert(exists(f));
                auto cotets = get(f).tet_keys();
                assert((get(f).is_boundary() && cotets.size() == 1) || (!get(f).is_boundary() && cotets.size() == 2));
                assert(find(cotets.begin(), cotets.end(), tit.key()) != cotets.end());
                for (auto f2 : faces) {
                    assert(f == f2 || get_edge(f, f2).is_valid());
                }

                // Check edges:
                auto edges = get(f).edge_keys();
                assert(edges.size() == 3);
                for (auto e : edges)
                {
                    assert(exists(e));
                    auto cofaces = get(e).face_keys();
                    assert(find(cofaces.begin(), cofaces.end(), f) != cofaces.end());
                    for (auto e2 : edges) {
                        assert(e == e2 || get_node(e, e2).is_valid());
                    }

                    // Check nodes:
                    auto nodes = get(e).node_keys();
                    assert(nodes.size() == 2);
                    for (auto n : nodes)
                    {
                        assert(exists(n));
                        auto coedges = get(n).edge_keys();
                        assert(find(coedges.begin(), coedges.end(), e) != coedges.end());
                    }
                }

            }

            assert(get_edges(tit.key()).size() == 6);
            assert(get_nodes(tit.key()).size() == 4);
        }
        cout << "PASSED" << endl;

        cout << "Testing for inverted tetrahedra: ";
        for (auto & tit : tetrahedra())
        {
            assert(!is_inverted(tit.key()));
        }
        cout << "PASSED" << endl;
        if (!skip_boundary_check){
            cout << "Testing for corrupted interface or boundary: ";
            for (auto & eit : edges())
            {
                int boundary = 0;
                int interface = 0;
                for (auto f : eit.face_keys()) {
                    if(get(f).is_boundary())
                    {
                        boundary++;
                    }
                    if(get(f).is_interface())
                    {
                        interface++;
                    }
                }
                assert((eit.is_interface() && interface >= 2) || (!eit.is_interface() && interface == 0)); // Check that the interface is not corrupted
                assert((eit.is_boundary() && boundary == 2) || (!eit.is_boundary() && boundary == 0)); // Check that the boundary is not corrupted
            }
            cout << "PASSED" << endl;
        }
    }

    long ISMesh::add_gc_listener(function<void(const GarbageCollectDeletions&)> fn) {
        static long globalId = 0;
        globalId++;
        m_gc_listeners[globalId] = fn;
        return globalId;
    }

    bool ISMesh::remove_gc_listener(long id) {
        return m_gc_listeners.erase(id)>0;
    }

    long ISMesh::add_label_listener(std::function<void(const TetrahedronKey &tid, unsigned int oldValue)> fn) {
        static long globalId = 0;
        globalId++;
        m_set_label_listeners[globalId] = fn;
        return globalId;
    }

    bool ISMesh::remove_label_listener(long id) {
        return m_set_label_listeners.erase(id) > 0;
    }

    long ISMesh::add_split_listener(std::function<void(const NodeKey &nid_new, const NodeKey &nid1, const NodeKey &nid2)> fn) {
        static long globalId = 0;
        globalId++;
        m_split_listeners[globalId] = fn;
        return globalId;
    }

    bool ISMesh::remove_split_listener(long id) {
        return m_split_listeners.erase(id) > 0;
    }

    long ISMesh::add_collapse_listener(std::function<void(const NodeKey &nid, const NodeKey &nid_removed, double weight)> fn) {
        static long globalId = 0;
        globalId++;
        m_collapse_listeners[globalId] = fn;
        return globalId;
    }

    bool ISMesh::remove_collapse_listener(long id) {
        return m_collapse_listeners.erase(id) > 0;
    }

    unsigned int ISMesh::get_max_node_key() const { return (unsigned int)m_node_kernel.capacity(); }

    unsigned int ISMesh::get_max_edge_key() const { return (unsigned int)m_edge_kernel.capacity(); }

    unsigned int ISMesh::get_max_face_key() const { return (unsigned int)m_face_kernel.capacity(); }

    unsigned int ISMesh::get_max_tet_key() const { return (unsigned int)m_tetrahedron_kernel.capacity(); }

    std::shared_ptr<Geometry> ISMesh::get_subdomain() {
        return subdomain;
    }

    void ISMesh::clear_subdomain() {
        set_subdomain({});
    }

    void ISMesh::set_subdomain(std::shared_ptr<Geometry> subdomain) {
        this->subdomain = subdomain;

        m_node_kernel.revert_excluded();
        m_edge_kernel.revert_excluded();
        m_face_kernel.revert_excluded();
        m_tetrahedron_kernel.revert_excluded();

        if (this->subdomain){
            Query query(this);
            auto m_nodes = query.nodes(subdomain.get());
            auto m_edges = query.edges(m_nodes);
            auto m_faces = query.faces(m_edges);
            auto m_tetrahedra = query.tetrahedra(m_faces);
            query.filter_subset(m_nodes, m_edges, m_faces, m_tetrahedra);


            m_node_kernel.exclude_using_include_set(m_nodes);
            m_edge_kernel.exclude_using_include_set(m_edges);
            m_face_kernel.exclude_using_include_set(m_faces);
            m_tetrahedron_kernel.exclude_using_include_set(m_tetrahedra);
        }
    }

    vec3 ISMesh::get_barycenter(const SimplexSet<NodeKey>& nids, bool interface) {
        vec3 avg_pos(0.);
        int i = 0;
        for (auto n : nids)
        {
            if (!interface || get(n).is_interface())
            {
                avg_pos += get(n).get_pos();
                i++;
            }
        }

        assert(i != 0);

        return avg_pos / static_cast<double>(i);
    }

    int ISMesh::get_number_of_threads() const {
        return m_number_of_threads;
    }

    void ISMesh::set_number_of_threads(int m_number_of_threads) {
        ISMesh::m_number_of_threads = m_number_of_threads;
    }
}
