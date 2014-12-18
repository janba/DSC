#include "simplex.h"

#include <utility>

namespace is_mesh
{

    Node::Node(ISMesh *owner) : Simplex<Key, EdgeKey>(owner) {

    }


    Node::Node(ISMesh *owner,vec3 _p) : p(_p), p_new(_p), Simplex<Key, EdgeKey>(owner) {

    }

    Node::Node(Node&& other) : p (std::move(other.p)), p_new(std::move(other.p_new)),
                               flags(std::move(other.flags)),Simplex<Key, EdgeKey>(std::move(other)) {
    }


    Node &Node::operator=(Node&& other) {
        if (this != &other){
            p = other.p;
            p_new = std::move(other.p_new);
            flags = std::move(other.flags);
            ((Simplex<Key, EdgeKey>*)this)->operator=(std::move(other));
        }
        return *this;
    }

    const SimplexSet<EdgeKey> &Node::edge_keys() const {
        return get_co_boundary();
    }

    Edge::Edge(ISMesh *owner) : Simplex<NodeKey, FaceKey>(owner) {

    }

    Edge::Edge(Edge&& other) : flags(std::move(other.flags)), Simplex<NodeKey, FaceKey>(std::move(other)) {

    }

    Edge &Edge::operator=(Edge&& other) {
        if (this != &other){
            flags = std::move(other.flags);
            ((Simplex<NodeKey, FaceKey>*)this)->operator=(std::move(other));
        }
        return *this;
    }

    const SimplexSet<NodeKey> &Edge::node_keys()  const {
        return get_boundary();
    }

    const SimplexSet<FaceKey> &Edge::face_keys()  const {
        return get_co_boundary();
    }

    Face::Face(ISMesh *owner) : Simplex<EdgeKey, TetrahedronKey>(owner) {

    }

    Face::Face(Face&& other) : flags(std::move(other.flags)), Simplex<EdgeKey, TetrahedronKey>(std::move(other)) {}

    Face &Face::operator=(Face&& other) {
        if (this != &other){
            flags = std::move(other.flags);
            ((Simplex<EdgeKey, TetrahedronKey>*)this)->operator=(std::move(other));
        }
        return *this;
    }

    const SimplexSet<EdgeKey> &Face::edge_keys()  const {
        return get_boundary();
    }

    const SimplexSet<TetrahedronKey> &Face::tet_keys()  const {
        return get_co_boundary();
    }

    Tetrahedron::Tetrahedron(ISMesh *owner) : Simplex<FaceKey, Key>(owner) {

    }

    Tetrahedron::Tetrahedron(Tetrahedron&& other)
            : l(other.l), Simplex<FaceKey, Key>(std::move(other)) {}

    Tetrahedron &Tetrahedron::operator=(Tetrahedron&& other) {
        if (this != &other){
            l = other.l;
            ((Simplex<FaceKey, Key>*)this)->operator=(std::move(other));
        }
        return *this;
    }

    const SimplexSet<FaceKey> &Tetrahedron::face_keys()  const {
        return get_boundary();
    }
}