#include "simplex.h"

namespace is_mesh
{

    Node::Node(ISMesh *owner) : Simplex<Key, EdgeKey>(owner) {

    }


    Node::Node(ISMesh *owner,const NodeAttributes & t) : NodeAttributes(t), Simplex<Key, EdgeKey>(owner) {

    }

    Node::Node(Node&& other) : NodeAttributes(std::move(other)), Simplex<Key, EdgeKey>(std::move(other)) {}


    Node &Node::operator=(Node&& other) {
        if (this != &other){
            ((NodeAttributes*)this)->operator=(std::move(other));
            ((Simplex<Key, EdgeKey>*)this)->operator=(std::move(other));
        }
        return *this;
    }

    const SimplexSet<EdgeKey> &Node::edge_keys() const {
        return get_co_boundary();
    }

    Edge::Edge(ISMesh *owner) : Simplex<NodeKey, FaceKey>(owner) {

    }

    Edge::Edge(ISMesh *owner,const EdgeAttributes & t) : EdgeAttributes(t), Simplex<NodeKey, FaceKey>(owner) {

    }

    Edge::Edge(Edge&& other) : EdgeAttributes(std::move(other)), Simplex<NodeKey, FaceKey>(std::move(other)) {

    }

    Edge &Edge::operator=(Edge&& other) {
        if (this != &other){
            ((EdgeAttributes*)this)->operator=(std::move(other));
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

    Face::Face(ISMesh *owner, const FaceAttributes & t) : FaceAttributes(t), Simplex<EdgeKey, TetrahedronKey>(owner) {

    }

    Face::Face(Face&& other) : FaceAttributes(std::move(other)), Simplex<EdgeKey, TetrahedronKey>(std::move(other)) {}

    Face &Face::operator=(Face&& other) {
        if (this != &other){
            ((FaceAttributes*)this)->operator=(std::move(other));
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

    Tetrahedron::Tetrahedron(ISMesh *owner,const TetAttributes & t) : TetAttributes(t), Simplex<FaceKey, Key>(owner) {

    }

    Tetrahedron::Tetrahedron(Tetrahedron&& other)
            : TetAttributes(std::move(other)), Simplex<FaceKey, Key>(std::move(other)) {}

    Tetrahedron &Tetrahedron::operator=(Tetrahedron&& other) {
        if (this != &other){
            ((TetAttributes*)this)->operator=(std::move(other));
            ((Simplex<FaceKey, Key>*)this)->operator=(std::move(other));
        }
        return *this;
    }

    const SimplexSet<FaceKey> &Tetrahedron::face_keys()  const {
        return get_boundary();
    }
}