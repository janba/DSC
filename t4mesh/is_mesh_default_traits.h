#ifndef PROGRESSIVE_MESH_INCIDENCE_MESH_TRAITS_H
#define PROGRESSIVE_MESH_INCIDENCE_MESH_TRAITS_H

namespace OpenTissue
{
  namespace is_mesh
  {
    struct default_node_traits
    {
      float x,y,z;
      default_node_traits() {}
      default_node_traits(float x, float y, float z) : x(x), y(y), z(z) {}
      default_node_traits(const default_node_traits& t) { x = t.x; y = t.y; z = t.z; }
      void set (const default_node_traits& t) { x = t.x; y = t.y; z = t.z; }
      void set_coords(float const & i, float const & j, float const & k)
      {
        x = i; y = j; z = k;
      }
    };

    struct default_edge_traits {
      default_edge_traits(){}
      default_edge_traits(const default_edge_traits&) {}
      void set(default_edge_traits){};
    };

    struct default_face_traits {
      default_face_traits(){}
      default_face_traits(const default_face_traits&){}
      void set(default_face_traits){};
    };

    struct default_tetrahedron_traits {
      default_tetrahedron_traits(){}
      default_tetrahedron_traits(const default_tetrahedron_traits&){}
      void set(default_tetrahedron_traits){};
    };
  }
} //namespace OpenTissue
#endif //PROGRESSIVE_MESH_INCIDENCE_MESH_TRAITS_H
