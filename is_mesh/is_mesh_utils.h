#ifndef PROGRESSIVE_MESH_INCIDENCE_SIMPLICIAL_UTILS_H
#define PROGRESSIVE_MESH_INCIDENCE_SIMPLICIAL_UTILS_H

// Disclaimer / copyrights / stuff

#include <vector>

namespace OpenTissue
{
  namespace is_mesh
  {
    namespace util
    {
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Simplex Tags and typebinding traits
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
      //tags
      struct node_simplex_tag{};
      struct edge_simplex_tag{};
      struct face_simplex_tag{};
      struct tetrahedron_simplex_tag{};

      template<typename mesh_type, int dim>
      struct simplex_traits{};

      template<typename mesh_type>
      struct simplex_traits<mesh_type, 0>
      {
        typedef typename mesh_type::node_type             simplex_type;
        typedef typename mesh_type::edge_type             co_boundary_type;
        typedef typename mesh_type::node_iterator         iterator;
        typedef typename mesh_type::edge_iterator         co_boundary_iterator;
        typedef typename mesh_type::node_key_type         key_type;
        typedef          node_simplex_tag                 simplex_tag;
      };

      template<typename mesh_type>
      struct simplex_traits<mesh_type, 1>
      {
        typedef typename mesh_type::edge_type             simplex_type;
        typedef typename mesh_type::node_type             boundary_type;
        typedef typename mesh_type::face_type             co_boundary_type;
        typedef typename mesh_type::edge_iterator         iterator;
        typedef typename mesh_type::node_iterator         boundary_iterator;
        typedef typename mesh_type::face_iterator         co_boundary_iterator;
        typedef typename mesh_type::edge_key_type         key_type;
        typedef          edge_simplex_tag                 simplex_tag;
      };

      template<typename mesh_type>
      struct simplex_traits<mesh_type, 2>
      {
        typedef typename mesh_type::face_type             simplex_type;
        typedef typename mesh_type::edge_type             boundary_type;
        typedef typename mesh_type::tetrahedron_type      co_boundary_type;
        typedef typename mesh_type::face_iterator         iterator;
        typedef typename mesh_type::edge_iterator         boundary_iterator;
        typedef typename mesh_type::tetrahedron_iterator  co_boundary_iterator;
        typedef typename mesh_type::face_key_type         key_type;
        typedef          face_simplex_tag                 simplex_tag;
      };

      template<typename mesh_type>
      struct simplex_traits<mesh_type, 3>
      {
        typedef typename mesh_type::tetrahedron_type      simplex_type;
        typedef typename mesh_type::face_type             boundary_type;
        typedef typename mesh_type::tetrahedron_iterator  iterator;
        typedef typename mesh_type::face_iterator         boundary_iterator;
        typedef typename mesh_type::tetrahedron_key_type  key_type;
        typedef          tetrahedron_simplex_tag          simplex_tag;
      };
    } //namespace util
  } //namespace is_mesh
} //namespace OpenTissue
#endif //PROGRESSIVE_MESH_INCIDENCE_SIMPLICIAL_UTILS_H
