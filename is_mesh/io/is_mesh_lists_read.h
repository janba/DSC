#ifndef PROGRESSIVE_MESH_INCIDENCE_SIMPLICIAL_LISTS_READ_H
#define PROGRESSIVE_MESH_INCIDENCE_SIMPLICIAL_LISTS_READ_H
//
//

//disbaling
// 4996 : strncpy deprecated
#pragma warning(disable: 4996)

#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <algorithm>
#include <utility>
#include <string>

#include <is_mesh/io/is_mesh_xml_read.h>
//#include <OpenTissue/math/is_number.h>

namespace OpenTissue
{
  namespace is_mesh
  {

    /**
    * Read XML Method.
    *
    * @param filename      The path and filename of the tetgen file to be
    *                      read (without extensions).
    * @param mesh          The mesh which the file data is read into.
    *
    * @return              A boolean indicating success or failure.
    */
    template<typename mesh_type>
    bool lists_read(const std::string & vl_filename,
                    const std::string & tl_filename,
                    mesh_type & mesh)
    {
      typedef          util::edge_key           edge_key;
      typedef          util::face_key           face_key;
      typedef          std::map<edge_key, int>       edge_map_type;
      typedef          std::map<face_key, int>       face_map_type;
      typedef typename edge_map_type::iterator  edge_iterator;
      typedef typename face_map_type::iterator  face_iterator;
      typedef typename mesh_type::node_type     node_type;

      //assign fout maps to keep track of created indices
      edge_map_type edge_map;
      face_map_type face_map;

      mesh.clear();

      std::ifstream vl_stream(vl_filename.c_str());
      if (!vl_stream)
      {
        std::cerr << "file not found" << std::endl;
        return false;
      }

      std::ifstream tl_stream(tl_filename.c_str());
      if (!tl_stream)
      {
        std::cerr << "file not found" << std::endl;
        return false;
      }

      int cnt_nodes = 0;
      while (vl_stream.good() && !vl_stream.eof())
      {
        float x, y, z;
        vl_stream >> x >> y >> z;
        mesh.insert_node();
        mesh.find_node(cnt_nodes).set(typename mesh_type::node_traits(x,y,z));
        ++cnt_nodes;
      }

      while (tl_stream.good() && !tl_stream.eof())
      {
        int idx[4];
        tl_stream >> idx[0] >> idx[1] >> idx[2] >> idx[3];
        idx[0] -= 1;
        idx[1] -= 1;
        idx[2] -= 1;
        idx[3] -= 1;
        //std::cout << idx[0] << " " << idx[1] << " " << idx[2] << " " << idx[3] << std::endl;
        //we now have the four indicies for the tetrahedra.. create it from the bottom up.

        //now get or create the 6 edges
        int edges[6];
        edges[0] = util::create_edge(idx[0], idx[1], edge_map, mesh); //edge 01
        edges[1] = util::create_edge(idx[0], idx[2], edge_map, mesh); //edge 02
        edges[2] = util::create_edge(idx[0], idx[3], edge_map, mesh); //edge 03
        edges[3] = util::create_edge(idx[1], idx[2], edge_map, mesh); //edge 12
        edges[4] = util::create_edge(idx[1], idx[3], edge_map, mesh); //edge 13
        edges[5] = util::create_edge(idx[2], idx[3], edge_map, mesh); //edge 23

        int faces[4];
        faces[0] = util::create_face(edges[3], edges[5], edges[4], face_map, mesh); //12-23-31
        faces[1] = util::create_face(edges[1], edges[5], edges[2], face_map, mesh); //02-23-30
        faces[2] = util::create_face(edges[0], edges[4], edges[2], face_map, mesh); //01-13-30
        faces[3] = util::create_face(edges[0], edges[3], edges[1], face_map, mesh); //01-12-20

        mesh.insert_tetrahedron( faces[0], faces[1], faces[2], faces[3] );
      }

      vl_stream.close();
      tl_stream.close();

      //compress mesh
      mesh.compress_all();
      return true;
    }

    template<typename mesh_type>
    bool vectors_read(std::vector<double> & points,
                      std::vector<int>    & tets,
                      mesh_type & mesh)
    {
      typedef          util::edge_key           edge_key;
      typedef          util::face_key           face_key;
      typedef          std::map<edge_key, int>  edge_map_type;
      typedef          std::map<face_key, int>  face_map_type;
      typedef typename edge_map_type::iterator  edge_iterator;
      typedef typename face_map_type::iterator  face_iterator;
      typedef typename mesh_type::node_type     node_type;

      edge_map_type edge_map;
      face_map_type face_map;

      mesh.clear();

      int cnt_nodes = 0;
      for (unsigned int i = 0; i < points.size()/3; ++i)
      {
        float x, y, z;
        x = points[3*i];
        y = points[3*i+1];
        z = points[3*i+2];
        mesh.insert_node();
        mesh.find_node(cnt_nodes).set(typename mesh_type::node_traits(x,y,z));
        ++cnt_nodes;
      }

      for (unsigned int j = 0; j < tets.size()/4; ++j)
      {
        int idx[4];
        idx[0] = tets[4*j];
        idx[1] = tets[4*j+1];
        idx[2] = tets[4*j+2];
        idx[3] = tets[4*j+3];

        int edges[6];
        edges[0] = util::create_edge(idx[0], idx[1], edge_map, mesh); //edge 01
        edges[1] = util::create_edge(idx[0], idx[2], edge_map, mesh); //edge 02
        edges[2] = util::create_edge(idx[0], idx[3], edge_map, mesh); //edge 03
        edges[3] = util::create_edge(idx[1], idx[2], edge_map, mesh); //edge 12
        edges[4] = util::create_edge(idx[1], idx[3], edge_map, mesh); //edge 13
        edges[5] = util::create_edge(idx[2], idx[3], edge_map, mesh); //edge 23

        int faces[4];
        faces[0] = util::create_face(edges[3], edges[5], edges[4], face_map, mesh); //12-23-31
        faces[1] = util::create_face(edges[1], edges[5], edges[2], face_map, mesh); //02-23-30
        faces[2] = util::create_face(edges[0], edges[4], edges[2], face_map, mesh); //01-13-30
        faces[3] = util::create_face(edges[0], edges[3], edges[1], face_map, mesh); //01-12-20

        mesh.insert_tetrahedron( faces[0], faces[1], faces[2], faces[3] );
      }


      mesh.compress_all();
      return true;
    }
  } // namespace is_mesh
} //namespace OpenTissue
#endif //PROGRESSIVE_MESH_INCIDENCE_SIMPLICIAL_XML_READ_H
