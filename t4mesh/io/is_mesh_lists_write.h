#ifndef OPENTISSUE_T4MESH_IO_T4MESH_LISTS_WRITE_H
#define OPENTISSUE_T4MESH_IO_T4MESH_LISTS_WRITE_H
//
// OpenTissue, A toolbox for physical based simulation and animation.
// Copyright (C) 2007 Department of Computer Science, University of Copenhagen
//
#include <is_mesh/configuration.h>

#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>

//THIS FILE DONT WORK!!//

namespace OpenTissue
{
  namespace is_mesh
  {
    namespace util 
    {

      template<typename mesh_type>
      bool make_lists(std::string const & vl_file,
                      std::string const & tl_file,
                      mesh_type & mesh)
      {
        typedef typename mesh_type::node_iterator node_iterator;
        typedef typename mesh_type::node_key_type node_key_type;
        typedef typename mesh_type::node_iterator const_node_iterator;
        typedef typename mesh_type::tetrahedron_iterator tetrahedron_iterator;
        typedef typename mesh_type::tetrahedron_iterator const_tetrahedron_iterator;
        typedef typename mesh_type::simplex_set_type simplex_set_type;
        typedef typename simplex_set_type::node_set_iterator node_set_iterator;

        typedef typename mesh_type::simplex_set_type simplex_set_type;

        std::ofstream vl_stream(vl_file.c_str());
        std::ofstream tl_stream(tl_file.c_str());

        //the nodes might not be numbered from 1 to n, so label them all
        int i = 1;
        const_node_iterator ne = mesh.nodes_end();
        for(const_node_iterator nodes = mesh.nodes_begin(); nodes != ne; ++nodes) 
        {
          nodes->set_label(i);
          ++i;
        }

        const_tetrahedron_iterator te = mesh.tetrahedra_end();
        for(const_tetrahedron_iterator tetrahedron = mesh.tetrahedra_begin(); tetrahedron != te; ++tetrahedron)
        {
          //get data
          std::vector<node_key_type> verts(4);
          mesh.vertices(tetrahedron.key(), verts);

          tl_stream << mesh.find_node(verts[1]).get_label() << " "
                    << mesh.find_node(verts[0]).get_label() << " "
                    << mesh.find_node(verts[2]).get_label() << " "
                    << mesh.find_node(verts[3]).get_label() << std::endl;
        }

        ne = mesh.nodes_end();
        for(const_node_iterator node_itr = mesh.nodes_begin(); node_itr != ne; ++node_itr) 
        {
          vl_stream << node_itr->v[0] << " " << node_itr->v[1] << " " << node_itr->v[2] << std::endl;
        }

        ne = mesh.nodes_end();
        for(const_node_iterator nodes = mesh.nodes_begin(); nodes != ne; ++nodes) 
        {
          nodes->reset_label();
        }

        return true;
      } 
    } //namespace util

  } // namespace t4mesh
} // namespace OpenTissue

//OPENTISSUE_T4MESH_IO_T4MESH_XML_WRITE_H
#endif
