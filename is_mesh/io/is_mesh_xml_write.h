#ifndef OPENTISSUE_T4MESH_IO_T4MESH_XML_WRITE_H
#define OPENTISSUE_T4MESH_IO_T4MESH_XML_WRITE_H
//
// OpenTissue, A toolbox for physical based simulation and animation.
// Copyright (C) 2007 Department of Computer Science, University of Copenhagen
//
#include <is_mesh/configuration.h>

#include <TinyXML/tinyxml.h>

#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>

namespace OpenTissue
{
  namespace is_mesh
  {
    namespace util 
    {

      /**
       * XML Document creation routine. Useful lower­level routine for when
       * you want to keep a reference to the created doc.
       *
       * @param doc          Reference to XML doc to be written.
       * @param mesh         Mesh to be serialized.
       * @param points       Node coordinates will be read from this container.
       *
       * @return             If succesfully written to the specified file then the
       *                     return value is true otherwise it is false.
       *
       * @example
       * @code
       *
       * template<typename skin_type>
       * void skin_xml_embed_extra(skin_type const & skin, TiXmlDocument & doc)
       * {
       *   // iterate over 'POINT' nodes and add normals
       *   TiXmlElement * point = TiXmlHandle(&doc).FirstChild("T4MESH").FirstChild("POINT").Element();
       *   for( int i = 0; point; point=point->NextSiblingElement("POINT"), ++i )
       *   {
       *     skin_type::const_node_iterator n = skin.const_node(i);
       * 
       *     std::stringstream normal;
       *     normal << n->m_original_normal;
       * 
       *     point->SetAttribute("normal", normal.str());
       *   }
       * }
       * 
       * // write the tetrahedral mesh
       * mesh_type skin;
       * TiXmlDocument doc;
       * t4mesh::xml_make_doc(doc, skin, default_point_container<mesh_type>(&skin));
       * skin_xml_embed_extra(skin, doc);
       * doc.SaveFile("skin.xml");
       *
       * @endcode
       */
      template<typename mesh_type>
      bool xml_make_doc(TiXmlDocument & doc, mesh_type & mesh)
      {
        typedef typename mesh_type::node_iterator node_iterator;
        typedef typename mesh_type::node_key_type node_key_type;
        typedef typename mesh_type::node_iterator const_node_iterator;
        typedef typename mesh_type::tetrahedron_iterator tetrahedron_iterator;
        typedef typename mesh_type::tetrahedron_iterator const_tetrahedron_iterator;
        typedef typename mesh_type::simplex_set_type simplex_set_type;
        typedef typename simplex_set_type::node_set_iterator node_set_iterator;

        typedef typename mesh_type::simplex_set_type simplex_set_type;

        TiXmlDeclaration * decl = new TiXmlDeclaration( "1.0", "", "" );
        TiXmlElement * meshelem = new TiXmlElement( "T4MESH" );
        meshelem->SetAttribute( "nodes", static_cast<unsigned int>(mesh.size_nodes()) );

        doc.LinkEndChild(decl);
        doc.LinkEndChild(meshelem);

        //the nodes might not be numbered from 0 to n, so label them all
        int i = 0;
        const_node_iterator ne = mesh.nodes_end();
        for(const_node_iterator nodes = mesh.nodes_begin(); nodes != ne; ++nodes) 
        {
          nodes->set_label(i);
          ++i;
        }

        const_tetrahedron_iterator te = mesh.tetrahedra_end();
        for(const_tetrahedron_iterator tetrahedron=mesh.tetrahedra_begin();tetrahedron!=te;++tetrahedron)
        {
          //get data
          std::vector<node_key_type> verts(4);
          mesh.vertices(tetrahedron.key(), verts);

          TiXmlElement * elem = new TiXmlElement( "TETRAHEDRON" );

          elem->SetAttribute( "i", mesh.find_node(verts[0]).get_label() );
          elem->SetAttribute( "j", mesh.find_node(verts[1]).get_label() );
          elem->SetAttribute( "k", mesh.find_node(verts[2]).get_label() );
          elem->SetAttribute( "m", mesh.find_node(verts[3]).get_label() );

          meshelem->LinkEndChild( elem );
        }

        ne = mesh.nodes_end();
        for(const_node_iterator node_itr = mesh.nodes_begin(); node_itr != ne; ++node_itr) 
        {
          TiXmlElement * elem = new TiXmlElement( "POINT" );

          std::stringstream s;
          s.precision(15);
          s << "[" << node_itr->v[0] << "," << node_itr->v[1] << "," << node_itr->v[2] << "]";
          elem->SetAttribute( "idx", node_itr->get_label() );
          elem->SetAttribute( "coord", s.str().c_str() );

          meshelem->LinkEndChild( elem );
        }


        //reset the labels
        ne = mesh.nodes_end();
        for(const_node_iterator nodes = mesh.nodes_begin(); nodes != ne; ++nodes) 
        {
          nodes->reset_label();
        }

        return true;
      } 
    } //namespace util

      /**
       * Write t4mesh to XML file. 
       *
       * @param filename     File to write to.
       * @param mesh         Mesh to be serialized.
       * @param points       Node coordinates will be read from this container.
       *
       * @return             If succesfully written to the specified file then the
       *                     return value is true otherwise it is false.
       */
      template<typename mesh_type>
      bool xml_write(std::string const & filename, mesh_type const & mesh)
      {
        // build document
        TiXmlDocument doc;
        if (!OpenTissue::is_mesh::util::xml_make_doc(doc, mesh))
          return false;

        // write the document
        #ifdef TIXML_USE_STL
        doc.SaveFile(filename);
        #else
        doc.SaveFile(filename.c_str());
        #endif

        return true;
      }
  } // namespace t4mesh
} // namespace OpenTissue

//OPENTISSUE_T4MESH_IO_T4MESH_XML_WRITE_H
#endif
