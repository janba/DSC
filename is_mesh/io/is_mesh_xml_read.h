#ifndef PROGRESSIVE_MESH_INCIDENCE_SIMPLICIAL_XML_READ_H
#define PROGRESSIVE_MESH_INCIDENCE_SIMPLICIAL_XML_READ_H
//
//

//disbaling
// 4996 : strncpy deprecated
#pragma warning(disable: 4996)

#include <TinyXML/tinyxml.h>

#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <algorithm>
#include <utility>
#include <string>

namespace OpenTissue
{
  namespace is_mesh
  {
    namespace util
    {
      struct edge_key {
        int k1, k2;
        edge_key(int i, int j) : k1(i), k2(j) {}
        bool operator<(const edge_key& k) const
        {
          return k1 < k.k1 || (k1 == k.k1 && k2 < k.k2);
        }
      };
      struct face_key
      {
        int k1, k2, k3;
        face_key(int i, int j, int k) : k1(i), k2(j), k3(k){}
        bool operator<(const face_key& k) const
        {
  //        return k1 < k.k1 || (k1 == k.k1 && k2 < k.k2 || (k2 == k.k2 && k3 < k.k3)) ;
          if (k1 < k.k1) return true;
          if (k1 == k.k1 && k2 < k.k2) return true;
          if (k1 == k.k1 && k2 == k.k2 && k3 < k.k3) return true;
          return false;
        }
      };
      struct tetrahedron_key { int k1, k2, k3, k4; };

      template<typename map_type, typename mesh_type>
      inline int create_edge(int i, int j, map_type& edge_map, mesh_type& mesh)
      {
        int a,b;
        if (i <= j) { a = i; b = j; }
        else { a = j; b = i; }
        edge_key key(a, b);
        typename map_type::iterator it = edge_map.find(key);
        if (it == edge_map.end())
        {
          int n = mesh.insert_edge(i, j); //non-sorted
          it = edge_map.insert(std::pair<edge_key,int>(key, n)).first;
        }
        return it->second;
      }

      template<typename map_type, typename mesh_type>
      inline int create_face(int i, int j, int k, map_type& face_map, mesh_type& mesh)
      {
        int a[3] = {i, j, k};
        std::sort(a, a+3);
        face_key key(a[0], a[1], a[2]);
        typename map_type::iterator it = face_map.find(key); //lookup in sorted order
        if (it == face_map.end())
        {
          int a = mesh.insert_face(i, j, k); //create in supplied order
          it = face_map.insert(std::pair<face_key,int>(key, a)).first;
        }
        return it->second;
      }

    }

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
    bool xml_read(const std::string & filename,mesh_type & mesh)
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

  #ifdef TIXML_USE_STL
      TiXmlDocument xml_document(filename);
  #else
      TiXmlDocument xml_document(filename.c_str());
  #endif

      if(!xml_document.LoadFile())
      {
        std::cerr << "file not found" << std::endl;
        return false;
      }
      TiXmlHandle document_handle( &xml_document );

      TiXmlElement * xml_mesh = document_handle.FirstChild( "T4MESH" ).Element();
      assert(xml_mesh || !"Oh no, could not find a T4MESH tag?");

      int cnt_nodes = 0;
      if(xml_mesh->Attribute("nodes"))
      {
        std::istringstream str_stream(xml_mesh->Attribute("nodes"));
        str_stream >> cnt_nodes;
      }
      assert(cnt_nodes>0 || !"Node count was not positive");

      //create all the nodes - asumes that indices will be set correctly
      for(int i=0;i< cnt_nodes;++i)
        mesh.insert_node();

      TiXmlElement * xml_tetrahedron = document_handle.FirstChild( "T4MESH" ).FirstChild( "TETRAHEDRON" ).Element();
      for( ; xml_tetrahedron; xml_tetrahedron=xml_tetrahedron->NextSiblingElement("TETRAHEDRON") )
      {
        int idx[4];
        for(int i=0; i<4; ++i) idx[i] = -1;

        if(!xml_tetrahedron->Attribute("i"))
        {
          std::cerr << "t4mesh::xml_read(): Missing i index on tetrahedron" << std::endl;
          return false;
        }
        else
        {
          std::istringstream str_stream(xml_tetrahedron->Attribute("i"));
          str_stream >> idx[0];
        }
        if(!xml_tetrahedron->Attribute("j"))
        {
          std::cerr << "t4mesh::xml_read(): Missing j index on tetrahedron" << std::endl;
          return false;
        }
        else
        {
          std::istringstream str_stream(xml_tetrahedron->Attribute("j"));
          str_stream >> idx[1];
        }
        if(!xml_tetrahedron->Attribute("k"))
        {
          std::cerr << "t4mesh::xml_read(): Missing k index on tetrahedron" << std::endl;
          return false;
        }
        else
        {
          std::istringstream str_stream(xml_tetrahedron->Attribute("k"));
          str_stream >> idx[2];
        }
        if(!xml_tetrahedron->Attribute("m"))
        {
          std::cerr << "t4mesh::xml_read(): Missing m index on tetrahedron" << std::endl;
          return false;
        }
        else
        {
          std::istringstream str_stream(xml_tetrahedron->Attribute("m"));
          str_stream >> idx[3];
        }
        if((idx[0]<0 ||idx[1]<0 ||idx[2]<0 ||idx[3]<0)||(idx[0]>=cnt_nodes ||idx[1]>=cnt_nodes ||idx[2]>=cnt_nodes ||idx[3]>=cnt_nodes))
        {
          std::cerr << "t4mesh::xml_read(): Illegal node index on tetrahedron" << std::endl;
          return false;
        }
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
  //* Not loading attributes just yet... make sure everything else works first
      TiXmlElement * xml_point = document_handle.FirstChild( "T4MESH" ).FirstChild( "POINT" ).Element();
      for( ; xml_point; xml_point=xml_point->NextSiblingElement("POINT") )
      {
        int idx;

        if(!xml_point->Attribute("idx"))
        {
          std::cerr << "t4mesh::xml_read(): Missing index on point" << std::endl;
          return false;
        }
        else
        {
          std::istringstream str_stream(xml_point->Attribute("idx"));
          str_stream >> idx;
        }
        if(idx<0 || idx>=cnt_nodes)
        {
          std::cerr << "t4mesh::xml_read(): Illegal index on point" << std::endl;
          return false;
        }

        //vector3_type value;
        float x,y,z;

        if(!xml_point->Attribute("coord"))
        {
          std::cerr << "t4mesh::xml_read(): Missing coord on point" << std::endl;
          return false;
        }
        else
        {
          std::string s(xml_point->Attribute("coord"));
          int pos = static_cast<int>(s.find('['));
          s.erase(pos,1);
          pos = static_cast<int>(s.find(']'));
          s.erase(pos,1);
          while ((pos = static_cast<int>(s.find(','))) >= 0)
              s.replace(pos,1," ");
          std::istringstream str_stream(s);
          str_stream >> x >> y >> z;
          mesh.find_node(idx).set(mesh_type::node_traits(x,y,z));
        }

        /*assert(is_number(value(0)) || !"First coordinate was not a number");
        assert(is_finite(value(0)) || !"First coordinate was not finite");
        assert(is_number(value(1)) || !"Second coordinate was not a number");
        assert(is_finite(value(1)) || !"Second coordinate was not finite");
        assert(is_number(value(2)) || !"Third coordinate was not a number");
        assert(is_finite(value(2)) || !"Third coordinate was not finite");*/

        /*points[idx](0) = value(0);
        points[idx](1) = value(1);
        points[idx](2) = value(2);*/
      }
 // */
      xml_document.Clear();
      //compress mesh
      mesh.compress_all();
      return true;
    }
  } // namespace is_mesh
} //namespace OpenTissue
#endif //PROGRESSIVE_MESH_INCIDENCE_SIMPLICIAL_XML_READ_H
