//
//  Deformabel Simplicial Complex (DSC) method
//  Copyright (C) 2013  Technical University of Denmark
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  See licence.txt for a copy of the GNU General Public License.
//
//#include "draw.h"
//
//#ifdef WIN32
//#include <GL/glew.h>
//#include <GL/glut.h>
//#include <GLGraphics/SOIL.h>
//#else
//#include <GEL/GL/glew.h>
//#include <GLUT/glut.h>
//#include <GEL/GLGraphics/SOIL.h>
//#endif
//
//
//void Painter::save_painting(int width, int height, std::string folder, int time_step)
//{
//    std::ostringstream s;
//    if (folder.length() == 0) {
//        s << "scr";
//    }
//    else {
//        s << folder << "/scr";
//    }
//    
//    if (time_step >= 0)
//    {
//        s << std::string(Util::concat4digits("_", time_step));
//    }
//    s << ".png";
//    glPixelStorei(GL_PACK_ALIGNMENT, 1);
//    int success = SOIL_save_screenshot(s.str().c_str(), SOIL_SAVE_TYPE_PNG, 0, 0, width, height);
//    if(!success)
//    {
//        std::cout << "ERROR: Failed to take screen shot: " << s.str().c_str() << std::endl;
//        return;
//    }
//}
//
//
//void Painter::begin()
//{
//    glClearColor(BACKGROUND_COLOR[0], BACKGROUND_COLOR[1], BACKGROUND_COLOR[2],0);
//    glClear(GL_COLOR_BUFFER_BIT);
//}
//
//void Painter::end()
//{
//    glFinish();
//    glutSwapBuffers();
//}
//
//
//void Painter::draw_complex(const SimplicialComplex *complex)
//{
//    draw_domain(complex->get_design_domain());
//    std::vector<Domain*> objects = complex->get_object_domains();
////    for (auto obj = objects.begin(); obj != objects.end(); obj++)
////    {
////        draw_domain(*obj, DARK_GRAY);
////    }
//    draw_faces(complex);
//    draw_edges(complex);
//    draw_vertices(complex);
//}
//
//void Painter::draw_domain(const Domain *domain, CGLA::Vec3d color)
//{
//    std::vector<CGLA::Vec2d> corners = domain->get_corners();
//    glColor3dv(&color[0]);
//    CGLA::Vec3d cor;
//    glBegin(GL_POLYGON);
//    for (int j = 0; j < corners.size(); j++)
//    {
//        cor = CGLA::Vec3d(corners[j][0], corners[j][1], 0.);
//        glVertex3dv(&cor[0]);
//    }
//    glEnd();
//}
//
//void Painter::draw_vertices(const SimplicialComplex *complex)
//{
//    HMesh::VertexAttributeVector<CGLA::Vec3d> colors = complex->get_vertex_colors();
//    glPointSize(std::max(std::floor(POINT_SIZE*complex->get_avg_edge_length()), 1.));
//	glBegin(GL_POINTS);
//    CGLA::Vec3d p;
//	for(HMesh::VertexIDIterator vi = complex->vertices_begin(); vi != complex->vertices_end(); ++vi)
//    {
//        p = CGLA::Vec3d(complex->get_pos(*vi)[0], complex->get_pos(*vi)[1], 0.);
//        glColor3dv(&colors[*vi][0]);
//        glVertex3dv(&p[0]);
//    }
//	glEnd();
//}
//
//void Painter::draw_interface(const SimplicialComplex *complex, CGLA::Vec3d color)
//{
//    glPointSize(std::max(std::floor(POINT_SIZE*complex->get_avg_edge_length()), 1.));
//	glBegin(GL_POINTS);
//    CGLA::Vec3d p;
//    glColor3dv(&color[0]);
//	for(HMesh::VertexIDIterator vi = complex->vertices_begin(); vi != complex->vertices_end(); ++vi)
//    {
//        if (complex->is_movable(*vi)) {
//            CGLA::Vec3d temp(complex->get_pos_new(*vi)[0], complex->get_pos_new(*vi)[1], 0.);
//            glVertex3dv(&temp[0]);
//        }
//    }
//	glEnd();
//    glLineWidth(std::max(std::floor(LINE_WIDTH*complex->get_avg_edge_length()), 1.));
//    CGLA::Vec3d p1, p2;
//	glBegin(GL_LINES);
//    for(HMesh::HalfEdgeIDIterator hei = complex->halfedges_begin(); hei != complex->halfedges_end(); ++hei)
//    {
//        HMesh::Walker hew = complex->walker(*hei);
//        if (complex->is_movable(hew.halfedge()) && (complex->is_movable(hew.vertex()) || complex->is_movable(hew.opp().vertex())))
//        {
//            p1 = CGLA::Vec3d(complex->get_pos_new(hew.vertex())[0], complex->get_pos_new(hew.vertex())[1], 0.);
//            p2 = CGLA::Vec3d(complex->get_pos_new(hew.opp().vertex())[0], complex->get_pos_new(hew.opp().vertex())[1], 0.);
//            glVertex3dv(&p1[0]);
//            glVertex3dv(&p2[0]);
//        }
//    }
//	glEnd();
//}
//
//void Painter::draw_arrows(const SimplicialComplex *complex, const HMesh::VertexAttributeVector<CGLA::Vec2d> &arrows, CGLA::Vec3d color)
//{
//    glColor3dv(&color[0]);
//    glLineWidth(std::max(std::floor(LINE_WIDTH*complex->get_avg_edge_length()), 1.));
//    CGLA::Vec3d arrow, a_hat, p;
//    for(HMesh::VertexIDIterator vi = complex->vertices_begin(); vi != complex->vertices_end(); ++vi)
//    {
//        arrow = CGLA::Vec3d(arrows[*vi][0], arrows[*vi][1], 0.f);
//        if(arrow.length() > EPSILON)
//        {
//            a_hat = CGLA::Vec3d(-arrow[1], arrow[0], 0.f);
//            p = CGLA::Vec3d(complex->get_pos(*vi)[0], complex->get_pos(*vi)[1], 0.);
//#ifdef DEBUG
//            if (complex->is_movable(*vi)) {
//                p = CGLA::Vec3d(complex->get_pos_new(*vi)[0], complex->get_pos_new(*vi)[1], 0.);
//            }
//#endif
//            glBegin(GL_LINES);
//            glVertex3dv(&p[0]);
//            glVertex3dv(&(p + 0.7*arrow)[0]);
//            glEnd();
//            
//            glBegin(GL_POLYGON);
//            glVertex3dv(&(p + arrow)[0]);
//            glVertex3dv(&(p + 0.6*arrow + 0.13*a_hat)[0]);
//            glVertex3dv(&(p + 0.6*arrow - 0.13*a_hat)[0]);
//            glEnd();
//        }
//    }
//}
//
//
//void Painter::draw_lines(const SimplicialComplex *complex, const HMesh::VertexAttributeVector<CGLA::Vec2d> &lines, CGLA::Vec3d color)
//{
//    glColor3dv(&color[0]);
//    glLineWidth(std::max(std::floor(LINE_WIDTH*complex->get_avg_edge_length()), 1.));
//    CGLA::Vec3d line, p;
//    for(HMesh::VertexIDIterator vi = complex->vertices_begin(); vi != complex->vertices_end(); ++vi)
//    {
//        line = CGLA::Vec3d(lines[*vi][0], lines[*vi][1], 0.f);
//        if(line.length() > EPSILON)
//        {
//            p = CGLA::Vec3d(complex->get_pos(*vi)[0], complex->get_pos(*vi)[1], 0.);
//            
//            glBegin(GL_LINES);
//            glVertex3dv(&p[0]);
//            glVertex3dv(&(p + line)[0]);
//            glEnd();
//        }
//    }
//}
//
//void Painter::draw_edges(const SimplicialComplex *complex)
//{
//    HMesh::HalfEdgeAttributeVector<CGLA::Vec3d> colors = complex->get_edge_colors();
//    glLineWidth(std::max(std::floor(LINE_WIDTH*complex->get_avg_edge_length()), 1.));
//    CGLA::Vec3d p1, p2;
//	glBegin(GL_LINES);
//	for(HMesh::HalfEdgeIDIterator hei = complex->halfedges_begin(); hei != complex->halfedges_end(); ++hei)
//    {
//        glColor3dv(&colors[*hei][0]);
//        
//        HMesh::Walker hew = complex->walker(*hei);
//        p1 = CGLA::Vec3d(complex->get_pos(hew.vertex())[0], complex->get_pos(hew.vertex())[1], 0.);
//        p2 = CGLA::Vec3d(complex->get_pos(hew.opp().vertex())[0], complex->get_pos(hew.opp().vertex())[1], 0.);
//        glVertex3dv(&p1[0]);
//        glVertex3dv(&p2[0]);
//    }
//	glEnd();
//}
//
//void Painter::draw_faces(const SimplicialComplex *complex)
//{
//    HMesh::FaceAttributeVector<CGLA::Vec3d> colors = complex->get_face_colors();
//    draw_faces(complex, colors);
//}
//
//void Painter::draw_faces(const SimplicialComplex *complex, const HMesh::FaceAttributeVector<CGLA::Vec3d> &colors)
//{
//	for(HMesh::FaceIDIterator fi = complex->faces_begin(); fi != complex->faces_end(); ++fi)
//    {
//        if(colors[*fi] != INVISIBLE)
//        {
//            glColor3dv(&colors[*fi][0]);
//            glBegin(GL_POLYGON);
//            for (HMesh::Walker hew = complex->walker(*fi); !hew.full_circle(); hew = hew.circulate_face_cw())
//            {
//                glVertex3dv(&complex->get_pos(hew.vertex())[0]);
//            }
//            glEnd();
//        }
//    }
//}
//
