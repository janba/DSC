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

#pragma once

#include <bitset>
#include "util.h"

namespace is_mesh {
    
    class NodeAttributes
    {
        vec3 p;
        vec3 p_new;
        std::bitset<3> flags;
        
    public:

        NodeAttributes()
        {
        }

        NodeAttributes(vec3 _p) : p(_p), p_new(_p)
        {
            
        }

        NodeAttributes(const NodeAttributes& other)
                :p{other.p}, p_new{other.p_new}, flags{other.flags}
        {}

        NodeAttributes(NodeAttributes&& other)
        :p{other.p}, p_new{other.p_new}, flags{other.flags}
        {}

        NodeAttributes& operator=(NodeAttributes&& other){
            if (this != &other){
                std::swap(p, other.p);
                std::swap(p_new, other.p_new);
                std::swap(flags, other.flags);
            }
            return *this;
        }
        
        /**
         * Returns the position of the node.
         */
        const vec3& get_pos() const
        {
            return p;
        }
        
        /**
         * Returns the destination of the node.
         */
        const vec3& get_destination() const
        {
            return p_new;
        }
        
        void set_pos(const vec3& p_)
        {
            p = p_;
        }
        
        void set_destination(const vec3& p_)
        {
            p_new = p_;
        }
        
        bool is_crossing() const
        {
            return flags[2];
        }
        
        bool is_boundary() const
        {
            return flags[1];
        }
        
        bool is_interface() const
        {
            return flags[0];
        }
    private:
        void set_crossing(bool b)
        {
            flags[2] = b;
        }
        
        void set_boundary(bool b)
        {
            flags[1] = b;
        }
        
        void set_interface(bool b)
        {
            flags[0] = b;
        }
        friend class ISMesh;
    };
    
    class EdgeAttributes
    {
        std::bitset<3> flags;
        
    public:
        EdgeAttributes() {}

        EdgeAttributes(const EdgeAttributes& other)
                :flags(other.flags)
        {
        }

        EdgeAttributes(EdgeAttributes&& other)
        :flags(other.flags)
        {
        }

        EdgeAttributes& operator=(EdgeAttributes&& other){
            if (this != &other){
                std::swap(flags,other.flags);
            }
            return *this;
        }

        bool is_crossing()
        {
            return flags[2];
        }
        
        bool is_boundary()
        {
            return flags[1];
        }
        
        bool is_interface()
        {
            return flags[0];
        }
    private:
        void set_crossing(bool b)
        {
            flags[2] = b;
        }
        
        void set_boundary(bool b)
        {
            flags[1] = b;
        }
        
        void set_interface(bool b)
        {
            flags[0] = b;
        }
        friend class ISMesh;
    };
    
    class FaceAttributes
    {
        std::bitset<2> flags;
        
    public:
        FaceAttributes() {}

        FaceAttributes(const FaceAttributes& other)
                :flags(other.flags)
        {

        }

        FaceAttributes(FaceAttributes&& other)
        :flags(other.flags)
        {

        }

        FaceAttributes& operator=(FaceAttributes&& other){
            if (this != &other)
            {
                std::swap(flags, other.flags);
            }
            return *this;
        }
        
        bool is_boundary()
        {
            return flags[1];
        }
        
        bool is_interface()
        {
            return flags[0];
        }
    private:
        void set_boundary(bool b)
        {
            flags[1] = b;
        }
        
        void set_interface(bool b)
        {
            flags[0] = b;
        }
        friend class ISMesh;
    };
    
    class TetAttributes
    {
        unsigned int l = 0;
        
    public:
        TetAttributes() {}

        TetAttributes(const TetAttributes& other)
                :l(other.l)
        {

        }


        TetAttributes(TetAttributes&& other)
        :l(other.l)
        {

        }

        TetAttributes& operator=(TetAttributes&& other){
            if (this != &other)
            {
                l = other.l;
            }
            return *this;
        }

        
        int label()
        {
            return l;
        }
        
        void label(unsigned int _label)
        {
            l = _label;
        }
        
    };
    
}
