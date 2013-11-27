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
        
        NodeAttributes(vec3 _p) : p(_p), p_new(_p)
        {
            
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
    };
    
    class EdgeAttributes
    {
        std::bitset<3> flags;
        
    public:
        EdgeAttributes() {}
        
        
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
    };
    
    class FaceAttributes
    {
        std::bitset<2> flags;
        
    public:
        FaceAttributes() {}
        
        bool is_boundary()
        {
            return flags[1];
        }
        
        bool is_interface()
        {
            return flags[0];
        }
        
        void set_boundary(bool b)
        {
            flags[1] = b;
        }
        
        void set_interface(bool b)
        {
            flags[0] = b;
        }
    };
    
    class TetAttributes
    {
        unsigned int l = 0;
        
    public:
        TetAttributes() {}
        
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