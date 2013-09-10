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

template<typename MT>
class default_node_traits
{
    typedef typename MT::real_type      T;
    typedef typename MT::vector3_type   V;
    
    V p;
    V p_new;
    std::bitset<3> flags;
    
public:
    
    default_node_traits() : p(0.,0.,0.), p_new(0.,0.,0.) {}
    default_node_traits(T const & x, T const & y, T const & z) : p(x,y,z), p_new(x,y,z) {}
    
    V get_pos()
    {
        return p;
    }
    
    V get_destination()
    {
        return p_new;
    }
    
    void set(default_node_traits const & t)
    {
	    p = t.p;
        p_new = t.p_new;
		flags = t.flags;
    }
    
    void set_pos(V p_)
    {
        p = p_;
    }
    
    void set_destination(V p_)
    {
        p_new = p_;
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

class default_edge_traits
{
    std::bitset<3> flags;
    
public:
    default_edge_traits() {}
    
    
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

class default_face_traits
{
    std::bitset<2> flags;
    
public:
    default_face_traits() {}
        
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

class default_tetrahedron_traits
{
    unsigned int l;
    
public:
    default_tetrahedron_traits() : l(0) {}
    
    int label()
    {
        return l;
    }
    
    void label(unsigned int _label)
    {
        l = _label;
    }

};
