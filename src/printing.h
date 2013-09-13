//
//  printing.h
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

#include <string.h>

namespace DSC {
    
    inline void print_out(std::string text)
    {
        std::cout << text << std::endl;
    }
    
    template <class VT>
    inline void print(VT const & a)
    {
        std::cout << "x " << a[0] << " y " << a[1] << " z " << a[2] << std::endl;
    }
    
    template <class VT>
    inline void print(std::vector<VT> const & a)
    {
        for(unsigned int i = 0; i < a.size(); i++)
            print(a[i]);
    }
    
    template <class VT>
    inline void print_diff(std::vector<VT> const & a, std::vector<VT> const & b)
    {
        if(length(a - b) > EPSILON)
        {
            std::cout << "Difference:" << std::endl;
            print(a);
            print(b);
        }
    }
    
    template <class VT>
    inline void print_diff(real const & a, real const & b)
    {
        if(abs(a - b) > EPSILON)
        {
            std::cout << "Difference:" << std::endl;
            std::cout << "a: " << a << std::endl;
            std::cout << "b: " << b << std::endl;
        }
    }
    
}
