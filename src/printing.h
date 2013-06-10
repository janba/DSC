//
//  printing.h
//  DSC
//
//  Created by Asger Nyman Christiansen on 6/10/12.
//  Copyright (c) 2012 DTU Informatics. All rights reserved.
//

#ifndef PRINTING_H
#define PRINTING_H

#include <string.h>


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
inline void print_diff(double const & a, double const & b)
{
	if(abs(a - b) > EPSILON)
	{
        std::cout << "Difference:" << std::endl;
		std::cout << "a: " << a << std::endl;
		std::cout << "b: " << b << std::endl;
	}
}

#endif
