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

#ifndef OTB_TYPES_H
#define OTB_TYPES_H


//-------------------------------------------------

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/operation_sparse.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>

namespace ublas = boost::numeric::ublas;

#include <OpenTissue/core/math/math_basic_types.h>
#include <OpenTissue/core/math/math_is_number.h>

//-------------------------------------------------

template <typename T>
class OpenTissueBoostTypes
{
public:
    
    typedef T                    real_type;
    
protected:
    
    typedef OpenTissue::math::BasicMathTypes<T,int>  OT;
    
public:
    
    typedef typename OT::vector3_type    vector3_type;
    typedef typename OT::matrix3x3_type  matrix3x3_type;
    
    typedef T vector4_type[4];
    typedef T matrix4x4_type[4][4];
    
public:
    
    typedef ublas::vector<T>                vectorN_type;
    typedef ublas::compressed_matrix<T>     matrixNxN_type;
    
public:
    
    static real_type dot(vector3_type const & v1, vector3_type const & v2)
    {
        return OpenTissue::math::dot(v1,v2);
    }
    
    static vector3_type cross(vector3_type const & v1, vector3_type const & v2)
    {
        return OpenTissue::math::cross(v1, v2);
    }
    
    static real_type length(vector3_type const & v)
    {
        return OpenTissue::math::length(v);
    }
    
    static bool is_nan(double const & t)
    {
        return !(is_number(t));
    }
    
};

// OPENTISSUEBOOST_TYPES_H
#endif
