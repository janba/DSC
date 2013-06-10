#ifndef GEL_TYPES_H
#define GEL_TYPES_H

/** 2013 Mark Viinblad Jensen **/

#include <CGLA/Vec3d.h>
#include <CGLA/Vec4d.h>
#include <CGLA/Mat4x4d.h>
#include <CGLA/Mat3x3d.h>
#include <CGLA/Mat2x3d.h>
#include <CGLA/eigensolution.h>

#include <LinAlg/Matrix.h>
#include <LinAlg/LapackFunc.h>

class GELTypes //GELTypes
{
public:
    
    typedef double            real_type;
    
public:
    
    class CMatrixWrapper : public LinAlg::CMatrix
    {
    public:
        
        typedef LinAlg::CMatrix base_type;
        
    public:
        
        CMatrixWrapper(size_t const & M, size_t const & N)
        : base_type(M,N)
        {}
        
        CMatrixWrapper(base_type const & A)
        : base_type(A)
        {}
        
    public:
        
        real_type const & operator()(int const & i, int const & j) const
        {
            return base_type::get(i, j);
        }
        
        real_type & operator()(int const & i, int const & j)
        {
            return (*this)[i][j];
        }
        
        
        CMatrixWrapper & operator=(CMatrixWrapper const & rhs)
        {
            if (&rhs != this)
            {
                base_type * base = this;
                base->operator=(rhs);
                //::base_type::operator=
            }
            return *this;
        }
        
    };
    
    static CMatrixWrapper transpose(CMatrixWrapper const & A)
    {
        return CMatrixWrapper( A.Transposed() );
    }
    
    class CVectorWrapper : public LinAlg::CVector
    {
    public:
        
        typedef LinAlg::CVector base_type;
        
        
    public:
        
        CVectorWrapper(size_t const & N)
        : base_type(N)
        {}
        
        CVectorWrapper(base_type const & v)
        : base_type(v)
        {}
        
    public:
        
        real_type const & operator()(int const & i) const
        {
            return base_type::get(i);
        }
        
        real_type & operator()(int const & i)
        {
            return (*this)[i];
        }
        
    };
    
public:
    
    typedef CGLA::Vec3d       vector3_type;
    typedef CGLA::Vec4d       vector4_type;
    typedef CVectorWrapper    vectorN_type;
    typedef CGLA::Mat4x4d     matrix4x4_type;
    typedef CGLA::Mat3x3d     matrix3x3_type;
    typedef CGLA::Mat2x3d     matrix2x3_type;
    typedef CMatrixWrapper    matrixMxN_type;
    typedef CGLA::Axis        axis_type;
    
public:
    
    static real_type dot(vector3_type const & v1, vector3_type const & v2)
    {
        return CGLA::dot(v1, v2);
    }
    
    static vector3_type cross(vector3_type const & v1, vector3_type const & v2)
    {
        return CGLA::cross(v1, v2);
    }
    
    static real_type length(vector3_type const & v)
    {
        return CGLA::length(v);
    }
    
    static bool is_nan(real_type const & t)
    {
        return CGLA::isnan(t);
    }
    
    static real_type determinant(matrix3x3_type const & m)
    {
        return CGLA::determinant(m);
    }
    
    static real_type determinant(matrix4x4_type const & m)
    {
        return CGLA::determinant(m);
    }
    
    
    static matrix4x4_type transpose(matrix4x4_type const & m)
    {
        return CGLA::transpose(m);
    }
    
    static matrix4x4_type invert(matrix4x4_type const & m)
    {
        return CGLA::invert(m);
    }
    
    
    static real_type sqr_length(vector3_type const & v)
    {
        return CGLA::sqr_length(v);
    }
    
    static vector3_type normalize(vector3_type const & v)
    {
        return CGLA::normalize(v);
    }
    
    static void orthogonal(vector3_type const & v1, vector3_type & v2, vector3_type & v3)
    {
        CGLA::orthogonal(v1, v2, v3);
    }
    
    static void eigen(matrix3x3_type AA, matrix3x3_type Q, matrix3x3_type L)
    {
        CGLA::power_eigensolution(AA, Q, L);
    }
    
    static matrix3x3_type get_rotation_matrix(axis_type const & axis, real_type const & angle)
    {
        return CGLA::rotation_Mat3x3d(axis, angle);
    }
    
    static axis_type get_z_axis()
    {
        return CGLA::ZAXIS;
    }
    
    static vector3_type mul(matrix3x3_type const & m, vector3_type const & v)
    {
        return m*v;
    }
    
    static vectorN_type solve(matrixMxN_type const & A, vectorN_type const & b)
    {
        vectorN_type x = LinAlg::LinearLSSolve(A,b);
        return x;
    }
    
};

// GEL_TYPES_H
#endif
