/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * @brief Abstract vector class
 * ----------------------------------------------------------------------- */

/** @file ArithVec.h
 * @brief Abstract vector class
 */

#ifndef __CGLA_ARITHVEC_H__
#define __CGLA_ARITHVEC_H__

#include <iostream>
#include "CGLA.h"

#include <numeric>

namespace CGLA 
{

  /** \brief Template representing generic arithmetic vectors.  

			The three parameters to the template are

      T - the scalar type (i.e. float, int, double etc.)

      V - the name of the vector type. This template is always (and
      only) used as ancestor of concrete types, and the name of the
      class _inheriting_ _from_ this class is used as the V argument.

      N - The final argument is the dimension N. For instance, N=3 for a
      3D vector.

      This class template contains all functions that are assumed to be
      the same for any arithmetic vector - regardless of dimension or
      the type of scalars used for coordinates.

      The template contains no virtual functions which is important
      since they add overhead.
  */

  template <class T, class V, unsigned int N> 
    class ArithVec
    {

    protected:

      /// The actual contents of the vector.
      T data[N];

    protected:

      //----------------------------------------------------------------------
      // Constructors
      //----------------------------------------------------------------------

      /// Construct uninitialized vector
      ArithVec() 
	{
	}

      /// Construct a vector where all coordinates are identical
      explicit ArithVec(T _a)
	{
	  std::fill_n(data, N, _a);
	}

      /// Construct a 2D vector 
      ArithVec(T _a, T _b)
	{
	  assert(N==2);
	  data[0] = _a;
	  data[1] = _b;
	}

      /// Construct a 3D vector
      ArithVec(T _a, T _b, T _c)
	{
	  assert(N==3);
	  data[0] = _a;
	  data[1] = _b;
	  data[2] = _c;
	}

      /// Construct a 4D vector
      ArithVec(T _a, T _b, T _c, T _d)
	{
	  assert(N==4);
	  data[0] = _a;
	  data[1] = _b;
	  data[2] = _c;
	  data[3] = _d;
	}


    public:

      /// For convenience we define a more meaningful name for the scalar type
      typedef T ScalarType;
	
      /// A more meaningful name for vector type
      typedef V VectorType;

      /// Return dimension of vector
      static unsigned int get_dim() {return N;}
	
      /// Set all coordinates of a 2D vector.
      void set(T _a, T _b)
	{
	  assert(N==2);
	  data[0] = _a;
	  data[1] = _b;
	}

      /// Set all coordinates of a 3D vector.
      void set(T _a, T _b, T _c)
	{
	  assert(N==3);
	  data[0] = _a;
	  data[1] = _b;
	  data[2] = _c;
	}

      /// Set all coordinates of a 4D vector.
      void set(T _a, T _b, T _c, T _d)
	{
	  assert(N==4);
	  data[0] = _a;
	  data[1] = _b;
	  data[2] = _c;
	  data[3] = _d;
	}

  /// Const index operator
  const T& operator [] ( unsigned int i ) const
	{
	  assert(i<N);
	  return data[i];
	}

  /// Non-const index operator
  T& operator [] ( unsigned int i ) 
	{
	  assert(i<N);
	  return data[i];
	}

	/// Const index operator
  const T& operator () ( unsigned int i ) const
	{
	  assert(i<N);
	  return data[i];
	}

  /// Non-const index operator
  T& operator () ( unsigned int i ) 
	{
	  assert(i<N);
	  return data[i];
	}


      /** Get a pointer to first element in data array.
	  This function may be useful when interfacing with some other API 
	  such as OpenGL (TM) */
      T* get() {return &data[0];}

      /** Get a const pointer to first element in data array.
	  This function may be useful when interfacing with some other API 
	  such as OpenGL (TM). */
      const T* get() const {return &data[0];}

      //----------------------------------------------------------------------
      // Comparison operators
      //----------------------------------------------------------------------

      /// Equality operator
      bool operator==(const V& v) const 
	{
	  return std::inner_product(data, &data[N], v.get(), true,
				    std::logical_and<bool>(), std::equal_to<T>());
	}

#define for_all_i(expr) for(unsigned int i=0;i<N;i++) {expr;}
      /// Equality wrt scalar. True if all coords are equal to scalar
      bool operator==(T k) const 
	{ 
	  for_all_i(if (data[i] != k) return false)
	    return true;
	}
#undef for_all_i  

      /// Inequality operator
      bool operator!=(const V& v) const 
	{
	  return std::inner_product(data, &data[N], v.get(), false,
				    std::logical_or<bool>(), std::not_equal_to<T>());
	}

      /// Inequality wrt scalar. True if any coord not equal to scalar
      bool operator!=(T k) const 
	{ 
	  return !(*this==k);
	}


      //----------------------------------------------------------------------
      // Comparison functions ... of geometric significance 
      //----------------------------------------------------------------------

      /** Compare all coordinates against other vector. ( < )
	  Similar to testing whether we are on one side of three planes. */
      bool  all_l  (const V& v) const
	{
	  return std::inner_product(data, &data[N], v.get(), true, 
				    std::logical_and<bool>(), std::less<T>());
	}

      /** Compare all coordinates against other vector. ( <= )
	  Similar to testing whether we are on one side of three planes. */
      bool  all_le (const V& v) const
	{
	  return std::inner_product(data, &data[N], v.get(), true, 
				    std::logical_and<bool>(), std::less_equal<T>());
	}

      /** Compare all coordinates against other vector. ( > )
	  Similar to testing whether we are on one side of three planes. */
      bool  all_g  (const V& v) const
	{
	  return std::inner_product(data, &data[N], v.get(), true, 
				    std::logical_and<bool>(), std::greater<T>());
	}

      /** Compare all coordinates against other vector. ( >= )
	  Similar to testing whether we are on one side of three planes. */
      bool  all_ge (const V& v) const
	{
	  return std::inner_product(data, &data[N], v.get(), true, 
				    std::logical_and<bool>(), std::greater_equal<T>());
	}


      //----------------------------------------------------------------------
      // Assignment operators
      //----------------------------------------------------------------------

      /// Assignment multiplication with scalar.
      const V& operator *=(T k) 
	{ 
	  std::transform(data, &data[N], data, std::bind2nd(std::multiplies<T>(), k));
	  return static_cast<const V&>(*this);
	}

      /// Assignment division with scalar.
      const V& operator /=(T k)
	{ 
	  std::transform(data, &data[N], data, std::bind2nd(std::divides<T>(), k));
	  return static_cast<const V&>(*this);
	}

      /// Assignment addition with scalar. Adds scalar to each coordinate.
      const V& operator +=(T k) 
	{
	  std::transform(data, &data[N], data, std::bind2nd(std::plus<T>(), k));
	  return  static_cast<const V&>(*this);
	}

      /// Assignment subtraction with scalar. Subtracts scalar from each coord.
      const V& operator -=(T k) 
	{ 
	  std::transform(data, &data[N], data, std::bind2nd(std::minus<T>(), k));
	  return  static_cast<const V&>(*this);
	}

      /** Assignment multiplication with vector. 
	  Multiply each coord independently. */
      const V& operator *=(const V& v) 
	{ 
	  std::transform(data, &data[N], v.get(), data, std::multiplies<T>());
	  return  static_cast<const V&>(*this);
	}

      /// Assigment division with vector. Each coord divided independently.
      const V& operator /=(const V& v)
	{
	  std::transform(data, &data[N], v.get(), data, std::divides<T>());
	  return  static_cast<const V&>(*this);
	}

      /// Assignmment addition with vector.
      const V& operator +=(const V& v) 
	{
	  std::transform(data, &data[N], v.get(), data, std::plus<T>());
	  return  static_cast<const V&>(*this);
	}
		
      /// Assignment subtraction with vector.
      const V& operator -=(const V& v) 
	{ 
	  std::transform(data, &data[N], v.get(), data, std::minus<T>());
	  return  static_cast<const V&>(*this);
	}


      //----------------------------------------------------------------------
      // Unary operators on vectors
      //----------------------------------------------------------------------

      /// Negate vector.
      const V operator - () const
	{
	  V v_new;
	  std::transform(data, &data[N], v_new.get(), std::negate<T>());
	  return v_new;
	}

      //----------------------------------------------------------------------
      // Binary operators on vectors
      //----------------------------------------------------------------------

      /** Multiply vector with vector. Each coord multiplied independently
	  Do not confuse this operation with dot product. */
      const V operator * (const V& v1) const
	{
	  V v_new;
	  std::transform(data, &data[N], v1.get(), v_new.get(), std::multiplies<T>());
	  return v_new;
	}

      /// Add two vectors
      const V operator + (const V& v1) const
	{
	  V v_new;
	  std::transform(data, &data[N], v1.get(), v_new.get(), std::plus<T>());
	  return v_new;
	}

      /// Subtract two vectors. 
      const V operator - (const V& v1) const
	{
	  V v_new;
	  std::transform(data, &data[N], v1.get(), v_new.get(), std::minus<T>());
	  return v_new;
	}

      /// Divide two vectors. Each coord separately
      const V operator / (const V& v1) const
	{
	  V v_new;
	  std::transform(data, &data[N], v1.get(), v_new.get(), std::divides<T>());
	  return v_new;
	}

      //----------------------------------------------------------------------
      // Binary operators on vector and scalar
      //----------------------------------------------------------------------

      /// Multiply scalar onto vector.
      const V operator * (T k) const
	{
	  V v_new;
	  std::transform(data, &data[N], v_new.get(), std::bind2nd(std::multiplies<T>(), k));
	  return v_new;
	}
  

      /// Divide vector by scalar.
      const V operator / (T k) const
	{
	  V v_new;
	  std::transform(data, &data[N], v_new.get(), std::bind2nd(std::divides<T>(), k));
	  return v_new;      
	}


      /// Return the smallest coordinate of the vector
      const T min_coord() const 
	{
	  return *std::min_element(data, &data[N]);
	}

      /// Return the largest coordinate of the vector
      const T max_coord() const
	{
	  return *std::max_element(data, &data[N]);
	}

    };

    template <class T, class V, unsigned int N> 
    inline std::ostream& operator<<(std::ostream&os, const ArithVec<T,V,N>& v)
    {
      os << "[ ";
      for(unsigned int i=0;i<N;i++) os << v[i] << " ";
      os << "]";
      return os;
    }

  /// Get from operator for ArithVec descendants. 
    template <class T,class V, unsigned int N>
    inline std::istream& operator>>(std::istream&is, ArithVec<T,V,N>& v)
    {
        is >> std::ws;
        if (is.peek() == '[')
            is.ignore();
        is >> std::ws;
        for (int c=0; c<N; ++c)
        {
            is >> v(c) >> std::ws;
        }
        if (is.peek() == ']')
            is.ignore();
        return is;
    }


  /** Dot product for two vectors. The `*' operator is 
      reserved for coordinatewise	multiplication of vectors. */
  template <class T,class V, unsigned int N>
    inline T dot(const ArithVec<T,V,N>& v0, const ArithVec<T,V,N>& v1)
    {
      return std::inner_product(v0.get(), v0.get() + N, v1.get(), T(0));
    }

  /** Compute the sqr length by taking dot product of vector with itself. */
  template <class T,class V, unsigned int N>
    inline T sqr_length(const ArithVec<T,V,N>& v)
    {
      return dot(v,v);
    }

  /** Multiply double onto vector. This operator handles the case 
      where the vector is on the righ side of the `*'.
	 
      \note It seems to be optimal to put the binary operators inside the
      ArithVec class template, but the operator functions whose 
      left operand is _not_ a vector cannot be inside, hence they
      are here.
      We need three operators for scalar * vector although they are
      identical, because, if we use a separate template argument for
      the left operand, it will match any type. If we use just T as 
      type for the left operand hoping that other built-in types will
      be automatically converted, we will be disappointed. It seems that 
      a float * ArithVec<float,Vec3f,3> function is not found if the
      left operand is really a double.
  */

  template<class T, class V, unsigned int N>
    inline const V operator * (double k, const ArithVec<T,V,N>& v) 
    {
      return v * k;
    }

  /** Multiply float onto vector. See the note in the documentation
      regarding multiplication of a double onto a vector. */
  template<class T, class V, unsigned int N>
    inline const V operator * (float k, const ArithVec<T,V,N>& v) 
    {
      return v * k;
    }

  /** Multiply unsigned int onto vector. See the note in the documentation
      regarding multiplication of a double onto a vector. */
  template<class T, class V, unsigned int N>
    inline const V operator * (int k, const ArithVec<T,V,N>& v) 
    {
      return v * k;
    }

  /** Returns the vector containing for each coordinate the smallest
      value from two vectors. */
  template <class T,class V, unsigned int N>
    inline V v_min(const ArithVec<T,V,N>& v0, const ArithVec<T,V,N>& v1)
    {
      V v;
      std::transform(v0.get(), v0.get() + N, v1.get(), v.get(), std::ptr_fun(s_min<T>));
      return v;
    }

  /** Returns the vector containing for each coordinate the largest 
      value from two vectors. */
  template <class T,class V, unsigned int N>
    inline V v_max(const ArithVec<T,V,N>& v0, const ArithVec<T,V,N>& v1)
    {
      V v;
      std::transform(v0.get(), v0.get() + N, v1.get(), v.get(), std::ptr_fun(s_max<T>));
      return v;
    }


}

#endif
