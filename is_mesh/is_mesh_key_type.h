#pragma once


namespace is_mesh
{
    namespace util
    {
        struct CustomKey
        {
            //this is where the magic happens...
            unsigned int key;
            
            //default constructor - sets key to 0
            CustomKey() : key() { }
            
            //conversion from int
            CustomKey(unsigned int k) : key(k) { }
            
            //conversion to int
            operator unsigned int() { return key; }
            
            friend inline bool operator==(CustomKey    const & a, CustomKey    const & b)   { return a.key == b.key; }
            friend inline bool operator==(CustomKey          & a, CustomKey          & b)   { return a.key == b.key; }
            friend inline bool operator==(unsigned int const & k, CustomKey    const & b)   { return   k   == b.key; }
            friend inline bool operator==(CustomKey    const & a, unsigned int const & k)   { return a.key ==   k;   }
            friend inline bool operator!=(CustomKey    const & a, CustomKey    const & b)   { return a.key != b.key; }
            friend inline bool operator!=(CustomKey          & a, CustomKey          & b)   { return a.key != b.key; }
            friend inline bool operator!=(unsigned int const & k, CustomKey    const & b)   { return   k   != b.key; }
            friend inline bool operator!=(CustomKey    const & a, unsigned int const & k)   { return a.key !=   k;   }
            friend inline bool operator< (CustomKey    const & a, CustomKey    const & b)   { return a.key <  b.key; }
            friend inline bool operator< (CustomKey          & a, CustomKey          & b)   { return a.key <  b.key; }
            friend inline bool operator< (unsigned int const & k, CustomKey    const & b)   { return   k   <  b.key; }
            friend inline bool operator< (CustomKey    const & a, unsigned int const & k)   { return a.key <    k;   }
            
            friend std::ostream& operator<< (std::ostream & os, CustomKey const & a) { return (os << a.key); }
            friend std::istream& operator>> (std::istream & is, CustomKey       & a) { return (is >> a.key); }
        };
    }
    
    struct NodeKey : public util::CustomKey
    {
        NodeKey() : CustomKey() {}
        NodeKey(unsigned int k) : CustomKey(k) {}
        static const int dim = 0;
    };
    struct EdgeKey : public util::CustomKey
    {
        EdgeKey() : CustomKey() {}
        EdgeKey(unsigned int k) : CustomKey(k) {}
        static const int dim = 1;
    };
    struct FaceKey : public util::CustomKey
    {
        FaceKey() : CustomKey() {}
        FaceKey(unsigned int k) : CustomKey(k) {}
        static const int dim = 2;
    };
    struct TetrahedronKey : public util::CustomKey
    {
        TetrahedronKey() : CustomKey() {}
        TetrahedronKey(unsigned int k) : CustomKey(k) {}
        static const int dim = 3;
    };
    
    template<int n>struct key_traits
    {
        typedef unsigned int    key_type;
    };
    
    template<>
    struct key_traits<0>
    {
        typedef NodeKey         key_type;
    };
    template<>
    struct key_traits<1>
    {
        typedef EdgeKey         key_type;
    };
    template<>
    struct key_traits<2>
    {
        typedef FaceKey         key_type;
    };
    template<>
    struct key_traits<3>
    {
        typedef TetrahedronKey   key_type;
    };
    
}
