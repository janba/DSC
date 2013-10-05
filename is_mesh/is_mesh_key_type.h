#pragma once


namespace is_mesh
{
    struct Key
    {
        //this is where the magic happens...
        unsigned int key;
        
        //default constructor - sets key to 0
        Key() : key() { }
        
        //conversion from int
        Key(unsigned int k) : key(k) { }
        
        //conversion to int
        operator unsigned int() { return key; }
        
        friend inline bool operator==(Key    const & a, Key    const & b)   { return a.key == b.key; }
        friend inline bool operator==(Key          & a, Key          & b)   { return a.key == b.key; }
        friend inline bool operator==(unsigned int const & k, Key    const & b)   { return   k   == b.key; }
        friend inline bool operator==(Key    const & a, unsigned int const & k)   { return a.key ==   k;   }
        friend inline bool operator!=(Key    const & a, Key    const & b)   { return a.key != b.key; }
        friend inline bool operator!=(Key          & a, Key          & b)   { return a.key != b.key; }
        friend inline bool operator!=(unsigned int const & k, Key    const & b)   { return   k   != b.key; }
        friend inline bool operator!=(Key    const & a, unsigned int const & k)   { return a.key !=   k;   }
        friend inline bool operator< (Key    const & a, Key    const & b)   { return a.key <  b.key; }
        friend inline bool operator< (Key          & a, Key          & b)   { return a.key <  b.key; }
        friend inline bool operator< (unsigned int const & k, Key    const & b)   { return   k   <  b.key; }
        friend inline bool operator< (Key    const & a, unsigned int const & k)   { return a.key <    k;   }
        
        friend std::ostream& operator<< (std::ostream & os, Key const & a) { return (os << a.key); }
        friend std::istream& operator>> (std::istream & is, Key       & a) { return (is >> a.key); }
    };
    
    struct NodeKey : public Key
    {
        NodeKey() : Key() {}
        NodeKey(unsigned int k) : Key(k) {}
        static const int dim = 0;
    };
    
    struct EdgeKey : public Key
    {
        EdgeKey() : Key() {}
        EdgeKey(unsigned int k) : Key(k) {}
        static const int dim = 1;
    };
    
    struct FaceKey : public Key
    {
        FaceKey() : Key() {}
        FaceKey(unsigned int k) : Key(k) {}
        static const int dim = 2;
    };
    
    struct TetrahedronKey : public Key
    {
        TetrahedronKey() : Key() {}
        TetrahedronKey(unsigned int k) : Key(k) {}
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
