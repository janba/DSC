#pragma once


namespace is_mesh
{
    class Key
    {
    protected:
        unsigned int key;
        unsigned int dim;
        
        Key() : Key(-1,-1)
        {
            
        }
        
        Key(unsigned int _key, unsigned int _dim) : key(_key), dim(_dim)
        {
            
        }
        
    public:
        bool is_valid() const
        {
            return dim != static_cast<unsigned int>(-1);
        }
        
        const int get_dim() const
        {
            return dim;
        }
        
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
    
    class NodeKey : public Key
    {
    public:
        NodeKey() : Key() {}
        NodeKey(unsigned int k) : Key(k, 0) {}
    };
    
    class EdgeKey : public Key
    {
    public:
        EdgeKey() : Key() {}
        EdgeKey(unsigned int k) : Key(k, 1) {}
    };
    
    class FaceKey : public Key
    {
    public:
        FaceKey() : Key() {}
        FaceKey(unsigned int k) : Key(k, 2) {}
    };
    
    class TetrahedronKey : public Key
    {
    public:
        TetrahedronKey() : Key() {}
        TetrahedronKey(unsigned int k) : Key(k, 3) {}
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
