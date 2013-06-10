#ifndef ATTRIBUTES_H
#define ATTRIBUTES_H

/** 2013 Marek K. Misztal and Mark Viinblad Jensen **/

template<typename MT>
class default_node_traits
{
private:
    typedef typename MT::real_type      T;
    typedef typename MT::vector3_type   V;
    
    V p;
    V p_new;
    unsigned int flags;
    
public:
    
    default_node_traits() : p(0.,0.,0.), p_new(0.,0.,0.), flags(0) {}
    default_node_traits(T const & x, T const & y, T const & z) : p(x,y,z), p_new(x,y,z), flags(0) {}
    
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
    
    bool is_locked()    { return (flags%2 == 1); }
    bool is_boundary()  { return ((flags>>1)%2 == 1); }
    bool is_interface()    { return ((flags>>3)%2 == 1); }
    bool is_stable()    { return ((flags>>4)%2 == 1); }
    
    void set_locked(bool b)
    {
		unsigned int mask = 1;
		if (b) flags |= mask;
		else flags &= ~mask;
    }
    
    void set_boundary(bool b)
    {
		unsigned int mask = 2;
		if (b) flags |= mask;
		else flags &= ~mask;
    }
    
    void set_interface(bool b)
    {
		unsigned int mask = 8;
		if (b) flags |= mask;
		else flags &= ~mask;
    }
    
    void set_stable(bool b)
    {
		unsigned int mask = 16;
		if (b) flags |= mask;
		else flags &= ~mask;
    }
};

template <typename MT>
class default_edge_traits
{
public:
    
    typedef typename MT::real_type      T;
    typedef typename MT::vector3_type   V;
    
public:
    unsigned int flags;
    
    default_edge_traits() : flags(0) {}
    default_edge_traits(const default_edge_traits& t) : flags(t.flags) {}
    
    bool is_locked()    { return (flags%2 == 1); }
    bool is_boundary()  { return ((flags>>1)%2 == 1); }
    bool is_interface()    { return ((flags>>3)%2 == 1); }
    
    void set_locked(bool b)
    {
		unsigned int mask = 1;
		if (b) flags |= mask;
		else flags &= ~mask;
    }
    
    void set_boundary(bool b)
    {
		unsigned int mask = 2;
		if (b) flags |= mask;
		else flags &= ~mask;
    }
    
    void set_interface(bool b)
    {
		unsigned int mask = 8;
		if (b) flags |= mask;
		else flags &= ~mask;
    }
};

template <typename MT>
class default_face_traits
{
public:
    
    typedef typename MT::real_type      T;
    typedef typename MT::vector3_type   V;
    
public:
    
    unsigned int flags;
    
    default_face_traits() : flags(0) {}
    default_face_traits(const default_face_traits& t) : flags(t.flags) {}
    
    bool is_locked()    { return (flags%2 == 1); }
    bool is_boundary()  { return ((flags>>1)%2 == 1); }
    bool is_processed() { return ((flags>>2)%2 == 1); }
    bool is_interface() { return ((flags>>3)%2 == 1); }
    bool is_error() { return ((flags>>4)%2 == 1); }
    
    void set_error(bool b)
    {
		unsigned int mask = 16;
		if (b) flags |= mask;
		else flags &= ~mask;
    }
    
    void set_locked(bool b)
    {
		unsigned int mask = 1;
		if (b) flags |= mask;
		else flags &= ~mask;
    }
    
    void set_boundary(bool b)
    {
		unsigned int mask = 2;
		if (b) flags |= mask;
		else flags &= ~mask;
    }
    
    void set_processed(bool b)
    {
		unsigned int mask = 4;
		if (b) flags |= mask;
		else flags &= ~mask;
    }
    
    void set_interface(bool b)
    {
		unsigned int mask = 8;
		if (b) flags |= mask;
		else flags &= ~mask;
    }
};

template <typename MT>
class default_tetrahedron_traits
{
public:
    
    typedef typename MT::real_type      T;
    typedef typename MT::vector3_type   V;
    
public:
    unsigned int label;
    
    default_tetrahedron_traits() : label(0) {}
    default_tetrahedron_traits(const default_tetrahedron_traits& t) : label(t.label) {}
    
};

//ATTRIBUTES_H
#endif