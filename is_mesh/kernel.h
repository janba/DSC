
#pragma once

#include <cassert>
#include <iostream>

#include <is_mesh/kernel_iterator.h>

namespace is_mesh
{
    
    /**
     * namespace that defines auxiliary data structures and helper methods.
     */
    namespace util
    {
        
        /**
         * Kernel Element struct
         * A struct to hold the actual data along with auxiliary informations, like
         * the state of af memory cell, and list pointers.
         */
        template<
        typename value_t_
        , typename key_t_
        >
        struct kernel_element
        {
            typedef value_t_            value_type;
            typedef key_t_              key_type;
            
            enum state_type { VALID, MARKED, EMPTY };
            
            kernel_element() : value() { }
            
            kernel_element(kernel_element const & other)
            {
                value = other.value;
                key   = other.key;
                next  = other.next;
                prev  = other.prev;
            }
            
            value_type              value;
            key_type                key;
            state_type              state;  //one of
            key_type                next;
            key_type                prev;
        };
    }
    
    /**
     * Memory Kernel developed for the DSC project.
     * The kernel uses the supplied allocator to allocate memory for data structures,
     * like the IS mesh used in DSC. The kernel uses array-based allocation. The kernel supports
     * undo operations. Each cell in the kernel uses an excess of 12 bytes, which is used to
     * support fast iterators through the kernel and the undo functionality.
     *
     * @param value_type The type of the elements that is to be stored in the kernel. The value_type
     *        must have the typedef type_traits.
     * @param key_type The type of the keys used in the kernel. Should be an integer type.
     * @param allocator_type This type must conform to the C++ standard allocator concept. The
     *        default argument is the default STL allocator.
     */
    template<typename value_t_, typename key_t_, typename allocator_t_ = std::allocator<util::kernel_element<value_t_, key_t_> > >
    class kernel
    {
    public:
        typedef          value_t_                                       value_type;
        typedef          key_t_                                         key_type;
        typedef          allocator_t_                                   allocator_type;
        typedef          kernel<value_type, key_type, allocator_type>   kernel_type;
        typedef typename allocator_type::template rebind<value_type>    rebind_type;
        typedef typename rebind_type::other                             value_allocator;
        typedef          util::kernel_element<value_type, key_type>     kernel_element;
        typedef          kernel_iterator<kernel_type>                   iterator;
        typedef          iterator const                                 const_iterator;
        typedef typename allocator_type::size_type                      size_type;
        typedef typename allocator_type::reference                      reference;
        typedef typename allocator_type::const_reference                const_reference;
        typedef          value_type*                                    pointer;
        
        friend class kernel_iterator<kernel_type>;
        
    private:
        typedef          kernel_element*                                element_pointer;
        typedef typename value_type::type_traits                        type_traits;
        typedef typename kernel_element::state_type                     state_type;
        
        allocator_type        m_alloc;               //the allocator
        kernel_element*       m_mem;                 //the allocated memory
        key_type              m_first;               //Holds the first allocated element in the collection
        key_type              m_last;                //Holds the last allocated element in the collection
        key_type              m_first_marked;        //Points to the first element (or last) that is marked for deletion
        key_type              m_last_marked;         //Points to the first element (or last) that is marked for deletion
        key_type              m_first_empty;         //points to the first empty cell (the next to be allocated to).
        key_type              m_last_empty;          //points to the first empty cell (the next to be allocated to).
        key_type              m_past_the_end;        //used to mark the end of iterations
        
        
        size_type             m_size;                //How many elements are in the collection
        size_type             m_shadow_size;         //How many elements are really allocated (along with elements marked for deletion)
        size_type             m_capacity;            //How many elements can currently be allocated without expansion
        
        size_type             m_initial_size;        //the size by wich we grow
        
    private:
        /**
         * Converts an indirect reference to the direct memory reference currently allocated by the cell.
         *
         * @param k       The index or key of the cell that is to be looked up.
         *
         * @return        A reference to the memory occupied by the cell.
         */
        kernel_element& lookup(key_type k)
        {
            //asume key_type is integer type
            assert(k >= 0 || !"looked up with negative element");
            assert((int)k < m_capacity || !"k out of range");
            return m_mem[k];
        }
        
        /**
         * Finds the next free cell in the kernel. If neccesary new memory is allocated.
         *
         * @return The handle to the new cell.
         */
        key_type get_next_free_cell()
        {
            key_type k;
            if (m_shadow_size == m_capacity || m_first_empty == m_past_the_end) grow(); //are we full?
            //look for fragments in mem
            k = m_first_empty;
            assert((int)k < m_capacity || !"Get next free cell extendended the allocated memory.");
            return k;
        }
        
        /**
         * Links a cell at the end of a double-linked list.
         *
         * @param cur       The cell that is changed.
         * @param front     The head of the list that should add the cell.
         * @param end       The tail of the list that should add the cell.
         */
        void link_at_end(kernel_element& cur, key_type& front, key_type& end)
        {
            //pointer arithmetics
            if (front == m_past_the_end)
            {
                //first' element in collection
                front = cur.key;
                cur.prev = m_past_the_end;
            } else {
                lookup(end).next = cur.key;
                cur.prev = lookup(end).key;
            }
            end       = cur.key;
            cur.next  = m_past_the_end;
        }
        
        /**
         * Links a cell at the head of a double-linked list.
         *
         * @param cur       The cell that is changed.
         * @param front     The head of the list that should add the cell.
         * @param end       The tail of the list that should add the cell.
         */
        void link_at_head(kernel_element& cur, key_type& front, key_type& end)
        {
            //pointer arithmetics
            if (front == m_past_the_end)
            {
                //first' element in collection
                end = cur.key;
                cur.next = m_past_the_end;
            } else {
                //lookup(front).next = cur.key;
                lookup(front).prev = cur.key;
                cur.next = front;
            }
            front     = cur.key;
            cur.prev  = m_past_the_end;
        }
        
        /**
         * Unlinks a cell from a double-linked list.
         *
         * @param cur       The cell that is changed.
         * @param front     The head of the list from which the cell is removed.
         * @param end       The tail of the list from which the cell is removed.
         */
        void unlink(const kernel_element& cur, key_type& front, key_type& end)
        {
            if (cur.key != front)
            {
                lookup(cur.prev).next = cur.next;
            } else {
                front = cur.next;
            }
            if (cur.key != end)
            {
                lookup(cur.next).prev = cur.prev;
            } else {
                end = cur.prev;
            }
        }
        
        /**
         * Unlinks a cell from a double-linked list without updating the front and the end pointers.
         * This is a helper function for undo() method, and shouldn't be used elsewhere.
         *
         * @param cur       The cell that is changed.
         */
        void unsafe_unlink(const kernel_element& cur)
        {
            if (cur.prev != m_past_the_end)
            {
                lookup(cur.prev).next = cur.next;
            }
            if (cur.next != m_past_the_end)
            {
                lookup(cur.next).prev = cur.prev;
            }
        }
        
        /**
         * Links a cell into a linked list - makes the assumption that the cell knows where
         * it is to be inserted into the list, that is that the cells next and prev indicies
         * should be valid.
         *
         * @param cur       The cell that is relinked.
         * @param front     The head of the list to which the cell is re-inserted.
         * @param end       The tail of the list to which the cell is re-inserted.
         */
        void relink(kernel_element& cur, key_type& front, key_type& end)
        {
            if (cur.next != m_past_the_end)
                lookup(cur.next).prev = cur.key;
            else
                end = cur.key;
            if (cur.prev != m_past_the_end)
                lookup(cur.prev).next = cur.key;
            else
                front = cur.key;
        }
        
        /**
         * Doubles the size of the allocated memory for the kernel
         */
        void grow()
        {
            resize(m_capacity << 1); //double size
        }
        
        /**
         * Changes the size of the kernel - either grow or shrink.
         * All allocated memory is copied.
         *
         * @param new_size        The new size of the kernel.
         */
        void resize(size_type new_size)
        {
            //realloc mem
            kernel_element* new_mem = m_alloc.allocate(new_size);        //allocate new space
            memcpy(new_mem, m_mem, m_capacity * sizeof(kernel_element));  //copy data to new mem location
            m_alloc.deallocate(m_mem, m_capacity);                         //deallocate old mem
            //update values
            for (unsigned int i = static_cast<unsigned int>(m_capacity); i < static_cast<unsigned int>(new_size); ++i)
            {
                new_mem[i].state = kernel_element::EMPTY;
                new_mem[i].key = i;
            }
            
            if (m_first_empty == m_past_the_end)
            {
                m_first_empty = static_cast<unsigned int>(m_capacity);
                new_mem[m_capacity].prev = m_past_the_end;
            }
            else
            {
                assert (m_last_empty != m_past_the_end);
                new_mem[m_capacity].prev = m_last_empty;
                new_mem[m_last_empty].next = static_cast<int>(m_capacity);
            }
            
            m_last_empty = static_cast<unsigned int>(new_size) - 1;
            new_mem[m_capacity].next   = static_cast<unsigned int>(m_capacity) + 1;
            new_mem[m_last_empty].prev = m_last_empty-1;
            new_mem[m_last_empty].next = m_past_the_end;
            
            for (unsigned int i = static_cast<unsigned int>(m_capacity) + 1; i < static_cast<unsigned int>(new_size) - 1; ++i)
            {
                new_mem[i].prev = i-1;
                new_mem[i].next = i+1;
            }
            
            m_capacity = new_size;
            m_mem = new_mem;
        }
        
    public:
        
        /**
         * Default constructor for the kernel. This will initialize all lists and memory used by the kernel.
         * The initialsize should always be greater than 0.
         *
         * @param size The initial size of the kernel. Has a default value.
         */
        kernel(size_type size =64) : m_initial_size(size)
        { //default constructor
            assert(size > 0 || !"Cannot construct kernels with no initial memory");
            //allocate some uninitialized memory
            m_mem = m_alloc.allocate(m_initial_size);
            //allocate past the end
            m_past_the_end = static_cast<unsigned int>(-1);
            //initilize points
            m_first        = m_past_the_end;
            m_last         = m_past_the_end;
            m_first_marked = m_past_the_end;
            m_last_marked  = m_past_the_end;
            m_first_empty  = m_past_the_end;
            m_last_empty   = m_past_the_end;
            
            m_capacity     = m_initial_size;
            m_size         = 0;
            m_shadow_size  = 0;
            
            for (unsigned int i = 0; i < m_capacity; ++i)
            {
                m_mem[i].state = kernel_element::EMPTY;
                m_mem[i].key = i;
            }
            
            m_first_empty = 0;
            m_last_empty = static_cast<unsigned int>(m_capacity)-1;
            m_mem[m_first_empty].prev = m_past_the_end;
            m_mem[m_first_empty].next = m_first_empty+1;
            m_mem[m_last_empty].prev  = m_last_empty-1;
            m_mem[m_last_empty].next  = m_past_the_end;
            
            for (unsigned int i = 1; i < m_capacity-1; ++i)
            {
                m_mem[i].prev = i-1;
                m_mem[i].next = i+1;
            }
        }
        
        /**
         * Kernel destructor, frees allocated memory by the kernel.
         */
        ~kernel()
        { //destructor
            //call the deallocate of all managed objects
            for(size_type i=0; i<m_capacity; ++i)
            {
                if (m_mem[i].state == kernel_element::VALID || m_mem[i].state == kernel_element::MARKED)
                    m_alloc.destroy(&m_mem[i]);
            }
            //deallocate kernel
            m_alloc.deallocate(m_mem, m_capacity);
        }
        
        /**
         * The size of the kernel. That is the number of valid elements in the kernel.
         */
        size_type size()     { return m_size; }
        
        /**
         * Returns the maximum size the kernel can be. This is dertemined by the allocator.
         */
        size_type max_size() { return m_alloc.max_size(); }
        
        /**
         * Returns the size of the allocated memory of the kernel.
         */
        size_type capacity() { return m_capacity; }
        
        /**
         * Returns a boolean value indicating if the size is zero.
         */
        bool      empty()    { return m_size == 0; }
        
        /**
         * Creates (or inserts) a new element into the kernel. From this point forward the
         * memory mangement of the element is controlled by the kernel.
         *
         * @param k     A key indicating the handle the new element should have. If the cell coresponding
         *              to the cell is not available an assertion is thrown.
         * @param d     The type traits of the element to be inserted.
         *
         * @return      An iterator pointing to the element.
         */
        const_iterator create(key_type k, const type_traits& d)
        //const_iterator create(key_type k, type_traits & d)
        {
            //element_pointer cur = get_next_free_cell();
            kernel_element& cur = lookup(k);
            cur.key = k;
            assert(cur.state != kernel_element::VALID || !"Cannot create new element, duplicate key.");
            assert(cur.state != kernel_element::MARKED || !"Attempted to overwrite a marked element.");
            /*  if (cur.state == kernel_element::MARKED)
             {
             unlink(cur, m_first_marked, m_last_marked);
             } else */ {
                 //only setting value if mem was empty (or unused)
                 //we copy the memory instead of assinging, as assigning would create a temp object.
                 value_allocator(m_alloc).construct( &cur.value , value_type(d) );
             }
            if (cur.state == kernel_element::EMPTY)
            {
                unlink(cur, m_first_empty, m_last_empty);
            }
            cur.state = kernel_element::VALID;
            ++m_size;
            ++m_shadow_size;
            link_at_end(cur, m_first, m_last);
            return iterator(this, cur.key);
        }
        
        /**
         * Creates (or inserts) a new element into the kernel. From this point forward the
         * memory mangement of the element is controlled by the kernel and all access to the
         * element should be through an iterator.
         * Default type traits is assinged as value for the element.
         *
         * @param k     A key indicating the handle the new element should have. If the cell coresponding
         *              to the cell is not available an assertion is thrown.
         * @return      An iterator pointing to the element.
         */
        const_iterator create(key_type k)
        {
            return create(k, type_traits());
        }
        
        /**
         * Creates (or inserts) a new element into the kernel. From this point forward the
         * memory mangement of the element is controlled by the kernel and all access to the
         * element should be through an iterator.
         * The next free cell is chosen for the element.
         *
         * @param d     The type traits of the element to be inserted.
         *
         * @return      An iterator pointing to the element.
         */
        const_iterator create(type_traits d)
        {
            key_type key = get_next_free_cell();
            return create(key, d);
        }
        
        /**
         * Creates (or inserts) a new element into the kernel. From this point forward the
         * memory mangement of the element is controlled by the kernel and all access to the
         * element should be through an iterator.
         * Default type traits is assinged as value for the element.
         * The next free cell is chosen for the element.
         *
         * @return      An iterator pointing to the element.
         */
        const_iterator create()
        {
            return create(type_traits());
        }
        
        /**
         * Returns the one-past-the-end iterator to the kernel.
         */
        const_iterator end() const
        {
            return iterator(this, m_past_the_end);
        }
        
        /**
         * Returns an iterator to the first element of the kernel.
         */
        const_iterator begin() const
        {
            return iterator(this, m_first);
        }
        
        /**
         * Deletes an element in the kernel. The element is not actually deleted, it is merely
         * marked for deletion. If no element is found in the cell nothing is performed.
         *
         * @param k   The key or handle to the element to be deleted.
         */
        void erase(key_type const & k)
        {
            kernel_element& p = lookup(k);
            assert (p.state == kernel_element::VALID || !"Attempted to remove a non-valid element!");
            if (p.state != kernel_element::VALID)
            {
                //No element with that key.
                return;
            }
            //relink neighbors
            unlink(p, m_first, m_last);
            p.state = kernel_element::MARKED;
            --m_size;
            //put in marked list
            link_at_end(p, m_first_marked, m_last_marked);
        }
        
        /**
         * Deletes an element in the kernel. The element is not actually deleted, it is merely
         * marked for deletion. If no element is found in the cell nothing is performed.
         *
         * @param it   An iterator pointing to an element in the kernel.
         */
        void erase(const iterator& it)
        {
            erase(it.key());
        }
        
        /**
         * Permanently deletes all elements in the kernel and resets all internal lists and stacks.
         */
        void clear()
        {
            m_alloc.deallocate(m_mem, m_capacity);
            m_mem = m_alloc.allocate(m_capacity);
            //reset pointers
            m_first        = m_past_the_end;
            m_last         = m_past_the_end;
            m_first_marked = m_past_the_end;
            m_last_marked  = m_past_the_end;
            
            for (unsigned int i = 0; i < m_capacity; ++i)
            {
                m_mem[i].state = kernel_element::EMPTY;
                m_mem[i].key = i;
            }
            
            m_first_empty = 0;
            m_last_empty = static_cast<unsigned int>(m_capacity)-1;
            m_mem[m_first_empty].prev = m_past_the_end;
            m_mem[m_first_empty].next = m_first_empty+1;
            m_mem[m_last_empty].prev  = m_last_empty-1;
            m_mem[m_last_empty].next  = m_past_the_end;
            
            for (unsigned int i = 1; i < m_capacity-1; ++i)
            {
                m_mem[i].prev = i-1;
                m_mem[i].next = i+1;
            }
            
            m_size         = 0;
            m_shadow_size  = 0;
        }
        
        /**
         * Converts a key or handle to an iterator pointing to the same cell.
         *
         * @param k     The handle of the cell.
         */
        iterator find_iterator(key_type const & k)
        {
            //we dont just return iterator(this, k) as this is not defensive enough, we need to return valid values.
            kernel_element& tmp = lookup(k);
            if (tmp.state == kernel_element::VALID && tmp.key == k)
                return iterator(this, k);
            else
                return end();
        }
        
        /**
         * Returns a managed object. Beware of deallocating or other memory
         * handlings of the returned object, as this might lead to undefined
         * behavior in the kernel.
         *
         * @param k     The handle to the object.
         */
        value_type & find(key_type const & k)
        {
            kernel_element& tmp = lookup(k);
            assert(tmp.state == kernel_element::VALID);
            assert(tmp.key == k);
            return tmp.value;
        }
        
        /**
         * Returns a managed object.
         *
         * @param k     The handle to the object.
         */
        value_type const & find_const(key_type const & k)
        {
            return find(k);
        }
        
        /**
         * Returns the status of the cell given its key.
         *
         * @param k     The handle to the object.
         * @returns     True if the object is a valid element, false if it is marked for deletion or k refers to an empty cell.
         */
        bool is_valid(key_type const & k)
        {
            kernel_element& tmp = lookup(k);
            if (tmp.state == kernel_element::VALID) return true;
            return false;
        }
        
        /**
         * Commits all the changes in the kernel, and permanently removes all the marked elements.
         */
        void commit_all()
        {
            key_type m = m_first_marked;
            while (m != m_past_the_end)
            {
                kernel_element* p = &m_mem[m_first_marked];
                key_type n = p->next;
                
                unlink(*p, m_first_marked, m_last_marked);
                link_at_head(*p, m_first_empty, m_last_empty);
                
                m_alloc.destroy(p);
                
                p->state = kernel_element::EMPTY;
                p->key = m;
                
                m = n;
            }
            
            m_shadow_size = m_size;
        }
        
        /**
         * Reorders all the lists so that cache trashing will be at a minimum.
         */
        void reorder_lists()
        {
            //first reset list pointers
            m_first        = m_past_the_end;
            m_last         = m_past_the_end;
            m_first_marked = m_past_the_end;
            m_last_marked  = m_past_the_end;
            m_first_empty  = m_past_the_end;
            m_last_empty   = m_past_the_end;
            
            //now rebuild lists in incremental order of key
            for(size_type i = 0; i < m_capacity; ++i) 
            {
                kernel_element & p = lookup(static_cast<int>(i));
                p.next = m_past_the_end; //always set this just to be sure to end the lists - likely to be overwritten
                if ( p.state == kernel_element::VALID) 
                {
                    link_at_end(p, m_first, m_last);
                } 
                else 
                { 
                    if ( p.state == kernel_element::MARKED )
                    {
                        std::cout << "something's wrong" << std::endl;
                        //should've been handled by commit
                        m_alloc.destroy(&p);
                    }
                    p.state = kernel_element::EMPTY;
                    link_at_end(p, m_first_empty, m_last_empty);
                }
            }
            m_shadow_size = m_size;
        }
        
        /**
         * Commits all changes and permanently deletes all marked elements.
         * The garbage collect routine performs two functions. First it commits all changes on the undo 
         * stack and clears all undo marks. Secondly the routine cleans up the kernel and reorders all
         * lists in the kernel. 
         * The operation runs in O(n) - where n is m_capacity or the size of allocated memory (not 
         * efficient, but cleans lists).
         */
        void garbage_collect() 
        {
            commit_all();
            reorder_lists();
        }
    };
    
}
