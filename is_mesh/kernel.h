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

#include <cassert>
#include <iostream>
#include <vector>
#include <set>

#include "kernel_iterator.h"

namespace is_mesh
{

    class ISMesh;

    /**
     * namespace that defines auxiliary data structures and helper methods.
     */
    namespace util
    {
        
        /**
         * Kernel Element struct
         * A struct to hold the actual data along with the state of the memory cell.
         */
        template< typename key_t_, typename value_t_ >
        struct kernel_element
        {
            using value_type = value_t_;

            enum state_type { VALID, MARKED, EMPTY, EXCLUDED };
            
            kernel_element() : value(nullptr), state{kernel_element::EMPTY} { }
            
            kernel_element(const kernel_element& ke) = delete;

            kernel_element(kernel_element&& ke) : value{std::move(ke.value)}, state{std::move(ke.state)} {
            }

            kernel_element&& operator=(kernel_element&& other) {
                if (this != &other) {
                    value = std::move(other.value);
                    state = other.state;
                }
                return *this;
            }
            
            value_type              value;
            state_type              state;
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
     */
    template<typename key_type, typename value_type>
    class kernel
    {
    public:
        using kernel_type = kernel<key_type, value_type>;
        using kernel_key_type = key_type;
        using kernel_element = util::kernel_element<key_type, value_type>;
        using iterator = kernel_iterator<kernel_type>;
        using const_iterator = iterator const;

        friend class kernel_iterator<kernel_type>;

    private:
        std::vector<kernel_element> m_data;
        std::vector<key_type> m_data_freelist;
        std::vector<key_type> m_data_marked_for_deletion;
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
//          assume key_type is integer type
            assert(k >= 0 || !"looked up with negative element");
            assert((int)k < m_data.size() || !"k out of range");
            return m_data[k];
        }
        
        /**
         * Finds the next free cell in the kernel. If necessary new memory is allocated.
         *
         * @return The handle to the new cell.
         */
        kernel_element& get_next_free_cell(unsigned int &key)
        {
            if (m_data_freelist.size()==0){
                key = (unsigned int)m_data.size();
                m_data.emplace_back();
                return m_data.back();
            } else {
                key = (unsigned int)m_data_freelist.back();
                m_data_freelist.pop_back();
                return m_data[key];
            }
        }
    public:
        
        /**
         * Default constructor for the kernel. This will initialize all lists and memory used by the kernel.
         * The initial should always be greater than 0.
         *
         * @param reservedSize The initial capacity of the kernel.
         */
        kernel(size_t reservedSize =64)
        {
            m_data.reserve(reservedSize);
        }
        
        /**
         * Kernel destructor, frees allocated memory by the kernel.
         */
        ~kernel()
        {
        }
        
        /**
         * The size of the kernel. That is the number of valid elements in the kernel.
         */
        size_t size() const     { return m_data.size() - m_data_freelist.size(); }


        /**
         * The capacity of the kernel. Note not all elements may in use.
         */
        size_t capacity() const     { return m_data.size(); }

        /**
         * Returns a boolean value indicating if the size is zero.
         */
        bool      empty()    { return size() == 0; }

        /**
         * Creates (or inserts) a new element into the kernel. From this point forward the
         * memory management of the element is controlled by the kernel.
         *
         * @param attributes     The type traits of the element to be inserted.
         *
         * @return      An iterator pointing to the element.
         */
        template<typename... Values>
        const_iterator create(ISMesh* isMesh, Values... values)
        {
            unsigned int key;
            kernel_element& cur = get_next_free_cell(key);

            assert(cur.state != kernel_element::VALID || !"Cannot create new element, duplicate key.");
            assert(cur.state != kernel_element::MARKED || !"Attempted to overwrite a marked element.");
            assert(cur.state != kernel_element::EXCLUDED || !"Attempted to overwrite an excluded element.");

            cur.value = value_type{isMesh, values...};
            cur.state = kernel_element::VALID;
            return iterator(this, key);
        }
        
        /**
         * Returns the one-past-the-end iterator to the kernel.
         */
        const_iterator end() const
        {
            return iterator(this, key_type{(unsigned int)m_data.size()});
        }
        
        /**
         * Returns an iterator to the first element of the kernel.
         */
        const_iterator begin() const
        {
            unsigned int i = 0;
            // find first valid element (if any)
            for (;i<m_data.size();i++){
                if (m_data[i].state == kernel_element::VALID){
                    break;
                }
            }
            return iterator(this, key_type{i});
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
            assert (p.state == kernel_element::VALID || p.state == kernel_element::EXCLUDED || !"Attempted to remove a non-valid element!");
            if (p.state != kernel_element::VALID && p.state != kernel_element::EXCLUDED )
            {
                //No element with that key.
                return;
            }
            p.state = kernel_element::MARKED;
            m_data_marked_for_deletion.push_back(k);
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
            m_data.clear();
            m_data_freelist.clear();
            m_data_marked_for_deletion.clear();
        }
        
        /**
         * Converts a key or handle to an iterator pointing to the same cell.
         *
         * @param k     The handle of the cell.
         */
        iterator find_iterator(key_type const & k)
        {
            if (is_valid(k))
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
            assert(tmp.state == kernel_element::VALID || tmp.state == kernel_element::EXCLUDED);
            return tmp.value;
        }

        /**
        *  Returns pointer to the underlying array serving as element storage.
        */
        char * data(){
            return (char*)m_data.data();
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
            if (tmp.state == kernel_element::VALID || tmp.state == kernel_element::EXCLUDED ) return true;
            return false;
        }
        
        /**
         * Commits all the changes in the kernel, and permanently removes all the marked elements.
         */
        std::vector<key_type> commit_all()
        {
            std::vector<key_type> deletedKeys = m_data_marked_for_deletion;
            for (auto key : m_data_marked_for_deletion){
                auto & p = m_data[key];

                //m_alloc.destroy(&p);  // needed?
                p.state = kernel_element::EMPTY;
                m_data_freelist.push_back(key);
            }
            m_data_marked_for_deletion.clear();
            return deletedKeys;
        }

        void revert_excluded(){
            for (auto &e : m_data){
                if (e.state == kernel_element::state_type::EXCLUDED){
                    e.state = kernel_element::state_type::VALID;
                }
            }
        }

        void exclude_using_include_set(std::set<key_type> include_set){
            for (unsigned int i=0;i<m_data.size();i++){
                auto& e= m_data[i];
                key_type k{i};
                if (e.state == kernel_element::state_type::VALID){
                    if (include_set.find(k) == include_set.end()) {
                        e.state = kernel_element::state_type::EXCLUDED;
                    }
                }
            };
        }
        
        /**
         * Commits all changes and permanently deletes all marked elements.
         * The garbage collect routine performs one functions. First it commits all changes on the undo
         * stack and clears all undo marks.
         * The operation runs in O(n) - where n is m_capacity or the size of allocated memory (not 
         * efficient, but cleans lists).
         */
        std::vector<key_type> garbage_collect()
        {
            return commit_all();
        }
    };
}
