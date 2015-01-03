
#pragma once

namespace is_mesh
{
    template <class kernel_t_>
    class kernel_iterator_value {
        typedef typename kernel_t_::kernel_element element_type;
        typedef typename element_type::value_type value_type;
        typedef typename element_type::key_type key_type;
        key_type     m_key;
        kernel_t_*     m_kernel;
    public:
        kernel_iterator_value(key_type m_key, kernel_t_ *m_kernel) : m_key(m_key), m_kernel(m_kernel) {
        }

        key_type key() const {
            return m_key;
        }

        value_type * get()
        {
            assert(m_kernel->lookup(key_type{m_key}).state == element_type::VALID);
            return &m_kernel->m_data[m_key].value;
        }

        /**
        * The member access operator.
        *
        * @return Pointer to the element contained within the kernel cell.
        */
        value_type * operator->()
        {
            assert(m_kernel->lookup(key_type{m_key}).state == element_type::VALID);
            return &m_kernel->m_data[m_key].value;
        }
    };

    /**
     * An iterator class used by the is_mesh::kernel.
     * The kernel wraps a kernel handle and uses indirect access to the kernel.
     * Iterators are valid through out their existence - only if an element pointed
     * to by an iterator is it invalidated, but if the delete is undone the iterator
     * will be valid again after the undo operation.
     *
     * The garbage collect routine in the kernel might change the ordering of
     * elements in the kernel and provide different results is one iterates through
     * the kernel before and after a garbage collect.
     * The undo operation is guaranteed to keep the ordering before and after a mark/undo
     * pair of operations.
     *
     * @see kernel.
     */
    template <class kernel_t_>
    class kernel_iterator
    {
    private:
        typedef typename kernel_t_::kernel_element element_type;
    public:
        typedef          kernel_t_ kernel_type;
        typedef          kernel_iterator<kernel_type> iterator;
        typedef typename kernel_type::kernel_element kernel_element;
        typedef typename element_type::value_type value_type;
        typedef typename element_type::key_type key_type;
        
    private:
        std::vector<key_type> *m_keys;
        unsigned int     m_key;
        kernel_type*     m_kernel;
        
    public:
        /**
         * The only constructor.
         * Creates a kernel iterator. Should only be created from the kernel.
         */
        kernel_iterator(kernel_type const * const kernel, unsigned int const & key) : m_key(key), m_keys(nullptr)
        {
            m_kernel = const_cast<kernel_type*>(kernel);
        }

        /**
        * The only constructor.
        * Creates a kernel iterator. Should only be created from the kernel.
        */
        kernel_iterator(kernel_type const * const kernel, std::vector<key_type> *keys) : m_key(0), m_keys(keys)
        {
            m_kernel = const_cast<kernel_type*>(kernel);
        }
        
    public:
        /**
         * Converts the iterator to a handle or key.
         */
        key_type     key()    const {
            if (m_keys){
                return m_keys->at(m_key);
            }
            return key_type{m_key};
        }
        /**
         * Returns a pointer to the kernel that the iterator is bound to.
         */
        kernel_type* kernel() const { return m_kernel; }
        
        /**
         * Compares two iterators for equality.
         */
        friend bool operator==(const iterator& i, const iterator& j)
        {
            return (i.m_key == j.m_key) && (i.m_kernel == j.m_kernel) && (i.m_keys == j.m_keys);
        }
        /**
         * Compares two iterators for inequality.
         */
        friend bool operator!=(const iterator& i, const iterator& j)
        {
            return (i.m_key != j.m_key) || (i.m_kernel != j.m_kernel) || (i.m_keys != j.m_keys);
        }
        
        /**
         * The member access operator.
         *
         * @return Pointer to the element contained within the kernel cell.
         */
        value_type* operator->()
        {
            if (m_keys){
                return &m_kernel->m_data[(int)m_keys->at(m_key)].value;
            }
            assert(m_kernel->lookup(key_type{m_key}).state == element_type::VALID);
            return &m_kernel->m_data[m_key].value;
        }
        
        /**
         * The dereference operator.
         *
         * @return The element that is contained within the kernel cell.
         */
        kernel_iterator_value<kernel_t_> operator*()
        {
            if (m_keys){
                return kernel_iterator_value<kernel_t_>(m_keys->at(m_key), m_kernel);
            }
            return kernel_iterator_value<kernel_t_>(key_type{m_key}, m_kernel);
        }
        
        /**
         * The pre-increment operator.
         *
         * @return The next iterator in the kernel sequence.
         */
        iterator& operator++()
        {
            if (m_keys){
                m_key++;
                if (m_key>= m_keys->size()){
                    m_key = m_kernel->m_data.size();
                }
                return *this;
            }
            do {
                m_key++;
            } while (m_key < m_kernel->m_data.size() && m_kernel->m_data[m_key].state != element_type::VALID);
            return *this;
        }
        
        /**
         * The post-increment operator.
         *
         * @return The iterator
         */
        iterator operator++(int)
        {
            iterator tmp = *this;
            ++*this;
            return tmp;
        }
    };
}
