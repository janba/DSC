
#pragma once

namespace is_mesh
{
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
        using element_type = typename kernel_t_::kernel_element;
    public:
        using kernel_type = kernel_t_;
        using iterator = kernel_iterator<kernel_type>;
        using kernel_element = typename kernel_type::kernel_element;
        using value_type = typename element_type::value_type;
        using key_type = typename kernel_t_::kernel_key_type;
        
    private:
        unsigned int     m_key;
        kernel_type*     m_kernel;
        
    public:
        /**
         * The only constructor.
         * Creates a kernel iterator. Should only be created from the kernel.
         */
        kernel_iterator(kernel_type const * const kernel, unsigned int const & key) : m_key(key)
        {
            m_kernel = const_cast<kernel_type*>(kernel);
        }

        
    public:
        /**
         * Converts the iterator to a handle or key.
         */
        key_type     key()    const {
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
            return (i.m_key == j.m_key) && (i.m_kernel == j.m_kernel);
        }
        /**
         * Compares two iterators for inequality.
         */
        friend bool operator!=(const iterator& i, const iterator& j)
        {
            return (i.m_key != j.m_key) || (i.m_kernel != j.m_kernel);
        }
        
        /**
         * The member access operator.
         *
         * @return Pointer to the element contained within the kernel cell.
         */
        value_type* operator->()
        {
            assert(m_kernel->lookup(key_type{m_key}).state == element_type::VALID);
            return &m_kernel->m_data[m_key].value;
        }
        
        /**
         * The dereference operator.
         *
         * @return The element that is contained within the kernel cell.
         */
        value_type& operator*()
        {
            assert(m_kernel->lookup(key_type{m_key}).state == element_type::VALID);
            return m_kernel->m_data[m_key].value;
        }
        
        /**
         * The pre-increment operator.
         *
         * @return The next iterator in the kernel sequence.
         */
        iterator& operator++()
        {
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
