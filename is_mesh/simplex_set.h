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

#include <vector>
#include <cassert>
#include "key.h"

namespace is_mesh
{
    
    template<typename key_type>
    class SimplexSet
    {
        std::vector<key_type> set;
        
    public:
        
        SimplexSet() : set()
        {
            
        }
        
        SimplexSet(std::initializer_list<key_type> il) : set(il)
        {
            
        }
        
        SimplexSet(const SimplexSet& ss) : set(ss.set)
        {
            
        }
        
        SimplexSet& operator=(const SimplexSet& ss)
        {
            set = std::vector<key_type>(ss.set);
            return *this;
        }
        
        SimplexSet(SimplexSet&& ss) : set(std::move(ss.set))
        {
            
        }
        
        SimplexSet& operator=(SimplexSet&& ss)
        {
            if (this != &ss){
                set = std::move(ss.set);
            }
            return *this;
        }
        
        ~SimplexSet()
        {
            
        }

        typename std::vector<key_type>::const_iterator begin() const
        {
            return set.begin();
        }
        
        typename std::vector<key_type>::const_iterator end() const
        {
            return set.end();
        }
        
        unsigned int size() const
        {
            return static_cast<unsigned int>(set.size());
        }
        
        const key_type& front() const
        {
            assert(set.size() > 0);
            return set.front();
        }
        
        const key_type& back() const
        {
            assert(set.size() > 0);
            return set.back();
        }
        
        const key_type& operator[](unsigned int i) const
        {
            assert(size() > i);
            return set[i];
        }
        
        bool contains(const key_type& k) const
        {
            return std::find(set.begin(), set.end(), k) != end();
        }
        
        int index(const key_type& k) const
        {
            for (int i = 0; i < set.size(); i++) {
                if(set[i] == k)
                {
                    return i;
                }
            }
            return -1;
        }
        
        void push_front(const key_type& k)
        {
            assert(!contains(k));
            set.insert(set.begin(), k);
        }
        
        void push_back(const key_type& k)
        {
            assert(!contains(k));
            set.push_back(k);
        }
        
        void push_back(key_type&& k)
        {
            assert(!contains(k));
            set.push_back(std::move(k));
        }
        
        void swap(unsigned int i = 0, unsigned int j = 1)
        {
            assert(size() > i);
            assert(size() > j);
            std::swap(set[i], set[j]);
        }
        
        SimplexSet<key_type>& operator+=(const SimplexSet<key_type>& ss)
        {
            for (const key_type& k : ss) {
                *this += k;
            }
            return *this;
        }
        
        SimplexSet<key_type>& operator+=(SimplexSet<key_type>&& ss)
        {
            for (key_type& k : ss.set) {
                *this += std::move(k);
            }
            return *this;
        }
        
        SimplexSet<key_type>& operator+=(const key_type& key)
        {
            if(!contains(key))
            {
                set.push_back(key);
            }
            return *this;
        }
        
        SimplexSet<key_type>& operator+=(key_type&& key)
        {
            if(!contains(key))
            {
                set.push_back(std::move(key));
            }
            return *this;
        }
        
        SimplexSet<key_type>& operator-=(const SimplexSet<key_type>& set)
        {
            for (auto &k : set) {
                *this -= k;
            }
            return *this;
        }
        
        SimplexSet<key_type>& operator-=(const key_type& key)
        {
            auto iter = std::find(begin(), end(), key);
            if(iter != end())
            {
                set.erase(iter);
            }
            return *this;
        }

        bool operator==(const key_type& key)
        {
            if (set.size() != key.set.size()){
                return false;
            }
            for (auto k : set){
                if (!key.contains(k)){
                    return false;
                }
            }

            return true;
        }
    };
    
    template<typename key_type>
    bool operator==(const SimplexSet<key_type>& A, const SimplexSet<key_type>& B)
    {
        if(A.size() == B.size())
        {
            for (auto k : A)
            {
                if(!B.contains(k))
                {
                    return false;
                }
            }
            return true;
        }
        return false;
    }
    
    /**
     *  Returns the union of the two sets A and B.
     */
    template<typename key_type>
    SimplexSet<key_type> operator+(const SimplexSet<key_type>& A, const SimplexSet<key_type>& B)
    {
        SimplexSet<key_type> C = A;
        return C += B;
    }
    
    /**
     *  Returns the union of the two sets A and B.
     */
    template<typename key_type>
    SimplexSet<key_type>&& operator+(SimplexSet<key_type>&& A, const SimplexSet<key_type>& B)
    {
        return std::move(A += B);
    }
    
    /**
     *  Returns the union of the two sets A and B.
     */
    template<typename key_type>
    SimplexSet<key_type>&& operator+(SimplexSet<key_type>&& A, SimplexSet<key_type>&& B)
    {
        return std::move(A) + B;
    }

    /**
     *  Returns the union of the two sets A and B.
     */
    template<typename key_type>
    SimplexSet<key_type>&& operator+(const SimplexSet<key_type>& A, SimplexSet<key_type>&& B)
    {
        return std::move(B) + A;
    }
    
    /**
     *  Returns set A without the elements in set B.
     */
    template<typename key_type>
    SimplexSet<key_type> operator-(const SimplexSet<key_type>& A, const SimplexSet<key_type>& B)
    {
        SimplexSet<key_type> C = A;
        return C -= B;
    }
    
    /**
     *  Returns set A without the elements in set B.
     */
    template<typename key_type>
    SimplexSet<key_type>&& operator-(SimplexSet<key_type>&& A, const SimplexSet<key_type>& B)
    {
        return std::move(A -= B);
    }
    
    /**
     *  Returns set A without the element key.
     */
    template<typename key_type>
    SimplexSet<key_type> operator-(const SimplexSet<key_type>& A, const key_type& key)
    {
        SimplexSet<key_type> B = {key};
        return A - B;
    }
    
    /**
     *  Returns set A without the element key.
     */
    template<typename key_type>
    SimplexSet<key_type>&& operator-(SimplexSet<key_type>&& A, const key_type& key)
    {
        return std::move(A -= key);
    }
    
    /**
     *  Returns the intersection of sets A and B.
     */
    template<typename key_type>
    SimplexSet<key_type> operator&(const SimplexSet<key_type>& A, const SimplexSet<key_type>& B)
    {
        SimplexSet<key_type> C = A;
        return C - (A - B);
    }
    
    /**
     *  Returns the intersection of sets A and B.
     */
    template<typename key_type>
    SimplexSet<key_type>&& operator&(SimplexSet<key_type>&& A, const SimplexSet<key_type>& B)
    {
        return std::move(std::move(A) - (A - B));
    }
}
