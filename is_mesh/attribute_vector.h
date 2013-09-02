//
//  attribute_vector.h
//  DSC
//
//  Created by Asger Nyman Christiansen on 8/20/13.
//  Copyright (c) 2013 Asger Nyman Christiansen. All rights reserved.
//

#ifndef DSC_attribute_vector_h
#define DSC_attribute_vector_h

#include <cassert>
#include <vector>
#include <map>

/** Abstract class for attribute vectors. This class is used for storing all attributes
 associated with any type of simplicial complex entity - be it vertex, edge, face or tetrahedron. Also the position attribute
 is stored in an AttributeVector class.
 */
template<typename ITEM, typename ITEMID>
class AttributeVector
{
protected:
    /// Construct from optional size and item
    AttributeVector(size_t _size = 0, ITEM item = ITEM())
    {
        
    }
    
public:
    /// const reference to item given by ID
    const ITEM& operator [](ITEMID id) const
    {
        return items.at(id);
    }
    
    /// reference to item given by ID
    ITEM& operator [](ITEMID id)
    {
        return items[id];
    }
    
    /// number of attribute items in vector
    size_t size() const
    {
        return items.size();
    }
    
    /// clear the vector
    void clear()
    {
        items.clear();
    }
    
    /// clenup unused items from the vector deleted
    void cleanup(const std::vector<ITEMID>& deleted)
    {
        for(auto it : deleted)
        {
            items.erase(*it);
        }
    }
    
private:
    std::map<ITEMID, ITEM> items;
};

#endif
