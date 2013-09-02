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

/** Abstract class for HMesh entity attribute vectors. This class is used for storing all attributes
 associated with any type of mesh entity - be it Vertex, HalfEdge, or Face. Also the position attribute
 is stored in an AttributeVector class.
 */
template<typename ITEM, typename ITEMID>
class AttributeVector
{
protected:
    /// Construct from optional size and item (size should be identical to number of entities in associated container
    AttributeVector(size_t _size = 0, ITEM item = ITEM()) : items(_size, item)
    {
        
    }
    
public:
    /// const reference to item given by ID
    const ITEM& operator [](ITEMID id) const
    {
        assert(id < items.size());
        return items[id];
    }
    
    /// reference to item given by ID
    ITEM& operator [](ITEMID id)
    {
        if(id >= items.size())
        {
            items.resize(id + 1);
        }
        return items[id];
    }
    
    /// resize the vector (may be necessary if associated container size grows)
    void resize(size_t _size, ITEM item = ITEM())
    {
        items.resize(_size, item);
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
    
    /// clenup unused items from the vector, given by remap from associated container
    void cleanup(const std::map<ITEMID, ITEMID>& map)
    {
        std::map<ITEMID, ITEM> new_items(map.size());
        for(auto it : map)
        {
            assert(it->second.index < map.size());
            new_items[it->second.index] = items[it->first.index];
        }
        std::swap(items, new_items);
    }
    
private:
    std::map<ITEMID, ITEM> items;
};

#endif
