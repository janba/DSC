//
//  cache.hpp
//  interface_tracking
//
//  Created by Tuan Nguyen Trung on 11/7/16.
//  Copyright © 2016 Tuan Nguyen Trung. All rights reserved.
//

#ifndef cache_hpp
#define cache_hpp

#include <stdio.h>
#include <vector>
#include <map>
#include <iostream>
#include <queue>

#define BUFFER_SIZE 1.2 // Buffer of cache preserved for adding entities

#define DSC_CACHE

enum key_type{
    key_vertex_type = 0,
    key_edge_type = 1,
    key_face_type = 2,
    key_tet_type = 3
};

class cache_item_base
{
protected:
    key_type _key_type;
public:
//    cache_item_base(key_type type){type = _key_type;};
//    ~cache_item_base(){};
    
    const key_type &  get_key_type(){return _key_type;}
    
//    virtual void release(){};
    virtual void * get(size_t idx){return nullptr;};
    virtual void set(void * newdata, size_t idx){};
    virtual void mark_dirty(size_t idx){};
    virtual void resize(size_t new_size){};
};

template<typename type>
class cache_item : public cache_item_base
{
    std::vector<type *> _data;
    
public:
    cache_item(key_type ktype, size_t length){
        _key_type = ktype;
        _data = std::vector<type *> (length, nullptr);
    }
    
    ~cache_item(){
        for(auto & d : _data){
            if(d){
                delete d;
                d = nullptr;
            }
        }
    }
    
    virtual void * get(size_t idx){
        assert(idx < _data.size());
        return _data[idx];
    }
    
    virtual void set(void * newdata, size_t idx){
        assert(idx < _data.size());
        
        if(_data[idx])
        {
            delete _data[idx];
            _data[idx] = nullptr;
        }
        _data[idx] = (type*)newdata;
    };
    
    virtual void mark_dirty(size_t idx){
        if(_data[idx]){
            delete _data[idx];
            _data[idx] = nullptr;
        }
    }
    
    virtual void resize(size_t new_size){
        for (size_t i = new_size; i < _data.size(); i++)
        {
            if(_data[i])
                delete _data[i];
        }
        
        _data.resize(new_size);
    }
};

class cache
{
    std::map<size_t, cache_item_base*> _data;
    
    size_t no_element[4];
public:
    // Get cache
    template<typename type, key_type ktype>
    type* get_cache(size_t cache_id, size_t item_id, std::function<type*(size_t)>compute)
    {
        // create if not exist?
        auto map_cache = _data.find(cache_id);
        if (map_cache == _data.end())
        {
            cache_item<type> * new_array = new cache_item<type>(ktype, no_element[ktype]);
            _data.insert(std::make_pair(cache_id, (cache_item_base*)new_array));
        }
        
        // Compute if not cached
        auto cached_array = _data.find(cache_id)->second;
        if(!cached_array->get(item_id))
        {
            cached_array->set((void*)compute(item_id), item_id);
        }
        
        return (type*)cached_array->get(item_id);
    }
    
    // Update cache size
    void resize(key_type ktype, size_t new_size)
    {
        no_element[ktype] = new_size*BUFFER_SIZE;
        for (auto & c : _data)
        {
            if (c.second->get_key_type() == ktype)
            {
                c.second->resize(no_element[ktype]);
            }
        }
    }
    
    void mark_dirty(key_type ktype, size_t item_id)
    {
        for(auto it = _data.begin(); it != _data.end(); it++)
        {
            auto cache_array = it->second;
            if (cache_array->get_key_type() == ktype)
            {
                cache_array->mark_dirty(item_id);
            }
        }
    }

    
    static cache * get_instance(){
        static cache instance;
        return &instance;
    }
    
private:
    cache(){}
    ~cache(){
        for(auto it = _data.begin(); it != _data.end(); it++)
        {
            delete it->second;
        }
        _data.clear();
    }
};



#endif /* cache_hpp */
