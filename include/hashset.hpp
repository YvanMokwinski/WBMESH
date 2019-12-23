#pragma once
#include <unordered_set>

template <typename K> class hashset : public std::unordered_set<K>
{
public:
  using it_t = typename std::unordered_set< K>::iterator;
  
  inline bool has(K key_) 
  {
    it_t it = this->find(key_);
    return it!=this->end();
  };
  
  
  inline bool add(K key_)
  {
    auto p = this->insert(key_);
    return p.second;
  };
  
  inline bool remove(K key_)
  {
    it_t it 	= this->find(key_);
    if (it != this->end())
      {
	this->erase(it);
	return true;
      }
    else
      {
	return false;
      }
  };
  
};
