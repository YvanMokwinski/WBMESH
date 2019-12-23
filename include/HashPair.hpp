#pragma once
template <typename _type_pair,typename _object> class HashPair : public std::map< std::pair<_type_pair,_type_pair>, _object>
{
public:
  typedef typename std::pair<_type_pair,_type_pair> MapPairKey;
  typedef typename std::map< MapPairKey, _object>::iterator MapPairIterator;

  inline bool Get(const _type_pair&  		i1_,
		  const _type_pair&  		i2_) 
  {
    const _type_pair  i1 	= Min<_type_pair>(i1_,i2_);
    const _type_pair  i2 	= Max<_type_pair>(i1_,i2_);
    MapPairIterator it 		= this->find(MapPairKey(i1,i2));
    return (it!=this->end());
  };
  

  inline bool Add(const _type_pair& 		i1_,
		  const _type_pair& 		i2_,
		  const _object& 		obj_)
  {
    const _type_pair  i1 	= Min<_type_pair>(i1_,i2_);
    const _type_pair  i2 	= Max<_type_pair>(i1_,i2_);

    std::pair<MapPairIterator,bool> insertBack 	= this->insert(std::pair< MapPairKey, _type_pair >(MapPairKey(i1,i2),obj_));
    return insertBack.second;
  };

  
  inline bool Remove(const _type_pair&	i1_,
		     const _type_pair&	i2_)
  {
    const _type_pair  i1 	= Min<_type_pair>(i1_,i2_);
    const _type_pair  i2 	= Max<_type_pair>(i1_,i2_);
    MapPairIterator 	it 	= this->find(MapPairKey(i1,i2));
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

  
  void print()
  {
    for (MapPairIterator it = this->begin();it!=this->end();++it)
      {
	MapPairKey key = it->first;
	std::cout << "record edge " << key.first << " " << key.second << " " <<  it->second << std::endl;
      }

  };  

};
