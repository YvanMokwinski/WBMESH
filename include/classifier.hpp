#pragma once
#include "hashset.hpp"

class classifier : public hashset<int_t>
{
public: 
  template <typename criterion_t>
  void push(int_t edge, criterion_t criterion)
  {
    if (criterion(edge))
      {		
	this->add(edge);
      }
  };
  
  void pop(int_t edge)
  {
    this->remove(edge);      
  };
  
};
