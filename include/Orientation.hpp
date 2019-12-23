#pragma once


struct Orientation
{
  typedef enum
    {
      BACKWARD=0,
      FORWARD,
      NEUTRAL,
    } Kind;    
 
    
  template <typename _float_type> inline static Kind Create(const _float_type& value_)
  {
    return (value_<Math<_float_type>::Zero) ? BACKWARD : ( (value_>Math<_float_type>::Zero) ? FORWARD : NEUTRAL  );
  };   

    
  inline static Kind Opposite(const Kind&orientation_)
  {
    switch(orientation_)
      {
      case BACKWARD:
	{
	  return FORWARD;
	}
      case FORWARD:
	{
	  return BACKWARD;
	}
      case NEUTRAL:
	{
	  return NEUTRAL;
	}
      }
    std::cerr << "invalid switch" << std::endl;
    exit(1);
    return BACKWARD;
  };
    
    
};

namespace std
{
  std::ostream& operator<< (std::ostream &out_, 
			    const Orientation::Kind &value_)    
  {
    switch(value_)
      {
      case Orientation::BACKWARD:
	{
	  out_<<"BACKWARD";
	  return out_;
	}
      case Orientation::FORWARD:
	{
	  out_<<"FORWARD";
	  return out_;
	}
      case Orientation::NEUTRAL:
	{
	  out_<<"NEUTRAL";
	  return out_;
	}
      }
    std::cerr << "invalid switch" << std::endl;
    exit(1);

    return out_;      
  };
};
