#pragma once
template <typename tObject> class HashCavityTable
{
  tObject *m_objects;
  int_t * m_i;
  int_t * m_j;
  int_t * m_link;
  int_t * m_first;
  int_t m_n;
  int_t m_mx;
#define STOP_VALUE -10000000
public:    HashCavityTable()
  {
    //  std::cout << "hello"  << std::endl;
    m_mx = 64000*2;
    m_objects = (tObject*)malloc(sizeof(tObject)*m_mx);
    m_i = (int_t*)malloc(sizeof(int_t)*m_mx);
    m_j = (int_t*)malloc(sizeof(int_t)*m_mx);
    m_link = (int_t*)malloc(sizeof(int_t)*m_mx);
    m_first = (int_t*)malloc(sizeof(int_t)*2048);
    for (int i=0;i<2048;++i)
      {
	m_first[i] = STOP_VALUE;
      }
  };

public:
  inline void clear()
  {
    m_n = 0;
    for (int i=0;i<2048;++i)
      {
	m_first[i] = STOP_VALUE;
      }
    for (int i=0;i<m_n;++i)
      {
	m_link[i] = STOP_VALUE;
      }
  };
  inline tObject Get(const int_t  		i1_,
		     const int_t  		i2_) const 
  {
    int_t  i1 = std::min(i1_,i2_);
    int_t  i2 = std::max(i1_,i2_);      
    int_t  hash = (i1%2048);

    if (m_first[hash] != STOP_VALUE)
      {
	for (int_t k = m_first[hash];k != STOP_VALUE;k = m_link[-k-1])
	  {	  
	    int_t i = m_i[-k-1];
	    int_t j = m_j[-k-1];
	    //	      std::cout << "comparing  " << i << " " << j << " " << i1 << " " << i2 << std::endl;
	    if ( ((i==i1)&&(j==i2)) ||
		 ((j==i1)&&(i==i2)) )
	      {
		//  std::cout << "found " << std::endl;
		return m_objects[-k-1];
	      }
	  }
	return 0;
      }
    else
      {
	return 0;
      }
  };


  void print()
  {
#if 0
    std::cout << "----"  << std::endl;
    for (int i=0;i<m_n;++i)
      {
	std::cout << "m_link[" << i << "] " << m_link[i] << std::endl;
      }
    for (int i=0;i<m_n;++i)
      {
	std::cout << "m_i[" << i << "] " << m_i[i] << std::endl;
      }
    for (int i=0;i<m_n;++i)
      {
	std::cout << "m_j[" << i << "] " << m_j[i] << std::endl;
      }
    for (int i=0;i<32;++i)
      {
	std::cout << "m_first[" << i << "] " << m_first[i] << std::endl;
      }
#endif      
  }
    
  inline bool Add(const int_t 		i1_,
		  const int_t 		i2_,
		  tObject 			object_)
  {
    int_t  i1 = std::min(i1_,i2_);
    int_t  i2 = std::max(i1_,i2_);      
    int_t  hash = (i1%2048);
    //      std::cout << "hash " << i1 << " " << i2 << " " << hash << std::endl;
    if (m_first[hash] != STOP_VALUE)
      {
	for (int_t k = m_first[hash];k != STOP_VALUE;k = m_link[-k-1])
	  {	  
	    // std::cout << "k " <<k << std::endl;
	    int_t i = m_i[-k-1];
	    int_t j = m_j[-k-1];
	    if ( ((i==i1)&&(j==i2)) ||
		 ((j==i1)&&(i==i2)) )
	      {
		return false;
	      }
	  }
	m_i[m_n] = i1;
	m_j[m_n] = i2;
	m_objects[m_n] = object_;
	m_link[m_n] = m_first[hash];
	m_first[hash] = -(++m_n);
	//	  std::cout<<"rokrokorkrorko" << std::endl;
	print();
	return true;
      }
    else
      {
	m_i[m_n] = i1;
	m_j[m_n] = i2;
	m_objects[m_n] = object_;
	m_link[m_n] = m_first[hash];
	m_first[hash] = -(++m_n);
	//	  std::cout<<"rokrokorkrorko" << std::endl;
	print();
	return true;
      }
  };
  
  inline bool Remove(const int_t	i1_,
		     const int_t	i2_)
  {
    std::cout << "remove " << std::endl;
    exit(1);
    int_t  i1 = std::min(i1_,i2_);
    int_t  i2 = std::max(i1_,i2_);      
    int_t  hash = (i1%2048);

    if (m_first[hash] != STOP_VALUE)
      {
	int_t previous = m_first[hash];
	for (int_t k = m_first[hash];k != STOP_VALUE;k = m_link[-k-1])
	  {	  
	    int_t i = m_i[-k-1];
	    int_t j = m_j[-k-1];
	    if ( ((i==i1)&&(j==i2)) ||
		 ((j==i1)&&(i==i2)) )
	      {
		if (k!=previous)
		  {
		    m_link[-previous-1] = m_link[-k-1];
		  }
		else
		  {
		    m_first[hash] = m_link[-k-1];
		  }
		return true;
	      }
	    previous = k;
	  }
	return false;
      }
    else
      {
	return false;
      }
  };

};
