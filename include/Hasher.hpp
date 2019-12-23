template <typename int_t> class Hasher
{
  
private: int_t m_size;
private: int_t*__restrict__ m_link;

public: static constexpr const int_t s_default_hash = - std::numeric_limits<int_t>::max();
  
public: inline Hasher(int_t size_) noexcept
  : m_size(size_),
    m_link(new int_t[size_])
  {
    this->Reset();
  };
  
public: inline ~Hasher() noexcept
  {
    if (nullptr != this->m_link)
      {
	delete[] this->m_link;
	this->m_link = nullptr;
      }
  };
  
public: inline int_t CalculateMaxSize(int_t*		array_)
  {
    int_t mx = 0;    
    for (int_t i=this->m_size;i>0;--i)
      {
	int_t j = i-1;
	int_t lmx = 0;
	for (int_t h = this->m_link[j];h != s_default_hash;h = array_[-(h+1)])
	  {
	    ++lmx;
	  }
	mx = (mx < lmx) ? lmx : mx;
      }
    
    return mx;
  };

  
public: inline const int_t&operator[](int_t hash_value_) const noexcept
  {
    return this->m_link[hash_value_];
  };
  
public: inline int_t&operator[](int_t hash_value_) noexcept
  {
    return this->m_link[hash_value_];
  };  
  
public: inline void Reset() noexcept
  {
    for (int_t i=0;i<this->m_size;++i)
      {
	this->m_link[i] = s_default_hash;
      }
  };
  
};
