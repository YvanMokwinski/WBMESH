#pragma once

struct triangle_t
{

private:  int_t m_edges[3] {};
private:  int   m_ways[3] {};
private:  int_t m_id {};

public:  inline triangle_t() noexcept;  
public:  inline triangle_t(int_t id_,
			   int_t edges_[3],
			   int ways_[3]) noexcept;
  
public:  inline int_t 	id	() const noexcept;
public:  inline int_t 	edge	(unsigned int i) const noexcept;
public:  inline int_t 	way	(unsigned int i) const noexcept;
public:  inline void 	id	(int_t id) noexcept;
public:  inline void 	edge	(unsigned int i,
				 int_t j )  noexcept;
public:  inline void 	way	(unsigned int i,
				 int w)  noexcept;  

public:  inline int  	get_local_edge_index(int_t edge) const noexcept;  
};

inline triangle_t::triangle_t() noexcept
{
};
  
inline triangle_t::triangle_t(int_t id_,
			      int_t edges_[3],
			      int ways_[3]) noexcept
  : m_id(id_)
{
  this->m_edges[0] = edges_[0];
  this->m_edges[1] = edges_[1];
  this->m_edges[2] = edges_[2];

  this->m_ways[0] = ways_[0];
  this->m_ways[1] = ways_[1];
  this->m_ways[2] = ways_[2];
};
  
inline int_t 	triangle_t::id() const noexcept { return this->m_id; };
inline int_t 	triangle_t::edge(unsigned int i) const noexcept { return this->m_edges[i]; };
inline int_t 	triangle_t::way(unsigned int i) const noexcept { return this->m_ways[i]; };

inline void 	triangle_t::id(int_t id) noexcept { this->m_id = id; };
inline void 	triangle_t::edge(unsigned int i,int_t j )  noexcept { this->m_edges[i] = j; };
inline void 	triangle_t::way(unsigned int i,int w)  noexcept { this->m_ways[i] = w; };
inline int 	triangle_t::get_local_edge_index(int_t edge) const noexcept
{
  for (int i=0;i<3;++i)
    {
      if (edge == this->m_edges[i])
	{
	  return i;
	}
    }
  return -1;
};
  
