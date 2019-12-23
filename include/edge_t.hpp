#pragma once
//!
//! @brief Struct storing an edge topology.
//!
struct edge_t
{
  //!
  //! @brief Empty constructor.
  //!  
public: inline edge_t() noexcept;

  //!
  //! @brief Constructor.
  //!  
public: inline edge_t(int_t i,
		      int_t j,
		      int_t id_) noexcept;
  
  //!
  //! @brief Constructor.
  //!  
public: inline edge_t(int_t nodes_[2],
		      int_t triangles_[2],
		      int_t id_) noexcept;
  
  //!
  //! @brief Returns is boundary.
  //!  
public: inline bool		is_boundary() 		const 	noexcept;
  //!
  //! @brief Returns is non flipable.
  //!  
public: inline bool		is_nonflipable() 	const 	noexcept ;
  //!
  //! @brief Returns the starting node.
  //!  
public: inline int_t 		start		(int way_ = 1) 		const noexcept;
  //!
  //! @brief Returns the stopping node.
  //!  
public: inline int_t		stop		(int way_ = 1)		const noexcept;
  //!
  //! @brief Returns the adjacent triangle.
  //!  
public: inline int_t		triangle	(unsigned int lindex_) 	const noexcept;
  //!
  //! @brief Returns the other adjacent triangle.
  //!  
public: inline int_t 		otriangle	(int_t triangle_)  	const noexcept;
  
  //!
  //! @brief Set is non flipable.
  //!  
public: inline void		is_nonflipable	(bool value)	noexcept;
  //!
  //! @brief Set is boundary.
  //!  
public: inline void 		is_boundary	(bool value) 	noexcept ;
  //!
  //! @brief Set the adjacent triangle.
  //!  
public: inline void		triangle	(unsigned int 	lindex_,
						 int_t 		triangle_) noexcept;

  //!
  //! @brief Deref the adjacent triangle.
  //!  
public: inline void		dereftriangle	(int_t triangle_) noexcept;  
  //!
  //! @brief Reference  adjacent triangle.
  //!  
public: inline void 		reftriangle	(int_t triangle_) noexcept;

  //!
  //! @brief Reference  adjacent triangle.
  //!  
public: inline void SetStart	(int_t point_) noexcept;

  //!
  //! @brief Reference  adjacent triangle.
  //!  
public: inline void SetEnd	(int_t point_) noexcept;

private:  int_t m_id {};
private:  int_t m_s  {};
private:  int_t m_e  {};
private:  int_t m_lt {};
private:  int_t m_rt {};
private:  bool 	m_isBoundary  {};
private:  bool	m_isUnflipable{};

};




inline edge_t::edge_t() noexcept {};
inline edge_t::edge_t(int_t i,
		      int_t j,
		      int_t id_) noexcept
  : m_s(i),m_e(j),m_lt(0),m_rt(0),m_id(id_)
{
};

inline edge_t::edge_t(int_t nodes_[2],
		      int_t triangles_[2],
		      int_t id_) noexcept
  : edge_t(nodes_[0],nodes_[1],id_)
{
  this->m_lt = triangles_[0];
  this->m_rt = triangles_[1];
};

inline bool	edge_t::is_boundary() 		const 	noexcept 	{ return this->m_isBoundary; };  
inline bool	edge_t::is_nonflipable() 	const 	noexcept 	{ return this->m_isUnflipable; };
inline int_t	edge_t::start		(int way_) 		const noexcept  { return (way_ > 0) ? this->m_s : this->m_e;};  
inline int_t	edge_t::stop		(int way_)		const noexcept  { return (way_ > 0) ? this->m_e : this->m_s;  };
inline int_t	edge_t::triangle	(unsigned int lindex_) 	const noexcept { return lindex_ == 0 ? this->m_lt : this->m_rt; };  
inline int_t	edge_t::otriangle	(int_t triangle_)  	const noexcept { return (this->m_lt == triangle_) ? this->m_rt : this->m_lt; };
inline void	edge_t::is_nonflipable	(bool value)	noexcept 	{ this->m_isUnflipable = value; };  
inline void	edge_t::is_boundary	(bool value) 	noexcept  	{ this->m_isBoundary = value; };
inline void	edge_t::triangle	(unsigned int 	lindex_, int_t	triangle_) noexcept
{
  if (lindex_==0)
    {
      this->m_lt = triangle_;
    }
  else
    {
      this->m_rt = triangle_;
    }
};

//!
//! @brief Deref the adjacent triangle.
//!  
inline void	edge_t::dereftriangle	(int_t triangle_) noexcept
  {
    if (this->m_lt == triangle_)
      {
	this->m_lt = this->m_rt;
	this->m_rt = 0;
      }
    else
      {
	this->m_rt = 0;
      }
  };
  
  //!
  //! @brief Reference  adjacent triangle.
  //!  
inline void edge_t::reftriangle(int_t triangle_) noexcept
  {
    if (this->m_lt)
      {
	this->m_rt = triangle_;
      }
    else
      {
	this->m_lt = triangle_;
      }
  };
inline void edge_t::SetStart	(int_t point_) noexcept	{ this->m_s = point_; };
inline void edge_t::SetEnd	(int_t point_) noexcept	{ this->m_e = point_; };

