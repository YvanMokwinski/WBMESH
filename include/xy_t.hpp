#pragma once

struct xy_t
{
  double x;
  double y;
  int_t  m_id;
  double m_metric[3]{};  

public: inline xy_t(double x_,double y_,int_t id_) noexcept;
public: inline xy_t() noexcept;
public: inline int_t id() const noexcept;
public: inline void  id(int_t id) noexcept;
};

inline xy_t::xy_t(double x_,double y_,int_t id_) noexcept
  : x(x_),
    y(y_),
    m_id(id_) { };
inline xy_t::xy_t() noexcept { };
inline int_t xy_t::id() const noexcept   { return  this->m_id; };
inline void  xy_t::id(int_t i) noexcept  { this->m_id = i; };
