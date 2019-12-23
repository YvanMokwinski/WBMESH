#pragma once
#include "curve_t.hpp"

struct delaunay_cad_t
{
private:
  int_t   m_npoints{0};
  int_t   m_nedges{0};
  int_t * m_e2n{nullptr};
  xy_t *  m_points{nullptr};
public:
  inline delaunay_cad_t() noexcept
  {
  };

  inline delaunay_cad_t(int_t npoints_,int_t nedges_) noexcept
    : m_npoints(npoints_),
      m_nedges(nedges_)
  {
    this->m_e2n 	= new int_t[3*(nedges_+1)];
    this->m_points 	= new xy_t[npoints_+1];    
  };

 static int comp2(const void * a,const void * b)
{
  const int_t * a_ = (const int_t * )a;
  const int_t * b_ = (const int_t * )b;
  if (a_[2] < b_[2]) return -1;
  else
    if (a_[2] > b_[2]) return 1;
  else
    return 0;
}

  inline void analysis()
  {
    qsort(m_e2n+3,m_nedges,3*sizeof(int_t), comp2);
  };
  
  inline void init(int_t npoints_,int_t nedges_) noexcept
  {
    this->m_npoints=npoints_;
    this->m_nedges=nedges_;
    this->m_e2n = new int_t[3*(nedges_+1)];
    this->m_points = new xy_t[npoints_+1];    
  };
  
  inline ~delaunay_cad_t() noexcept
  {
    if (nullptr != this->m_e2n)
      {
	delete [] this->m_e2n;
      }
  };
  
  inline void box(double pmin[],double pmax[]) const
  {
    double xmin = 1e+30,ymin=1e+30,xmax=0.0,ymax=0.0;
    for (int_t i=1;i<=m_npoints;++i)
      {
	double x = m_points[i].x;
	double y = m_points[i].y;
	xmin = std::min(xmin,x);
	ymin = std::min(ymin,y);
	xmax = std::max(xmax,x);
	ymax = std::max(ymax,y);
      }
    pmin[0] = xmin - (xmax-xmin)/8.0;
    pmin[1] = ymin  - (xmax-xmin)/8.0;
    pmax[0] = xmax + (xmax-xmin)/8.0;
    pmax[1] = ymax+ (xmax-xmin)/8.0;
  };
  
  inline int_t node(int_t edge_,unsigned int i) const noexcept { return this->m_e2n[3*edge_+i];}
  inline void  node(int_t edge_,unsigned int i,int_t node_) noexcept { this->m_e2n[3*edge_+i] = node_;}
  
  inline void get_point(int_t node_,double coo[]) const noexcept
  {
    coo[0] = this->m_points[node_].x;
    coo[1] = this->m_points[node_].y;
  };

  inline void set_point(int_t node_,const double coo[]) noexcept
  {
    this->m_points[node_].x = coo[0];
    this->m_points[node_].y = coo[1];
  };
  

  inline int_t 		npoints() const noexcept {return this->m_npoints;}
  inline int_t 		nedges() const noexcept {return this->m_nedges;}
};
