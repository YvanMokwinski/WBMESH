#pragma once

struct metric_t
{
  double m[3]{};
public: metric_t(const double*m_) : m{m_[0],m_[1],m_[2]}
  {
  };
  
public:metric_t(double a,double b,double c) : m{a,b,c}
  {
  };
  
public: inline double size(double hx,double hy) const noexcept
  {
    return hx*hx*m[0] + 2.0 * m[1]*hx*hy + m[2]*hy*hy;
  };
  
public: inline double size(double px,
			   double py,
			   double qx,
			   double qy) const 
  {
    return this->size(qx-px,qy-py);
  };
  
  
public: inline static void interpolate(double t,const metric_t * metric_start,const metric_t*metric_end,metric_t*metric) 
  {
    auto mx0 = (1.0 - t) * 1.0/sqrt(metric_start->m[0]) + t * 1.0/sqrt(metric_end->m[0]);
    auto mx1 = (1.0 - t) * 1.0/sqrt(metric_start->m[2]) + t * 1.0/sqrt(metric_end->m[2]);
    metric->m[0]  = 1.0 / (mx0*mx0);
    metric->m[1]  = 0.0;
    metric->m[2]  = 1.0 / (mx1*mx1);
  };
  
};
