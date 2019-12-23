#include <iostream>
#include <sys/time.h>
#include <string>
#include <math.h>
#include <stdlib.h>

struct EllapsedTime
{
 protected:
  struct timeval tval_start;    
  struct timeval tval_stop;    
 public:
  inline void Start()
  {
    gettimeofday (&this->tval_start, NULL);
  };
  inline void Stop()
  {
    gettimeofday (&this->tval_stop, NULL);
  };

  template <typename _float_type> inline _float_type Time() const
  {
    unsigned long int seconds   = this->tval_stop.tv_sec - this->tval_start.tv_sec;
    unsigned long int microseconds; 
    if (this->tval_stop.tv_usec<this->tval_start.tv_usec)
      {
	seconds-=1;
	microseconds = ((1000000 + this->tval_stop.tv_usec) - this->tval_start.tv_usec);
      }
    else
      {
	microseconds = this->tval_stop.tv_usec - this->tval_start.tv_usec;
      }
    return (_float_type(seconds)) + (_float_type(microseconds))*(_float_type(1.0e-6));
  };
  
  inline std::string ToString() const
  {
    unsigned long int seconds   = this->tval_stop.tv_sec - this->tval_start.tv_sec;
    unsigned long int microseconds; 
    if (this->tval_stop.tv_usec<this->tval_start.tv_usec)
      {
	seconds-=1;
	microseconds = ((1000000 + this->tval_stop.tv_usec) - this->tval_start.tv_usec);
      }
    else
      {
	microseconds = this->tval_stop.tv_usec - this->tval_start.tv_usec;
      }
    char str[64];
    sprintf(str,"%lum%2.5fs",seconds/60,(double(seconds%60))+(double(microseconds))*double(1.0e-6));
    return std::string(str);
  };

};

namespace std
{
  std::ostream& operator<< (std::ostream &out_, 
			    const EllapsedTime& ellapsedTime_)    
    {
      out_ << ellapsedTime_.ToString();
      return out_;
    }
};


using int_t = long long int;

#include "delaunay_calculator_t.hpp"


void usage(const char * name)
{
  std::cerr << name << " <input_file> -t <float> -o <output_file> [-b] [-r] [-p] [-d] [-v]" << std::endl;
  std::cerr << " [-b] benchmark" << std::endl;
  std::cerr << " [-p] optimize" << std::endl;
  std::cerr << " [-r] remove bounding box" << std::endl;
  std::cerr << " [-d] debug" << std::endl;
  std::cerr << " [-v] verbose" << std::endl;
}



//
//
// s(x,y) = 0.5*(ax^2 + 2bxy + cy^2)
// ds(x,y)/dx = ax'x + x'by
// ds(x,y)/dy = bxy' + cy'y
//
//
// J0(x,y) = s(x-x0,y-y0) - s(x-x1,y-y1) = 0
// J1(x,y) = s(x-x0,y-y0) - s(x-x2,y-y2) = 0
//
// dJ0/dx = a(x-x0) + b(y-y0) - a(x-x1) - b(y-y1) = a(x1 - x0) + b(y1 - y0)
// dJ0/dy = b(x-x0) + c(y-y0) - b(x-x1) - c(y-y1) = b(x1 - x0) + c(y1 - y0)
//
// dJ1/dx = a(x-x0) + b(y-y0) - a(x-x2) - b(y-y2) = a(x2 - x0) + b(y2 - y0)
// dJ1/dy = b(x-x0) + c(y-y0) - b(x-x2) - c(y-y2) = b(x2 - x0) + c(y2 - y0)
//
//

#if 0
double edge_squared_length(double hx,
			   double hy,
			   double a,
			   double b,
			   double c)
{
  return hx*hx*a + 2.0 * b*hx*hy + c*hy*hy;
}

double delaunay_criterion(double x,
			  double y,
			  double x0,
			  double y0,
			  double x1,
			  double y1,
			  double x2,
			  double y2,
			  double a,
			  double b,
			  double c,
			  double *ccx,
			  double *ccy,
			  double *r)
{
  double x01 = x1 - x0;
  double y01 = y1 - y0;
  double x02 = x2 - x0;
  double y02 = y2 - y0;

  double a00 = a*x01 + b*y01;
  double a01 = b*x01 + c*y01;
  double a10 = a*x02 + b*y02;
  double a11 = b*x02 + c*y02;

  double x_0 = x - x0;
  double y_0 = y - y0;
  double x_1 = x - x1;
  double y_1 = y - y1;
  double x_2 = x - x2;
  double y_2 = y - y2;
  double s0  = (a*x_0*x_0 + 2.0*b*x_0*y_0 + c*y_0*y_0)/2.0;
  double s1  = (a*x_1*x_1 + 2.0*b*x_1*y_1 + c*y_1*y_1)/2.0;
  double s2  = (a*x_2*x_2 + 2.0*b*x_2*y_2 + c*y_2*y_2)/2.0;
  double b0  = -(s0 - s1);
  double b1  = -(s0 - s2);

  //
  // Solve the linear system.
  //
  double d  = a00 * a11 - a10 * a01;
  double hx = (b0 * a11 - b1 * a01)/d;
  double hy = (a00 * b1 - a10 * b0)/d;

  ccx[0] = x + hx;
  ccx[1] = y + hy;
  r[0] = edge_squared_length(x0-ccx[0],y0-ccx[1],a,b,c);
  r[1] = edge_squared_length(hx,hy,a,b,c);
  return r[1] / r[0];
}

double edge_squared_length(const double * xy0,
			   const double * xy1,
			   double metric[])
{
  double vx  = xy1[0] - xy0[0];
  double vy  = xy1[1] - xy0[1];
  double a00 = metric[0];
  double a01 = metric[1];
  double a11 = metric[2];
  return vx * (vx * a00 + vy * a01) + vy * (vx * a01 + vy * a11);
}

double edge_squared_length(const xy_t * x0,
			   const xy_t * x1,
			   double metric[])
{
  return edge_squared_length(&x0->x,&x1->x,metric);
}



int comp(const void * a,const void * b)
{
  const long int * a_ = (const long int * )a;
  const long int * b_ = (const long int * )b;
  if (a_[2] < b_[2]) return -1;
  else
    if (a_[2] > b_[2]) return 1;
  else
    return 0;
}

#endif


bool has_option(const int argc_,
		char ** argv_,
		const char * opt)
{
  for (int i=1;i<argc_;++i)
    {
      if (!strcmp(opt,argv_[i]))
	{
	  return true;
	}
    }
  return false;
}

int main(const int argc_,
	 char ** argv_)
{
  if (argc_ < 6)
    {
      usage(argv_[0]);
      return 1;
    }
  
  const char * ifilename = argv_[1];
  const char * stol = argv_[2];
  double target_size     = atof(argv_[3]);
  const char * sofilename = argv_[4];
  const char * ofilename = argv_[5];
  if (strcmp(stol,"-t"))
    {
      usage(argv_[0]);
      return 1;
    }

  if (strcmp(sofilename,"-o"))
    {
      usage(argv_[0]);
      return 1;
    }

  bool bench 		= has_option(argc_, argv_,"-b");
  bool remove_box 	= has_option(argc_, argv_,"-r");
  bool optimize 	= has_option(argc_, argv_,"-p");
  bool debug 		= has_option(argc_, argv_,"-d");
  bool verbose 		= has_option(argc_, argv_,"-v");

  if (verbose)
    {
      std::cout << "// " << __DATE__ << std::endl;
      std::cout << "// " << argv_[0];
      for (int i=1;i<argc_;++i)
	{
	  std::cout << " " << argv_[i];
	}
      std::cout << std::endl;
    }
  
  delaunay_calculator_t*  delaunay_calculator = new delaunay_calculator_t(ifilename);

  //
  // Set the debug mode.
  //
  delaunay_calculator->debug_mode(debug);
  
  //
  // INITIALIZE THE BOX.
  //
  delaunay_calculator->build_box_mesh();

  if (bench)
    {
      //
      // INSERT CAD POINTS.
      //
      delaunay_calculator->insert_random();
      delaunay_calculator->output(ofilename);      
      return 0;
    }

  delaunay_calculator->insert_cad_points();
  
  //
  // RETRIEVE CAD BOUNDARY.
  //
  delaunay_calculator->retrieve_boundary();

  //
  // REMOVE BOX.
  //
  if (remove_box)
    {
      delaunay_calculator->remove_box(false == remove_box);
    }
  
  //
  // INITIALIZE THE METRIC.
  //
  delaunay_calculator->init_metric(target_size);
  
  //
  // OPTIMIZE THE METRIC.
  //
  if (optimize)
    {
      delaunay_calculator->optimize();
    }
  
  //
  // OUTPUT
  //
  delaunay_calculator->output(ofilename);

  delete delaunay_calculator;
  return 0;
}

