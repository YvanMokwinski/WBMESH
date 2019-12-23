# include <ctime>
# include <iomanip>
# include <iostream>

using namespace std;

# include "hilbert_curve.hpp"

//****************************************************************************80

void d2xy ( int m, int d, int &x, int &y )

//****************************************************************************80
//
//  Purpose:
//
//    D2XY converts a 1D Hilbert coordinate to a 2D Cartesian coordinate.
//
//  Modified:
//
//    24 December 2015
//
//  Parameters:
//
//    Input, int M, the index of the Hilbert curve.
//    The number of cells is N=2^M.
//    0 < M.
//
//    Input, int D, the Hilbert coordinate of the cell.
//    0 <= D < N * N.
//
//    Output, int &X, &Y, the Cartesian coordinates of the cell.
//    0 <= X, Y < N.
//
{
  int n;
  int rx;
  int ry;
  int s;
  int t = d;

  n = i4_power ( 2, m );

  x = 0;
  y = 0;
  for ( s = 1; s < n; s = s * 2 )
  {
    rx = ( 1 & ( t / 2 ) );
    ry = ( 1 & ( t ^ rx ) );
    rot ( s, x, y, rx, ry );
    x = x + s * rx;
    y = y + s * ry;
    t = t / 4;
  }
  return;
}
//****************************************************************************80

int i4_power ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_POWER returns the value of I^J.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, J, the base and the power.  J should be nonnegative.
//
//    Output, int I4_POWER, the value of I^J.
//
{
  int k;
  int value;

  if ( j < 0 )
  {
    if ( i == 1 )
    {
      value = 1;
    }
    else if ( i == 0 )
    {
      cerr << "\n";
      cerr << "I4_POWER - Fatal error!\n";
      cerr << "  I^J requested, with I = 0 and J negative.\n";
      exit ( 1 );
    }
    else
    {
      value = 0;
    }
  }
  else if ( j == 0 )
  {
    if ( i == 0 )
    {
      cerr << "\n";
      cerr << "I4_POWER - Fatal error!\n";
      cerr << "  I^J requested, with I = 0 and J = 0.\n";
      exit ( 1 );
    }
    else
    {
      value = 1;
    }
  }
  else if ( j == 1 )
  {
    value = i;
  }
  else
  {
    value = 1;
    for ( k = 1; k <= j; k++ )
    {
      value = value * i;
    }
  }
  return value;
}
//****************************************************************************80


uint64_t i4_power ( uint64_t i, uint64_t j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_POWER returns the value of I^J.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, J, the base and the power.  J should be nonnegative.
//
//    Output, int I4_POWER, the value of I^J.
//
{
  uint64_t  k;
  uint64_t value;

//  if ( j < 0 )
//  {
//    if ( i == 1 )
//    {
//      value = 1;
//    }
//    else if ( i == 0 )
//    {
//      cerr << "\n";
//      cerr << "I4_POWER - Fatal error!\n";
//      cerr << "  I^J requested, with I = 0 and J negative.\n";
//      exit ( 1 );
//    }
//    else
//    {
//      value = 0;
//    }
//  }
//  else
  if ( j == 0 )
  {
    if ( i == 0 )
    {
      cerr << "\n";
      cerr << "I4_POWER - Fatal error!\n";
      cerr << "  I^J requested, with I = 0 and J = 0.\n";
      exit ( 1 );
    }
    else
    {
      value = 1;
    }
  }
  else if ( j == 1 )
  {
    value = i;
  }
  else
  {
    value = 1;
    for ( k = 1; k <= j; k++ )
    {
      value = value * i;
      //      std::cout << "val " << value  << std::endl;
    }
  }

  return value;
}
//****************************************************************************80

void rot ( int n, int &x, int &y, int rx, int ry ) 

//****************************************************************************80
//
//  Purpose:
//
//    ROT rotates and flips a quadrant appropriately
//
//  Modified:
//
//    24 December 2015
//
//  Parameters:
//
//    Input, int N, the length of a side of the square.  N must be a power of 2.
//
//    Input/output, int &X, &Y, the input and output coordinates of a point.
//
//    Input, int RX, RY, ???
//
{
  int t;

  if ( ry == 0 )
  {
//
//  Reflect.
//
    if ( rx == 1 )
    {
      x = n - 1 - x;
      y = n - 1 - y;
    }
//
//  Flip.
//
    t = x;
    x = y;
    y = t;
  }
  return;
}
//****************************************************************************80


void rot ( uint64_t n, uint64_t &x, uint64_t &y, uint64_t rx, uint64_t ry ) 

//****************************************************************************80
//
//  Purpose:
//
//    ROT rotates and flips a quadrant appropriately
//
//  Modified:
//
//    24 December 2015
//
//  Parameters:
//
//    Input, int N, the length of a side of the square.  N must be a power of 2.
//
//    Input/output, int &X, &Y, the input and output coordinates of a point.
//
//    Input, int RX, RY, ???
//
{
  uint64_t t;

  if ( ry == 0 )
  {
//
//  Reflect.
//
    if ( rx == 1 )
    {
      x = n - 1 - x;
      y = n - 1 - y;
    }
//
//  Flip.
//
    t = x;
    x = y;
    y = t;
  }
  return;
}


void timestamp ( )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    31 May 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct std::tm *tm_ptr;
  size_t len;
  std::time_t now;

  now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );

  len = std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

  std::cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
//****************************************************************************80

int xy2d ( int m, int x, int y )

//****************************************************************************80
//
//  Purpose:
//
//    XY2D converts a 2D Cartesian coordinate to a 1D Hilbert coordinate.
//
//  Discussion:
//
//    It is assumed that a square has been divided into an NxN array of cells,
//    where N is a power of 2.
//
//    Cell (0,0) is in the lower left corner, and (N-1,N-1) in the upper 
//    right corner.
//
//  Modified:
//
//    24 December 2015
//
//  Parameters:
//
//    Input, int M, the index of the Hilbert curve.
//    The number of cells is N=2^M.
//    0 < M.
//
//    Input, int X, Y, the Cartesian coordinates of a cell.
//    0 <= X, Y < N.
//
//    Output, int XY2D, the Hilbert coordinate of the cell.
//    0 <= D < N * N.
//
{
  int d = 0;
  int n;
  int rx;
  int ry;
  int s;

  n = i4_power ( 2, m );
  //  std::cout << "hello " << m << " " << n << " " << n / 2 << std::endl;
  for ( s = n / 2; s > 0; s = s / 2 )
  {
    rx = ( x & s ) > 0;
    ry = ( y & s ) > 0;
    d = d + s * s * ( ( 3 * rx ) ^ ry );
    //    std::cout << "hello " << d << std::endl;
    rot ( s, x, y, rx, ry );
  }
  return d;
}


uint64_t convert ( uint64_t m, uint64_t x, uint64_t y )
{
  uint64_t d = 0;
  uint64_t rx;
  uint64_t ry;
  uint64_t s;
  
  uint64_t n = i4_power ( ((uint64_t)2), m );
  //  std::cout << "hello " << m << " " << n << " " << n / 2 << std::endl;
  for ( s = n / 2; s > 0; s = s / 2 )
  {
    rx = ( x & s ) > 0;
    ry = ( y & s ) > 0;
    d = d + s * s * ( ( 3 * rx ) ^ ry );
    //    std::cout << "hello " << d << std::endl;
    rot ( s, x, y, rx, ry );
  }
  return d;
}



double convert ( uint64_t m, double x, double y )
{
  uint64_t * px = (uint64_t*)&x;
  uint64_t * py = (uint64_t*)&y;
  *px -= ((uint64_t)9223372036854775808);
  *py -= ((uint64_t)9223372036854775808);
  uint64_t n = i4_power ( ((uint64_t)2), m );
  std::cout << *px << " " << *py << " " << n << std::endl;
  
  uint64_t h =  convert(m,*px,*py);
  //  h += ((uint64_t)9223372036854775808);
  std::cout << "h " << h << std::endl;
  double * hh = (double*)&h;
  std::cout << "h " << *hh <<  std::endl;
  return *hh;
#if 0
  uint64_t d = 0;
  uint64_t rx;
  uint64_t ry;
  uint64_t s;
  
  uint64_t n = i4_power ( ((uint64_t)2), m );
  //  std::cout << "hello " << m << " " << n << " " << n / 2 << std::endl;
  for ( s = n / 2; s > 0; s = s / 2 )
  {
    rx = ( x & s ) > 0;
    ry = ( y & s ) > 0;
    d = d + s * s * ( ( 3 * rx ) ^ ry );
    //    std::cout << "hello " << d << std::endl;
    rot ( s, x, y, rx, ry );
  }
  return d;
#endif
}

//  uint64_t, and subtract = 9223372036854775808 and use the result as an unsigned 64-bit inte




static void coordinates2d(long long int				N_, 
			  double				box_[],
			  double 				vertex_[],
			  long long int                         vertexoff_,
			  unsigned long long int *  	cod2_,
			  long long int				codoff_)
{
  static const unsigned long long m = 1LL<<62;
  static const int BitTab[2] 	= {1,2};
  static const int GeoCod[4]	= {1,2,0,3};
  static const int HilCod[4][4] = {{0,3,2,1}, {0,1,2,3}, {0,1,2,3}, {2,1,0,3}};
  static const double len 	= 4.611686018427387904e+18;
  double box[4],dbl;
  unsigned long long int IntCrd[2],cod;
  int rot[4];
  box[0] 	= box_[0];
  box[1] 	= box_[1];
  box[2] 	= len / (box_[2] - box_[0]);
  box[3] 	= len / (box_[3] - box_[1]);
  double loc[2];
  { int i;
    for(i=0; i<N_; i++)
      {
	loc[0] = vertex_[vertexoff_*i+0];
	loc[1] = vertex_[vertexoff_*i+1];
	/* Convert double precision coordinates to integers */
	dbl 	= (loc[0] - box[0]) * box[0+2];
	IntCrd[0] = (unsigned long long int)dbl;
	dbl 	= (loc[1] - box[1]) * box[1+2];
	IntCrd[1] = (unsigned long long int)dbl;
	/* Binary hilbert renumbering loop */
	cod = 0;
	rot[0] = GeoCod[0];
	rot[1] = GeoCod[1];
	rot[2] = GeoCod[2];
	rot[3] = GeoCod[3];
	{ int b;
	  for(b=0;b<31;b++)
	    {
	      int GeoWrd = 0;

	      if(IntCrd[0] & m)
		GeoWrd |= BitTab[0];
	      IntCrd[0] = IntCrd[0]<<1;

	      if(IntCrd[1] & m)
		GeoWrd |= BitTab[1];
	      IntCrd[1] = IntCrd[1]<<1;

	      const int NewWrd = rot[ GeoWrd ];
	      cod = cod<<2 | NewWrd;
	      rot[0] = HilCod[ NewWrd ][ rot[0] ];
	      rot[1] = HilCod[ NewWrd ][ rot[1] ];
	      rot[2] = HilCod[ NewWrd ][ rot[2] ];
	      rot[3] = HilCod[ NewWrd ][ rot[3] ];
	    } }
	cod2_[codoff_*i+0] 	= cod;
	cod2_[codoff_*i+1] 	= i;
      } }

}

