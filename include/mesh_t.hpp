#pragma once

#include <math.h>
#include <float.h>
#include <string.h>
#include <stack>

#include <sstream>
#include <limits>
#include <map>
#include <set>
#include <algorithm>

#include "xy_t.hpp"  
#include "edge_t.hpp"
#include "triangle_t.hpp"

#include "metric_t.hpp"
#include "delaunay_cad_t.hpp"
#include "InputFile/Medit.h"
#include "classifier.hpp"
#include "HashPair.hpp"
#include "Hasher.hpp"
#include "Orientation.hpp"
#include "HashCavityTable.hpp"
#include "CavityOperator.hpp"


#if 0

//! @brief Set of mathematical static functions. 
//! @tparam _float_type The real type to use.
template <typename _float_type> struct Math
{
  //! @brief Compute the squared root.
  //! @param x_ The value from which we compute the squared root.
  //! @return The squared root.
  static constexpr _float_type Sqrt(const _float_type&x_) noexcept;
  //! @brief Compute the absolute value.
  //! @param x_ The value from which we compute the absolue value.
  //! @return The absolute value.
  static constexpr _float_type Abs(const _float_type&x_) noexcept;
  //! @brief Compute the power x^y.
  //! @param x_ The base of the power.
  //! @param y_ The exponent of the power.
  //! @return The power x^y.
  static constexpr _float_type Pow(const _float_type&x_,const _float_type&y_) noexcept;
  //! @brief Compute the maximum between two real.
  //! @param a_ The first operand.
  //! @param b_ The second operand.
  //! @return MAX(a,b)
  static constexpr _float_type Max(const _float_type&a_,const _float_type&b_) noexcept;
  //! @brief Compute the maximum between two real.
  //! @param a_ The first operand.
  //! @param b_ The second operand.
  //! @return MAX(a,b)
  static constexpr _float_type Min(const _float_type&a_,const _float_type&b_) noexcept;
#if 0
  //! @brief Get the maximum value of the _float_type.
  //! @return The maximum value of the _float_type.
  static constexpr _float_type MaxValue() noexcept;
  //! @brief Get the minimum value of the _float_type.
  //! @return The minimum value of the _float_type.
  static constexpr _float_type MinValue() noexcept;
#endif
  //! @brief Set the value to zero if the absolute value is less than a tolerance.
  //! @param The value.
  //! @param The tolerance.
  //! @return The filtered value.
  static constexpr _float_type Filter(const _float_type&value_,const _float_type&tolerance_) noexcept;
};

//! @copydoc Math
template <> struct Math<double>
{
  static constexpr unsigned int NumDigits = 15;
  static constexpr double Zero = 0.0;
  static constexpr double MachineEpsilon = 2.22044604925031308e-16;
  static constexpr double machineEpsilon = 2.22044604925031308e-16;
  //! @copydoc Math::Pow(const _float_type&x_,const _float_type&y_)
  static constexpr double Pow(const double&x_,const double&y_) noexcept
  {
    return pow(x_,y_);
  };
  //! @copydoc Math::Sqrt(const _float_type&x_)
  static constexpr double Sqrt(const double&x_) noexcept
  {
    return  sqrt(x_);
  };
  //! @copydoc Math::Abs(const _float_type&x_)
  static constexpr double Abs(const double&x_) noexcept
  {
    return  (x_<0.0) ? -x_ : x_;
  };
#if 1
  //! @copydoc Math::MaxValue()
  static constexpr double MaxValue() noexcept { return DBL_MAX; };
  //! @copydoc Math::MinValue()
  static constexpr double MinValue() noexcept { return DBL_MIN; };
#endif
  //! @copydoc Math::Max(const _float_type&,const _float_type&)
  static constexpr double Max(const double&a_,const double&b_) noexcept
  {
    return (a_>=b_) ? a_ : b_;
  };
  //! @copydoc Math::Min(const _float_type&,const _float_type&)
  static constexpr double Min(const double&a_,const double&b_) noexcept
  {
    return (a_<=b_) ? a_ : b_;
  };

  //! @copydoc Math::Filter(const _float_type&,const _float_type&)
  static constexpr double Filter(const double&value_,const double&tolerance_) noexcept
  {
    return (value_ < -tolerance_) ? value_ : ( (value_ > tolerance_) ? value_ : 0.0  );
  };

};
//constexpr double Math<double>::MachineEpsilon;
//constexpr double Math<double>::machineEpsilon;


//! @copydoc Math
template <> struct Math<float>
{
  static constexpr unsigned int NumDigits = 7;
  static constexpr float Zero = 0.0f;
  static constexpr float MachineEpsilon = 1.19209290e-07f;
  //! @copydoc Math::Pow(const _float_type&x_,const _float_type&y_)
  static constexpr float Pow(const float&x_,const float&y_) noexcept
  {
    return powf(x_,y_);
  };
  //! @copydoc Math::Sqrt(const _float_type&x_)
  static constexpr float Sqrt(const float&x_) noexcept
  {
    return sqrtf(x_);
  };
  //! @copydoc Math::Abs(const _float_type&x_)
  static constexpr float Abs(const float&x_) noexcept
  {
    return  (x_<0.0) ? -x_ : x_;
  };
  //! @copydoc Math::MaxValue()
#if 0
  static constexpr float MaxValue() noexcept { return FLT_MAX; };
  //! @copydoc Math::MinValue()
  static constexpr float MinValue() noexcept { return FLT_MIN; };
#endif
  //! @copydoc Math::Max(const _float_type&,const _float_type&)
  static constexpr float Max(const float&a_,const float&b_) noexcept
  {
    return (a_>=b_) ? a_ : b_;
  };
  //! @copydoc Math::Min(const _float_type&,const _float_type&)
  static constexpr float Min(const float&a_,const float&b_) noexcept
  {
    return (a_<=b_) ? a_ : b_;
  };

  //! @copydoc Math::Filter(const _float_type&,const _float_type&)
  static constexpr float Filter(const float &value_,const float&tolerance_) noexcept
  {
    return (value_ < -tolerance_) ? value_ : ( (value_ > tolerance_) ? value_ : 0.0f  );
  };

};
//constexpr float Math<float>::MachineEpsilon;

//! @copydoc Math
template <> struct Math<long double>
{
  static constexpr unsigned int NumDigits = 17;
  static constexpr long double Zero = 0.0L;
  static constexpr long double MachineEpsilon = 1.08420217248550443401e-19L;
  //! @copydoc Math::Pow(const _float_type&x_,const _float_type&y_)
  static constexpr long double Pow(const long double&x_,const long double&y_) noexcept
  {
    return powl(x_,y_);
  };
  //! @copydoc Math::Sqrt(const _float_type&x_)
  static constexpr long double Sqrt(const long double&x_) noexcept
  {
    return sqrtl(x_);
  };
  //! @copydoc Math::Abs(const _float_type&x_)
  static constexpr long double Abs(const long double&x_) noexcept
  {
    return  (x_<0.0) ? -x_ : x_;
  };
#if 0
  //! @copydoc Math::MaxValue()
  static constexpr long double MaxValue() noexcept { return LDBL_MAX; };
  //! @copydoc Math::MinValue()
  static constexpr long double MinValue() noexcept { return LDBL_MIN; };
#endif
  //! @copydoc Math::Max(const _float_type&,const _float_type&)
  static constexpr long double Max(const long double&a_,const long double&b_) noexcept
  {
    return (a_>=b_) ? a_ : b_;
  };
  //! @copydoc Math::Min(const _float_type&,const _float_type&)
  static constexpr long double Min(const long double&a_,const long double&b_) noexcept
  {
    return (a_<=b_) ? a_ : b_;
  };

  //! @copydoc Math::Filter(const _float_type&,const _float_type&)
  static constexpr long double Filter(const long double &value_,const long double&tolerance_) noexcept
  {
    return (value_ < -tolerance_) ? value_ : ( (value_ > tolerance_) ? value_ : 0.0L  );
  };

};
//constexpr long double Math<long double>::MachineEpsilon;
#endif


struct mesh_t
{

private: char m_debug_filebasename[256];
public: inline const char* debug_filebasename() const noexcept
  {
    return m_debug_filebasename;
  };
  
public: inline void debug_filebasename(const char * n)  noexcept
  {
    strcpy(m_debug_filebasename,n);
  };
  
private: bool m_debug_mode {};  
public: inline void debug_mode(bool debug_mode) noexcept
  {
    this->m_debug_mode = debug_mode;
  };
public: inline bool debug_mode() const noexcept
  {
    return this->m_debug_mode;
  };

public:
  inline void optimize()
  {
    //
    // Analyze the boundary.
    //
    while_any< edge_split_boundary<edge_select_boundary_saturation> >();

#if 0
    int_t triangle = 1;
    double dx = 1.0/sqrt(xy[1].m_metric[0]);
    double dy = 1.0/sqrt(xy[1].m_metric[2]);
    int nx = (xy[2].x-xy[1].x) / dx;
    int ny = (xy[3].y-xy[1].y) / dy; 
    std::cout << nx << " " << ny << std::endl;
    for (int i=0;i<nx-1;++i)
      {
	for (int j=0;j<ny-1;++j)
	  {
	    int_t node = create_node(xy[1].x + (i+1)*dx, xy[1].y + (j+1)*dy);
	    xy[node].m_metric[0] = xy[1].m_metric[0];
	    xy[node].m_metric[1] = xy[1].m_metric[1];
	    xy[node].m_metric[2] = xy[1].m_metric[2];
	    triangle = find_triangle(node,triangle);
	    if (triangle)
	      {		
		this->InsertPoint(triangle,node,false);
	      }
	    else
	      {
    double coo[6];
    double x = xy[node].x;
    double y = xy[node].y;

    for (int_t i=1;i<=m_ntriangles;++i)
      {
	int_t seed = i;
    	triangle_coordinates(seed,coo);
	
	double x0 = coo[0];
	double y0 = coo[1];
	double x1 = coo[2];
	double y1 = coo[3];
	double x2 = coo[4];
	double y2 = coo[5];
	//	std::cout << "seed " << seed << std::endl;
	//	next = false;
	Orientation::Kind oo[3];
	oo[0] = ComputeOrientation(x0,y0,x1,y1,x,y);
	oo[1] = ComputeOrientation(x1,y1,x2,y2,x,y);
	oo[2] = ComputeOrientation(x2,y2,x0,y0,x,y);
	//	std::cout << "seed " << oo[0] <<  " " << oo[1] << " " << oo[2] << std::endl;

	if (Orientation::BACKWARD!=oo[0] && Orientation::BACKWARD!=oo[1] && Orientation::BACKWARD!=oo[2])
	  {
	    //	    std::cout << "yo " << std::endl;
	    triangle = seed;
	    break;
	  }
      }

	    if (triangle)
	      {		
		this->InsertPoint(triangle,node,false);
	      }
	    else
	      {
		std::cout << "triangle not found" << std::endl;
		triangle = 1;
	      }
	      }
	  }	
      }
#endif
    
#if 0
    //
    //
    //
    for (int_t i=1;i<=this->m_nedges;++i)
      {
	if (edge_split_boundary_to_interior<edge_select_boundary>::criterion(this, i))
	  {	
	    edge_split_boundary_to_interior<edge_select_boundary>::apply(this, i);
	  }
	
      }
#endif
    //    while_any<  >();
    
    //
    // Analyze interior points.
    //
    while_any< edge_split<edge_select_interior> >();

    //
    // Optimize the node positions
    //
    //    optimize_node_positions();    
  };
  
  inline void optimize_interior()
  {
    while_any< edge_split<edge_select_interior> >();
  };
  
  inline void optimize_boundary()
  {
    while_any< edge_split_boundary<edge_select_boundary_saturation> >();
  };

  

  //  using FaceEdge = typename std::pair<int_t,int>;
  

  inline int_t markboundaryedge(int_t edge_,
				int_t a,
				int_t b)
  {
    auto edge = get_edge_ptr(edge_);
    if ((edge->start() == a) && (edge->stop() == b))
      {
	edge->is_boundary(true);
	edge->is_nonflipable(true);
	return edge_;
      }
    else if ((edge->start() == b) && (edge->stop() == a))
      {
	edge->is_boundary(true);
	edge->is_nonflipable(true);
	return edge_;
      }
    return 0;
  };

  
  inline void RebuildWeakBoundary(int_t ep0, int_t ep1, int_t point_,CavityOperator& cavityOperator_,int_t* newedge_)
  {
#if 1
    cavityOperator_.reset_faces();
    cavityOperator_.clear_interior_edges();// m_listCavityInteriorEdges.clear();
    cavityOperator_.clear_hash();
    int ith=0;
    
    //    typename std::list<FaceEdge>::iterator it = cavityOperator_.m_listCavityBoundaryEdges.begin();
    FaceEdge 	triangleEdge 	= cavityOperator_.boundary_edge(0);
    int_t edge		= triangleEdge.first;
    int_t is_boundary;
    int 		way 		= triangleEdge.second;
    auto edge_ptr = get_edge_ptr(edge);
    if (way)
      {
	int_t startingEdge 	= edge_ptr->start();
	int_t endingEdge 	= edge_ptr->stop();
	//	tEdge* edgep 		= (way>0) ? new tEdge(point_,startingEdge) 	: new tEdge(point_,endingEdge);
	//	tEdge* edgen 		= (way>0) ? new tEdge(endingEdge,point_) 	: new tEdge(startingEdge,point_);
	int_t edgep	= (way>0) ? create_edge(point_,startingEdge) 	: create_edge(point_,endingEdge);
	int_t edgen	= (way>0) ? create_edge(endingEdge,point_) 	: create_edge(startingEdge,point_);
#if 0
	std::cout << " edge interior " << edgep << " " << ", m_nedges = " << m_nedges << std::endl;
	std::cout << get_edge_ptr(edgep)->start() << " " << get_edge_ptr(edgep)->stop() << " " << std::endl;
	std::cout << " edge interior " << edgen << " " << ", m_nedges = " << m_nedges << std::endl;
	std::cout << get_edge_ptr(edgen)->start() << " " << get_edge_ptr(edgen)->stop() << " " << std::endl;
#endif	    
	//	    list.push_front(edgep);
	//	    list.push_front(edgen);	    
	//	      std::cout << "##################11111111111" << std::endl;
	is_boundary = markboundaryedge(edgep,point_,ep0);
	if (is_boundary)
	  {
	    newedge_[0]  = edgep;
	  }
	is_boundary = markboundaryedge(edgep,point_,ep1);
	if (is_boundary)
	  {
	    newedge_[1]  = edgep;
	  }
	    
	//	      std::cout << "##################22222222222" << std::endl;
	    
	is_boundary = markboundaryedge(edgen,point_,ep0);
	if (is_boundary)
	  {
	    newedge_[0]  = edgen;
	  }
	is_boundary = markboundaryedge(edgen,point_,ep1);
	if (is_boundary)
	  {
	    newedge_[1]  = edgen;
	  }
	//	      std::cout << "##################33333" << std::endl;
	    
	cavityOperator_.add_edge(m_edges[edgep].start(),
				 m_edges[edgep].stop(),
				 edgep);
	cavityOperator_.add_edge(m_edges[edgen].start(),
					 m_edges[edgen].stop(),
					 edgen);

	cavityOperator_.add_interior_edge(edgep);
	cavityOperator_.add_interior_edge(edgen);
	//	cavityOperator_.add_interior_edge(edgep);
	//	cavityOperator_.add_interior_edge(edgen);

	int_t face = create_triangle(edgep,1,edge,(way>0) ? 2 : -2,edgen,3,__LINE__);
	cavityOperator_.AddFace(face);
	//	cavityOperator_.AddFace(face);	
	//	cavityOperator_.AddFace(new tFace(edgep,1,edge,(way>0) ? 2 : -2,edgen,3));
      } 

    for (int_t i=1;i<cavityOperator_.nboundary_edges();++i)
      //     for (it++; it != cavityOperator_.m_listCavityBoundaryEdges.end(); it++)
      {
	//	FaceEdge 	triangleEdge 	= *it;
	FaceEdge 	triangleEdge 	= cavityOperator_.boundary_edge(i);
	int_t 		edge		= triangleEdge.first;
	//    std::cout << " edge boundary " << edge << " " << ", m_nedges = " << m_nedges << std::endl;
	//    std::cout << get_edge_ptr(edge)->start() << " " << get_edge_ptr(edge)->stop() << " " << std::endl;
	edge_ptr = get_edge_ptr(edge);
	int 		way 		= (triangleEdge.second>0) ? 1 : -1;
	int_t		p0 		= edge_ptr->start(triangleEdge.second);
	int_t		p1 		= edge_ptr->stop(triangleEdge.second);

	int_t 		edgep		= cavityOperator_.try_get_edge(point_,p0);
	int_t 		edgen		= cavityOperator_.try_get_edge(point_,p1);	    
	const int 	wayp 		= (edgep) ? -3 : 3;
	const int 	wayn 		= (edgen) ? -2 : 2;
	if (!edgep)
	  {
	    edgep = create_edge(point_,p0);
	    cavityOperator_.add_edge(point_,
				     p0,
				     edgep);
	    
#if 0
	    std::cout << " edge interior " << edgep << " " << ", m_nedges = " << m_nedges << std::endl;
	    std::cout << get_edge_ptr(edgep)->start() << " " << get_edge_ptr(edgep)->stop() << " " << std::endl;
#endif
	
	    is_boundary = markboundaryedge(edgep,point_,ep0);
	    if (is_boundary)
	      {
		newedge_[0]  = edgep;
	      }

	    is_boundary = markboundaryedge(edgep,point_,ep1);
	    if (is_boundary)
	      {
		newedge_[1]  = edgep;
	      }

	    
	    cavityOperator_.add_interior_edge(edgep);
	    
	    //	    cavityOperator_.add_interior_edge(edgep);
	    //		  list.push_front(edgep);

	  }	    
	if (!edgen)
	  {
	    edgen = create_edge(p1,point_);
	    cavityOperator_.add_edge(p1,
				     point_,
				     edgen);
		  
	    is_boundary = markboundaryedge(edgen,point_,ep0);
	    if (is_boundary)
	      {
		newedge_[0]  = edgen;
	      }
	    is_boundary = markboundaryedge(edgen,point_,ep1);
	    if (is_boundary)
	      {
		newedge_[1]  = edgen;
	      }
	    
#if 0
	    std::cout << " edge interior " << edgen << " " << ", m_nedges = " << m_nedges << std::endl;
	    std::cout << get_edge_ptr(edgen)->start() << " " << get_edge_ptr(edgen)->stop() << " " << std::endl;
#endif
	    
	    cavityOperator_.add_interior_edge(edgen);

	    //	    cavityOperator_.add_interior_edge(edgen);
	    //		  list.push_front(edgen);
		  
	  }
	int_t face = create_triangle(edge,way,
				     edgen,wayn,
				     edgep,wayp,__LINE__);
	cavityOperator_.AddFace(face);
	//	cavityOperator_.AddFace(face);	
      }
#endif
  }


	
  inline void RebuildBoundary(int_t edge_, int_t ep0, int_t ep1, int_t point_,CavityOperator& cavityOperator_,int_t* newedge_)
  {
	
    cavityOperator_.reset_faces();
    cavityOperator_.clear_interior_edges();
    cavityOperator_.clear_hash();

    weak_remove_edge(edge_);
    
    int_t first_edge  = create_edge(ep1,point_);
    int_t second_edge = create_edge(point_,ep0);
    cavityOperator_.add_edge(ep1, point_,first_edge);
    cavityOperator_.add_edge(point_,ep0, second_edge);
    
    get_edge_ptr(first_edge)->is_boundary(true);
    get_edge_ptr(first_edge)->is_nonflipable(true);
    get_edge_ptr(second_edge)->is_boundary(true);
    get_edge_ptr(second_edge)->is_nonflipable(true);
    newedge_[0] = first_edge;
    newedge_[1] = second_edge;
    
    for (int_t i=0;i<cavityOperator_.nboundary_edges();++i)
      // for (auto fe :  cavityOperator_.m_listCavityBoundaryEdges)
      {
	auto fe  = cavityOperator_.boundary_edge(i);
	int_t edge		= fe.first;
	int_t is_boundary;
	int   way 		= fe.second;
	if (edge == edge_)
	  {
	    
	  }
	else
	  {
	    //	    
	    auto edge_ptr = get_edge_ptr(edge);
	    int_t	p0 		= edge_ptr->start(way);
	    int_t	p1 		= edge_ptr->stop(way);
	    //	    std::cout << "ooooooooooooooooo "  << p0 << " " << p1 << std::endl;
	    
	    int_t 	edgep		= cavityOperator_.try_get_edge(point_,p0);
	    int_t 	edgen		= cavityOperator_.try_get_edge(point_,p1);

	    //	    std::cout << edgep << " eeeeeeeeeee ? " << edgen << std::endl;
	    int 	wayp 		= (edgep) ? -3 : 3;
	    int 	wayn 		= (edgen) ? -2 : 2;

	    bool ff = false;
	    bool ff0 = false;
	    if (!edgep)
	      {
		ff = true;
		edgep	= create_edge(point_,p0);
#if 0
		if (way>0)
		  {
		    std::cout << "create edgep " << point_ << " " << p0 << std::endl;
		  }
		else
		  {
		    std::cout << "create edgep " << point_ << " " << p1 << std::endl;
		  }
#endif
	      }
	    
	    if (!edgen)
	      {
		ff0 = true;
		edgen	= create_edge(p1,point_);
#if 0
		if (way>0)
		  {
		    std::cout << "create edgen " << p1 << " " << point_ << std::endl;
		  }
		else
		  {
		    std::cout << "create edgen " << p0 << " " << point_ << std::endl;

		  }
#endif
	      }
	    
	    if (ff)
	      {
		cavityOperator_.add_edge(m_edges[edgep].start(wayp),
					 m_edges[edgep].stop(wayp),
					 edgep);
	      }
	    if (ff0)
	      {
		cavityOperator_.add_edge(m_edges[edgen].start(wayn),
					 m_edges[edgen].stop(wayn),	 
					 edgen);
	      }

#if 0
	    std::cout << "oooooooooooooooooDONE "  << ep0 << " " << ep1 << std::endl;
	    std::cout << "oooooooooooooooooDONE "  << p0 << " " << p1 << std::endl;
	    
	    std::cout << "oooooooooooooooooDONE "  << get_edge_ptr(edge)->start(way) << " " << get_edge_ptr(edge)->stop(way) << std::endl;
	    std::cout << "oooooooooooooooooDONE "  << get_edge_ptr(edgen)->start(wayn) << " " << get_edge_ptr(edgen)->stop(wayn) << std::endl;
	    std::cout << "oooooooooooooooooDONE "  << get_edge_ptr(edgep)->start(wayp) << " " << get_edge_ptr(edgep)->stop(wayp) << std::endl;
#endif	    
	    int_t face = create_triangle(edge,way,
					 edgen,wayn,
					 edgep,wayp,__LINE__);

#if 0	    
	    if (get_edge_ptr(edge)->stop(way) != get_edge_ptr(edgen)->start(wayn))
	      {

		std::cout << "hello0" << std::endl;
		exit(1);
	      }
	    if (get_edge_ptr(edgen)->stop(wayn) != get_edge_ptr(edgep)->start(wayp))
	      {

		std::cout << "hello1" << std::endl;
		exit(1);
	      }
	    if (get_edge_ptr(edgep)->stop(wayp) != get_edge_ptr(edge)->start(way))
	      {

		std::cout << "hello2" << std::endl;
		exit(1);
	      }
#endif
	    //	    int_t nodes[3];
	    //	    triangle_to_nodes(triangles[face],nodes);
	    //	    std::cout << "connec " << nodes[0] << " " << nodes[1] << " " << nodes[2] << std::endl;
#if 0
	    if (nodes[0]==nodes[1] || nodes[1] == nodes[2])
	      {
		std::cout << "hello" << std::endl;
		exit(1);
	      }
	    else
	      {

	      }
#endif
	    
	    cavityOperator_.AddFace(face);  
	  }
      }
    
#if 1
    if (m_edges[newedge_[0]].start() == ep1 || m_edges[newedge_[0]].stop() == ep1)
      {
	int_t tmp = newedge_[0];
	newedge_[0] = newedge_[1];
	newedge_[1] = tmp;
      }
#endif
  }




  double edge_squared_length(double hx,
			     double hy,
			     double a,
			     double b,
			     double c)
  {
    return hx*hx*a + 2.0 * b*hx*hy + c*hy*hy;
  }

  double edge_squared_length(double hx,
			     double hy)
  {
    return hx*hx+hy*hy;
  }

  double edge_squared_length(double px,
			     double py,
			     double qx,
			     double qy)
  {
    return edge_squared_length(qx-px,qy-py);
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
    ccy[0] = y + hy;
    r[0] = edge_squared_length(x0-ccx[0],y0-ccx[1],a,b,c);
    r[1] = edge_squared_length(hx,hy,a,b,c);
  
    std::cout << "ratio " << r[1] / r[0] << std::endl;
    return r[1] / r[0];
  }



  inline bool InsertPointWeakBoundary(int_t point_,
				      int_t& edge_,
				      bool iii,
				      int_t*newedge_) 
  {
    auto edge_ptr = get_edge_ptr(edge_);
    auto t0       = edge_ptr->triangle(0);
    auto t1       = edge_ptr->triangle(1);
    
    //    std::cerr << "INSERT POINT WEAK BOUNDARY triangles " << t0 << " " << t1 << "edges " << edge_ << std::endl;
    if (!t0 || !t1)
      {
	std::cerr << "ffffffffffffexit " << t0 << " " << t1 << std::endl;
	exit(1);
	//	return InsertPointBoundary( (t0) ? t0 : t1, point_, edge_);
      }
    
    //    int jj = get_triangle_ptr(t0)->get_local_edge_index(edge_);
    //    auto p0 = edge_ptr->start(get_triangle_ptr(t0)->way(jj));
    //    auto p1 = edge_ptr->stop(get_triangle_ptr(t0)->way(jj));
    auto p0 = edge_ptr->start();
    auto p1 = edge_ptr->stop();
    
    //    std::cout << "edge!!!!!!!!!!!!!!!!! nodes " << p0 << ", " << p1 <<  ", triangle " << edge_ptr->triangle(0) << std::endl;
    bool success = m_cavityOperator.Compute(edge_ptr->triangle(0),
					    point_,
					    [this] (int_t face,
						    int_t point)
					    {
					      int_t nodes[3];
					      triangle_to_nodes(face,nodes);
#if 1
					      //						std::cout << "testing triangle " << face << std::endl;
					      return (Orientation::BACKWARD != 
						      ComputeCircleOrientation(xy[nodes[0]].x,
									       xy[nodes[0]].y,
									       xy[nodes[1]].x,
									       xy[nodes[1]].y,
									       xy[nodes[2]].x,
									       xy[nodes[2]].y,
									       xy[point].x,
									       xy[point].y));
#endif
					    },
					    
					    [this,edge_](int_t edge)
					    {
					      // std::cout << " cccc " << edge << std::endl;
					      //						return !this->get_edge_ptr(edge)->triangle(1) && edge != edge_;				
					      return m_edges[edge].is_boundary() && edge != edge_;
					      //  return edge->is_boundary();
					    },
					    					    
					    [this](int_t face)
					    {
					      return get_triangle_ptr(face);
					    },
					    
					    [this](int_t edge)
					    {
					      return get_edge_ptr(edge);
					    }
					    
					    );
    
    if (success)
      {
	//
	// Remove the cavity.
	//
	//	    std::cout << "BEFORE REBUILD " << ntriangles << std::endl;
	*this -= m_cavityOperator;

	RebuildWeakBoundary(p0, p1, point_, m_cavityOperator, newedge_);

	print_mesh(m_debug_filebasename);
      }
    else
      {
	std::cout << "unable to compute the cavity" << std::endl;
	return false;
      }
    
    return success;
  };

  inline bool InsertPointBoundary(int_t point_,
				  int_t& edge_,
				  bool iii,
				  int_t*newedge_) 
  {
    auto edge_ptr = get_edge_ptr(edge_);
    auto t0       = edge_ptr->triangle(0);
    auto t1       = edge_ptr->triangle(1);
    
    //    std::cerr << "INSERT POINT WEAK BOUNDARY triangles " << t0 << " " << t1 << "edges " << edge_ << std::endl;
    if (t1)
      {
	std::cerr << "ffffffffffffexit " << t0 << " " << t1 << std::endl;
	exit(1);
	//	return InsertPointBoundary( (t0) ? t0 : t1, point_, edge_);
      }
    
    //    int jj = get_triangle_ptr(t0)->get_local_edge_index(edge_);
    //    auto p0 = edge_ptr->start(get_triangle_ptr(t0)->way(jj));
    //    auto p1 = edge_ptr->stop(get_triangle_ptr(t0)->way(jj));


    int i = get_triangle_ptr(t0)->get_local_edge_index(edge_);    
    auto p0 = edge_ptr->start(get_triangle_ptr(t0)->way(i));
    auto p1 = edge_ptr->stop(get_triangle_ptr(t0)->way(i));
    
    //    std::cout << "edge!!!!!!!!!!!!!!!!! nodes " << p0 << ", " << p1 <<  ", triangle " << edge_ptr->triangle(0) << std::endl;
    bool success = m_cavityOperator.Compute(t0,
					    point_,
					    [this] (int_t face,
						    int_t point)
					    {
					      int_t nodes[3];
					      triangle_to_nodes(face,nodes);
#if 1
					      //						std::cout << "testing triangle " << face << std::endl;
					      return (Orientation::BACKWARD != 
						      ComputeCircleOrientation(xy[nodes[0]].x,
									       xy[nodes[0]].y,
									       xy[nodes[1]].x,
									       xy[nodes[1]].y,
									       xy[nodes[2]].x,
									       xy[nodes[2]].y,
									       xy[point].x,
									       xy[point].y));
#endif
					    },
					    
					    [this,edge_](int_t edge)
					    {
					      // std::cout << " cccc " << edge << std::endl;
					      //						return !this->get_edge_ptr(edge)->triangle(1) && edge != edge_;				
					      return m_edges[edge].is_boundary();// && edge != edge_;
					      //  return edge->is_boundary();
					    },
					    					    					    
					    [this](int_t face)
					    {
					      return get_triangle_ptr(face);
					    },
					    
					    [this](int_t edge)
					    {
					      return get_edge_ptr(edge);
					    });
    
    if (success)
      {
	//
	// Remove the cavity.
	//
	//	    std::cout << "BEFORE REBUILD " << ntriangles << std::endl;
	*this -= m_cavityOperator;

	RebuildBoundary(edge_,p0, p1, point_, m_cavityOperator, newedge_);
	this->check(__LINE__);

	print_mesh();
      }
    else
      {
	std::cout << "unable to compute the cavity" << std::endl;
	return false;
      }
    
    return success;
  };



  inline int_t edge_saturation_size(int_t edge)
  {
    const edge_t* edge_ptr = get_edge_ptr(edge);
    auto node0_ptr = get_node_ptr(edge_ptr->start());
    auto node1_ptr = get_node_ptr(edge_ptr->stop());
    auto metric0   = &node0_ptr->m_metric[0];
    auto metric1   = &node1_ptr->m_metric[0];
    auto coo0      = &node0_ptr->x;
    auto coo1      = &node1_ptr->x;
    double 
      x0 = coo0[0],
      y0 = coo0[1],
      
      x1 = coo1[0],
      y1 = coo1[1];

    double
      ex = x1 - x0,
      ey = y1 - y0;
    
    //
    //           h0^-2 0       ex
    //  ex ey    0     h1^-2   ey
    //
    double s0  = 1.0/sqrt(metric0[0]);
    double s1  = 1.0/sqrt(metric1[0]);

    double l   = sqrt( (x0-x1)*(x0-x1)+(y0-y1)*(y0-y1));
    
    double neg = (s0 > s1) ? -1 : 1;
    
    double tmp = s0;
    s0 = std::min(tmp,s1);
    s1 = std::max(tmp,s1);
    return  floor((2.0 * l) / (s0 + s1) - 1); // + d;
  };
  
  inline double edge_saturation(int_t edge,
				int_t n,
				double xpts[],
				double ypts[],
				double tpts[])
  {
      
    auto edge_ptr  = get_edge_ptr(edge);
    auto node0_ptr = get_node_ptr(edge_ptr->start());
    auto node1_ptr = get_node_ptr(edge_ptr->stop());
    auto metric0   = &node0_ptr->m_metric[0];
    auto metric1   = &node1_ptr->m_metric[0];
    auto coo0      = &node0_ptr->x;
    auto coo1      = &node1_ptr->x;
    
    double 
      x0 = coo0[0],
      y0 = coo0[1],
      
      x1 = coo1[0],
      y1 = coo1[1];
    double s0  = 1.0/sqrt(metric0[0]);
    double s1  = 1.0/sqrt(metric1[0]);
    double l   = sqrt( (x0-x1)*(x0-x1)+(y0-y1)*(y0-y1));
    
    double neg = (s0 > s1) ? -1 : 1;
    if (n==1)
      {
	xpts[0] = (x0+x1)*0.5;
	ypts[0] = (y0+y1)*0.5;
	tpts[0] = 0.5;
	return neg;
      }
    //    std::cout << std::endl << "neg " << neg << std::endl;
    double tmp = s0;
    s0 = std::min(tmp,s1);
    s1 = std::max(tmp,s1);
    auto r =  2.0 *(l - n * s0 - s1) / ( n*(n+1) -2 );    
    double sum = 0.0;
    double a = s0 + r;
    //      sum=a;
    double t = (s0 + r);
    for (int i=1;i<n;++i)
      {
	//	std::cout << "t " << t/l << std::endl;
	double s = t / l;
	if (neg < 0)
	  {
	    s = 1.0 - t / l;
	  }
	double x = (1.0 - s) * x0 + s * x1;
	double y = (1.0 - s) * y0 + s * y1;
	xpts[i-1] = x;
	ypts[i-1] = y;
	tpts[i-1] = s;
	a+=r;
	t+=a;
      }
    
    {
      //      std::cout << "t " << t/l << std::endl;
      double s = t / l;
      if (neg < 0)
	{
	  s = 1.0 - t / l;
	}
      double x = (1.0 - s) * x0 + s * x1;
      double y = (1.0 - s) * y0 + s * y1;
      tpts[n-1] = s;
      xpts[n-1] = x;
      ypts[n-1] = y;
#if 0
      std::cout << "s" <<s << std::endl;
      std::cout << "l" <<l << std::endl;
      std::cout << "t" <<t << std::endl;
      std::cout << "x0 " << x0 << std::endl;
      std::cout << "y0 " << y0 << std::endl;
      std::cout << "x1 " << x1 << std::endl;
      std::cout << "y1 " << y1 << std::endl;
#endif
    }
    
#if 0
    std::cout << "X1 " << x1 << " " << y1 << std::endl;
    a = (s1-r);
    t+=a;
    {
      double s = t / l;
      if (neg < 0)
	{
	  s = 1.0 - t / l;
	}
      
      double x = (1.0 - s) * x0 + s * x1;
      double y = (1.0 - s) * y0 + s * y1;
      std::cout << "FINAL " << x << " " << y << std::endl;
      
    }
    std::cout << "t " << t/l << std::endl;
    exit(1);
#endif
    return neg;    
  }

  
  inline int_t opposite_node(int_t face,
			     int_t edge) const noexcept
  {
    auto face_ptr = get_triangle_ptr(face);
    auto edge_ptr = get_edge_ptr(edge);
    auto e = get_edge_ptr((edge == face_ptr->edge(0) ) ? face_ptr->edge(1) : face_ptr->edge(0));
    int_t startEdge 	= e->start();
    int_t endEdge 	= e->stop();
    return ( ((startEdge == edge_ptr->start()) || (startEdge == edge_ptr->stop())) ) ? endEdge : startEdge;
  };

  
public: inline  void FlipEdge(int_t edgeid)
  {    
    auto edge_ = this->get_edge_ptr(edgeid);
    int_t faceid0 = edge_->triangle(0);
    int_t faceid1 = edge_->triangle(1);
    triangle_t * face0 = (faceid0 > 0) ? &triangles[faceid0] : nullptr;
    triangle_t * face1 = (faceid1 > 0) ? &triangles[faceid1] : nullptr;
    /* REMOVE OLD TRIANGLES FROM ADT  */
    const int signedEdgeWayInFace0 = (edgeid==face0->edge(0)) ? face0->way(0)  : ( (edgeid==face0->edge(1)) ? face0->way(1) : face0->way(2));

    int_t faceNegativeId 	= (signedEdgeWayInFace0>0) ? faceid1 : faceid0;
    triangle_t * faceNegative = &triangles[faceNegativeId];
    
    int_t facePositiveId 	= (signedEdgeWayInFace0>0) ? faceid0 : faceid1;
    triangle_t * facePositive = &triangles[facePositiveId];
    
    const int 	edgeIndexInFacePositive	= (edgeid==facePositive->edge(0)) ? 0 : ( (edgeid==facePositive->edge(1)) ? 1 : 2);
    const int 	edgeIndexInFaceNegative	= (edgeid==faceNegative->edge(0)) ? 0 : ( (edgeid==faceNegative->edge(1)) ? 1 : 2);
    auto P_oppositeToFacePositive 		= opposite_node(facePositiveId,edgeid);
    auto P_oppositeToFaceNegative 		= opposite_node(faceNegativeId,edgeid);

    //
    // REMOVE ITSELF FROM THE LIST OF EDGES STORED IN EACH VERTEX
    //
    //      edge_->start()->RemoveEdge(edge_);
    //      edge_->stop()->RemoveEdge(edge_);

    //
    // REDEFINE THE EDGE
    //
    edge_->SetStart(P_oppositeToFaceNegative);
    edge_->SetEnd(P_oppositeToFacePositive);
      
    //
    // INSERT ITSELF IN THE LIST OF EDGES STORED IN EACH VERTEX
    //
    //      edge_->start()->AddEdge(edge_);
    //      edge_->stop()->AddEdge(edge_);

    auto  	nextEdgeInFacePositive		= facePositive->edge	( (edgeIndexInFacePositive+1)%3 );
    auto  	previousEdgeInFacePositive	= facePositive->edge	( (edgeIndexInFacePositive+2)%3 );
    const auto nextEdgeWayInFacePositive	= facePositive->way	( (edgeIndexInFacePositive+1)%3 );
    const auto previousEdgeWayInFacePositive	= facePositive->way	( (edgeIndexInFacePositive+2)%3 );
    auto  	nextEdgeInFaceNegative		= faceNegative->edge	( (edgeIndexInFaceNegative+1)%3 );
    auto  	previousEdgeInFaceNegative	= faceNegative->edge	( (edgeIndexInFaceNegative+2)%3 );
    const auto nextEdgeWayInFaceNegative	= faceNegative->way	( (edgeIndexInFaceNegative+1)%3 );
    const auto previousEdgeWayInFaceNegative	= faceNegative->way	( (edgeIndexInFaceNegative+2)%3 );
      
    //
    // REDEFINE THE FIRST TRIANGLE
    // 
    facePositive->edge(0,previousEdgeInFacePositive);
    facePositive->edge(1,nextEdgeInFaceNegative);
    facePositive->edge(2,edgeid);
    facePositive->way(0, (previousEdgeWayInFacePositive>0) ? 1 : -1);
    facePositive->way(1, (nextEdgeWayInFaceNegative>0) ? 2 : -2 );
    facePositive->way(2, 3);      

    //
    // REDEFINE THE SECOND TRIANGLE
    // 
    faceNegative->edge(0,previousEdgeInFaceNegative);
    faceNegative->edge(1,nextEdgeInFacePositive);
    faceNegative->edge(2,edgeid);
      
    faceNegative->way(0, (previousEdgeWayInFaceNegative>0) ? 1 : -1);
    faceNegative->way(1, (nextEdgeWayInFacePositive>0) ? 2 : -2);
    faceNegative->way(2, -3);
      
    // 
    // REDEFINE ADJACENCY OF THE BOUNDARY EDGES
    // FOR THE FLIPPED EDGE, WE DON'T NEED, THE TRIANGLE ARE THE SAME
    // 	      
    m_edges[nextEdgeInFacePositive].triangle( (facePositiveId == m_edges[nextEdgeInFacePositive].triangle(0)) ? 0 : 1,
					      faceNegativeId);		
    m_edges[nextEdgeInFaceNegative].triangle( (faceNegativeId == m_edges[nextEdgeInFaceNegative].triangle(0)) ? 0 : 1,
						 facePositiveId);
  };


  
  inline bool is_able_to_flip(int_t edge) const
  {
    auto edge_ptr 	= this->get_edge_ptr(edge);
    auto t0 		= edge_ptr->triangle(0);
    auto t1 		= edge_ptr->triangle(1);
    if (t0 && t1)
      {
	auto p0 		= edge_ptr->start();
	auto p1 		= edge_ptr->stop();
	auto p_t0 		= opposite_node(t0,edge);
	auto p_t1 		= opposite_node(t1,edge);
	  
	auto side0		= ComputeOrientation(xy[p_t0].x,
						     xy[p_t0].y,
						     xy[p_t1].x,
						     xy[p_t1].y,
						     xy[p0].x,
						     xy[p0].y);
	
	auto side1		= ComputeOrientation(xy[p_t0].x,
						     xy[p_t0].y,
						     xy[p_t1].x,
						     xy[p_t1].y,
						     xy[p1].x,
						     xy[p1].y);
	
	return (( (side0==Orientation::FORWARD)  && (side1==Orientation::BACKWARD) )
		|| 
		( (side0==Orientation::BACKWARD) && (side1==Orientation::FORWARD) ));	
      }
    else
      {	
	std::cout << "unable to flip edge" << std::endl;
	exit(1);
	return false;
      }
  };
  

  int_t boundary_marker(const bool&keepBoundingBox_)
  {
    int_t edgeBox = 0;
    for (int_t edge = 1;edge <= m_nedges;++edge)
      {
	auto edge_ptr   = this->get_edge_ptr(edge);	
	auto p0 	= edge_ptr->start();
	auto p1 	= edge_ptr->stop();
	if (false == edge_ptr->is_boundary())
	  {
	    if ( this->m_hashEdges.Get(p0,p1) )
	      {
		edge_ptr->is_boundary(true);
		if (!this->m_hashEdges.Remove(p0,p1))
		  {
		    printf("boundarmarker no deleted\n");
		    exit(1);
		  }
		edge_ptr->is_nonflipable(true);
	      }
	  }
	
	if ( (p0 <= 4) && (p1 <= 4) )
	  {
	    edgeBox = edge;
	    edge_ptr->is_nonflipable(true);	 
	    edge_ptr->is_boundary(keepBoundingBox_);
	  }
      }
    
    return edgeBox;	
  };



  inline bool DoIntersect(int_t a,int_t b,int_t p,int_t q) const
  {
    return DoIntersect(xy[a].x,
		       xy[a].y,
		       xy[b].x,
		       xy[b].y,
		       xy[p].x,
		       xy[p].y,
		       xy[q].x,
		       xy[q].y);
  }
  inline bool DoIntersect(double ax,
			  double ay,
			  double bx,
			  double by,
			  double px,
			  double py,
			  double qx,
			  double qy) const
  {      
    const Orientation::Kind side_a = ComputeOrientation(px,py,qx,qy,ax,ay);
    const Orientation::Kind side_b = ComputeOrientation(px,py,qx,qy,bx,by);
    bool ans = false;
    switch(side_a)
      {
      case  Orientation::FORWARD:
	{
	  switch(side_b)
	    {
	    case  Orientation::FORWARD:
	      {
		ans = false;
		break;
	      }
	    case  Orientation::BACKWARD:
	      {
		ans = true;
		break;
	      }
	    case  Orientation::NEUTRAL:
	      {
		ans = false;
		break;
	      }
	    }
	  break;
	}
      case  Orientation::BACKWARD:
	{
	  switch(side_b)
	    {
	    case  Orientation::FORWARD:
	      {
		ans = true;
		break;
	      }
	    case  Orientation::BACKWARD:
	      {
		ans = false;
		break;
	      }
	    case  Orientation::NEUTRAL:
	      {
		ans = false;
		break;
	      }
	    }
	  break;
	}
      case  Orientation::NEUTRAL:
	{
	  switch(side_b)
	    {
	    case  Orientation::FORWARD:
	      {
		ans = false;
		break;
	      }
	    case  Orientation::BACKWARD:
	      {
		ans = false;
		break;
	      }
	    case  Orientation::NEUTRAL:
	      {
		ans = false;
		break;
	      }
	    }
	  break;
	}
      }

    if (ans == false)
      {
	return ans;
      }
    else
      {
	  
	const Orientation::Kind side_p = ComputeOrientation(ax,ay,
							    bx,by,
							    px,py);
	  
	const Orientation::Kind side_q = ComputeOrientation(ax,ay,
							    bx,by,
							    qx,qy);
	  
	switch(side_p)
	  {
	  case  Orientation::FORWARD:
	    {
	      switch(side_q)
		{
		case  Orientation::FORWARD:
		  {
		    ans = false;
		    break;
		  }
		case  Orientation::BACKWARD:
		  {
		    ans = true;
		    break;
		  }
		case  Orientation::NEUTRAL:
		  {
		    ans = false;
		    break;
		  }
		}
	      break;
	    }
	  case  Orientation::BACKWARD:
	    {
	      switch(side_q)
		{
		case  Orientation::FORWARD:
		  {
		    ans = true;
		    break;
		  }
		case  Orientation::BACKWARD:
		  {
		    ans = false;
		    break;
		  }
		case  Orientation::NEUTRAL:
		  {
		    ans = false;
		    break;
		  }
		}
	      break;
	    }
	  case  Orientation::NEUTRAL:
	    {
	      switch(side_q)
		{
		case  Orientation::FORWARD:
		  {
		    ans = false;
		    break;
		  }
		case  Orientation::BACKWARD:
		  {
		    ans = false;
		    break;
		  }
		case  Orientation::NEUTRAL:
		  {
		    ans = false;
		    break;
		  }
		}
	      break;
	    }
	  }
	return ans;
      }
  };

  inline bool try_edge_flip(int_t edge)
  {
    auto edge_ptr 		= this->get_edge_ptr(edge);
    auto isNotABoundaryEdge 	= false == edge_ptr->is_boundary();
    auto isAFlipableEdge 	= is_able_to_flip(edge);    
    if ( isNotABoundaryEdge && isAFlipableEdge )
      {  	
	FlipEdge(edge);	  
	return true;
      }
    return false;
  };


  
  inline unsigned int force_edge(int_t a,int_t b)
  {
    //
    // FIND ALGORITHM
    //      
    std::set<int_t> edgesToFlip;    
    for (int_t triangle = 1;triangle <= m_ntriangles;++triangle)
      {
	auto triangle_ptr = get_triangle_ptr(triangle);
	for (unsigned int edgeIndex=0;edgeIndex<3;++edgeIndex)
	  {
	    int_t edge = triangle_ptr->edge(edgeIndex);
	    auto edge_ptr = this->get_edge_ptr(edge);
	    if ( (false == m_edges[edge].is_boundary()) && DoIntersect(a,
								      b,
								      edge_ptr->start(),
								      edge_ptr->stop() ) )
	      {
		edgesToFlip.insert(edge);
	      }
	  }	
      }

      
    unsigned int numFlips = 0;
    bool hasBeenFlip = true;
    while (hasBeenFlip)
      {
	hasBeenFlip = false;
	for (typename std::set<int_t>::iterator it = edgesToFlip.begin(); it != edgesToFlip.end(); it++)
	  {
	    int_t edge = *it;
	    auto edge_ptr = this->get_edge_ptr(edge);
	    hasBeenFlip = try_edge_flip(edge);
	    if (hasBeenFlip)
	      {
		print_mesh();
		int_t startingEdge 	= edge_ptr->start();
		int_t endingEdge 	= edge_ptr->stop();	
		++numFlips;		  
		if (this->m_hashEdges.Get(startingEdge, endingEdge))
		  {
		    edge_ptr->is_boundary(true);
		    edge_ptr->is_nonflipable(true);
		  }
		
		if (false == DoIntersect(a,
					 b,
					 startingEdge,
					 endingEdge))
		  {
		    edgesToFlip.erase(it);
		    break;
		  }
	      }	      
	  }
	
	for (typename std::set<int_t>::iterator it = edgesToFlip.begin(); it != edgesToFlip.end(); it++)
	  {
	    int_t edge	= *it;
	    auto edge_ptr = this->get_edge_ptr(edge);
	    hasBeenFlip = try_edge_flip(edge);
	    if (hasBeenFlip)
	      {
		int_t startingEdge 	= edge_ptr->start();
		int_t  endingEdge 	= edge_ptr->stop();	
		++numFlips;		  
		if (this->m_hashEdges.Get(startingEdge,endingEdge))
		  {
		    edge_ptr->is_boundary(true);
		    edge_ptr->is_nonflipable(true);
		  }
		
		if (false == DoIntersect(a,
					 b,
					 startingEdge,
					 endingEdge))
		  {
		    edgesToFlip.erase(it);
		    break;
		  }
	      }	      
	  }
      }
    //    std::cout << "numFlips " << numFlips << std::endl;
    return numFlips;
  };
  
  inline unsigned int ForceEdges()
  {
    unsigned int numFlips = 0;
    for (typename HashPair<unsigned int,unsigned int>::MapPairIterator it = this->m_hashEdges.begin();
	 it != this->m_hashEdges.end();
	 it++)
      {
	typename HashPair<unsigned int,unsigned int>::MapPairKey key = it->first;
	
	numFlips += force_edge(key.first,
			       key.second);
      }
    return numFlips;
  };


  
  int coloring(int * colors,int_t seed_edge)
  {
#if 0
    for (int i=1;i<=m_nedges;++i)
      {
	auto e = get_edge_ptr(i);
	if (e->is_boundary())
	  {
	    std::cout << "ffffffffffffffffffffffffff " << e->start() << " " << e->stop() << std::endl;
	  }
      }
#endif
    for (int i=1;i<=m_ntriangles;++i)
      {
	colors[i]=0;
      }
    //      std::cout << "coloring" << std::endl;
    auto seed_triangle = get_edge_ptr(seed_edge)->triangle(0);

    auto triangle = seed_triangle;
    bool all = false;
    int ncolors = 1;
    std::list<int_t> lst;
    do
      {
	all = false;	
	colors[triangle] = ncolors;
	for (int i=0;i<3;++i)
	  {
	    auto e = get_edge_ptr(get_triangle_ptr(triangle)->edge(i));
	    //	  std::cout << e->is_boundary() << std::endl;
	      
	    if (!e->is_boundary())
	      {
		int_t nei_triangle = e->otriangle(triangle);
		if (nei_triangle && !colors[nei_triangle])
		  {
		    lst.push_front(nei_triangle);
		  }
	      }
	    else
	      {
		//		  std::cout << "ffffffffffffffffffffffffff " << e->start() << " " << e->stop() << std::endl;
	      }
	  }
	if (!lst.empty())
	  {
	    all = true;
	    triangle = lst.front();
	    lst.pop_front();
	  }
      } while (all);

    triangle = 0;      
    for (int_t i=1;i<=m_ntriangles;++i)
      {
	if (!colors[i])
	  {
	    triangle = i;
	    break;
	  }
      }
      
    while (triangle)
      {
	all = false;
	++ncolors;
	lst.clear();
	do
	  {
	    all = false;	
	    colors[triangle] = ncolors;
	    for (int i=0;i<3;++i)
	      {
		auto e = get_edge_ptr(get_triangle_ptr(triangle)->edge(i));
		//	  std::cout << e->is_boundary() << std::endl;
		  
		if (!e->is_boundary())
		  {
		    int_t nei_triangle = e->otriangle(triangle);
		    if (nei_triangle && !colors[nei_triangle])
		      {
			lst.push_front(nei_triangle);
		      }
		  }
		else
		  {
		    //		  std::cout << "ffffffffffffffffffffffffff " << e->start() << " " << e->stop() << std::endl;
		  }
	      }
	    if (!lst.empty())
	      {
		all = true;
		triangle = lst.front();
		lst.pop_front();
	      }
	  } while (all);



	triangle = 0;      
	for (int_t i=1;i<=m_ntriangles;++i)
	  {
	    if (!colors[i])
	      {
		triangle = i;
		break;
	      }
	  }
      }
      

      
#if 0
    for (int i=1;i<=m_ntriangles;++i)
      {
	std::cout << "color[" << i << "] = "  << colors[i] << std::endl;
      }
#endif     
      
#if 0
    bool all = false;
    do
      {
	all = false;
	bool assigned = false;
	do
	  {
	    assigned = false;
	    for (int_t triangle = 1;triangle <= m_ntriangles;++triangle)
	      {
		if (!colors[triangle])
		  {
		    bool same_color = true;
		    for (int i=0;i<3;++i)
		      {
			auto e = get_edge_ptr(get_triangle_ptr(triangle)->edge(i));
			if (!e->is_boundary())
			  {
			    int_t nei_triangle = e->otriangle(triangle);
			    if (nei_triangle && colors[nei_triangle]>0)
			      {
				same_color = false;
				break;
			      }
			  }
		      }
		    if (same_color)
		      {
			colors[triangle] = ncolors;
			assigned = true;
		      }
		  }
	      }
	    
	  } while (assigned);

#if 0
	for (int_t triangle = 1;triangle <= m_ntriangles;++triangle)
	  std::cout << "color[" << triangle << "] = " << colors[triangle] << std::endl;
	for (int_t triangle = 1;triangle <= m_ntriangles;++triangle)
	  {
	    if (!colors[triangle])
	      {
		std::cout << "change color " << std::endl;
		++ncolors;
		all = true;
		break;
	      }
	  }
#endif
	break;
      } while (all);
#endif
    
    //      std::cout << "coloring done " << ncolors << std::endl;
    return ncolors;

  }

  int_t boundary_marker(const delaunay_cad_t& cad_)
  {
    for (int_t cad_edge=1;cad_edge<=cad_.nedges();++cad_edge)
      {
	this->m_hashEdges.Add(cad_.node(cad_edge,0) + 4,
			      cad_.node(cad_edge,1) + 4,
			      cad_edge);
      }      
    return boundary_marker(false);
  };
  
  inline void retrieve_boundary(const delaunay_cad_t& cad_)
  {
    int_t edgeBox 	        = this->boundary_marker(cad_);  
    unsigned int numFlips 	= this->ForceEdges();
    int * colors                = new int[m_ntriangles+1];
    int ncolors                 = coloring(colors,edgeBox);
    for (int_t i = 1;i <= this->m_ntriangles;++i)
      {
	this->get_triangle_ptr(i)->id(colors[i]);
      }
    delete[]colors;
    
    print_mesh();      
  };
    
  inline void remove_box(const bool&keepBoundingBox_ = false)
  {
    int_t color = keepBoundingBox_ ? 2 : 1;
    for (int_t i =1;i<=m_ntriangles;++i)
      {
	int_t icolor = this->get_triangle_ptr(i)->id();
	if (icolor == color)
	  {
	    //
	    // Remove triangle
	    //
	    this->weak_remove_triangle(i);
	  }
      }
      
    while (false == m_stack.empty())
      {
	auto a = m_stack.top();
	m_stack.pop();
	remove_face(a);	  
      }
    
    //
    // 
    //
    //
    for (int_t edge=1;edge<=m_nedges;++edge)
      {
	if (!m_edges[edge].triangle(0))
	  {
	    weak_remove_edge(edge);
	  }
      }
      

    while (false == m_stack_edges.empty())
      {
	auto a = m_stack_edges.top();
	m_stack_edges.pop();
	remove_edge(a);
      }

    print_mesh();
    
  };

  inline unsigned int RemoveBigBox(const bool&keepBoundingBox_)
  {

    int_t edgeBox 	        = this->boundary_marker(keepBoundingBox_);  
    unsigned int numFlips 	= this->ForceEdges();
    int * colors = new int[m_ntriangles+1];
    int ncolors = coloring(colors,edgeBox);
    for (int_t i =1;i<=m_ntriangles;++i)
      {
	this->get_triangle_ptr(i)->id(colors[i]);
      }
    
    delete[]colors;
    print_mesh();
    return numFlips;
  };

  static inline double ComputeSignedArea(double x0,double y0,
					 double x1,double y1,
					 double x2,double y2) 
  {
    const auto  signed_area = (x0*(y1-y2) - x1*(y0-y2) + x2*(y0-y1))*double(0.5);
    return Math<double>::Filter(signed_area,2.22044604925031308e-16);
  };

  static inline Orientation::Kind ComputeOrientation(double x0,double y0,
						     double x1,double y1,
						     double x2,double y2) 
  {
    return Orientation::Create(ComputeSignedArea(x0,y0,x1,y1,x2,y2));
  };
  
  //! @brief Determine if a point belongs to the triangle defined by three points.
  //! @param p_ The point.
  //! @param a_ The first point of the triangle.
  //! @param b_ The second point of the triangle.
  //! @param c_ The third point of the triangle.
  //! @return true if the point is inside or on the boundary, false otherwise.
  static inline bool is_inside(double x,double y,
			       double x0,double y0,
			       double x1,double y1,
			       double x2,double y2) 
  {
    if (Orientation::BACKWARD==ComputeOrientation(x0,y0,x1,y1,x,y)) return false;
    if (Orientation::BACKWARD==ComputeOrientation(x1,y1,x2,y2,x,y)) return false;
    if (Orientation::BACKWARD==ComputeOrientation(x2,y2,x0,y0,x,y)) return false;
    return true;
  };      

  void triangle_coordinates(int_t id,double coo[])
  {
    int_t e0 = triangles[id].edge(0);
    int_t w0 = triangles[id].way(0);

    coo[0] = xy[m_edges[e0].start(w0)].x;
    coo[1] = xy[m_edges[e0].start(w0)].y;

    coo[2] = xy[m_edges[e0].stop(w0)].x;
    coo[3] = xy[m_edges[e0].stop(w0)].y;

    int_t e1 = triangles[id].edge(1);
    int_t w1 = triangles[id].way(1);
    
    coo[4] = xy[m_edges[e1].stop(w1)].x;
    coo[5] = xy[m_edges[e1].stop(w1)].y;    
  };

  int_t find_triangle(int_t p,int_t seed)
  {
    check(__LINE__);
    double coo[6];
    double x = xy[p].x;
    double y = xy[p].y;
#if 0
    for (int_t i=1;i<=m_ntriangles;++i)
      {
	seed = i;
    	triangle_coordinates(seed,coo);
	
	double x0 = coo[0];
	double y0 = coo[1];
	double x1 = coo[2];
	double y1 = coo[3];
	double x2 = coo[4];
	double y2 = coo[5];
	//	std::cout << "seed " << seed << std::endl;
	//	next = false;
	Orientation::Kind oo[3];
	oo[0] = ComputeOrientation(x0,y0,x1,y1,x,y);
	oo[1] = ComputeOrientation(x1,y1,x2,y2,x,y);
	oo[2] = ComputeOrientation(x2,y2,x0,y0,x,y);
	//	std::cout << "seed " << oo[0] <<  " " << oo[1] << " " << oo[2] << std::endl;

	if (Orientation::BACKWARD!=oo[0] && Orientation::BACKWARD!=oo[1] && Orientation::BACKWARD!=oo[2])
	  {
	    //	    std::cout << "yo " << std::endl;
	    return seed;
	  }
      }
    return 0;
#else
    bool next;

    do
      {
	//	std::cout << "try find "  << seed << std::endl;
    	triangle_coordinates(seed,coo);
	
	double x0 = coo[0];
	double y0 = coo[1];
	double x1 = coo[2];
	double y1 = coo[3];
	double x2 = coo[4];
	double y2 = coo[5];
	//	std::cout << "seed " << seed << std::endl;
	next = false;
	Orientation::Kind oo[3];
	oo[0] = ComputeOrientation(x0,y0,x1,y1,x,y);
	oo[1] = ComputeOrientation(x1,y1,x2,y2,x,y);
	oo[2] = ComputeOrientation(x2,y2,x0,y0,x,y);
	//	std::cout << "seed " << oo[0] <<  " " << oo[1] << " " << oo[2] << std::endl;

	if (Orientation::BACKWARD!=oo[0] && Orientation::BACKWARD!=oo[1] && Orientation::BACKWARD!=oo[2])
	  {
	    //	    std::cout << "yo " << std::endl;
	    return seed;
	  }
	else
	  {
	    int_t nei[3];
	    nei[0] = get_edge_ptr(get_triangle_ptr(seed)->edge(0))->otriangle(seed);
	    nei[1] = get_edge_ptr(get_triangle_ptr(seed)->edge(1))->otriangle(seed);
	    nei[2] = get_edge_ptr(get_triangle_ptr(seed)->edge(2))->otriangle(seed);

	    
	    
	    for (int i=0;i<3;++i)
	      {
		//		std::cout << " " << oo[i] << " " << nei[i] << std::endl;
		if ( (oo[i] == Orientation::BACKWARD) && nei[i])
		  {
		    seed = nei[i];
		    next = true;
		    break;
		  }
	      }


	    

	    
	  }
	//	if (Orientation::BACKWARD==) return false;
	//	if (Orientation::BACKWARD==ComputeOrientation(x1,y1,x2,y2,x,y)) return false;
	//	if (Orientation::BACKWARD==ComputeOrientation(x2,y2,x0,y0,x,y)) return false;
	
      } while (next);
    //        std::cout << "done "  << std::endl;

    return 0;
#endif
#if 0
    for (int j = 1;j <= m_ntriangles;++j)
      {
	triangle_coordinates(j,coo);
	if (is_inside(x,y,
		      coo[0],coo[1],
		      coo[2],coo[3],
		      coo[4],coo[5]))
	  {
	    return j;
	  }	
      }
#endif
    return 0;
  };


  int_t find_triangle(int_t p)
  {
    return find_triangle(p,1);
  };

  //!
  //! Constructor.
  //!
  mesh_t()
  {
    sprintf(m_debug_filebasename,"initial");
    mxntriangles = 4096 * 200;
    mxnedges     = 4096 * 200;
    mxnnodes     = 4096 * 200;
    m_nnodes= 0;
    m_nedges=0;
    m_ntriangles=0;

    m_edges 	= (edge_t*)	calloc( mxnedges,
					sizeof(edge_t) );

    xy 		= (xy_t*)	calloc( mxnnodes,
					sizeof(xy_t) );
    
    triangles 	= (triangle_t*)	calloc( mxntriangles,
					sizeof(triangle_t) );
    
  };
protected:
public:
  

  inline void debug_view()
  {
#if 0
    std::cout << "debug view " << std::endl;
    for (int_t edge = 1; edge <= m_nedges;++edge)
      {
	std::cout << "neighbors edge " << edge << " : " << m_edges[edge].triangle(0) << " " << m_edges[edge].triangle(1) << std::endl;
      }
#endif
  };
  
  void triangle_to_nodes(int_t triangle,int_t nodes[3]) const noexcept
  {
    auto triangle_ptr 	= get_triangle_ptr(triangle);
    auto edge0_ptr 	= get_edge_ptr(triangle_ptr->edge(0));
    auto edge1_ptr 	= get_edge_ptr(triangle_ptr->edge(1));
    nodes[0] 		= edge1_ptr->stop(triangle_ptr->way(1));
    auto way0  		= triangle_ptr->way(0);
    nodes[1] 		= edge0_ptr->start(way0);
    nodes[2] 		= edge0_ptr->stop(way0);
  };
  
public:




  

  //
  // Check conformity.
  //
  inline void check(int line)
  {
    return;
    //
    // out of bounds triangles ? 
    // 
    for (int_t edge = 1;edge <= m_nedges;++edge)
      {
	for (unsigned int i=0;i<2;++i)
	  {
	    auto t = m_edges[edge].triangle(i);
	    if (t < 0 || t > m_ntriangles)
	      {
		std::cerr << "check(line " << line << ") invalid edge" << std::endl;
		std::cerr << "m_ntriangles : " << m_ntriangles << std::endl;
		print_cnc();
		exit(1);
	      }
	    if (t>0)
	      {
		if (-1 == triangles[t].get_local_edge_index(edge))
		  {
		    std::cerr << "check(line " << line << ") invalid edge" << std::endl;
		    std::cerr << "m_ntriangles : " << m_ntriangles << std::endl;
		    print_cnc();
		    exit(1);
		  }
	      }
	  }
      }
    
    //
    // out of bounds edges ? 
    // 
    for (int_t triangle = 1;triangle <= m_ntriangles;++triangle)
      {
	for (unsigned int i=0;i<3;++i)
	  {
	    int_t edge = triangles[triangle].edge(i);
	    if (edge > m_nedges || edge == 0)
	      {
		std::cerr << "check(line " << line << ") out of bound edge triangle["<< triangle << "].edge["  << i << "] = " <<  edge << std::endl;
		print_cnc();
		exit(1);
	      }

	    //
	    // 
	    //
	    if (m_edges[edge].triangle(0)!=triangle && m_edges[edge].triangle(1) != triangle)
	      {
		std::cerr << "check(line " << line << ") invalid edge-triangle["<< triangle << "].edge["  << i << "] = " <<  edge << std::endl;
		print_cnc();
		exit(1);
	      }
	  }
      }



    
#if 0
    for (int_t edge = 1;edge <= m_nedges;++edge)
      {
	for (unsigned int i=0;i<2;++i)
	  {
	    auto t = m_edges[edge].triangle(i);
	    if (t < 0 || t > m_ntriangles)
	      {
		std::cerr << "check(line " << line << ") invalid edge" << std::endl;
		std::cerr << "m_ntriangles : " << m_ntriangles << std::endl;
		print_cnc();
		exit(1);
	      }
	  }
      }
    
    
    for (int_t triangleid = 1;triangleid <= m_ntriangles;++triangleid)
      {
	
	for (unsigned int i=0;i<3;++i)
	  {
	    
	    int_t edge = triangles[triangleid].edge(i);
	    if (edge > m_nedges || edge == 0)
	      {
		std::cerr << "check(line " << line << ") out of bound edge triangle["<< triangleid << "].edge["  << i << "] = " <<  edge << std::endl;
		print_cnc();
		exit(1);
	      }
#if 0
	    if (m_edges[edge].triangle(0) == triangleid)
	      {
		
	      }
	    else if (m_edges[edge].triangle(1) == triangleid)
	      {

	      }
	    else
	      {

		std::cerr << "triangle " << triangleid << std::endl;
		std::cerr << "triangle.edge[" << i << "] = " << edge << std::endl;
		std::cerr << "edge.triangle[0] = " << m_edges[edge].triangle(0) << std::endl;
		std::cerr << "edge.triangle[1] = " << m_edges[edge].triangle(1) << std::endl;
		std::cerr << "check(line " << line << ") invalid edge - triangle" << std::endl;
		print_cnc();
		//exit(1);
	      }
#endif
	  }
	
      }
#endif
  };

  
  inline void RemoveFaceFromLocalization(int_t face_){};

  
  inline void move_face(int_t face_from, int_t face_to)
  {
    if (face_from != face_to)
      {
	for (unsigned int k=0;k<3;++k)
	  {
	    update_edge_triangle(this->triangles[face_from].edge(k),
				 face_from,
				 face_to);
	  }
	this->triangles[face_to] = this->triangles[face_from];
      }
  };

  inline void move_edge(int_t edge_from, int_t edge_to)
  {
    if (edge_from != edge_to)
      {
	for (unsigned int k=0;k<2;++k)
	  {
	    int_t face_from = this->m_edges[edge_from].triangle(k);
	    update_triangle_edge(face_from, edge_from, edge_to);
	  }
	this->m_edges[edge_to] = this->m_edges[edge_from];
      }
  };

  inline void update_edge_triangle(int_t edge,
				   int_t old_triangle,
				   int_t new_triangle)
  {
    this->m_edges[edge].dereftriangle(old_triangle);
    this->m_edges[edge].reftriangle(new_triangle);
  }
  
  inline void update_triangle_edge(int_t face,
				   int_t old_edge,
				   int_t new_edge)
  {
    if (face)
      {
	for (unsigned int localEdgeIndex=0;localEdgeIndex<3;++localEdgeIndex)
	  {
	    if (this->triangles[face].edge(localEdgeIndex) == old_edge)
	      {
		this->triangles[face].edge(localEdgeIndex, new_edge);
	      }
	  }
      }
  };
  
  inline void remove_edge(int_t i)
  {
    move_edge(m_nedges, i);

#if 0
    //
    // remove edge.
    //
    m_classifier.remove(i);
    if (i != m_nedges && m_classifier.has(m_nedges))
      {
	m_classifier.add(i);
      }
#endif
      
    this->m_edges[m_nedges] = edge_t();
    --m_nedges;
  };

  inline void remove_face(int_t i)
  {
    for (int k=0;k<3;++k)
      {
	this->m_edges[triangles[i].edge(k)].dereftriangle(i);
      }
    
    move_face(m_ntriangles, i);
    this->triangles[m_ntriangles] = triangle_t();
    --m_ntriangles;
  };

  inline void remove_faces(int n,int_t faces[])
  {
      
    bool all = true;
    while (all)
      {
	all = false;
	for (int k=0;k<n;++k)
	  {
	    if (faces[k] > 0)
	      {
		int_t N = m_ntriangles;
		remove_face(faces[k]);
		if (faces[k] != N)
		  {
		    //		    std::cout << " " << std::endl;
		    for (int i=0;i<n;++i)
		      {
			if (faces[i] == N)
			  {
			    //			      std::cout << "???????????????????" << std::endl;
			    faces[i] = faces[k];
			    break;
			  }
		      }
		  }		
		faces[k] = 0;
		all = true;
		break;
	      }
	  }
      }
    //    print_cnc();
    //    exit(1);
  };
  

  inline void remove_facesold(int n,int_t faces[])
  {

    bool all = true;
    while (all)
      {
	all = false;
	for (int k=0;k<n;++k)
	  {
	    if (faces[k] > 0)
	      {
		int_t N = m_ntriangles;
		remove_face(faces[k]);
		if (faces[k] != N)
		  {
		    //		    std::cout << " " << std::endl;
		    for (int i=0;i<n;++i)
		      {
			if (faces[i] == N)
			  {
			    //			      std::cout << "???????????????????" << std::endl;
			    faces[i] = faces[k];
			    break;
			  }
		      }
		  }		
		faces[k] = 0;
		all = true;
		break;
	      }
	  }
      }
    //    print_cnc();
    //    exit(1);
  };
  
  inline void remove_edges(int n,int_t edges[],std::list<FaceEdge>& 		lst)
  {

    bool all = true;
    while (all)
      {
	all = false;
	for (int k=0;k<n;++k)
	  {
	    if (edges[k] > 0)
	      {
#if 0
		std::cout << " remove = " << edges[k] << std::endl;
		std::cout << get_edge_ptr(edges[k])->start() << " " << get_edge_ptr(edges[k])->stop() << " " << std::endl;
#endif
		int_t N = m_nedges;
		remove_edge(edges[k]);
		if (edges[k] != N)
		  {

		    for (auto& f : lst)
		      {
			auto e = f.first;
			auto w = f.second;
			if (e == N)
			  {
			    f.first = edges[k];
			  }
		      }
		    
		    for (int i=0;i<n;++i)
		      {
			if (edges[i] == N)
			  {
			    edges[i] = edges[k];
			    break;
			  }
		      }
		    
		  }		
		edges[k] = 0;
		all = true;
		break;
	      }
	  }
      }
  };

  
  inline void RemoveFace(int_t i,int_t id)
  {
    //      fprintf(stdout,"BEFORE REMOVE FACE\n");
    //      print_cnc();
    this->m_edges[triangles[i].edge(0)].dereftriangle(i);
    this->m_edges[triangles[i].edge(1)].dereftriangle(i);
    this->m_edges[triangles[i].edge(2)].dereftriangle(i);
    const unsigned int n = m_ntriangles--;
    if (i!=n)
      {
#if 1
	this->m_edges[triangles[n].edge(0)].dereftriangle(n);
	this->m_edges[triangles[n].edge(1)].dereftriangle(n);
	this->m_edges[triangles[n].edge(2)].dereftriangle(n);
#endif
	
	this->triangles[i] = this->triangles[n];
	this->triangles[n] = triangle_t();
	//	this->m_edges[triangles[i].edge(0)].dereftriangle(i);
	//	this->m_edges[triangles[i].edge(1)].dereftriangle(i);
	//	this->m_edges[triangles[i].edge(2)].dereftriangle(i);
#if 1
	this->m_edges[triangles[i].edge(0)].reftriangle(i);
	this->m_edges[triangles[i].edge(1)].reftriangle(i);
	this->m_edges[triangles[i].edge(2)].reftriangle(i);
#endif
      }
    else
      {
	this->triangles[n] = triangle_t();
      }

    
    //      fprintf(stdout,"AFTER REMOVE FACE\n");
    //      print_cnc();
  };

  

  inline void RemoveFaces(int_t n,
			  int_t faces[])
  {
    
    int_t face;
    //
    // Remove edge->triangle adjacencies.
    //
    for (int i=0;i<n;++i)
      {
	face = faces[i];    
	this->m_edges[triangles[face].edge(0)].dereftriangle(face);
	this->m_edges[triangles[face].edge(1)].dereftriangle(face);
	this->m_edges[triangles[face].edge(2)].dereftriangle(face);
      }
    
    for (int i=0;i<n;++i)
      {
	face = faces[i];
	if (face>0)
	  {
	    this->m_edges[triangles[m_ntriangles].edge(0)].dereftriangle(m_ntriangles);
	    this->m_edges[triangles[m_ntriangles].edge(1)].dereftriangle(m_ntriangles);
	    this->m_edges[triangles[m_ntriangles].edge(2)].dereftriangle(m_ntriangles);
	    if (face != m_ntriangles)
	      {
		this->m_edges[triangles[m_ntriangles].edge(0)].reftriangle(face);
		this->m_edges[triangles[m_ntriangles].edge(1)].reftriangle(face);
		this->m_edges[triangles[m_ntriangles].edge(2)].reftriangle(face);

		this->triangles[face] = this->triangles[m_ntriangles];

		for (int j=0;j<n;++j)
		  {
		    if (faces[j] == m_ntriangles)
		      {
			faces[j] = face;
			break;
		      }
		  }       	
	      }	    
	    this->triangles[m_ntriangles] = triangle_t();
	    --m_ntriangles;
	    faces[i] = 0;
	  }
	//	  std::cout << "intermeiate" << std::endl;
	//	  print_cnc();
      }      
  };

  inline void RemoveFace(int_t i)
  {
    RemoveFace(i,i);
  };


  inline void Rebuild(int_t point_,
		      CavityOperator& cavityOperator_)
  {      
    
    cavityOperator_.clear_interior_edges();
    cavityOperator_.clear_hash();
    cavityOperator_.reset_faces();
      
    //    typename std::list<FaceEdge>::iterator it = cavityOperator_.m_listCavityBoundaryEdges.begin();
    FaceEdge          triangleEdge 	= cavityOperator_.boundary_edge(0);
    int_t	        edge		= triangleEdge.first;
    int 		way		= triangleEdge.second;
    if (way)
      {
	int_t  startingEdge	= m_edges[edge].start();
	int_t  endingEdge 	= m_edges[edge].stop();

	int_t edgep	= (way>0) ? create_edge(point_,startingEdge) 	: create_edge(point_,endingEdge);
	int_t edgen	= (way>0) ? create_edge(endingEdge,point_) 	: create_edge(startingEdge,point_);

	int_t nodepstart = m_edges[edgep].start();
	int_t nodepend   = m_edges[edgep].stop();
	int_t nodenstart = m_edges[edgen].start();
	int_t nodenend   = m_edges[edgen].stop();
	
	cavityOperator_.add_edge(nodepstart,
				 nodepend,
				 edgep);
	
	cavityOperator_.add_edge(nodenstart,
				 nodenend,	 
				 edgen);
	
	cavityOperator_.add_interior_edge(edgep);
	cavityOperator_.add_interior_edge(edgen);

	int_t face = create_triangle(edgep,1,edge,(way>0) ? 2 : -2,edgen,3,__LINE__);	  

	cavityOperator_.AddFace(face);
      }
    for (int_t i=1;i<cavityOperator_.nboundary_edges();++i)
      //     for (it++; it != cavityOperator_.m_listCavityBoundaryEdges.end(); it++)
      {

	FaceEdge	triangleEdge 	= cavityOperator_.boundary_edge(i);//*it;

	int_t 		edge		= triangleEdge.first;
	int 		way 		= (triangleEdge.second>0) ? 1 : -1;

	int_t		p0 		= m_edges[edge].start(triangleEdge.second);
	int_t		p1 		= m_edges[edge].stop(triangleEdge.second);

	//	  std::cout << "CAVITY EDGE IDS " << p0 << " " << p1 << std::endl;	  
	//	  std::cout << "THE THREE POINTS " << point_ << " " << p0 << " " << p1 << std::endl; 
	int_t 		edgep		= cavityOperator_.try_get_edge(point_,p0);
	int_t 		edgen		= cavityOperator_.try_get_edge(point_,p1);
	const int 	wayp 		= (edgep) ? -3 : 3;
	const int 	wayn 		= (edgen) ? -2 : 2;
	if (!edgep)
	  {
	    cavityOperator_.add_edge(point_,
					     p0,
					     edgep = create_edge(point_,p0));
	    cavityOperator_.add_interior_edge(edgep);
	  }
	    
	if (!edgen)
	  {
	    cavityOperator_.add_edge(p1,
					     point_,
					     edgen = create_edge(p1, point_));
	    
	    cavityOperator_.add_interior_edge(edgen);
	  }
	int_t face = create_triangle(edge,way,
				     edgen,wayn,
				     edgep,wayp,__LINE__);
	cavityOperator_.AddFace(face);
      }
      
  };

  //!
  //! @Remove the cavity from the mesh.
  //!
  mesh_t& operator += (CavityOperator& cavityOperator_)
  {
    
    //
    // Rebuild the cavity.
    //

#if 0
    //
    // Add faces.
    //
    for (auto face : this->m_cavityOperator.m_faces)
      {
	//	this->AddFace(face);
      }
#endif    
  };
  
  //!
  //! @Remove the cavity from the mesh.
  //!
  mesh_t& operator -= (CavityOperator& cavityOperator_)
  {
    //      m_stack.clear();
    cavityOperator_.apply_faces([this](int_t f) { weak_remove_triangle(f);});
    //      for (auto id : cavityOperator_.m_faces)
    //	{
    //	  weak_remove_triangle(id);
    //	}
    //      
    cavityOperator_.reset_faces();

    for (int_t i=0;i<cavityOperator_.ninterior_edges();++i)
      //    for (auto id : cavityOperator_.m_listCavityInteriorEdges)
      {
	auto id = cavityOperator_.interior_edge(i);
	weak_remove_edge(id);
      }
    cavityOperator_.clear_interior_edges();
    //    cavityOperator_.m_listCavityInteriorEdges.clear();
    return *this;
  };

#if 0
  void remove_triangle_color(int_t id_)
  {
    m_stack.push(id_);
    for (int k=0;k<3;++k)
      {
	int_t other = this->m_edges[triangles[id_].edge(k)].otriangle(id_);
	this->m_edges[triangles[id_].edge(k)].dereftriangle(id_);
	if (triangles[other].id() == triangles[id_].id())
	  {
	    remove_edge();
	  }	  
      }
  };
#endif
    
    


  static inline Orientation::Kind ComputeCircleOrientation(double x0,
							   double y0,
							   double x1,
							   double y1,
							   double x2,
							   double y2,
							   double x,
							   double y)
  { 
    const double  det1 	= (x0-x)* ( (y1-y) * ( (x2*x2-x*x)+(y2*y2-y*y) ) - ( (x1*x1-x*x)+(y1*y1-y*y) )* (y2-y) );
    const double  det2 	= (x1-x)* ( (y0-y) * ( (x2*x2-x*x)+(y2*y2-y*y) ) - ( (x0*x0-x*x)+(y0*y0-y*y) )* (y2-y) );
    const double  det3 	= (x2-x)* ( (y0-y) * ( (x1*x1-x*x)+(y1*y1-y*y) ) - ( (x0*x0-x*x)+(y0*y0-y*y) )* (y1-y) );
    return Orientation::Create(det1-det2+det3);
  };

  static inline void Circumcenter(double ax_,
				  double ay_,
				  
				  double bx_,
				  double by_,
				  
				  double cx_,
				  double cy_,
				  double * px_,
				  double * py_)
  {
    double cax 	= ax_ - cx_;
    double cay 	= ay_ - cy_;    
    double cbx 	= bx_ - cx_;
    double cby 	= by_ - cy_;    
    double mcax = (ax_ + cx_) * 0.5;
    double mcay = (ay_ + cy_) * 0.5;    
    double mcbx = (bx_ + cx_) * 0.5;
    double mcby = (by_ + cy_) * 0.5;    
    double r 	= cax * cby - cbx * cay;    
    *px_ = ( ( cax * mcax + cay * mcay ) * cby - (cbx * mcbx + cby * mcby) * cay ) / r;      
    *py_ = ( ( cbx * mcbx + cby * mcby ) * cax - (cax * mcax + cay * mcay) * cbx ) / r;      
  };

  
  class edge_select_boundary
  {    
  public: edge_select_boundary() = delete;
  public: static inline bool select(mesh_t*m,int_t edge) noexcept
    {
      return m->get_edge_ptr(edge)->is_boundary();
    };
  };

  class edge_select_boundary_saturation
  {    
  public: edge_select_boundary_saturation() = delete;
  public: static inline bool select(mesh_t*m,int_t edge) noexcept
    {
      if (m->get_edge_ptr(edge)->is_boundary())
	{
	  int_t n = m->edge_saturation_size(edge);		
	  return (n>1);
	}
      return false;
    };
  };

  class edge_select_interior
  {
  public: edge_select_interior() = delete;
  public: static inline bool select(mesh_t*m,int_t edge) noexcept
    {
      if (false == m->get_edge_ptr(edge)->is_boundary())	
	{
	  int_t n = m->edge_saturation_size(edge);		
	  return (n>1);
	}
      return false;
    };
  };
  
  template<typename S>
  class edge_split
  {
  public: edge_split() = delete;    
  public: inline static void apply(mesh_t*m,int_t edge) noexcept
    {      
      if (m->optimize_interior(edge))
	{
	  for (int_t i=0;i<m->m_cavityOperator.ninterior_edges();++i)
	    {
	      auto e = m->m_cavityOperator.interior_edge(i);
	      m->classify_edge(e);
	    }
#if 0
	  for (auto e : m->m_cavityOperator.m_listCavityInteriorEdges)
	    {
	      m->classify_edge(e);
	    }
#endif
	}      
    };
    
  public: inline static bool criterion(mesh_t*m,int_t edge) noexcept
    {
      if (S::select(m,edge))
	{
	  int_t n = m->edge_saturation_size(edge);		
	  return (n>1);
	}
      return false;
    };
    
  public: inline static bool cavity_criterion(mesh_t*m,int_t edge) noexcept
    {
      return m->get_edge_ptr(edge)->is_boundary();
    };
  };

  
  template<typename S>
  class edge_split_boundary
  {
  public: edge_split_boundary() = delete;    
  public: inline static void apply(mesh_t*m,int_t edge) noexcept
    {      
      m->optimize_boundary(edge);
	  for (int_t i=0;i<m->m_cavityOperator.ninterior_edges();++i)
	    {
	      auto e = m->m_cavityOperator.interior_edge(i);
	      m->classify_edge(e);
	    }
#if 0
      for (auto e : m->m_cavityOperator.m_listCavityInteriorEdges)
	{
	  m->classify_edge(e);
	}
#endif
    };
  public: inline static bool criterion(mesh_t*m,int_t edge) noexcept
    {
      if (S::select(m,edge))
	{
	  int_t n = m->edge_saturation_size(edge);		
	  return (n>1);
	}
      return false;
    };
  public: inline static bool cavity_criterion(mesh_t*m,int_t edge) noexcept
    {
      return S::cavity_criterion(m,edge);
    };
  };
  
  template<typename S>
  class edge_split_boundary_to_interior
  {
  public: edge_split_boundary_to_interior() = delete;    
  public: inline static void apply(mesh_t*m,int_t edge) noexcept
    {
      m->optimize_boundary_to_interior(edge);
#if 0
      for (auto e : m->m_cavityOperator.m_listCavityInteriorEdges)
	{
	  m->classify_edge(e);
	}
#endif
      	  for (int_t i=0;i<m->m_cavityOperator.ninterior_edges();++i)
	    {
	      auto e = m->m_cavityOperator.interior_edge(i);
	      m->classify_edge(e);
	    }
#if 0
      	  for (int_t i=0;i<cavityOperator_.ninterior_edges();++i)
	    {
	      auto e = cavityOperator_.interior_edge(i);
	      m->classify_edge(e);
	    }
#endif
    };
    
  public: inline static bool criterion(mesh_t*m,int_t edge) noexcept
    {
      if (S::select(m,edge))
	{
	  return true;
	}
      return false;
    };
    
  public: inline static bool cavity_criterion(mesh_t*m,int_t edge) noexcept
    {
      return S::cavity_criterion(m,edge);
    };
  };

  
#if 0
  double xs[2048*4];
  double ys[2048*4];
  double ts[2048*4];
  int_t n = this->edge_saturation_size(edge);
  int_t triangle = 1;
  if (n>1)
    {
      double neg = this->edge_saturation(edge,n,xs,ys,ts);
      double m0 = this->xy[this->get_edge_ptr(edge)->start()].m_metric[0];
      double m1 = this->xy[this->get_edge_ptr(edge)->stop()].m_metric[0];	
      int_t edge_backup  = edge;	
      int_t edgep0 = this->get_edge_ptr(edge)->start(); 
      int_t edgep1 = this->get_edge_ptr(edge)->stop();
      int_t new_edge[2];
      //	std::cout << "n " << n << std::endl;
      for (int i=0;i<n;++i)
	{
	  //	    std::cout << xs[i] << " " << ys[i] << std::endl;
	  int_t node = this->create_node(xs[i],ys[i]);
	  this->xy[node].x = xs[i];
	  this->xy[node].y = ys[i];			
	  auto mx = (1.0 - ts[i]) * 1.0/sqrt(m0) + ts[i] * 1.0/sqrt(m1);
	  this->xy[node].m_metric[0]  = 1.0 / (mx*mx);
	  this->xy[node].m_metric[1]  = 0.0;
	  this->xy[node].m_metric[2]  = 1.0 / (mx*mx);
	  triangle     = find_triangle(node,triangle);
	  if (triangle)
	    {
	      this->InsertPoint(triangle,node);
	    }
	}
      //  ....
    }
#endif  


  


  classifier m_classifier;


  double size_edge(const edge_t * edge_ptr)
  {
    //
    // Evaluate the size.
    //
    int_t node_start = edge_ptr->start();	
    int_t node_end = edge_ptr->stop();
    const metric_t * metric_start = (const metric_t *)&this->xy[node_start].m_metric[0];
    const metric_t * metric_end   = (const metric_t *)&this->xy[node_end].m_metric[0];			
    double s0 = metric_start->size(this->xy[node_start].x,this->xy[node_start].y,this->xy[node_end].x,this->xy[node_end].y);
    double s1 = metric_end->size(this->xy[node_start].x,this->xy[node_start].y,this->xy[node_end].x,this->xy[node_end].y);    
    return s0 >= s1 ? s0 : s1;
  };
  
  bool unacceptable_size(int_t edge)
  {
    //    std::cout << "size" << size_edge(edge) << std::endl;
    //    return size_edge(edge) > 1.0;
    int_t n = this->edge_saturation_size(edge);		
    return (n>1);
  }
  
  void classify_edge(int_t edge)
  {
#if 1
    this->m_classifier.push(edge,			    
			    [this](int_t edge)
			    {
			      return unacceptable_size(edge);			
			    });
#endif
  };

  void declassify_edge(int_t edge)
  {
#if 1
    this->m_classifier.pop(edge);
#endif
  };
  
  template<typename op>
  inline void while_anyold()
  {
#if 0
    bool tt = true;
    while (tt)
      {
	tt = false;
	std::stack<int_t> stack;
	for (int_t edge=1;edge<=this->m_nedges;++edge)
	  {
	    if (op::criterion(this, edge))
	      {		
		stack.push(edge);		
	      }
	  }
	
	tt = !stack.empty();    
	while (false == stack.empty())
	  {
	    int_t edge = stack.top();
	    op::apply(this, edge);
	    stack.pop();
	  }
      }
#endif


#if 0
    
    m_classifier.clear();
    for (int_t edge=1;edge<=this->m_nedges;++edge)
      {
	classify_edge(edge);
      }
    
    std::cout << "presize classifier " << m_classifier.size() << std::endl;
    while (m_classifier.size() > 0)
      {
	std::cout << "####### presize classifier " << m_classifier.size() << std::endl;
	int_t edge = 0;
	edge = *m_classifier.begin();
	m_classifier.remove(edge);
	
	op::apply(this, edge);	  
      }
    
    //    exit(1);
#if 0
    //
    // Initialize the set of too long edges.
    //
    for (int_t edge=1;edge<=this->m_nedges;++edge)
      {
	if (op::criterion(this, edge))
	  {		
	    m_setOfTooLongInteriorEdges.add(edge);
	  }
      }

    std::cout << "size " << m_setOfTooLongInteriorEdges.size() << " " << m_nedges << std::endl;
    
    while (false == m_setOfTooLongInteriorEdges.empty())
      {
	int_t edge = 0;

	edge = *m_setOfTooLongInteriorEdges.begin();
	m_setOfTooLongInteriorEdges.remove(edge);
	
	op::apply(this, edge);	  
      }
#endif    
#endif

#if 0
    m_classifier.clear();
    
    for (int_t edge=1;edge<=this->m_nedges;++edge)
      {

	if (op::criterion(this, edge))
	  {
	    classify_edge(edge);
	  }
      }
    
    for (int i=0;i<30;++i)
      {


	// td::cout << "presize classifier " << m_classifier.size() << std::endl;
	while (m_classifier.size() > 0)
	  {

	    int_t edge = 0;
	    edge = *m_classifier.begin();
	    m_classifier.remove(edge);
	    if (op::criterion(this, edge))
	      {	
		op::apply(this, edge);
	      }
	  }

    
#if 1
	for (int_t edge=1;edge<=this->m_nedges;++edge)
	  {
	    if (op::criterion(this, edge))
	      {
		classify_edge(edge);
	      }
	  }
#endif
	if (m_classifier.size()==0)
	  {
	    break;
	  }
#if 0
	if (m_classifier.size()>0)
	  {
	    std::cout << "#######2 presize classifier " << m_classifier.size() << std::endl;
	    exit(1);
	  }
#endif

      }
#endif
    
    //    std::cout << "#######final presize classifier " << m_classifier.size() << std::endl;    



      



#if 1
    bool tt = true;
    while (tt)
      {
	tt = false;
	for (int_t edge=1;edge<=this->m_nedges;++edge)
	  {
	
	    if (op::criterion(this, edge))
	      {
		//		std::cout << "edge " << edge << " " << get_edge_ptr(edge)->is_boundary() << std::endl;
		//		std::cout << "edge " << edge << " " << get_edge_ptr(edge)->start() << " " << get_edge_ptr(edge)->stop() << std::endl;
		//		std::cout << ip  << std::endl;
		op::apply(this, edge);
		tt = true;
	      }
	    // std::cout << "#######final presize classifier " << m_classifier.size() << std::endl;    
	    //    else
	    //	      {
	    //	++edge;
	    //	      }
	  }
      }
#endif

    //    std::cout << "#######final presize classifier " << m_classifier.size() << std::endl;    
  };





  template<typename op>
  inline void while_any()
  {
    m_classifier.clear();    
    for (int_t edge=1;edge<=this->m_nedges;++edge)
      {
	if (op::criterion(this, edge))
	  {
	    classify_edge(edge);
	  }
      }
    
    for (int i=0;i<30;++i)
      {
	//	std::cout << "presize classifier " << m_classifier.size() << std::endl;
	while (m_classifier.size() > 0)
	  {

	    int_t edge = 0;
	    edge = *m_classifier.begin();
	    m_classifier.remove(edge);
	    if (op::criterion(this, edge))
	      {	
		op::apply(this, edge);
	      }
	  }
#if 0
	for (int_t edge=1;edge<=this->m_nedges;++edge)
	  {
	    if (op::criterion(this, edge))
	      {
		classify_edge(edge);
	      }
	  }
#endif
	if (m_classifier.size()==0)
	  {
	    break;
	  }
      }

    //    std::cout << "#######final presize classifier " << m_classifier.size() << std::endl;    

  };

  static void vec_mmm(int_t n,double * x)
  {
    double mn = 1.0e+30,mx=0.0, mean=0.0, var;
    for (int_t i=1;i<=n;++i)
      {
	mean += x[i];
	mn = std::min(mn,x[i]);
	mx = std::max(mx,x[i]);	
      }
    mean /= double(n);
    var = 0.0;
    for (int_t i=1;i<=n;++i)
      {
	var += (x[i]-mean)*(x[i]-mean);
      }
    var /= double(n);
    var /= double(n);
    std::cout << " min  " << mn << std::endl;
    std::cout << " max  " << mx << std::endl;
    std::cout << " mean " << mean << std::endl;
    std::cout << " var  " << var << std::endl;
  }
  
  inline void analysis() const
  {
    //
    // 
    //

    {

      double * lengths = new double[m_nedges+1];
      for (int_t edge=1;edge<=m_nedges;++edge)
	{
	  auto p = get_edge_ptr(edge)->start();
	  auto q = get_edge_ptr(edge)->stop();
	  double hx = xy[p].x - xy[q].x;
	  double hy = xy[p].y - xy[q].y;
	  lengths[edge] = sqrt( hx*hx+hy*hy );
      }
      std::cout << "EDGE LENGTHS" << std::endl;
      vec_mmm(m_nedges,lengths);
      delete[] lengths;
    }

    {


      double * areas = new double[m_ntriangles+1];
      int_t cat[3] = {0,0,0};
      for (int_t face=1;face<=m_ntriangles;++face)
	{
	  int_t nodes[3];
	  triangle_to_nodes(face,nodes);
	  cat[ComputeOrientation(xy[nodes[0]].x,
				 xy[nodes[0]].y,
				 xy[nodes[1]].x,
				 xy[nodes[1]].y,
				 xy[nodes[2]].x,
				 xy[nodes[2]].y)]+=1;
	  
	  areas[face] = ComputeSignedArea(xy[nodes[0]].x,
					  xy[nodes[0]].y,
					  xy[nodes[1]].x,
					  xy[nodes[1]].y,
					  xy[nodes[2]].x,
					  xy[nodes[2]].y);
	}
      
      std::cout << "BACKWARD " << cat[0] << " / " << m_ntriangles << std::endl;
      std::cout << "FORWARD  " << cat[1] << " / " << m_ntriangles << std::endl;
      std::cout << "NEUTRAL  " << cat[2] << " / " << m_ntriangles << std::endl;
      vec_mmm(m_ntriangles,areas);
      delete[] areas;
      
    }


    
  };
  
  inline void optimize_boundary(int_t edge)
  {
    double xs[2048*4];
    double ys[2048*4];
    double ts[2048*4];
    int_t n = this->edge_saturation_size(edge);
    if (n>1)
      {
	double neg = this->edge_saturation(edge,n,xs,ys,ts);

	int_t node_start = this->get_edge_ptr(edge)->start();	
	int_t node_end = this->get_edge_ptr(edge)->stop();
	
	
	double m0 = this->xy[this->get_edge_ptr(edge)->start()].m_metric[0];
	double m1 = this->xy[this->get_edge_ptr(edge)->stop()].m_metric[0];	
	int_t edge_backup  = edge;	
	int_t edgep0 = this->get_edge_ptr(edge)->start(); 
	int_t edgep1 = this->get_edge_ptr(edge)->stop();
	int_t new_edge[2];
	for (int i=0;i<n;++i)
	  {
	    int_t node = this->create_node(xs[i],ys[i]);
	    this->xy[node].x = xs[i];
	    this->xy[node].y = ys[i];
	    const metric_t * metric_start = (const metric_t *)&this->xy[node_start].m_metric[0];
	    const metric_t * metric_end   = (const metric_t *)&this->xy[node_end].m_metric[0];
	    

	    metric_t * metric   = (metric_t *)&this->xy[node].m_metric[0];
	    //	    std::cout << "interpolate a" << std::endl;
	    metric_t::interpolate(ts[i],
				  metric_start,
				  metric_end,
				  metric);
	    
	    int_t triangle     = m_edges[edge_backup].triangle(0);
	    {
	      new_edge[0]=0;
	      new_edge[1]=0;
#if 1
	      if (!m_edges[edge_backup].triangle(1))
		{		  
		  this->InsertPointBoundary(node,
					    edge_backup,
					    neg < 0,
					    new_edge);
	      
		  if (m_edges[new_edge[0]].start() == edgep1 || m_edges[new_edge[0]].stop() == edgep1)
		    {
		      int_t tmp = new_edge[0];
		      new_edge[0] = new_edge[1];
		      new_edge[1] = tmp;
		    }

		  if (neg>0)
		    {
		      edge_backup = new_edge[1];
		    
		      //		      edgep1 = node;
		    }
		  else
		    {
		      edge_backup = new_edge[0];		      
		    }
		  
		  //		  else
		  {
		    //		      edge_backup = new_edge[0];
		    //		      edgep0 = node;
		  }
		}
	      else
#endif
		{
		  this->InsertPointWeakBoundary(node,
						edge_backup,
						neg < 0,
						new_edge);

		  if (m_edges[new_edge[0]].start() == edgep1 || m_edges[new_edge[0]].stop() == edgep1)
		    {
		      int_t tmp = new_edge[0];
		      new_edge[0] = new_edge[1];
		      new_edge[1] = tmp;
		    }
	      
		  if (neg<0)
		    {
		      edge_backup = new_edge[0];
			      
		      edgep1 = node;
		    }
		  else
		    {
		      edge_backup = new_edge[1];
		      edgep0 = node;
		    }
		}

			  
	    }
	  }
	//  ....
      }
  };


  


  inline void edge_subdivide(const edge_t* 		edge_,
			     double * x,
			     double * y)
  {      
    static constexpr double s3 		= double(3.0);
    static constexpr double s1_8	= double(0.125);
    
#if 0
    const auto 	face0 			= edge_->triangle(0);
    const auto	face1 			= edge_->triangle(1);      
    const auto 	p0 			= get_triangle_ptr(face0)->GetOppositePoint(edge_);
    const auto 	p1 			= get_triangle_ptr(face1)->GetOppositePoint(edge_);
#endif
    const auto 	es 			= edge_->start();
    const auto 	ee 			= edge_->stop();
    
    //    x[0] =  (s3 * ( xy[es].x + xy[ee].x) + xy[p0].x + xy[p1].x)*s1_8;
    //    y[0] =  (s3 * ( xy[es].y + xy[ee].y) + xy[p0].y + xy[p1].y)*s1_8;
    x[0] =  (xy[es].x + xy[ee].x) * 0.5;
    y[0] =  (xy[es].y + xy[ee].y) * 0.5;
  };

  
  inline bool optimize_interior(int_t edge)
  {
    bool has_rejected = false;
    double xs[2048*4];
    double ys[2048*4];
    double ts[2048*4];
    
    int_t n = this->edge_saturation_size(edge);
    if (n >= 2048*4)
      {
	std::cerr << "invalid saturation size" << std::endl;
	exit(1);
      }
    
    int_t triangle = get_edge_ptr(edge)->triangle(0);

    //    std::cout << "optmize interior n " << n << std::endl;
    if (n>1)
      {
	double neg = this->edge_saturation(edge,n,xs,ys,ts);
	int_t new_edge[2];
	int_t node_start = this->get_edge_ptr(edge)->start();	
	int_t node_end = this->get_edge_ptr(edge)->stop();
	//	std::cout << "start " << xy[node_start].x << " " << xy[node_start].y << std::endl;
	//	std::cout << "end   " << xy[node_end].x << " " << xy[node_end].y << std::endl;
	//	std::cout << "n " << n << std::endl;
	for (int i=0;i<n;++i)
	  {
	    int_t node = this->create_node(xs[i],ys[i]);
	    //  std::cout << xs[i] << " " << ys[i] << " " << ts[i] << std::endl;

	    //
	    // pointers might have moved.
	    //
	    const metric_t * metric_start = (const metric_t *)&this->xy[node_start].m_metric[0];
	    const metric_t * metric_end   = (const metric_t *)&this->xy[node_end].m_metric[0];
	    
	    this->xy[node].x = xs[i];
	    this->xy[node].y = ys[i];
	    
	    
	    metric_t * metric   = (metric_t *)&this->xy[node].m_metric[0];
	    //	    std::cout << "interpolate b" << node << " " << m_nnodes << " " << " " << mxnnodes << " " << this->xy[node].x << " " << this->xy[node].y << std::endl;
	    
	    metric_t::interpolate(ts[i],
				  metric_start,
				  metric_end,
				  metric);
	    //   std::cout << "interpolate done" << std::endl;
	    
	    //	    std::cout <<"ch "  << triangle << std::endl;
	    triangle = find_triangle(node,triangle);
	    if (triangle)
	      {
		if (!this->InsertPoint(triangle,node,true))
		  {
		    //		    std::cout << "rejected " << node << " " << m_nnodes << std::endl;
		    --m_nnodes;
		    has_rejected = true;
		  }
	      }
	    else
	      {
		std::cout<< "triangle not found " << xy[node].x << " " << xy[node].y << " " <<triangle<< std::endl;
		    --m_nnodes;
		    has_rejected = true;
		    //		exit(1);
	      }
	  }
	//  ....
      }
    return !has_rejected;
  };


  
  inline bool optimize_boundary_to_interior(int_t edge)
  {
    
    bool has_rejected = false;

    int_t 	triangle 	= get_edge_ptr(edge)->triangle(0);
    int 	lindex		= get_triangle_ptr(triangle)->get_local_edge_index(edge);
    int_t 	node_start 	= get_edge_ptr(edge)->start(get_triangle_ptr(triangle)->way(lindex));
    int_t 	node_stop 	= get_edge_ptr(edge)->stop(get_triangle_ptr(triangle)->way(lindex));

    double x1 			= get_node_ptr(node_start)->x;
    double y1 			= get_node_ptr(node_start)->y;
    double x2 			= get_node_ptr(node_stop)->x;
    double y2 			= get_node_ptr(node_stop)->y;
    std::cout << " " << x1 << " " << y1 << std::endl;
    std::cout << " " << x2 << " " << y2 << std::endl;

    
    double hx  = x2-x1;
    double hy  = y2-y1;
    double tx  = -hy;
    double ty  = hx;
    //    double l   = sqrt(hx*hx*this->xy[node_start].m_metric[0] + hy*hy*this->xy[node_start].m_metric[2] +
    //		      2.0*hx*hy * this->xy[node_start].m_metric[1]);
    double l   = sqrt(hx*hx+ hy*hy);    
    int_t node = this->create_node(x2 + tx,y2 + ty);
    
    std::cout << " " << x2 + tx << " " << y2 + ty << std::endl;
    
    this->xy[node].m_metric[0]=this->xy[node_start].m_metric[0];
    this->xy[node].m_metric[1]=this->xy[node_start].m_metric[1];
    this->xy[node].m_metric[2]=this->xy[node_start].m_metric[2];

    triangle = find_triangle(node,triangle);
    if (triangle)
      {
	if (!this->InsertPoint(triangle,node,true))
	  {
	    //		    std::cout << "rejected " << node << " " << m_nnodes << std::endl;
	    --m_nnodes;
	    has_rejected = true;
	  }
      }
    else
      {
	std::cout<< "triangle not found " << xy[node].x << " " << xy[node].y << " " <<triangle<< std::endl;
	--m_nnodes;
	has_rejected = true;
	//	exit(1);
      }
    return !has_rejected;
  };


  inline void node_optimize_position_barycentric(int_t node,int n,int_t nodes[])
  {
    double x = 0.0;
    double y = 0.0;
    //    std::cout << "ball["<<node<<"] = {" << std::endl;
    for (int i=0;i<n;++i)
      {
	//	std::cout << " " << nodes[i];
	x+=get_node_ptr(nodes[i])->x;
	y+=get_node_ptr(nodes[i])->y;
      }
    //    std::cout << "}" <<  std::endl;
    x/=double(n);
    y/=double(n);
    //    exit(1);
    
    get_node_ptr(node)->x = x;
    get_node_ptr(node)->y = y;
    //    std::cout <<     "eeeeee " << std::endl << get_node_ptr(node)->x << " " << x << std::endl;
    //    std::cout <<     get_node_ptr(node)->y << " " << y << std::endl;
  };


  inline void node_optimize_position(int_t node,int n,int_t nodes[])
  {
    double x = get_node_ptr(node)->x;
    double y = get_node_ptr(node)->y;
    double xnew = x;
    double ynew = y;
    //    std::cout << "ball["<<node<<"] = {" << std::endl;
    for (int i=0;i<n;++i)
      {
	//	std::cout << " " << nodes[i];
	double xp = get_node_ptr(nodes[i])->x;
	double yp = get_node_ptr(nodes[i])->y;
	xp = x-xp;
	yp = y-yp;
	double d = sqrt( xp*xp + yp*yp);
	xp /= d;
	yp /= d;
	double f = -d;//0.2* (1.0-d*d*d*d) * exp(-d*d*d*d);
	xnew += (xp * f)/double(n);
	ynew += (yp * f)/double(n);
      }
    
    //    std::cout << "}" <<  std::endl;
    //    x/=double(n);
    //    y/=double(n);
    //    exit(1);
    
    get_node_ptr(node)->x = xnew;
    get_node_ptr(node)->y = ynew;
    //    std::cout <<     "eeeeee " << std::endl << get_node_ptr(node)->x << " " << x << std::endl;
    //    std::cout <<     get_node_ptr(node)->y << " " << y << std::endl;
  };

  
  inline int node_to_nodes(int_t edge_,int_t nodes_[],int_t*marker,int_t*node)
  {    
    auto edge 		= edge_;
    auto edge_ptr 	= get_edge_ptr(edge);
    auto face 		= edge_ptr->triangle(0);
    auto face_ptr 	= get_triangle_ptr(face);
    int lindex		= face_ptr->get_local_edge_index(edge);
    int_t node0 	= edge_ptr->start(face_ptr->way(lindex));
    int_t node1 	= edge_ptr->stop(face_ptr->way(lindex));
    node[0] = node0;
    if (marker[node0])
      {
	return 0;
      }
    marker[node0] = 1;
    
    int_t nodes[3];
    triangle_to_nodes(face,nodes);    
    int_t n = 0;
    //    nodes_[n++] = node1;

    //
    // Connected to opposite node
    //
    do
      {
	int index_node1=0;
	for (int i=1;i<3;++i)
	  {
	    if (nodes[i] == node1)
	      {
		index_node1 = i;
		break;
	      }
	  }

	edge = face_ptr->edge(index_node1);
	if (edge != edge_)
	  {
	    edge_ptr 	= get_edge_ptr(edge);
	    face 	= edge_ptr->otriangle(face);
	    if (face)
	      {
		triangle_to_nodes(face,nodes);	
		face_ptr    	= get_triangle_ptr(face);		
		lindex  	= face_ptr->get_local_edge_index(edge);
		nodes_[n++] 	= nodes[lindex];		
		node1 		= edge_ptr->stop(face_ptr->way(lindex));	    
	      }
	    else
	      {
		return 0;
	      }
	  }
      }
    while (edge != edge_);
    
    std::cout << "[" << node0 << "]";
    for (int i=0;i<n;++i)
      {
    	std::cout << " " << nodes_[i];
      }
    std::cout << std::endl;
    
    return n;
  };



  inline void optimize_node_positionsold()
  {
    int_t * marker = (int_t*)calloc(m_nnodes+1,sizeof(int_t));
    int_t nodes[64];
    for (int i=0;i<30;++i)
      {
	for (int j=0;j<=m_nnodes;++j)
	  marker[j]=0;
	
	for (int_t edge = 1;edge <= m_nedges;++edge)
	  {
	    int_t node;
	    int n = node_to_nodes(edge,nodes,marker,&node);
	    if (n>0)
	      {
		node_optimize_position(node,n,nodes);
		print_mesh();
	      }
	  }
	for (int j=0;j<=m_nnodes;++j)
	  {
	    if (marker[j]==0)
	      {
		std::cout << " " << j << std::endl;
	      }
	  }

      }
    free(marker);
  };

  
  inline void optimize_node_positions()
  {
    
    int_t nodes[64];
    int_t n = 0;
    for (int_t i=0;i<100;++i)
      {
    for (int_t node = 1;node <= m_nnodes;++node)
      {
	
	n = 0;
	for (int_t edge = 1;edge <= m_nedges;++edge)
	  {
	    auto p = get_edge_ptr(edge);
	    if (p->start() == node || p->stop() == node)
	      {
		if (p->triangle(1) && false == p->is_boundary())
		  {
		    nodes[n++] = p->start() == node ? p->stop() : p->start();
		  }
		else
		  {
		    n = 0;
		    break;
		  }
	      }
	  }
	if (n>0)
	  {
	    node_optimize_position(node,n,nodes);
	  }
      }
      }
#if 0    
    int_t * marker = (int_t*)calloc(m_nnodes+1,sizeof(int_t));
    int_t nodes[64];
    for (int i=0;i<30;++i)
      {
	for (int j=0;j<=m_nnodes;++j)
	  marker[j]=0;
	
	for (int_t edge = 1;edge <= m_nedges;++edge)
	  {
	    int_t node;
	    int n = node_to_nodes(edge,nodes,marker,&node);
	    if (n>0)
	      {
		node_optimize_position(node,n,nodes);
		print_mesh();
	      }
	  }
	for (int j=0;j<=m_nnodes;++j)
	  {
	    if (marker[j]==0)
	      {
		std::cout << " " << j << std::endl;
	      }
	  }

      }
    free(marker);
#endif
  };
  
  inline bool InsertPoint(int_t triangle,
			  int_t point_,
			  bool filter=false)
  {    
    bool success = m_cavityOperator.Compute(triangle,
					    point_,
					    [this] (int_t face,
						    int_t point)
					    {
					      int_t nodes[3];
					      triangle_to_nodes(face,nodes);
					      return (Orientation::BACKWARD != 
						      ComputeCircleOrientation(xy[nodes[0]].x,
									       xy[nodes[0]].y,
									       xy[nodes[1]].x,
									       xy[nodes[1]].y,
									       xy[nodes[2]].x,
									       xy[nodes[2]].y,
									       xy[point].x,
									       xy[point].y));
					    },
					      
					    [this](int_t edge)
					    {
					      // std::cout << " cccc " << edge << std::endl;
					      return m_edges[edge].is_boundary();
					      // return edge->is_boundary();
					    },
					    					    					    
					    [this](int_t face)
					    {
					      return get_triangle_ptr(face);
					    },
					    
					    [this](int_t edge)
					    {
					      return get_edge_ptr(edge);
					    });
      
    if (success)
      {

	//
	// We have:
	// - the list of the faces to remove
	// - the list of the interior edges
	// - the list of the boundary edges
	//

	//
	// Transfer the ids of the triangles and interior edges to delete.
	// They won't be deleted but replaced with other triangles.
	//

	//
	//
	//
	if (filter)
	  {
	    // m_listCavityBoundaryEdges.size()
	    if (m_cavityOperator.nboundary_edges() < 2)
	      {
		return false;	  
	      }
	    //	    std::cout << "filter" << std::endl;
	    for (int_t i=0;i<m_cavityOperator.nboundary_edges();++i)
	    //	    for(auto K : m_cavityOperator.m_listCavityBoundaryEdges)
	      {
		auto K = m_cavityOperator.boundary_edge(i);
		auto p = get_edge_ptr(K.first)->start(K.second);
		auto q = get_edge_ptr(K.first)->stop(K.second);
		const metric_t * metric = (const metric_t *)&this->xy[point_].m_metric[0];
		double lq = metric->size(xy[point_].x,xy[point_].y,xy[p].x,xy[p].y);
		double lp = metric->size(xy[point_].x,xy[point_].y,xy[q].x,xy[q].y);
		if (lq < 1.0 || lp  < 1.0)
		  {
		    //	    std::cout << lq << " xxxxxxxxxxxxxxxxxxxxxxxxxxxxx "  << lp << std::endl;
		    return false;
		  }
	      }
	  }
	
	*this -= m_cavityOperator;

	// 
	// 
	//
	Rebuild(point_,m_cavityOperator);
	this->check(__LINE__);
      
	//
	// The cavity operator contains the new ids of triangles and interior edges and the list of boundary edges.
	//	  
	print_mesh();
      }
    else
      {
	std::cout << "unable to compute the cavity" << std::endl;
	return false;
      }
    
    return true;
  };


  
public: inline void print_mesh(const char * filename_)
  {
    if (m_debug_mode)
      {
	this->print(filename_,ip++);
      }
  }

public: inline void print_mesh()
  {
    print_mesh(this->debug_filebasename());
  }

private:inline void print_cnc()
  {
    fprintf(stdout,"\nINFO  nfaces: %lld m_nedges: %lld\n",m_ntriangles,m_nedges);
    int_t nodeids[3];
    for (int_t i=1;i<=m_ntriangles;++i)
      {
	triangle_to_nodes(i,nodeids);
	fprintf(stdout,"%lld %lld %lld %lld\n",
		nodeids[0],
		nodeids[1],
		nodeids[2],
		i);
      }
      
    fprintf(stdout,"\nedges->nodes\n");
    for (int_t i=1;i<=m_nedges;++i)
      {
	fprintf(stdout,"edge.node[%lld]  = %lld %lld \n",
		i,m_edges[i].start(),m_edges[i].stop());
      }
      
    fprintf(stdout,"\ntriangle->edges\n");
    for (int_t i=1;i<=m_ntriangles;++i)
      {
	fprintf(stdout,"triangle.edge[%lld] =  %lld %lld %lld\n",
		i,
		triangles[i].edge(0),
		triangles[i].edge(1),
		triangles[i].edge(2));
      }

    fprintf(stdout,"\nedges->triangles\n");
    for (int_t i=1;i<=m_nedges;++i)
      {
	fprintf(stdout,"edge.triangle[%lld] =  %lld %lld\n",
		i,m_edges[i].triangle(0),m_edges[i].triangle(1));
      }

  }
  
public: inline void print(const char * filename_,int i)
  {

    {
      char ctmp[256];
      sprintf(ctmp,"%s.%d.mesh",filename_,i);
      FILE * f = fopen(ctmp,"w");
      fprintf(f,"MeshVersionFormatted\n1\nDimension\n2\nVertices\n%lld\n",m_nnodes);
      for (int_t i=1;i<=m_nnodes;++i)
	{
	  fprintf(f,"%e %e %lld\n",xy[i].x,xy[i].y,i);
	}
      fprintf(f,"Triangles\n%lld\n",m_ntriangles);
      int_t nodeids[3];
      int_t edges[3];
      for (int_t i=1;i<=m_ntriangles;++i)
	{
	  triangle_to_nodes(i,nodeids);
	  fprintf(f,"%lld %lld %lld %lld\n",
		  nodeids[0],
		  nodeids[1],
		  nodeids[2],
		  triangles[i].id());
	}
      fprintf(f,"End");
      fclose(f);

    }

    {
      char ctmp[256];
      sprintf(ctmp,"%s.%d.bb",filename_,i);
      FILE * f = fopen(ctmp,"w");
      fprintf(f,"2 3 %lld 2\n",m_nnodes);
      for (int_t node=1;node<=m_nnodes;++node)
	{
	  fprintf(f,"%e %e %e\n",xy[node].m_metric[0],xy[node].m_metric[1],xy[node].m_metric[2]);
	}
      fclose(f);
    }

      
  };
    

public: inline void CreateInitialMesh(xy_t low,
				      xy_t up)
  {

    int_t p1 = CreateNode(low.x,low.y,0);
    int_t p2 = CreateNode(up.x,low.y,0);
    int_t p3 = CreateNode(up.x,up.y,0);
    int_t p4 = CreateNode(low.x,up.y,0);

    int_t a1 = CreateEdge(p1,p2);
    int_t a2 = CreateEdge(p2,p3);
    int_t a3 = CreateEdge(p3,p4);
    int_t a4 = CreateEdge(p4,p1);
    int_t a5 = CreateEdge(p1,p3);

    int_t f1 = CreateTriangle(a1,1,
			      a2,2,
			      a5,-3,__LINE__);

    int_t f2 = CreateTriangle(a5,1,
			      a3,2,
			      a4,3,__LINE__);
  };

public: inline int_t nnodes() 		const { return this->m_nnodes;};
public: inline int_t nedges() 		const { return this->m_nedges;};
public: inline int_t ntriangles() 	const { return this->m_ntriangles;};
  
public: inline const xy_t * 		get_node_ptr	(int_t node) const 	noexcept { return &this->xy[node]; };
public: inline const edge_t * 		get_edge_ptr	(int_t edge) const 	noexcept { return &this->m_edges[edge]; };
public: inline const triangle_t * 	get_triangle_ptr(int_t triangle) const 	noexcept { return &this->triangles[triangle]; };

public: inline xy_t * 			get_node_ptr	(int_t node) 		noexcept { return &this->xy[node]; };
public: inline edge_t * 		get_edge_ptr	(int_t edge) 		noexcept { return &this->m_edges[edge]; };
public: inline triangle_t * 		get_triangle_ptr(int_t triangle) 	noexcept { return &this->triangles[triangle]; };
  
public: inline int_t create_node(double x,
				 double y,
				 int id = 0)
  {      
    if (++m_nnodes > mxnnodes)
      {
	mxnnodes = mxnnodes * 3;
	void * tmp = realloc(xy,sizeof(xy_t)*mxnnodes);
	xy = (xy_t*)tmp;
      }      
    xy[m_nnodes] = { x, y, id };      
    return m_nnodes;
  }
  
public: inline int_t CreateNode(double x,
				double y,
				int id = 0)
  { 
    if (++m_nnodes > mxnnodes)
      {
	mxnnodes = mxnnodes * 3;
	void * tmp = realloc(xy,sizeof(xy_t)*mxnnodes);
	xy = (xy_t*)tmp;
      }
    
    xy[m_nnodes] = { x, y, id };
    return m_nnodes;
  }
  
  inline int_t CreateEdge(int_t i,
			  int_t j,
			  int id = 0)
  {
    if (++m_nedges > mxnedges)
      {
	mxnedges = mxnedges * 3;
	void * tmp = realloc(m_edges,sizeof(edge_t)*mxnedges);
	m_edges = (edge_t*)tmp;
      }
      
    m_edges[m_nedges] = edge_t( i, j, id );
    return m_nedges;
  }

  
  inline int_t CreateTriangle(int_t e0, int w0,int_t e1,int w1,int_t e2,int w2,int line_)
  {
    int_t edges[3] = {e0,e1,e2};
    int ways[3] = {w0,w1,w2};
    
    if (++m_ntriangles > mxntriangles)
      {
	mxntriangles = mxntriangles * 3;
	void * tmp = realloc(triangles,sizeof(triangle_t)*mxntriangles);
	triangles = (triangle_t*)tmp;
      }
    triangles[m_ntriangles] = triangle_t(m_ntriangles, edges ,ways );

    for (int i=0;i<3;++i)
      {
	//	std::cout << "   insert nei edge[" << i << "] " << edges[i] << std::endl;
	if (m_edges[edges[i]].triangle(0) == 0)
	  {
	    //	    std::cout << " A" << std::endl;
	    m_edges[edges[i]].triangle(0,m_ntriangles);
	  }
	else if (m_edges[edges[i]].triangle(1) == 0)
	  {
	    //	    std::cout << " B" << std::endl;
	    m_edges[edges[i]].triangle(1,m_ntriangles);
	  }
	else
	  {
	    std::cerr << "edge corrupted "<< edges[i] << std::endl;
	    exit(1);
	  }
	
      }

    return m_ntriangles;
  }

protected:
  inline void replace_edge(int_t edge,
			  int_t i,
			  int_t j,
			  int   id = 0)
  {
    int_t nodes[2] = {i,j};
    int_t tr[2] = {0,0};
    m_edges[edge] = edge_t( i, j, id );
  };

  inline void replace_triangle(int_t face, int_t e0, int w0,int_t e1,int w1,int_t e2,int w2,int line_)
  {
    //      std::cout << "Replace triangle " << face << std::endl;
    int_t edges[3] = {e0,e1,e2};
    int ways[3] = {w0,w1,w2};
    triangles[face] = triangle_t(face, edges ,ways );      
    for (int i=0;i<3;++i)
      {
	//	std::cout << "   insert nei edge[" << i << "] " << edges[i] << std::endl;
	if (m_edges[edges[i]].triangle(0) == 0)
	  {
	    //	    std::cout << " A" << std::endl;
	    m_edges[edges[i]].triangle(0,face);
	  }
	else if (m_edges[edges[i]].triangle(1) == 0)
	  {
	    //	    std::cout << " B" << std::endl;
	    m_edges[edges[i]].triangle(1,face);
	  }
	else
	  {
	    std::cerr << "edge corrupted "<< edges[i] << std::endl;
	    exit(1);
	  }
      }
  }

  
  inline int_t create_edge(int_t i,int_t j, int id = 0) noexcept
  {
    int_t edge;
    if (m_stack_edges.empty())
      {
	edge = CreateEdge(i,j, id);
	
	classify_edge(edge);	
      }
    else
      {
	edge = m_stack_edges.top();
	m_stack_edges.pop();
	replace_edge(edge,i,j,id);
	
	classify_edge(edge);
      }

    return edge;
  };

  inline int_t create_triangle(int_t e0, int w0,int_t e1,int w1,int_t e2,int w2,int line_) noexcept
  {
    if (m_stack.empty())
      {
	return CreateTriangle(e0,w0,e1,w1,e2,w2,line_);
      }
    else
      {
	int_t face = m_stack.top();
	m_stack.pop();
	replace_triangle(face,e0,w0,e1,w1,e2,w2,line_);
	return face;
      }
  };

  inline void weak_remove_edge(int_t id_) noexcept
  {
    
    declassify_edge(id_);
    m_stack_edges.push(id_);
  };

  inline void weak_remove_triangle(int_t id_) noexcept
  {
    m_stack.push(id_);
    for (int k=0;k<3;++k)
      {
	this->m_edges[triangles[id_].edge(k)].dereftriangle(id_);
      }
  };

  
private:
  
  //!
  //! @brief Store edges.
  //!
  int_t 	m_nedges{0};
  int_t 	mxnedges{0};
  edge_t* 	m_edges;
  
  //!
  //! @brief Store points.
  //!
  int_t 	m_nnodes{0};
  int_t 	mxnnodes{0};
  xy_t*    	xy;

  //!
  //! @brief Store triangles.
  //!
  int_t 	m_ntriangles{0};
  int_t 	mxntriangles{0};
  triangle_t* 	triangles;


  
  HashPair<unsigned int,unsigned int>  	m_hashEdges;

  std::stack<int_t>  m_stack;
  std::stack<int_t>  m_stack_edges;

  int ip{0};

  CavityOperator m_cavityOperator;

};
  
	  
