#pragma once

#include "hilbert_curve.cpp"
#include "hilbert_curve.hpp"

#include "mesh_t.hpp"

void read_medit(const char * filename,
		delaunay_cad_t& cad_)
{  
  //  HashPair<unsigned int,unsigned int>  	m_hashEdges;
  InputFile::InputFileMedit medit(filename);
  const unsigned int numPoints	= medit.GetNumPoints();
  const unsigned int numEdges	= medit.GetNumEdges();
  cad_.init(numPoints,numEdges);
  
  medit.ReadPoints();
  double coo[2];
  for (unsigned int i=5;i<=numPoints+4;++i)
    {
      unsigned int cod;
      //      points[i] = new MeshGenerator2D::point_t();
      medit.ReadPoint(coo,&cod);
      cad_.set_point(i - 4, coo);      
    }
  
  //
  // Connectivity.
  //
  int_t * 	data 	= new int_t[5*numEdges];
  int_t * 	edges 	= data;
  int_t * 	edgesToEdges = data + 3;
  int_t 	edgesLd = 5;
  int_t 	edgesToEdgesLd = 5;
  //  int_t * edges = new int_t[3*numEdges];  

    medit.ReadEdgeToNode();
    unsigned int	edgeToNode[2];
    unsigned int 	cod[1];
    for (unsigned int i=0;i<numEdges;++i)
      {
	medit.GetEdgeToNode(edgeToNode,
			    cod);
	int_t edge = i+1;
	cad_.node(edge,0,edgeToNode[0]+1);
	cad_.node(edge,1,edgeToNode[1]+1);
	cad_.node(edge,2,cod[0]);
      }
#if 0
    qsort(edges, numEdges, sizeof(int_t)*5,
	  [](const void*a,const void*b){ return (*((int_t*)a) < *((int_t*)b)) ? -1 :((*((int_t*)a) > *((int_t*)b))?1:0 ); });
    
    for (unsigned int i=0;i<numEdges;++i)
      {
	std::cout << " " << edges[edgesLd*i+0];
	std::cout << " " << edges[edgesLd*i+1] ;
	std::cout << " " << edges[edgesLd*i+2] << std::endl;
      }
    
    unsigned int numCurves = 1;
    for (unsigned int i=1;i<numEdges;++i)
      {
	if (edges[edgesLd*i+2] != edges[edgesLd*(i-1)+2])
	  {
	    ++numCurves;
	  }
      }
    std::cout << "numCurves " << numCurves << std::endl;
    
    int_t* numEdgesInCurve = new int_t[numCurves];
    for (int i=0;i<numCurves;++i)
      {
	numEdgesInCurve[i] = 1;
      }
    
    numCurves = 1;
    for (unsigned int i=1;i<numEdges;++i)
      {
	if (edges[edgesLd*i+2] != edges[edgesLd*(i-1)+2])
	  {
	    ++numCurves;
	  }
	else
	  {
	    numEdgesInCurve[numCurves-1]++;
	  }
      }

    for (int i=0;i<numCurves;++i)
      {
	std::cout << "curve[" << i << "].num_edges = " << numEdgesInCurve[i] << std::endl;
      }

    int_t ** cnc_curves = new int_t*[numCurves];
    cnc_curves[0] = &edges[0];
    numCurves = 1;    
    for (unsigned int i=1;i<numEdges;++i)
      {
	if (edges[edgesLd*i+2] != edges[edgesLd*(i-1)+2])
	  {
	    cnc_curves[i] = &edges[edgesLd*i+0];
	    ++numCurves;
	  }
      }
    
    

    
    bool *curves_is_closed        = new bool[numCurves];
    int_t *curves_num_fixed_nodes  = new int_t[numCurves];

    for (unsigned int i=0;i<numCurves;++i)
      {
	std::cout << "curve["<<i<<"]"<< std::endl;
	int_t * c = cnc_curves[i];
	for (int k=0;k< numEdgesInCurve[i];++k)
	  {
	    std::cout << "    " << c[edgesLd*k+0] << " " << c[edgesLd*k+1] << std::endl;
	  }
      }
    
    for (unsigned int i=0;i<numCurves;++i)
      {
	std::cout << "curve["<<i<<"]"<< std::endl;
	int_t * c = cnc_curves[i];
	for (int k=0;k< numEdgesInCurve[i];++k)
	  {
	    std::cout << "    " << c[edgesLd*k+0] << " " << c[edgesLd*k+1] << std::endl;
	  }
      }
#endif
    

				 
#if 0
    {
    //
    // Mark the boundary.
    //
    medit.ReadEdgeToNode();
    for (unsigned int i=0;i<numEdges;++i)
      {
	unsigned int cod,c[2];
	medit.GetEdgeToNode(c,&cod);
	//	std::cout << " " <<c[0]+1+4 << " " << c[1]+1+4 <<  std::endl;
	fast.m_hashEdges.Add(c[0]+1+4,
			     c[1]+1+4,
			     i+1);
	
      }
    }
#endif
    
}


struct delaunay_calculator_t
{
public:
  //!
  //! @brief Constructor.
  //! @param filename_ The file name.
  //!
  delaunay_calculator_t(const char * filename_)
    : m_filename(filename_)
  {
    read_medit(filename_,this->m_cad);
    this->m_cad.analysis();
    m_mesh = new mesh_t();
  };

  void debug_mode(bool debug_mode)
  {
    this->m_debug_mode = debug_mode;
    m_mesh->debug_mode(this->m_debug_mode);
  };
  
  //!
  //! @brief Build the box.
  //!
  void build_box_mesh()
  {
    double pmin[2];
    double pmax[2];
    
    this->m_cad.box(pmin,pmax);
    this->m_mesh->CreateInitialMesh(xy_t(pmin[0],pmin[1],0),
				    xy_t(pmax[0],pmax[1],0) );
    this->m_mesh->print_mesh("initial");
  };

  //!
  //! @brief Retrieve the boundary.
  //!
  void retrieve_boundary()
  {
    if (this->m_cad.nedges() > 0)
      {
	this->m_mesh->retrieve_boundary(this->m_cad);
      }
  };

  //!
  //! @brief Optimize.
  //!
  void init_metric(double target_size)
  {
    //
    // Init a metric
    //
    const int_t nnodes = m_mesh->nnodes();
    for (int_t node = 1; node <= nnodes;++node)
      {
	auto node_ptr = m_mesh->get_node_ptr(node);
	node_ptr->m_metric[0] = 1.0 / target_size / target_size;
	node_ptr->m_metric[1] = 0.0;
	node_ptr->m_metric[2] = 1.0 / target_size / target_size * 64.0;
      }
  };

  //!
  //! @brief Optimize.
  //!
  void optimize()
  {
    this->m_mesh->optimize();
  };

  inline void optimize_boundary(double target_size)
  {
    this->m_mesh->optimize_boundary();
  };

  inline void optimize_interior(double target_size)
  {
    this->m_mesh->optimize_interior();
  };
  
  inline void output(const char * filename)
  {
    this->m_mesh->print(filename,0);
  }

  inline void analysis()
  {
    this->m_mesh->analysis();
  }

  void remove_box(bool keep)
  {
    this->m_mesh->remove_box(keep);
  };

  static int comp(const void * a,const void * b)
{
  const long int * a_ = (const long int * )a;
  const long int * b_ = (const long int * )b;
  if (a_[2] < b_[2]) return -1;
  else
    if (a_[2] > b_[2]) return 1;
  else
    return 0;
}

  void insert_random()
  {


  int N = 30000;
  //int N2  = 2*2*2*2 * 2*2*2*2;
  int N2 = 1000;
  int length = (N2)*(N2);
  double * pos = new double[length*4];
  srandom(0);
  int next=0;
  for (int i=0;i<N2;++i)
    {
      for (int j=0;j<N2;++j)
	{
	  pos[4*next+0] = double(random())/double(RAND_MAX);
	  pos[4*next+1] = double(random())/double(RAND_MAX);
	  ++next;
	}
    }
  double bb[4] = {0,0,1,1};
#if 1
  coordinates2d(length, 
		bb,
		pos,
		4,
		(unsigned long long int * )&pos[2],
		4);

  qsort(pos,length,4*sizeof(double),comp);
#endif
  EllapsedTime e; 
  e.Start();
  
  std::cout << "benchmark insert "<< length << " random points " << std::endl;
  int_t triangle = 1;
  for (int i=0;i<length;++i)
    {
    //    if (i%1000==0)
      //      std::cout << i << std::endl;
      //      std::cout << pos[4*i+0] << " " << pos[4*i+1] << " " << *((unsigned long long int * )&pos[4*i+2]) << std::endl;
      	int_t node 	= this->m_mesh->create_node(pos[4*i+0],pos[4*i+1]);
	 triangle 	= this->m_mesh->find_triangle(node,triangle);
	if (triangle)
	  {

	    auto node_ptr = m_mesh->get_node_ptr(node);

	    double x = node_ptr->x;
	    double y = node_ptr->y;
	    double h = 0.125/2.0 + 0.25 * std::abs(x);
	    node_ptr->m_metric[0] = 1.0;
	    node_ptr->m_metric[1] = 0.0;
	    node_ptr->m_metric[2] = 1.0;
	    this->m_mesh->InsertPoint(triangle,
				      node);
	}
      else
	{
	  fprintf(stderr,"triangle not found\n");
	  exit(1);
	}

    }
  e.Stop();
  std::cout << "elapsed " << e << std::endl;

  };

  void insert_cad_points()
  {
    double xy[2];
    int_t cad_npoints = this->m_cad.npoints();
    for (int_t cad_node=1;cad_node<=cad_npoints;++cad_node)
      {
	m_cad.get_point(cad_node,xy);	
	int_t node 	= this->m_mesh->create_node(xy[0],xy[1]);

	int_t triangle 	= this->m_mesh->find_triangle(node);
	if (triangle)
	  {

	    auto node_ptr = m_mesh->get_node_ptr(node);

	    double x = node_ptr->x;
	    double y = node_ptr->y;
	    double h = 0.125/2.0 + 0.25 * std::abs(x);
	    node_ptr->m_metric[0] = 1.0;
	    node_ptr->m_metric[1] = 0.0;
	    node_ptr->m_metric[2] = 1.0;
	    this->m_mesh->InsertPoint(triangle, node, false);
	}
      else
	{
	  fprintf(stderr,"triangle not found\n");
	  exit(1);
	}
      }
  };

private:
  const char * m_filename;
  delaunay_cad_t m_cad;
  mesh_t* m_mesh;
  bool m_debug_mode{};
};
