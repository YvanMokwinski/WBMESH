#ifndef __header_MeditReader_h__
#define __header_MeditReader_h__
extern "C"
{
#include "extern_medit.h"
}

//#include "Common/Math.h"

//#include <list>
#include "Math.h"
#include <list>
#if 0
#include "Common/Type.h"
#include "Common/Debug.hpp"
#include "Common/Assert.hpp"
#include "Common/Exceptions/ExceptionSwitch.hpp"
#include "Common/Mpi.h"



#define SIMLAB_DEPRECATED_MSG(text) __attribute__((deprecated(# text)))
#endif

template<typename _value_type> inline  _value_type Min(const _value_type&a_,const _value_type&b_)
{
  return (a_ < b_) ? a_ : b_;
}


template<typename _value_type> inline  _value_type Max(const _value_type&a_,const _value_type&b_)
{
  return (a_ < b_) ? b_ : a_;
}

#if 0
// #define MNS_MAX_THREAD 16
class FreeDeleteObject
{
 protected:
  std::list<void*> m_references;
 public:
  FreeDeleteObject(){};
  virtual ~FreeDeleteObject(){};

  bool HasNoReference() const
  {
    return m_references.empty();
  };
  
  void InsertReference(void * reference_)
  {
    m_references.push_back(reference_);
  };

  void RemoveReference(void * reference_)
  {
    m_references.remove(reference_);
    if (this->HasNoReference())
      {
	delete this;
      }
  };
  
};
#endif
#if 0
#ifndef NDEBUG  
#define CHECK_INTERFACE_IMPLEMENTATION(__interface_method_to_call__)	\
  {									\
    static bool call = false;						\
    if( call == true )							\
      {									\
	cerr << #__interface_method_to_call__ << " not implemented" << endl;throw("Interface method not implemented!");\
      }									\
    call = true;							\
    try {								\
      (__interface_method_to_call__);					\
      call = false;							\
    }									\
    catch ( ... )							\
      {									\
	call = false;							\
	throw;								\
      }									\
  }
#else
#define CHECK_INTERFACE_IMPLEMENTATION(__interface_method_to_call__)
#endif  
#endif

//#define CHECK_INTERFACE_IMPLEMENTATION(__interface_method_to_call__)


//#include "Common/System.h"
//#include "Interfaces/Geometry/CRTP_Point.hpp"

namespace InputFile
{
  class InputFileMedit
  {
  protected:
    int 			m_version;
    int 			m_inm;
    unsigned int 		m_topologyDimension;
    unsigned int		m_dimension;
    unsigned int		m_numTetrahedrons;
    unsigned int		m_numPyramids;
    unsigned int		m_numWedges;
    unsigned int		m_numHexahedrons;
    unsigned int		m_numQuadrilaterals;
    unsigned int		m_numTriangles;
    unsigned int		m_numEdges;
    unsigned int		m_numPoints;
  public:
    
    inline InputFileMedit(const std::string& filename_)
      {
	std::string filename(filename_);
	this->m_inm = GmfOpenMesh(filename.c_str(),GmfRead,&this->m_version,&this->m_dimension);
	if (!this->m_inm)
	  {
	    std::cerr << "unable to openfile" << filename_ << std::endl;
	    exit(1);
	  }
	
	this->m_numTriangles		= (unsigned int)GmfStatKwd(this->m_inm,GmfTriangles);
	this->m_numQuadrilaterals	= (unsigned int)GmfStatKwd(this->m_inm,GmfQuadrilaterals);
	this->m_numPyramids		= (unsigned int)GmfStatKwd(this->m_inm,GmfPyramids);
	this->m_numWedges 		= (unsigned int)GmfStatKwd(this->m_inm,GmfWedges);
	this->m_numHexahedrons 		= (unsigned int)GmfStatKwd(this->m_inm,GmfHexahedra);
	this->m_numTetrahedrons 	= (unsigned int)GmfStatKwd(this->m_inm,GmfTetrahedra);
	this->m_numEdges 		= (unsigned int)GmfStatKwd(this->m_inm,GmfEdges);
	this->m_numPoints		= (unsigned int)GmfStatKwd(this->m_inm,GmfVertices);
      };
  

    inline ~InputFileMedit()
    {
      if (this->m_inm)
	{
	  this->m_inm = GmfCloseMesh(this->m_inm);     
	  }
    };
    


    template <typename T> inline void ReadPoint(T * c,
						unsigned int cod_[])
      { 
#ifndef NDEBUG
	Debug::IsTrue(__TRACE__,point_.GetDimension()==this->m_dimension);
#endif
	double c0[3];
	int icod;
	switch(this->m_dimension)
	  {
	  case 3:
	    {
	      GmfGetLin(this->m_inm,GmfVertices,&c0[0],&c0[1],&c0[2],&icod);
	      cod_[0] = icod;	  
	      break;
	    }
	  case 2:
	    {
	      GmfGetLin(this->m_inm,GmfVertices,&c0[0],&c0[1],&icod);
	      cod_[0] = icod;
	      break;
	    }
	  default:
	    {
	      break;
	    }
	  }
	for (unsigned int dimIndex=0;dimIndex<this->m_dimension;++dimIndex)
	  {
	    c[dimIndex] = c0[dimIndex];
	  }
      };


    inline unsigned int 		GetDimension() const { return this->m_dimension; };
    inline unsigned int 		GetTopologyDimension() const { return this->m_topologyDimension; }
    inline unsigned int		GetNumEdges() const { return this->m_numEdges; };
    inline unsigned int		GetNumPoints() const { return this->m_numPoints; };

    inline void 			ReadEdgeToNode()
    {
      GmfGotoKwd(this->m_inm,GmfEdges); 
    };

    inline void 			ReadPoints()
    {
      GmfGotoKwd(this->m_inm,GmfVertices); 
    };

    inline void 			GetEdgeToNode(unsigned int	edgeToNode_[2],
						      unsigned int 	cod_[1])
    {
      int localcod;
      int ci[2];
      GmfGetLin(this->m_inm,GmfEdges,&ci[0],&ci[1],&localcod);
      edgeToNode_[0] 	= ci[0]-1;
      edgeToNode_[1] 	= ci[1]-1;
      cod_[0]		= localcod;
    };
  };
};
  
#endif
