#pragma once
struct FaceEdge
{
  int_t first;
  int second;
  FaceEdge(int_t f,int s) : first(f),second(s)
  {
      
  };
};

  class CavityOperator
  {
  private:

    //
    // Store faces.
    //
    int_t m_numFaces;
    int_t m_mxNumFaces;
    int_t * m_faces;

    //
    // Waiging faces.
    //
    int_t m_numWaitingFaces;
    int_t m_mxNumWaitingFaces;
    int_t * m_waitingFaces;

    //
    // Store boundary edges.
    //
    int_t m_numBoundaryEdges;
    int_t m_mxNumBoundaryEdges;
    FaceEdge * m_boundaryEdges;

    //
    // Store interior edges.
    //    
    int_t m_numInteriorEdges;
    int_t m_mxNumInteriorEdges;
    int_t * m_interiorEdges;

    HashCavityTable<int_t>		m_hashCavity;
    
  public:
    CavityOperator()
    {
      this->m_numFaces = 0;
      this->m_mxNumFaces = 64;
      m_faces = (int_t*)malloc(sizeof(int_t)*this->m_mxNumFaces);

      this->m_numBoundaryEdges = 0;
      this->m_mxNumBoundaryEdges = 64;
      m_boundaryEdges = (FaceEdge*)malloc(sizeof(FaceEdge)*this->m_mxNumBoundaryEdges);

      this->m_numInteriorEdges = 0;
      this->m_mxNumInteriorEdges = 64;
      m_interiorEdges = (int_t*)malloc(sizeof(int_t)*this->m_mxNumInteriorEdges);


      this->m_numWaitingFaces = 0;
      this->m_mxNumWaitingFaces = 64;
      m_waitingFaces = (int_t*)malloc(sizeof(int_t)*this->m_mxNumWaitingFaces);
    };

    inline int_t interior_edge(int_t i) const noexcept
    {
      return this->m_interiorEdges[i];
    };

    inline int_t ninterior_edges() const noexcept
    {
      return this->m_numInteriorEdges;
    };
    
    inline FaceEdge boundary_edge(int_t i) const noexcept
    {
      return this->m_boundaryEdges[i];
    };
    inline int_t nboundary_edges() const noexcept
    {
      return this->m_numBoundaryEdges;
    };
    
    inline void clear_boundary_edges() noexcept
    {
      this->m_numBoundaryEdges = 0;
    };
    inline void clear_interior_edges() noexcept
    {
      this->m_numInteriorEdges = 0;
    };

    inline void clear_hash(){this->m_hashCavity.clear();};

    inline int_t try_get_edge(int_t p,int_t q) const
    {
      return this->m_hashCavity.Get(p,q);
    };

    inline int_t add_edge(int_t p,int_t q,int_t edge) 
    {
      return this->m_hashCavity.Add(p,q,edge);
    };
    

    //    void reset_faces() {m_faces.clear();};
    inline void reset_faces() { m_numFaces = 0;};
    template<typename f>
    inline void apply_faces(f ff)
    {
      for (int_t i=0;i<this->m_numFaces;++i)
	{
	  ff(this->m_faces[i]);
	}
    };
#if 0
    inline void Set(mesh_t*p)
    {
      this->m = p;
    };
#endif

  public : 
    
    inline void Reset()
    {
      reset_faces();
      clear_boundary_edges();
      clear_interior_edges();
    };
    
    inline void add_interior_edge(int_t edge)
    {
      if (this->m_numInteriorEdges == this->m_mxNumInteriorEdges)
	{
	  this->m_mxNumInteriorEdges *=2;
	  std::cout << "WARNING REALLOC CAVITY INTERIOR EDGES" << std::endl;
	  this->m_interiorEdges = (int_t*)realloc(this->m_interiorEdges,sizeof(int_t)*this->m_mxNumInteriorEdges);
	}
      this->m_interiorEdges[this->m_numInteriorEdges++] = edge;
      // this->m_listCavityInteriorEdges.push_front(edge);
    };
      
    inline void AddBoundaryEdge(int_t edge,int way)
    {
      
      if (this->m_numBoundaryEdges == this->m_mxNumBoundaryEdges)
	{
	  this->m_mxNumBoundaryEdges *=2;
	  std::cout << "WARNING REALLOC CAVITY BOUNDARY EDGES" << std::endl;
	  this->m_boundaryEdges = (FaceEdge*)realloc(this->m_boundaryEdges,sizeof(FaceEdge)*this->m_mxNumBoundaryEdges);
	}
      this->m_boundaryEdges[this->m_numBoundaryEdges++] = FaceEdge(edge,way);
      
      //      this->m_listCavityBoundaryEdges.push_front(FaceEdge(edge,way));
      
    };


    //
    // Check if the face exists.
    //
    inline bool HasFace(int_t face_) const noexcept
    {
      for (int_t i=0;i<this->m_numFaces;++i)
	{
	  if (face_ == this->m_faces[i])
	    {
	      return true;
	    }
	}
      return false;
    };

    
    inline void  AddFace(int_t face_) noexcept
    {
      if (false == this->HasFace(face_))
	{
	  if (this->m_numFaces == this->m_mxNumFaces)
	    {
	      this->m_mxNumFaces*=2;
	      std::cout << "WARNING REALLOC CAVITY FACES" << std::endl;
	      m_faces = (int_t*)realloc(m_faces,sizeof(int_t)*this->m_mxNumFaces);
	    }
	  
	  this->m_faces[this->m_numFaces++] = face_;

	}
    };

    inline void add_waiting_face(int_t face) noexcept
    {
      if (this->m_numWaitingFaces == this->m_mxNumWaitingFaces)
	{
	  this->m_mxNumWaitingFaces*=2;
	  std::cout << "WARNING REALLOC CAVITY WAITING FACES" << std::endl;
	  m_waitingFaces = (int_t*)realloc(m_waitingFaces,sizeof(int_t)*this->m_mxNumWaitingFaces);
	}
      
      this->m_waitingFaces[this->m_numWaitingFaces++] = face;
    };
    
  public:
    template <typename criterion_t,typename func, typename getface,typename getedge>
    inline bool Compute(int_t 	face_,
			int_t 	point_,
			criterion_t 	criterion,
			func      	conditionCavityBoundaryEdge_,
			getface gf,
			getedge ge)
    {
      if (!face_)
	{
	  return false;
	}
      this->Reset();
	
      this->AddFace(face_);


      //      std::list<int_t> waiting;
      //      waiting.push_front(face_);
      m_numWaitingFaces = 0;
      add_waiting_face(face_);

      while (m_numWaitingFaces > 0)
      //      while (false == waiting.empty())//(waitingSize>0)
	{	  
	  //
	  // Pop the list.
	  //

	  face_ = m_waitingFaces[--m_numWaitingFaces];
	  
	  // face_ = waiting.front();	    
	  //	  waiting.pop_front();


	  
	  /** ON PARCOURS LES ARETES **/
	  for (unsigned int k=0;k<3;k++)
	    {	  
	      //	      auto 	edge 	= m->triangles[face_].edge(k);
	      //	      auto 	neiFace	= m->m_edges[edge].otriangle(face_);
	      auto 	edge 	= gf(face_)->edge(k);
	      auto 	neiFace	= ge(edge)->otriangle(face_);

	      //		std::cout << "neighbors edge " << edge << " : " << m->m_edges[edge].triangle(0) << " " << m->m_edges[edge].triangle(1) << std::endl;
	      if (!neiFace)
		{
		  //		    std::cout << "add boundary edge" << std::endl;
		  AddBoundaryEdge(edge,
				  gf(face_)->way(k));
				  //				  m->triangles[face_].way(k));
		}
	      else
		{
		  /* 
		     on check si l arete est contrainte 
		     dans le cas ou l'arete est contrainte on skip le triangle
		  */
		  //		    std::cout << "interior edge " << this->HasFace(neiFace) << std::endl;
		  if (false == this->HasFace(neiFace))
		    {
		      //			std::cout << "apply criterion " << std::endl;
		      if (criterion(neiFace,point_))
			{
			  if (conditionCavityBoundaryEdge_(edge))
			    {
			      //			      this->AddBoundaryEdge(edge,m->triangles[face_].way(k));
			      this->AddBoundaryEdge(edge,gf(face_)->way(k));
			    }
			  else
			    {
			      this->add_interior_edge(edge);				
			      this->AddFace(neiFace);
			      add_waiting_face(neiFace);
			      //   waiting.push_front(neiFace);
			    }
			}
		      else
			{
			  this->AddBoundaryEdge(edge,
						gf(face_)->way(k));
						//m->triangles[face_].way(k));
			}
		    }
		  else
		    {

		    }
		}/* end for */      
	    }
	}

      return true;
    }

  };
