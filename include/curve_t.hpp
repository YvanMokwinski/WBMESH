#pragma once

struct curve_t
{
  bool          m_is_closed{};
  bool          m_is_non_manifold{};
  int_t 	m_n{};
  int_t * 	m_e2n{};
  int_t 	m_e2n_ld{};
  int_t * 	m_e2e{};
  int_t 	m_e2e_ld{};
  curve_t(int_t n_)
  {
    m_n = n_;
    m_e2n = (int_t*)malloc(sizeof(int_t)*2*n_);
    m_e2e = (int_t*)malloc(sizeof(int_t)*2*n_);
    m_e2e_ld = 2;
    m_e2n_ld = 2;
  };

  void e2n(int_t i,int_t node0_,int_t node1_)
  {
    this->m_e2n[this->m_e2n_ld*i+0] = node0_;
    this->m_e2n[this->m_e2n_ld*i+1] = node1_;
  };

  //
  // compute edge 2 edge
  //
  void compute_e2e(int_t       nedges,
		   const int_t*e2n,
		   int_t       e2nld,
		   int_t*      e2e,
		   int_t       e2eld)
  {
#define N 64
#define LL -100000
    int_t * link = new int_t[64];
    for (int_t i=0;i<64;++i) link[i] = LL;
    for (int_t i=0;i<nedges;++i)
      {
	for (int k=0;k<2;++k)
	  {
	    int_t node = e2n[e2nld*i+k];
	    int_t h = node % N;
	    int_t l = link[h];
	    e2e[e2eld*i+k] = l;
	    link[h] = -(1 + e2eld*i+k);
	  }
      }
    delete[]link;

    for (int_t i=nedges-1;i>=0;--i)
      {
	for (int k=1;k>=0;--k)
	  {
	    auto node = e2n[i*e2nld+k];	  
	    auto l = e2e[i*e2eld+k];
	    if (l < 0)
	      {
		int_t origin_at = i*e2eld+k;
		e2e[origin_at] = 0;

		//
		//
		//
		int_t jj = i + 1;
		for (; l!=LL; l = e2e[-l-1])
		  {		  
		    auto at = -l-1;
		    auto neinode = e2n[ (at / e2eld) *e2nld+at % e2eld];
		    if (neinode==node)
		      {
			e2e[at] = jj;
			jj = (at / e2eld) + 1;
		      }
		  }
	      
	      
		int_t previous_at = origin_at;
		int_t next_at = -l-1;
		while (l!=LL)
		  {
		    auto at = next_at;
		    next_at = e2e[next_at];
		    l = e2e[-l-1];
		  
		    auto inei = at / e2eld;
		    auto lp = at % e2eld;
		    auto neinode = e2n[inei*e2nld+lp];
		    if (neinode == node)
		      {
			e2e[origin_at] = inei;
			origin_at = at;

			if ( (next_at != LL) && (previous_at != origin_at) )
			  {
			    e2e[previous_at] = next_at;
			  }
		      }
		    previous_at = at;
		  }
	      }	  
	  }
      }

#undef N
  };
  
  void analyze()
  {
    compute_e2e(m_n,
		m_e2n,
		m_e2n_ld,
		m_e2e,
		m_e2e_ld);
    int_t nb=0;
    int_t nb1=0;
    int_t nbi=0;
    for (int_t i=0;i<m_n;++i)
      {
	if (m_e2e[i*2+0] && !m_e2e[i*2+1])
	  {
	    nb++;
	  }
	else if (m_e2e[i*2+1] && !m_e2e[i*2+0])
	  {
	    nb++;
	  }
	else if (!m_e2e[i*2+1] && !m_e2e[i*2+0])
	  {
	    nb1++;
	  }
	else
	  {
	    nbi++;
	  }
      }
    std::cout << "nb " << nb << std::endl;
    std::cout << "nb1 " << nb1 << std::endl;
    std::cout << "nbi " << nbi << std::endl;
    this->m_is_closed = (nb>0);
    //  bool          m_is_non_manifold{};    
  };
  
};
