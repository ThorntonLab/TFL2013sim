#ifndef __SAMPLING_FUNCTIONS_HPP__
#define __SAMPLING_FUNCTIONS_HPP__

#include <vector>
#include <map>
#include <functional>
#include <algorithm>
#include <utility>
#include <string>
#include <cassert>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


template<typename iterator_type>
std::vector<unsigned> population_sfs( iterator_type beg,
				      iterator_type end,
				      const unsigned & N)
{
  std::vector<unsigned> psfs(N-1,0);
  while(beg != end)
    {
      if(beg->n >0 && beg->n < N) psfs[beg->n-1]++;
      beg++;
    }
  return psfs;
}

template< typename gamete_type,
	  typename vector_type_allocator,
	  template<typename,typename> class vector_type>
	  std::vector<unsigned> sample(gsl_rng * r,
				       const vector_type<gamete_type,vector_type_allocator > & gametes,
				       const unsigned & n, const unsigned & N)
{
  std::vector<double> freqs(gametes.size(),0);
  std::vector<unsigned> counts(gametes.size(),0);
  for(unsigned i=0;i<gametes.size();++i)
    {
      double f = double(gametes[i].n)/double(N);
      freqs[i]=f;
    }
  gsl_ran_multinomial(r,gametes.size(),n,&freqs[0],&counts[0]);
  return counts;
}

template< typename gamete_type,
	  typename vector_type_allocator,
	  template<typename,typename> class vector_type>
	  std::vector<unsigned> sample_sfs(gsl_rng * r, 
					   const vector_type<gamete_type,vector_type_allocator > & gametes,
					   const unsigned & n, const unsigned & N)
{
  std::vector<unsigned> counts = sample(r,gametes,n,N);
  std::map<double,unsigned> samplemuts;
  std::map<double,unsigned>::iterator itr;
  for(unsigned i=0;i<gametes.size();++i)
    {
      if(counts[i]>0)
	{
	  for(unsigned j=0;j<gametes[i].mutations.size();++j)
	    {
	      itr = samplemuts.find(gametes[i].mutations[j]->pos);
	      if( itr == samplemuts.end() )
		{
		  samplemuts[gametes[i].mutations[j]->pos] = counts[i];
		}
	      else
		{
		  itr->second += counts[i];
		}
	    }
	  for(unsigned j=0;j<gametes[i].smutations.size();++j)
	    {
	      itr = samplemuts.find(gametes[i].smutations[j]->pos);
	      if( itr == samplemuts.end() )
		{
		  samplemuts[gametes[i].smutations[j]->pos] = counts[i];
		}
	      else
		{
		  itr->second += counts[i];
		}
	    }
	}
    }
  std::vector<unsigned> samplesfs(n,0);
  for(itr=samplemuts.begin();itr!=samplemuts.end();++itr)
    {
      samplesfs[itr->second-1]++;
    }
  return samplesfs;
}

struct find_mut_pos : public std::binary_function< std::pair<double,std::string>, double, bool >
{
  inline bool operator()(const std::pair<double,std::string> & pds, const double & d) const
  {
    return (std::fabs(pds.first-d) <= std::numeric_limits<double>::epsilon() );
  }
};

struct sortpos : public std::binary_function< std::pair<double,std::string>,
					      std::pair<double,std::string>, bool >
{
  inline bool operator()( const std::pair<double,std::string> & lhs,
			  const std::pair<double,std::string> & rhs ) const
  {
    return (lhs.first < rhs.first);
  }
};

template< typename gamete_type,
	  typename vector_type_allocator,
	  template<typename,typename> class vector_type>
std::vector< std::pair<double, std::string> > 
ms_sample(gsl_rng * r,
	  const vector_type<gamete_type,vector_type_allocator > & gametes,
	  const unsigned & n, const unsigned & N,
	  bool remove_fixed = true)
{
  std::vector<unsigned> counts = sample(r,gametes,n,N);	
  std::vector< std::pair<double, std::string> > rv;
  std::vector< std::pair<double, std::string> >::iterator itr;

  unsigned individual = 0;
  for(unsigned i=0;i<counts.size();++i)
    {
      for(unsigned j=0;j<counts[i];++j,++individual)
	{
	  for(unsigned mut = 0 ; mut < gametes[i].mutations.size() ; ++mut)
	    {
	      double mutpos = gametes[i].mutations[mut]->pos;
	      itr = std::find_if(rv.begin(),rv.end(),std::bind2nd(find_mut_pos(),mutpos));
	      if( itr == rv.end() )
		{
		  rv.push_back( std::make_pair(mutpos,std::string(n,'0')) );
		  rv[rv.size()-1].second[individual] = '1';
		}
	      else
		{
		  itr->second[individual] = '1';
		}
	    }
	  for(unsigned mut = 0 ; mut < gametes[i].smutations.size() ; ++mut)
	    {
	      double mutpos = gametes[i].smutations[mut]->pos;
	      itr = std::find_if(rv.begin(),rv.end(),std::bind2nd(find_mut_pos(),mutpos));
	      if( itr == rv.end() )
		{
		  rv.push_back( std::make_pair(mutpos,std::string(n,'0')) );
		  rv[rv.size()-1].second[individual] = '1';
		}
	      else
		{
		  itr->second[individual] = '1';
		}
	    }
	}
    }
  assert(individual==n);
		
  if(remove_fixed&&!rv.empty())
    {
      //remove fixations, to make it really ms-like
      itr = rv.end()-1;
      while( itr >= rv.begin() )
	{
	  if(unsigned(std::count(itr->second.begin(),itr->second.end(),'1'))==n)
	    {
	      rv.erase(itr);
	      itr=rv.end();
	    }
	  itr--;
	}
    }
  if(!rv.empty())
    {
      std::sort(rv.begin(),rv.end(),sortpos());
    }
  return rv;
}

template< typename gamete_type,
	  typename vector_type_allocator,
	  template<typename,typename> class vector_type>
std::pair< std::vector< std::pair<double, std::string> > ,
	   std::vector< std::pair<double, std::string> > >
ms_sample_separate(gsl_rng * r,
		   const vector_type<gamete_type,vector_type_allocator > & gametes,
		   const unsigned & n, const unsigned & N,
		   bool remove_fixed = true)
{
  std::vector<unsigned> counts = sample(r,gametes,n,N);	
  std::vector< std::pair<double, std::string> > rvsel,rvneut;
  std::vector< std::pair<double, std::string> >::iterator itr;

  unsigned individual = 0;
  for(unsigned i=0;i<counts.size();++i)
    {
      for(unsigned j=0;j<counts[i];++j,++individual)
	{
	  for(unsigned mut = 0 ; mut < gametes[i].mutations.size() ; ++mut)
	    {
	      double mutpos = gametes[i].mutations[mut]->pos;
	      itr = std::find_if(rvneut.begin(),rvneut.end(),std::bind2nd(find_mut_pos(),mutpos));
	      if( itr == rvneut.end() )
		{
		  rvneut.push_back( std::make_pair(mutpos,std::string(n,'0')) );
		  rvneut[rvneut.size()-1].second[individual] = '1';
		}
	      else
		{
		  itr->second[individual] = '1';
		}
	    }
	  for(unsigned mut = 0 ; mut < gametes[i].smutations.size() ; ++mut)
	    {
	      double mutpos = gametes[i].smutations[mut]->pos;
	      itr = std::find_if(rvsel.begin(),rvsel.end(),std::bind2nd(find_mut_pos(),mutpos));
	      if( itr == rvsel.end() )
		{
		  rvsel.push_back( std::make_pair(mutpos,std::string(n,'0')) );
		  rvsel[rvsel.size()-1].second[individual] = '1';
		}
	      else
		{
		  itr->second[individual] = '1';
		}
	    }
	}
    }
  assert(individual==n);
		
  if(remove_fixed)
    {
      if(!rvneut.empty())
	{
	  //remove fixations, to make it really ms-like
	  itr = rvneut.end()-1;
	  while( itr >= rvneut.begin() )
	    {
	      if(unsigned(std::count(itr->second.begin(),itr->second.end(),'1'))==n)
		{
		  rvneut.erase(itr);
		  itr=rvneut.end();
		}
	      itr--;
	    }
	}
      if(!rvsel.empty())
	{
	  //remove fixations, to make it really ms-like
	  itr = rvsel.end()-1;
	  while( itr >= rvsel.begin() )
	    {
	      if(unsigned(std::count(itr->second.begin(),itr->second.end(),'1'))==n)
		{
		  rvsel.erase(itr);
		  itr=rvsel.end();
		}
	      itr--;
	    }
	}

    }
  if(!rvneut.empty())
    {
      std::sort(rvneut.begin(),rvneut.end(),sortpos());
    }
  if(!rvsel.empty())
    {
      std::sort(rvsel.begin(),rvsel.end(),sortpos());
    }
  return std::make_pair(rvneut,rvsel);
}
#endif 
