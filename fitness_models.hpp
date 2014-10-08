#ifndef _FITNESS_MODELS_HPP_
#define _FITNESS_MODELS_HPP_

#include <forward_types.hpp>
#include <fwd_functional.hpp>
#include <cassert>
#include <functional>
#include <cmath>
#ifndef NDEBUG
#include <iostream>
#include <utility>
#include <limits>
#endif

template<typename iterator_type >
double no_selection( const iterator_type & g1, const iterator_type &g2)
{
  return 1.;
}

template<typename iterator_type>
struct find_mut_iterator : public std::binary_function<std::pair<iterator_type,unsigned>,iterator_type,bool >
{
  inline bool operator()(const std::pair<iterator_type,unsigned> & p,const iterator_type & i) const
  {
    return p.first == i;
  }
};

struct multiplicative_diploid
{
  typedef double result_type;
  template< typename iterator_type >
  inline double operator()(const iterator_type & g1, const iterator_type & g2, const double scaling = 1.) const
  //assumes addresses of mutations are stored in order...
  {
    double fitness = 1.;
    typename iterator_type::value_type::mutation_list_type_iterator ib1,ib2;
    typename iterator_type::value_type::mutation_container::const_iterator b1=g1->smutations.begin(),
      e1=g1->smutations.end(),
      b2=g2->smutations.begin(),
      e2=g2->smutations.end();
 
#ifndef NDEBUG
    //This block is an explicit, although inefficient, way of calculating the fitness of this
    //genotype.  It is used here for debugging purposes, to compare to the results of the 
    //faster algorithm below
    std::vector< std::pair< typename  iterator_type::value_type::mutation_list_type_iterator, unsigned > > mcounts;
    typename std::vector< std::pair< typename  iterator_type::value_type::mutation_list_type_iterator, unsigned > >::iterator mci;
    for( ; b1 < e1 ; ++b1)
      {
	assert( !(*b1)->neutral );
	for( mci = mcounts.begin() ; mci != mcounts.end() ; ++mci )
	  {
	    if( mci->first->pos == (*b1)->pos )
	      {
		break;
	      }
	  }
	if (mci < mcounts.end())
	  {
	    mci->second++;
	  }
	else
	  {
	    mcounts.push_back( make_pair(*b1,1) );
	  }
      }
    for( ; b2 < e2 ; ++b2)
      {
	assert( !(*b2)->neutral );
	for( mci = mcounts.begin() ; mci != mcounts.end() ; ++mci )
	  {
	    if( mci->first->pos == (*b2)->pos )
	      {
		break;
	      }
	  }
	if (mci < mcounts.end())
	  {
	    mci->second++;
	  }
	else
	  {
	    mcounts.push_back( make_pair(*b2,1) );
	  }
      } 
    double fcheck = 1.;
    unsigned nhet=0,nhom=0;
    for( mci = mcounts.begin() ; mci !=mcounts.end() ; ++mci )
      {
	if( mci->second == 1 )
	  {
	    fcheck *= (1.+scaling*mci->first->s*mci->first->h);
	    ++nhet;
	  }
	else if( mci->second == 2)
	  {
	    fcheck *= (1.+scaling*mci->first->s);
	    ++nhom;
	  }
	else
	  {
	    abort();
	  }
      }
    b1=g1->smutations.begin();
    e1=g1->smutations.end();
    b2=g2->smutations.begin();
    e2=g2->smutations.end();
    unsigned nhet2=0,nhom2=0;
#endif

    //This is a much faster way to calculate fitnesses,
    //as it just compares addresses in memory, and 
    //does little in the way of dereferencing and storing
    for( ; b1 < e1 ; ++b1 )
      {
	bool found = false;
	ib1 = *b1;
	if(ib1->n) //if mutation has frequency > 0
	  {
	    for( ; b2 < e2 && !((*b2)->pos > (ib1)->pos); ++b2)					
	      {
		ib2 = *b2;
		if ( ib2 == ib1 )
		  //mutation at same position (note this is a comparison of iterators, not values)
		  {
		    assert(ib1->s == ib2->s); //just making sure
		    assert(ib1->pos == ib2->pos);
		    fitness *= (1. + scaling*(ib2)->s);
#ifndef NDEBUG
		    nhom2++;
#endif
		    found=true;
		  }
		else if( (ib2)->pos < (ib1)->pos ) 
		  //b2 points to a unique mutation that comes before b1
		  //(this is the only numeric comparison)
		  {
		    assert(ib2->pos != ib1->pos);
		    assert(ib2->pos < ib1->pos);
		    //fitness += std::log((1. + (ib2)->h*(ib2)->s));
		    fitness *= (1. + scaling*(ib2)->h*(ib2)->s);
#ifndef NDEBUG
		    nhet2++;
#endif
		  }
		/*
		  else 
		  //position of b2 is > b1, and b1 is heterozygous
		  {
		  //this never happens, b/c of outer if statement
		  abort();
		  assert(ib2->pos > ib1->pos);
		  fitness += std::log((1. + (ib1)->h*(ib1)->s));
		  #ifndef NDEBUG
		  nhet2++;
		  #endif
		  found=true;
		  break;
		  }
		*/
	      }
	    if(!found)
	      {
		//fitness += std::log((1. + (ib1)->h*(ib1)->s));	
		fitness *= (1. + scaling*(ib1)->h*(ib1)->s);	
#ifndef NDEBUG
		nhet2++;
#endif
	      }
	  }
      }
    for( ; b2 < e2 ; ++b2)
      {	
	ib2=*b2;
	//fitness += std::log((1. + (ib2->h*ib2->s)));	
	if(ib2->n)
	  {
	    fitness *= (1. + scaling*(ib2->h*ib2->s));	
#ifndef NDEBUG
	    nhet2++;
#endif
	  }
      }
  //fitness = std::exp(fitness);
    assert (g1->smutations.size()+g2->smutations.size() == (2*nhom+nhet));
    assert ( std::fabs(fcheck-fitness) <= 10.*std::numeric_limits<double>::epsilon() );
#ifndef NDEBUG
    if(nhet != nhet2 || nhom != nhom2) 
      {
	std::cerr << fcheck << ' ' << fitness << ' '
		  << nhet << ' ' << nhom << ' '
		  << nhet2 << ' ' << nhom2 << '\n';
	b1=g1->smutations.begin();
	e1=g1->smutations.end();
	b2=g2->smutations.begin();
	e2=g2->smutations.end();
	for(;b1<e1;++b1)
	  {
	    std::cerr << (*b1)->pos << ' ' << (*b1)->s << ' ';
	  }
	std::cerr << '\n';
	for(;b2<e2;++b2)
	  {
	    std::cerr << (*b2)->pos << ' ' << (*b2)->s << ' ';
	  }
	std::cerr << '\n';
	assert(nhet == nhet2 && nhom == nhom2);
      }
#endif
  return (fitness);
  //return std::max(fitness,0.);
  }
};

struct additive_diploid
{
  typedef double result_type;
  template< typename iterator_type >
  inline double operator()(const iterator_type & g1, const iterator_type & g2) const
  //assumes addresses of mutations are stored in order...
  {
    double sum_sh = 0.;
    typename iterator_type::mutation_list_type::const_iterator ib1,ib2;
    typename iterator_type::mutation_container::const_iterator b1=g1->smutations.begin(),
      e1=g1->smutations.end(),
      b2=g2->smutations.begin(),
      e2=g2->smutations.end();

    for( ; b1 < e1 ; ++b1 )
      {
	bool found = false;
	ib1 = *b1;
	if(ib1->n)
	  {
	    for( ; b2 < e2 && !((*b2)->pos > (ib1)->pos); ++b2)					
	      {
		ib2 = *b2;
		if(ib2->n)
		  {
		    if ( ib2 == ib1 )
		      //mutation at same position (note this is a comparison of iterators, not values)
		      {
			assert(ib1->s == ib2->s); //just making sure
			sum_sh += (ib2)->s;
			found=true;
		      }
		    else if( (ib2)->pos < (ib1)->pos ) 
		      //b2 points to a unique mutation that comes before b1
		      //(this is the only numeric comparison)
		      {
			sum_sh += (ib2)->h*(ib2)->s;
		      }
		    else 
		      //position of b2 is > b1, and b1 is heterozygous
		      {
			assert(ib2->pos > ib1->pos);
			sum_sh += (ib1)->h*(ib1)->s;
			found=true;
			break;
		      }
		  }
	      }
	    if(!found)
	      {
		sum_sh +=(ib1)->h*(ib1)->s;	
	      }
	  }
      }
    for( ; b2 < e2 ; ++b2)
      {	
	ib2=*b2;
	if(ib2->n)
	  {
	    sum_sh += (ib2)->h*(ib2)->s;	
	  }
      }
    return std::max(1.+sum_sh,0.);
  }
};
#endif /* _FITNESS_MODELS_HPP_ */
