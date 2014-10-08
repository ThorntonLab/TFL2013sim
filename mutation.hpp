#ifndef _MUTATION_HPP_
#define _MUTATION_HPP_

#include <forward_types.hpp>
#include <fwd_functional.hpp>
#include <algorithm>
#include <limits>
#include <cmath>
#include <boost/bind.hpp>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

template<typename mutation_type,
	 typename list_type >
inline mutation_type neutral_mutations_inf_sites_big(gsl_rng * r,list_type * mutations)
{
  //double pos = gsl_rng_uniform(r);
  double pos = double(::round(gsl_ran_flat(r,0.,std::numeric_limits<unsigned>::max())));
  while( std::find_if(mutations->begin(),mutations->end(),boost::bind(mutation_at_pos(),_1,pos)) != mutations->end())	
    {
      //pos = gsl_rng_uniform(r);
      pos = double(::round(gsl_ran_flat(r,0.,std::numeric_limits<unsigned>::max())));
    }
  return mutation_type(pos,0.,1.);
}

template<typename mutation_type,
	 typename list_type >
inline mutation_type neutral_mutations_inf_sites(gsl_rng * r,list_type * mutations)
{
  double pos = gsl_rng_uniform(r);
  while( std::find_if(mutations->begin(),mutations->end(),boost::bind(mutation_at_pos(),_1,pos)) != mutations->end())	
    {
      pos = gsl_rng_uniform(r);
    }
  return mutation_type(pos,0.,1);
}

template<typename mutation_type,
	 typename list_type>
mutation_type background_selection(gsl_rng * r,const double & plethal,list_type * mutations)
{
	double pos = gsl_rng_uniform(r);
	while( std::find_if(mutations->begin(),mutations->end(),boost::bind(mutation_at_pos(),_1,pos)) != mutations->end())	
	{
		pos = gsl_rng_uniform(r);
	}
	if( gsl_rng_uniform(r) <= plethal )
	{
		return mutation_type(pos,-1.,1,1.);
	}
	return mutation_type(pos,0.,1);
}

template<typename mutation_type,
	 typename list_type>
mutation_type background_selection(gsl_rng * r,const double & plethal,const double & lethal_pos,list_type * mutations)
{
	if( gsl_rng_uniform(r) <= plethal )
	{
		return mutation_type(lethal_pos,-1.,1,1.);
	}
	double pos = gsl_rng_uniform(r);
	while( std::find_if(mutations->begin(),mutations->end(),boost::bind(mutation_at_pos(),_1,pos)) != mutations->end())	
	{
		pos = gsl_rng_uniform(r);
	}
	return mutation_type(pos,0.,1);
}

template<typename mutation_type,
	 typename list_type>
mutation_type weak_selection(gsl_rng * r, const double & twoN,const double & sigma, const double & offset,list_type * mutations)
{
	double pos = gsl_rng_uniform(r);
	while( std::find_if(mutations->begin(),mutations->end(),boost::bind(mutation_at_pos(),_1,pos)) != mutations->end())	
	{
		pos = gsl_rng_uniform(r);
	}
	return mutation_type(pos,(gsl_ran_gaussian(r,sigma)+offset)/twoN,1);
}

template< typename gamete_type,
	  typename mutation_model,
	  typename gamete_insertion_policy,
	  typename mutation_insertion_policy,
	  typename vector_type_allocator,
	  typename list_type_allocator,
	  template<typename,typename> class vector_type,
	  template<typename,typename> class list_type>
	  unsigned mutate(gsl_rng * r, 
			  vector_type<gamete_type,vector_type_allocator > * gametes, 
			  list_type<typename gamete_type::mutation_type,list_type_allocator > * mutations, 
			  const double & mu,
			  const mutation_model & mmodel,
			  const gamete_insertion_policy & gpolicy, 
			  const mutation_insertion_policy & mpolicy);
#endif /* _MUTATION_HPP_ */
#include <mutation.tcc>
