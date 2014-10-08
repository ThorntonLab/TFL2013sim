//  -*- C++ -*- 
#ifndef _MUTATION_TCC_
#define _MUTATION_TCC_

#include <boost/static_assert.hpp>
#include <boost/type_traits/is_base_and_derived.hpp>
#include <boost/type_traits/is_same.hpp>
#include <gsl/gsl_randist.h>
#include <algorithm>
#include <numeric>
#include <iostream>


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
		const mutation_insertion_policy & mpolicy)
{
  BOOST_STATIC_ASSERT( (boost::is_base_and_derived<mutation_base,typename gamete_type::mutation_type>::value) );
  typedef gamete_base< typename gamete_type::mutation_type, 
    typename gamete_type::mutation_list_type > gamete_base_type;
  BOOST_STATIC_ASSERT( (boost::is_base_and_derived<gamete_base_type,
			gamete_type>::value) || (boost::is_same<gamete_base_type,gamete_type >::value) );
  BOOST_STATIC_ASSERT( (boost::is_same< list_type<typename gamete_type::mutation_type,list_type_allocator >,
			typename gamete_type::mutation_list_type >::value) );

  unsigned ncurrent_classes = gametes->size();
  typename vector_type<gamete_type,vector_type_allocator>::iterator ibeg;
  unsigned nm=0;
  for(unsigned i=0;i<ncurrent_classes;++i)
    {
      ibeg=(gametes->begin()+i);
      unsigned nmuts = gsl_ran_poisson(r,double(ibeg->n)*mu);
      //unsigned nmuts_e = std::min(nmuts,ibeg->n);
      //nm += nmuts_e;
      nm += nmuts;
      //ibeg->n -= nmuts_e;	

      std::vector<unsigned> nm(ibeg->n,0u);
      std::vector<double> pm(ibeg->n,1./double(ibeg->n));
      gsl_ran_multinomial(r,ibeg->n,nmuts,&pm[0],&nm[0]);
#ifndef NDEBUG
      if( !(std::accumulate(nm.begin(),nm.end(),0u) == nmuts ))
	{
	  std::cerr << std::accumulate(nm.begin(),nm.end(),0u) << ' ' << nmuts << '\n';
	}
      assert( std::accumulate(nm.begin(),nm.end(),0u) == nmuts );
#endif
      std::sort(nm.begin(),nm.end(),std::greater<unsigned>());
      for(std::vector<unsigned>::const_iterator itr = nm.begin() ; 
	  itr < nm.end() && *itr>0 ; ++itr )
	{
	  ibeg->n--;
	  gamete_type new_gamete( 1,ibeg->mutations,ibeg->smutations );
	  for(unsigned j=0;j<*itr;++j)
	    {
	      nmuts--;

	      //create a new mutant type to enter the population
	      typename gamete_type::mutation_type new_mutant = mmodel(mutations);
	      
	      //insert the new mutant type into the list of mutations, and record the position of insertion
	      typename gamete_type::mutation_list_type_iterator mitr = mpolicy(mutations,new_mutant);
	      if(mitr->neutral)
		new_gamete.mutations.push_back(mitr);
	      else
		new_gamete.smutations.push_back(mitr);
	    }
	  std::sort(new_gamete.mutations.begin(),new_gamete.mutations.end(),fake_less());
	  std::sort(new_gamete.smutations.begin(),new_gamete.smutations.end(),fake_less());
	  gpolicy(new_gamete,gametes,ncurrent_classes,mutations);
	  ibeg=(gametes->begin()+i);
	}
      assert(nmuts==0);

      /*
      //add new mutations
      unsigned nma=0;
      for( unsigned j = 0; j < nmuts_e; ++j,++nma )
	{
	  //the new gamete is a copy of the i-th gamete + any mutations
	  gamete_type new_gamete( 1,ibeg->mutations );

	  //create a new mutant type to enter the population
	  typename gamete_type::mutation_type new_mutant = mmodel(mutations);
				
	  //insert the new mutant type into the list of mutations, and record the position of insertion
	  typename gamete_type::mutation_list_type_iterator itr = mpolicy(mutations,new_mutant);
	  new_gamete.mutations.push_back(itr);
	  //"correction for large mutation rate?" -- seems to make very little difference
	  if(j+1 == nmuts_e)
	    {
	      for(unsigned k=0 ; k < nmuts-nmuts_e ; ++k,++nma)
		{
		  new_mutant=mmodel(mutations);
		  itr = mpolicy(mutations,new_mutant);
		  new_gamete.mutations.push_back(itr);
		}
	    }
	  std::sort(new_gamete.mutations.begin(),new_gamete.mutations.end(),fake_less());
	  gpolicy(new_gamete,gametes,ncurrent_classes,mutations);
	  ibeg=(gametes->begin()+i);
	}
      assert(nma==nmuts);
      */
    }
  return nm;
}

#endif /* _MUTATION_TCC_ */
