#ifndef _UTIL_HPP_
#define _UTIL_HPP_

#include <isbinary.hpp>
#include <forward_types.hpp>
#include <set>
#include <map>
#include <algorithm>
#include <cassert>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <boost/static_assert.hpp>
#include <boost/type_traits/is_base_and_derived.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <iostream>
#include <Sequence/SimData.hpp>

template< typename gamete_type,
	  typename vector_allocator_type,
	  typename mutation_type,
	  typename list_allocator_type,
	  template<typename,typename> class vector_type,
	  template<typename,typename> class list_type >
void valid_copy( const vector_type<gamete_type,vector_allocator_type> & gametes,
		 const list_type<mutation_type,list_allocator_type> & mutations,
		 vector_type<gamete_type,vector_allocator_type> & gametes_destination,
		 list_type<mutation_type,list_allocator_type> & mutations_destination )
/*!
  If you ever need to store (and later restore) the state of the population, a naive copy operation
  is not sufficient, because of all the pointers from the gametes container to elements
  of the mutations container.  Use this function instead.
*/
{
  BOOST_STATIC_ASSERT( (boost::is_same<typename gamete_type::mutation_type,mutation_type>::value) );
  typedef typename list_type<mutation_type,list_allocator_type>::iterator literator;  
  typedef typename list_type<mutation_type,list_allocator_type>::const_iterator cliterator;
  typedef typename gamete_type::mutation_container::const_iterator gciterator;
  gametes_destination.clear();
  mutations_destination.clear();
  //copying the mutations is trivial
  std::map<double,literator> mutlookup;
  for( cliterator i = mutations.begin();
       i!=mutations.end();++i)
    {
      literator j = mutations_destination.insert(mutations_destination.end(),*i);
      mutlookup[j->pos]=j;
    }
  /*
    std::copy(mutations.begin(),mutations.end(),
    std::back_inserter(mutations_destination));
  */
  for(unsigned i=0;i<gametes.size();++i)
    {
      //copy construct so that all public, etc., data
      //are properly initialized
      gamete_type new_gamete(gametes[i]);
      new_gamete.mutations.clear();
      new_gamete.smutations.clear();
      for(gciterator itr = gametes[i].mutations.begin() ; 
	  itr != gametes[i].mutations.end() ; ++itr)
	{
	  new_gamete.mutations.push_back( mutlookup[(*itr)->pos] );
	  /*
	    literator itr2 = std::find(mutations_destination.begin(),
				     mutations_destination.end(),
				     **itr);
				     assert(itr2 != mutations_destination.end());
				     new_gamete.mutations.push_back(itr2);
	  */
	}
      for(gciterator itr = gametes[i].smutations.begin() ; 
	  itr != gametes[i].smutations.end() ; ++itr)
	{
	  new_gamete.smutations.push_back( mutlookup[(*itr)->pos] );
	  /*
	  literator itr2 = std::find(mutations_destination.begin(),
				     mutations_destination.end(),
				     **itr);
	  assert(itr2 != mutations_destination.end());
	  new_gamete.smutations.push_back(itr2);
	  */
	}
      gametes_destination.push_back(new_gamete);
    }
}

template<typename mutationtype,
	 typename vector_type_allocator1,
	 typename vector_type_allocator2,
	 typename list_type_allocator,
	 template <typename,typename> class vector_type,
	 template <typename,typename> class list_type >
void remove_fixed_lost( list_type<mutationtype,list_type_allocator> * mutations, 
			vector_type<mutationtype,vector_type_allocator1> * fixations, 
			vector_type<unsigned,vector_type_allocator2> * fixation_times,
			const unsigned & generation,const unsigned & N)
{
  BOOST_STATIC_ASSERT( ( boost::is_base_and_derived<mutation_base,mutationtype>::value) );
  typename list_type<mutationtype,list_type_allocator>::iterator i = mutations->begin();
  
  while(i != mutations->end())
    {
      assert(i->n <= N);			
      i->checked=false;
      if(i->n==N )
	{
	  fixations->push_back(*i);
	  fixation_times->push_back(generation);
	}
      if( i->n == 0 || i->n == N )
	{
	  mutations->erase(i);
	  i=mutations->begin();
	}
      else
	{
	  ++i;
	}
    }
}

template<typename gamete_type>
void adjust_mutation_counts( gamete_type * g , const unsigned & n)
/*!
  \note Will need a specialization if types have other data that need updating
*/
{
  for(unsigned j=0;j< g->mutations.size();++j)
    {
      if( g->mutations[j]->checked==false)
	{
	  g->mutations[j]->n=n;
	  g->mutations[j]->checked=true;
	}
      else
	{
	  g->mutations[j]->n += n;
	}
    }
  for(unsigned j=0;j< g->smutations.size();++j)
    {
      if( g->smutations[j]->checked==false)
	{
	  g->smutations[j]->n=n;
	  g->smutations[j]->checked=true;
	}
      else
	{
	  g->smutations[j]->n += n;
	}
    }
}

template< typename gamete_type,
	  typename vector_type_allocator,
	  template <typename,typename> class vector_type >
	  void update_gamete_list( vector_type<gamete_type,vector_type_allocator > * gametes, 
				   const unsigned & N)
{
  int g = gametes->size()-1;
  typename vector_type<gamete_type,vector_type_allocator >::iterator gbeg;
  while(g>=0)
    {
      gbeg=gametes->begin()+g;
      assert(gbeg->n <= N);
      if( gbeg->n == 0 )
	{
	  gametes->erase(gbeg);
	}
      else
	{
	  int j = gbeg->mutations.size()-1;
	  while(j >= 0)
	    {
	      assert( gbeg->mutations[j]->n <= N );
	      if(gbeg->mutations[j]->n == N|| gbeg->mutations[j]->n == 0)
		{
		  gbeg->mutations.erase(gbeg->mutations.begin()+j);
		}
	      --j;
	    }
	  j = gbeg->smutations.size()-1;
	  while(j >= 0)
	    {
	      assert( gbeg->smutations[j]->n <= N );
	      if(gbeg->smutations[j]->n == N|| gbeg->smutations[j]->n == 0)
		{
		  gbeg->smutations.erase(gbeg->smutations.begin()+j);
		}
	      --j;
	    }
	}
      --g;
    }
  //go backwards through gamete list, and see if there are 
  //any non-unique gametes that we can remove.

  //    seems to have no effect, which is good...
  /*
  typename vector_type<gamete_type,vector_type_allocator >::iterator itr;
  g = gametes->size()-1;
  unsigned s = gametes->size();
  while( g >= 0 )
    {
      gbeg = (gametes->begin()+g);
      itr = std::find(gametes->begin(),gbeg,*gbeg);
      if( itr < gbeg )
	{
	  std::cerr << "erasing: " << (itr-gbeg) << ' ' << g << '\n';
	  itr->n += gbeg->n;
	  size_t diff = g;
	  gametes->erase(gbeg);
	  //gbeg = gametes->begin();
	  //gend = gametes->end()-1;//gbeg+diff-1;
	}
      --g;
    }
  if( s != gametes->size() )
    {
      std::cerr << s << ' ' << gametes->size() << '\n';
    }
  */    
}

template< typename floating_type,
	  typename iterator_type>
void multinomial_sample_gametes(gsl_rng * r,iterator_type gbegin, 
				const size_t & ngametes, 
				const floating_type * efreqs, 
				const unsigned & N)
{
  BOOST_STATIC_ASSERT( (boost::is_convertible<floating_type, double>::value) );
  unsigned n = N;
  double sum=1.;
  for(unsigned i=0;i<ngametes;++i)
    {
      if( (*(efreqs+i)) > 0. )
	{
	  (gbegin+i)->n = gsl_ran_binomial(r,*(efreqs+i)/sum,n);
	  sum -= *(efreqs+i);
	}
      else
	{
	  (gbegin+i)->n = 0;
	}
      n -= (gbegin+i)->n;
      adjust_mutation_counts( &*(gbegin+i), (gbegin+i)->n );
    }
  //assert(n==0);
}

#endif /* _UTIL_HPP_ */
