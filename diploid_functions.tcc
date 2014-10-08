//  -*- C++ -*- 
#ifndef _DIPLOID_FUNCTIONS_TCC_
#define _DIPLOID_FUNCTIONS_TCC_

#include <cassert>
#include <algorithm>
#include <iterator>
#include <functional>
#include <util.hpp>
#include <cmath>
#include <numeric>
//#include <vector> //added for approx WF
#ifndef NO_KAHAN_SUMS
#include <kahan.hpp>
#endif
#include <boost/static_assert.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/type_traits/is_base_and_derived.hpp>

#ifndef NDEBUG
#include <iostream>
#endif

/*
template< typename gamete_type,
	  typename diploid_fitness_function,
	  typename floating_type,
	  typename vector_type_allocator,
	  template<typename,typename> class vector_type>
inline void calc_mean_fitness( floating_type *wbar, //std::vector< floating_type> & expected_gamete_frequencies,
			       floating_type * expected_gamete_frequencies,
			       const vector_type<gamete_type,vector_type_allocator > * gametes, 
			       const unsigned & twoN, const diploid_fitness_function & ff)
{
  if(gametes->empty()){std::cerr << "gametes empty\n"; abort();}
  typename vector_type<gamete_type,vector_type_allocator >::const_iterator ibeg,jbeg;
  double pqw;

  for(unsigned i=0;i<gametes->size();++i)
    {
      ibeg=(gametes->begin()+i);
      double w= ff( *ibeg, *ibeg ) ;
      //double lnp = std::log(double(ibeg->n)/double(twoN));
      double p = double(ibeg->n)/double(twoN);
      //pqw = std::exp(lnp+lnp + std::log(w));
      pqw = p*p*w;
      *wbar += (w>0.) ? pqw : 0.;
      *(expected_gamete_frequencies+i) += pqw;
      for(unsigned j=i+1 ; j < gametes->size() ; ++j)
	{
	  jbeg=(gametes->begin()+j);
	  w = ff( *ibeg, *jbeg );
	  //double lnq  = std::log(double(jbeg->n)/double(twoN));
	  double q = double(jbeg->n)/double(twoN);
	  //pqw  = std::exp( lnp + lnq + std::log(w));
	  pqw=p*q*w;
	  *wbar += (w>0.) ? 2.*pqw : 0.;
	  *(expected_gamete_frequencies+i) += pqw;
	  *(expected_gamete_frequencies+j) += pqw;
	}
    }
}
*/
template< typename gamete_type,
	  typename diploid_fitness_function,
	  typename floating_type,
	  typename vector_type_allocator,
	  template<typename,typename> class vector_type>
inline void calc_mean_fitness( floating_type *wbar,
			       floating_type * expected_gamete_frequencies,
			       const vector_type<gamete_type,vector_type_allocator > * gametes, 
			       const unsigned & twoN, const diploid_fitness_function & ff,
			       const double & f)
{
#ifndef NDEBUG
  if(gametes->empty()){std::cerr << "gametes empty\n"; abort();}
#endif
  typename vector_type<gamete_type,vector_type_allocator >::const_iterator ibeg=gametes->begin(),jbeg;
  double pqw,p,q,w;
  unsigned j;
  for(unsigned i=0;i<gametes->size();++i)
    {
      if(ibeg->n)
	{
	  p = double(ibeg->n)/double(twoN);
	  if ( (w= ff( ibeg, ibeg )) >0. )
	    {
	      pqw = (p*p*(1.-f) + p*f)*w;
	      *wbar += pqw;
	      *(expected_gamete_frequencies+i) += pqw;
	    }
	  
	  j=i+1;
	  jbeg=ibeg+1;
	  for( ; j < gametes->size() ; ++j)
	    {
	      if(jbeg->n)
		{
		  q = double(jbeg->n)/double(twoN);
		  if((w = ff( ibeg, jbeg ))>0.)
		    {
		      pqw=p*q*(1.-f)*w;
		      *wbar += 2.*pqw;
		      *(expected_gamete_frequencies+i) += pqw;
		      *(expected_gamete_frequencies+j) += pqw;
		    }
		}
	      ++jbeg;
	    }
	}
      ++ibeg;
    }
}



template< typename gamete_type,
	  typename diploid_fitness_function,
	  typename floating_type,
	  typename vector_type_allocator,
	  template<typename,typename> class vector_type>
inline void calc_mean_fitness_monte_carlo( gsl_rng * r,
					   floating_type *wbar,
					   floating_type * expected_gamete_frequencies,
					   const vector_type<gamete_type,vector_type_allocator > * gametes, 
					   const unsigned & twoN, const diploid_fitness_function & ff,
					   const double & f)
{
   for( unsigned ind = 0 ; ind < twoN/2 ;  ) //individual-based sampling
    {
      typedef typename vector_type<gamete_type,vector_type_allocator >::const_iterator giter;
      unsigned gam1 = unsigned(gsl_ran_flat(r,1.,double(twoN+1))),
	gam2 = unsigned(gsl_ran_flat(r,1.,double(twoN+1))),
	sum1=0,sum2=0;
      giter i=gametes->begin(),j=i;
      do
	{
	  sum1 += i->n;
	  if(gam1 <= sum1)
	    {
	      break;
	    }
	  ++i;
	}
      while(1);
      
      //account for inbreeding coefficient f (Crow & Kimura, p. 65)
      if(f>0.)
	{
	  if( gsl_rng_uniform(r) <= f )
	    {
	      j=i;
	    }
	}
      else
	{
	  do
	    {
	      sum2 += j->n;
	      if(gam2 <= sum2)
		{
		  break;
		}
	      ++j;
	    }
	  while(1);
	}
      assert(i < gametes->end());
      assert(j < gametes->end());
      double w = ff(i,j);
      //accept individual prop. to w
      if ( w > 0. && gsl_rng_uniform(r) <= w)
      {
	*wbar += w;
	expected_gamete_frequencies[i - gametes->begin()] += 1./double(twoN);
	expected_gamete_frequencies[j - gametes->begin()] += 1./double(twoN);
	++ind;
      }
    }
  *wbar /= double(twoN/2);
}

template< typename gamete_type,
	  typename diploid_fitness_function,
	  typename floating_type,
	  typename vector_type_allocator,
	  template<typename,typename> class vector_type>
inline void calc_mean_fitness_monte_carlo( gsl_rng * r,
					   floating_type *wbar,
					   floating_type * expected_gamete_frequencies,
					   const vector_type<gamete_type,vector_type_allocator > * gametes, 
					   const unsigned & twoN_curr, const unsigned & twoN_next,
					   const diploid_fitness_function & ff,
					   const double & f)
{
   for( unsigned ind = 0 ; ind < twoN_next/2 ;  ) //individual-based sampling
    {
      typedef typename vector_type<gamete_type,vector_type_allocator >::const_iterator giter;
      unsigned gam1 = unsigned(gsl_ran_flat(r,1.,double(twoN_curr+1))),
	gam2 = unsigned(gsl_ran_flat(r,1.,double(twoN_curr+1))),
	sum1=0,sum2=0;
      giter i=gametes->begin(),j=i;
      do
	{
	  sum1 += i->n;
	  if(gam1 <= sum1)
	    {
	      break;
	    }
	  ++i;
	}
      while(1);
      
      //account for inbreeding coefficient f (Crow & Kimura, p. 65)
      if(f>0.)
	{
	  if( gsl_rng_uniform(r) <= f )
	    {
	      j=i;
	    }
	}
      else
	{
	  do
	    {
	      sum2 += j->n;
	      if(gam2 <= sum2)
		{
		  break;
		}
	      ++j;
	    }
	  while(1);
	}
      assert(i < gametes->end());
      assert(j < gametes->end());
      double w = ff(i,j);
      if ( w > 0. && gsl_rng_uniform(r) <= w)
	{
	  *wbar += w;
	  expected_gamete_frequencies[i - gametes->begin()] += 1./double(twoN_next);
	  expected_gamete_frequencies[j - gametes->begin()] += 1./double(twoN_next);
	  ++ind;
	}
    }
  *wbar /= double(twoN_next/2);
}

template< typename gamete_type,
	  typename vector_type_allocator,
	  template<typename,typename> class vector_type,
	  typename diploid_fitness_function >
double sample_diploid(gsl_rng * r, vector_type<gamete_type,vector_type_allocator > * gametes, 
		      const unsigned & twoN, const diploid_fitness_function & ff,const double & f = 0.)
  //hermaphroditic diploids
{
  BOOST_STATIC_ASSERT( (boost::is_base_and_derived<mutation_base,typename gamete_type::mutation_type>::value) );
  typedef gamete_base< typename gamete_type::mutation_type, typename gamete_type::mutation_list_type > gamete_base_type;
  BOOST_STATIC_ASSERT( (boost::is_base_and_derived<gamete_base_type,
			gamete_type>::value) || (boost::is_same<gamete_base_type,gamete_type >::value) );

#ifdef NO_KAHAN_SUMS
  double wbar = 0.;
  std::vector<double> expected_gamete_frequencies(gametes->size(),0.);
#else
  kahan_sum<double> wbar;
  std::vector<kahan_sum<double> > expected_gamete_frequencies(gametes->size(),kahan_sum<double>());
#endif


  calc_mean_fitness(&wbar,&expected_gamete_frequencies[0],gametes,twoN,ff,f);
  //normalize expected frequqncies
#ifdef NO_KAHAN_SUMS
  std::transform(expected_gamete_frequencies.begin(),expected_gamete_frequencies.end(),
	    expected_gamete_frequencies.begin(),std::bind2nd(std::divides<double>(),wbar));
#else
  std::transform(expected_gamete_frequencies.begin(),expected_gamete_frequencies.end(),
	    expected_gamete_frequencies.begin(),std::bind2nd(std::divides<kahan_sum<double> >(),wbar));						
#endif
  multinomial_sample_gametes(r, gametes->begin(),gametes->size(),&expected_gamete_frequencies[0],twoN);
    
  update_gamete_list(gametes,twoN);
  return wbar;
}

template< typename gamete_type,
	  typename vector_type_allocator,
	  template<typename,typename> class vector_type,
	  typename diploid_fitness_function >
double sample_diploid_mc(gsl_rng * r, vector_type<gamete_type,vector_type_allocator > * gametes, 
			 const unsigned & twoN, const diploid_fitness_function & ff,const double & f = 0.)
  //hermaphroditic diploids, individual-based sampling
{
  BOOST_STATIC_ASSERT( (boost::is_base_and_derived<mutation_base,typename gamete_type::mutation_type>::value) );
  typedef gamete_base< typename gamete_type::mutation_type, typename gamete_type::mutation_list_type > gamete_base_type;
  BOOST_STATIC_ASSERT( (boost::is_base_and_derived<gamete_base_type,
			gamete_type>::value) || (boost::is_same<gamete_base_type,gamete_type >::value) );

#ifdef NO_KAHAN_SUMS
  double wbar = 0.;
  std::vector<double> expected_gamete_frequencies(gametes->size(),0.);
#else
  kahan_sum<double> wbar;
  std::vector<kahan_sum<double> > expected_gamete_frequencies(gametes->size(),kahan_sum<double>());
#endif
  calc_mean_fitness_monte_carlo(r,&wbar,&expected_gamete_frequencies[0],gametes,twoN,ff,f);
  typename vector_type<gamete_type,vector_type_allocator >::iterator i = gametes->begin();
  for(unsigned j = 0 ; j < expected_gamete_frequencies.size() ; ++j)
    {
      i->n = unsigned( round ( double(twoN)*expected_gamete_frequencies[j]) );
      adjust_mutation_counts(&*i,i->n);
      ++i;
    }
  update_gamete_list(gametes,twoN);
  return wbar;
}



template< typename gamete_type,
	  typename vector_type_allocator,
	  template<typename,typename> class vector_type,
	  typename diploid_fitness_function >
double sample_diploid(gsl_rng * r, vector_type<gamete_type,vector_type_allocator > * gametes, 
		      const unsigned & twoN_curr, const unsigned & twoN_next,
		      const diploid_fitness_function & ff,
		      const double & f)
  //hermaphroditic diploids
{
  BOOST_STATIC_ASSERT( (boost::is_base_and_derived<mutation_base,typename gamete_type::mutation_type>::value) );
  typedef gamete_base< typename gamete_type::mutation_type, typename gamete_type::mutation_list_type > gamete_base_type;
  BOOST_STATIC_ASSERT( (boost::is_base_and_derived<gamete_base_type,
			gamete_type>::value) || (boost::is_same<gamete_base_type,gamete_type >::value) );
#ifdef NO_KAHAN_SUMS
  double wbar = 0.;
  std::vector<double> expected_gamete_frequencies(gametes->size(),0.);
#else
  kahan_sum<double> wbar;
  std::vector<kahan_sum<double> > expected_gamete_frequencies(gametes->size(),kahan_sum<double>());
#endif
  calc_mean_fitness(&wbar,&expected_gamete_frequencies[0],gametes,twoN_curr,ff,f);
  //normalize expected frequqncies
#ifdef NO_KAHAN_SUMS
  std::transform(expected_gamete_frequencies.begin(),expected_gamete_frequencies.end(),
	    expected_gamete_frequencies.begin(),std::bind2nd(std::divides<double>(),wbar));
#else
  std::transform(expected_gamete_frequencies.begin(),expected_gamete_frequencies.end(),
	    expected_gamete_frequencies.begin(),std::bind2nd(std::divides<kahan_sum<double> >(),wbar));						
#endif
  multinomial_sample_gametes(r, gametes->begin(),gametes->size(),&expected_gamete_frequencies[0],twoN_next);
  update_gamete_list(gametes,twoN_next);
  return wbar;
}

template< typename gamete_type,
	  typename vector_type_allocator,
	  template<typename,typename> class vector_type,
	  typename diploid_fitness_function >
double sample_diploid_mc(gsl_rng * r, vector_type<gamete_type,vector_type_allocator > * gametes, 
		      const unsigned & twoN_curr, const unsigned & twoN_next,
		      const diploid_fitness_function & ff,
		      const double & f)
  //hermaphroditic diploids
{
  BOOST_STATIC_ASSERT( (boost::is_base_and_derived<mutation_base,typename gamete_type::mutation_type>::value) );
  typedef gamete_base< typename gamete_type::mutation_type, typename gamete_type::mutation_list_type > gamete_base_type;
  BOOST_STATIC_ASSERT( (boost::is_base_and_derived<gamete_base_type,
			gamete_type>::value) || (boost::is_same<gamete_base_type,gamete_type >::value) );
#ifdef NO_KAHAN_SUMS
  double wbar = 0.;
  std::vector<double> expected_gamete_frequencies(gametes->size(),0.);
#else
  kahan_sum<double> wbar;
  std::vector<kahan_sum<double> > expected_gamete_frequencies(gametes->size(),kahan_sum<double>());
#endif
  calc_mean_fitness_monte_carlo(r,&wbar,&expected_gamete_frequencies[0],gametes,twoN_curr,twoN_next,ff,f);
  typename vector_type<gamete_type,vector_type_allocator >::iterator i = gametes->begin();
  for(unsigned j = 0 ; j < expected_gamete_frequencies.size() ; ++j)
    {
      i->n = unsigned( round ( double(twoN_next)*expected_gamete_frequencies[j]) );
      adjust_mutation_counts(&*i,i->n);
      ++i;
    }
  update_gamete_list(gametes,twoN_next);
  return wbar;
  return wbar;
}

template< typename gamete_type,
	  typename vector_type_allocator,
	  template<typename,typename> class vector_type>
void pick_gametes(gsl_rng * r,
		  vector_type<gamete_type,vector_type_allocator > * gametes,
		  typename vector_type<gamete_type,vector_type_allocator >::iterator * ibeg,
		  typename vector_type<gamete_type,vector_type_allocator >::iterator * jbeg,
		  const unsigned & max_elements)
/*!
  picks 2 non-identical gametes
*/
{
  BOOST_STATIC_ASSERT( (boost::is_base_and_derived<mutation_base,typename gamete_type::mutation_type>::value) );
  typedef gamete_base< typename gamete_type::mutation_type, typename gamete_type::mutation_list_type > gamete_base_type;
  BOOST_STATIC_ASSERT( (boost::is_base_and_derived<gamete_base_type,
			gamete_type>::value) || (boost::is_same<gamete_base_type,gamete_type >::value) );

  *ibeg=gametes->end(),*jbeg=gametes->end();
  //get current gamete counts
  unsigned sum=0;
  unsigned i=0;
  typedef typename vector_type<gamete_type,vector_type_allocator >::iterator vtype_iterator;
  vtype_iterator itr;
  for(itr=gametes->begin();itr<gametes->begin()+max_elements;++itr)
    {
      sum+=itr->n;
    }
  double p = gsl_rng_uniform(r);
  unsigned temp=0;
  unsigned ibeg_index=0;
  itr=gametes->begin();
  for( i=0;i<max_elements;++i,++itr)
    {
      temp += itr->n;
      if( p <= double(temp)/double(sum) )
	{
	  ibeg_index=i;
	  *ibeg = itr;
	  sum -= itr->n;
	  i=max_elements;
	}
    }
  assert(*ibeg < (gametes->begin()+max_elements));	
  p = gsl_rng_uniform(r);
  temp=0;
  itr=gametes->begin();
  for( i=0;i<max_elements;++i,++itr)
    {
      if(i!=ibeg_index)
	{
	  temp += itr->n;
	  if( p < double(temp)/double(sum))
	    {
	      *jbeg = itr;
	      i=max_elements;
	    }
	}
    }
  assert(*jbeg < (gametes->begin()+max_elements));	
  assert(*ibeg != *jbeg);
}



template< typename gamete_type,
	  typename vector_type_allocator,
 	  template<typename,typename> class vector_type,
	  typename map_function>
unsigned recombine(gsl_rng * r, vector_type<gamete_type,vector_type_allocator > * gametes,
		   const unsigned & twoN, const double & littler,const map_function & mf,
		   const double & f = 0)
{
  BOOST_STATIC_ASSERT( (boost::is_base_and_derived<mutation_base,typename gamete_type::mutation_type>::value) );
  typedef gamete_base< typename gamete_type::mutation_type, typename gamete_type::mutation_list_type > gamete_base_type;
  BOOST_STATIC_ASSERT( (boost::is_base_and_derived<gamete_base_type,
			gamete_type>::value) || (boost::is_same<gamete_base_type,gamete_type >::value) );

  typedef typename gamete_type::mutation_container mcont;
  typedef typename mcont::iterator mcont_iterator;
  typedef typename vector_type<gamete_type, vector_type_allocator>::iterator vtype_iterator;
  //"individual-based recombination"
  //1. determine # recombinants in whole pop
  //2. assign them randomly to pairs of gametes chosen from the gametes list
  //3. if they are the same gamete, ignore, else do recombination thing.

  //calc freq of homozygotes
  
  double fAA=0.;
  if(gametes->size() < twoN)
    {
      const double twoNsq = double(twoN)*double(twoN);
      for(vtype_iterator i=gametes->begin() ;i<gametes->end();++i)
	{
	  fAA += ( (double(i->n)*double(i->n)*(1.-f))/twoNsq + (double(i->n)/twoN)*f );
	}
    }
  unsigned NRECS = unsigned(gsl_ran_poisson(r,(1.-fAA)*double(twoN)*littler));
  unsigned doublecheck=NRECS;
#ifndef NDEBUG
  unsigned NRECS_DONE = 0;
#endif
  vtype_iterator ibeg,jbeg;
  const unsigned ncurrent_classes = gametes->size();
  while(NRECS > 0)
    {
      pick_gametes(r,gametes,&ibeg,&jbeg, ncurrent_classes);
      unsigned nm1=ibeg->mutations.size()+ibeg->smutations.size(),
	nm2=jbeg->mutations.size()+jbeg->smutations.size();
      //if one gametes carries 0 mutations, and the other carries 1,
      //it is impossible for recombination to generate any new types,
      //so we skip those cases...
      if(!(std::min(nm1,nm2)==0 && std::max(nm1,nm2)==1))
	{
	  ibeg->n--;
	  jbeg->n--;
	      
	  gamete_type new_gamete1(1, mcont(),mcont()),new_gamete2(new_gamete1);
	  double pos = mf();
	  //pointer arithmetic over a range of pointers.  apologies...
	  for( mcont_iterator itr = ibeg->mutations.begin() ;
	       itr != ibeg->mutations.end() ; ++itr )
	    {
	      if( (*itr)->pos <= pos )
		{
		  new_gamete1.mutations.push_back(*itr);
		}
	      else
		{
		  new_gamete2.mutations.push_back(*itr);
		}
	    }
	      
	  for( mcont_iterator itr = jbeg->mutations.begin() ; 
	       itr != jbeg->mutations.end() ; ++itr )
	    {
	      if( (*itr)->pos <= pos )
		{
		  new_gamete2.mutations.push_back(*itr);
		}
	      else
		{
		  new_gamete1.mutations.push_back(*itr);
		}
	    }
	      
	  for( mcont_iterator itr = ibeg->smutations.begin() ;
	       itr != ibeg->smutations.end() ; ++itr )
	    {
	      if( (*itr)->pos <= pos )
		{
		  new_gamete1.smutations.push_back(*itr);
		}
	      else
		{
		  new_gamete2.smutations.push_back(*itr);
		}
	    }
	      
	  for( mcont_iterator itr = jbeg->smutations.begin() ; 
	       itr != jbeg->smutations.end() ; ++itr )
	    {
	      if( (*itr)->pos <= pos )
		{
		  new_gamete2.smutations.push_back(*itr);
		}
	      else
		{
		  new_gamete1.smutations.push_back(*itr);
		}
	    }
	      
	      
	  std::sort(new_gamete1.mutations.begin(),new_gamete1.mutations.end(),fake_less());
	  std::sort(new_gamete1.smutations.begin(),new_gamete1.smutations.end(),fake_less());
	  std::sort(new_gamete2.mutations.begin(),new_gamete2.mutations.end(),fake_less());
	  std::sort(new_gamete2.smutations.begin(),new_gamete2.smutations.end(),fake_less());
	      
	  vtype_iterator itr = std::find(gametes->begin(),gametes->end(),new_gamete1);
	  if( itr == gametes->end() )
	    {
	      gametes->push_back(new_gamete1);
	    }
	  else
	    {
	      itr->n++;
	    }
	  itr=find(gametes->begin(),gametes->end(),new_gamete2);
	  if( itr == gametes->end() )
	    { 
	      gametes->push_back(new_gamete2);
	    }
	  else
	    {
	      itr->n++;
	    }	
	}
      NRECS--;
#ifndef NDEBUG
      NRECS_DONE++;
#endif
    }
  assert(NRECS_DONE==doublecheck);
#ifndef NDEBUG
  unsigned sum=0;
  for(ibeg=gametes->begin();ibeg!=gametes->end();++ibeg)
    {
      sum+=ibeg->n;
    }
  assert(sum==twoN);
#endif
  return doublecheck;
}

template< typename gamete_type,
	  typename vector_type_allocator,
	  template<typename,typename> class vector_type>
void pick_other_gamete(gsl_rng * r,
		       vector_type<gamete_type,vector_type_allocator > * gametes,
		       const unsigned & max_elements, const unsigned & i,
		       typename vector_type<gamete_type,vector_type_allocator >::iterator * jbeg)
{
  unsigned sum=0,temp=0;
  typedef typename vector_type<gamete_type,vector_type_allocator >::iterator vtype_iterator;
  vtype_iterator itr = gametes->begin();
  for(unsigned j=0;j<max_elements;++j,++itr)
    {
      if(j!=i)
	{
	  sum += itr->n;
	}
    }
  double p = gsl_rng_uniform(r);
  itr=gametes->begin();
  for(unsigned j=0;j<max_elements;++j,++itr)
    {
      if(j!=i)
	{
	  temp += itr->n;
	}
      if( p <= double(temp)/double(sum) )
	{
	  *jbeg = itr;
	  j=max_elements;
	}
    }
  assert( (*jbeg)->n > 0 );
}

template< typename gamete_type,
	  typename vector_type_allocator,
 	  template<typename,typename> class vector_type,
	  typename map_function>
unsigned recombine_multi(gsl_rng * r, vector_type<gamete_type,vector_type_allocator > * gametes,
			 const unsigned & twoN, const double & littler,const map_function & mf)
{
  BOOST_STATIC_ASSERT( (boost::is_base_and_derived<mutation_base,typename gamete_type::mutation_type>::value) );
  typedef gamete_base< typename gamete_type::mutation_type, typename gamete_type::mutation_list_type > gamete_base_type;
  BOOST_STATIC_ASSERT( (boost::is_base_and_derived<gamete_base_type,
			gamete_type>::value) || (boost::is_same<gamete_base_type,gamete_type >::value) );

  typedef typename gamete_type::mutation_container mcont;
  typedef typename mcont::iterator mcont_iterator;
  typedef typename vector_type<gamete_type, vector_type_allocator>::iterator vtype_iterator;

  unsigned NRECS=0;
  vtype_iterator ibeg,jbeg;
  const unsigned ncurrent_classes = gametes->size();
  unsigned sum = 0;
  //bool mflag=false;

  for(unsigned i=0 ; i < ncurrent_classes ; ++i)
    {
      //1. Figure out how many recombination events originate in this gamete class
      ibeg = (gametes->begin()+i);
      //double p = double(ibeg->n)/double(twoN);
      unsigned nr = gsl_ran_poisson(r,double(ibeg->n)*littler);
      NRECS += nr;
      //2. Figure out how many rec events occur on each gamete
      std::vector<unsigned> nrg(ibeg->n,0u);
      std::vector<double> pr(ibeg->n,1./double(ibeg->n));
      gsl_ran_multinomial(r,ibeg->n,nr,&pr[0],&nrg[0]);
      std::sort(nrg.begin(),nrg.end(),std::greater<unsigned>());
      //3. Perform the recombination events
      for(std::vector<unsigned>::const_iterator itr = nrg.begin() ;
	  itr < nrg.end() && *itr>0 ; ++itr )
	{
	  pick_other_gamete(r,gametes,ncurrent_classes,i,&jbeg);
	  assert(ibeg != jbeg);
	  assert(jbeg != (gametes->begin()+ncurrent_classes));
	  unsigned nm1=ibeg->mutations.size()+ibeg->smutations.size(),
	    nm2=jbeg->mutations.size()+jbeg->smutations.size();
	  //if one gametes carries 0 mutations, and the other carries 1,
	  //it is impossible for recombination to generate any new types,
	  //so we skip those cases...
	  if(!(std::min(nm1,nm2)==0 && std::max(nm1,nm2)==1))
	    {
	      ibeg->n--;
	      jbeg->n--;

	  gamete_type new_gamete1(1, mcont(),mcont()),new_gamete2(new_gamete1);
	  if( *itr == 1 )
	    {
	      double pos = mf();
	      //pointer arithmetic over a range of pointers.  apologies...
	      for( mcont_iterator mitr = ibeg->mutations.begin() ;
		   mitr != ibeg->mutations.end() ; ++mitr )
		{
		  if( (*mitr)->pos <= pos )
		    {
		      new_gamete1.mutations.push_back(*mitr);
		    }
		  else
		    {
		      new_gamete2.mutations.push_back(*mitr);
		    }
		}
	      
	      for( mcont_iterator mitr = jbeg->mutations.begin() ; 
		   mitr != jbeg->mutations.end() ; ++mitr )
		{
		  if( (*mitr)->pos <= pos )
		    {
		      new_gamete2.mutations.push_back(*mitr);
		    }
		  else
		    {
		      new_gamete1.mutations.push_back(*mitr);
		    }
		}
	      
	      for( mcont_iterator mitr = ibeg->smutations.begin() ;
		   mitr != ibeg->smutations.end() ; ++mitr )
		{
		  if( (*mitr)->pos <= pos )
		    {
		      new_gamete1.smutations.push_back(*mitr);
		    }
		  else
		    {
		      new_gamete2.smutations.push_back(*mitr);
		    }
		}
	      
	      for( mcont_iterator mitr = jbeg->smutations.begin() ; 
		   mitr != jbeg->smutations.end() ; ++mitr )
		{
		  if( (*mitr)->pos <= pos )
		    {
		      new_gamete2.smutations.push_back(*mitr);
		    }
		  else
		    {
		      new_gamete1.smutations.push_back(*mitr);
		    }
		}
	    }
	  else
	    {
	      assert(*itr>1);
	      //mflag=true;
	      std::cout << *itr << ' ';
	      //3a. get the positions of the recombination events
	      std::vector<double> positions(*itr);
	      for(unsigned j=0;j<*itr;++j)
		{
		  positions[j] = mf();
		}
	      std::sort(positions.begin(),positions.end());

	      //3b. create the recombinant gametes
	      bool odd = true;
	      std::vector<double>::const_iterator vd1=positions.begin(),vd2=positions.begin();
	      mcont_iterator neut_itr = ibeg->mutations.begin(),sel_itr = ibeg->smutations.begin();
	      for( ; vd1 < positions.end() ; ++vd1 )
		{
		  for( ; neut_itr < ibeg->mutations.end() && (*neut_itr)->pos <= *vd1 ; 
		       ++neut_itr )
		    {
		      if(odd)
			{
			  new_gamete1.mutations.push_back(*neut_itr);
			}
		      else
			{
			  new_gamete2.mutations.push_back(*neut_itr);
			}
		    }
		  for( ; sel_itr < ibeg->smutations.end() && (*sel_itr)->pos <= *vd1 ; 
		       ++sel_itr )
		    {
		      if(odd)
			{
			  new_gamete1.smutations.push_back(*sel_itr);
			}
		      else
			{
			  new_gamete2.smutations.push_back(*sel_itr);
			}
		    }
		  odd=!odd;
		}
	      for( ; neut_itr < ibeg->mutations.end() ; ++neut_itr )
		{
		  if(odd)
		    {
		      new_gamete1.mutations.push_back(*neut_itr);
		    }
		  else
		    {
		      new_gamete2.mutations.push_back(*neut_itr);
		    }
		}
	      for( ; sel_itr < ibeg->smutations.end() ; ++sel_itr)
		{
		  if(odd)
		    {
		      new_gamete1.smutations.push_back(*sel_itr);
		    }
		  else
		    {
		      new_gamete2.smutations.push_back(*sel_itr);
		    }
		}
	      neut_itr = jbeg->mutations.begin(),sel_itr = jbeg->smutations.begin();
	      odd = true;
	      for( ; vd2 < positions.end() ; ++vd2 )
		{
		  for( ; neut_itr < jbeg->mutations.end()&&(*neut_itr)->pos <= *vd2 ;
		       ++neut_itr )
		    {
		      if(odd)
			{
			  new_gamete2.mutations.push_back(*neut_itr);
			}
		      else
			{
			  new_gamete1.mutations.push_back(*neut_itr);
			}
		    }
		  for( ; sel_itr < jbeg->smutations.end()&&(*sel_itr)->pos <= *vd2 ; ++sel_itr )
		    {
		      if(odd)
			{
			  new_gamete2.smutations.push_back(*sel_itr);
			}
		      else
			{
			  new_gamete1.smutations.push_back(*sel_itr);
			}
		    }
		  odd=!odd;
		}
	      for( ; neut_itr < jbeg->mutations.end() ; ++neut_itr )
		{
		  if(odd)
		    {
		      new_gamete2.mutations.push_back(*neut_itr);
		    }
		  else
		    {
		      new_gamete1.mutations.push_back(*neut_itr);
		    }
		}
	      for( ; sel_itr < jbeg->smutations.end() ; ++sel_itr)
		{
		  if(odd)
		    {
		      new_gamete2.smutations.push_back(*sel_itr);
		    }
		  else
		    {
		      new_gamete1.smutations.push_back(*sel_itr);
		    }
		}
	      assert( neut_itr == jbeg->mutations.end() );
	      assert( sel_itr == jbeg->smutations.end() );
	    }
	  std::sort(new_gamete1.mutations.begin(),new_gamete1.mutations.end(),fake_less());
	  std::sort(new_gamete1.smutations.begin(),new_gamete1.smutations.end(),fake_less());
	  std::sort(new_gamete2.mutations.begin(),new_gamete2.mutations.end(),fake_less());
	  std::sort(new_gamete2.smutations.begin(),new_gamete2.smutations.end(),fake_less());
	  assert( (ibeg->smutations.size()+ibeg->mutations.size() +
		   jbeg->smutations.size()+jbeg->mutations.size()) ==
		  (new_gamete1.mutations.size()+new_gamete1.smutations.size() +
		   new_gamete2.mutations.size()+new_gamete2.smutations.size()) );
	  vtype_iterator itr = std::find(gametes->begin(),gametes->end(),new_gamete1);
	  if( itr == gametes->end() )
	    {
	      gametes->push_back(new_gamete1);
	    }
	  else
	    {
	      itr->n++;
	    }
	  itr=find(gametes->begin(),gametes->end(),new_gamete2);
	  if( itr == gametes->end() )
	    { 
	      gametes->push_back(new_gamete2);
	    }
	  else
	    {
	      itr->n++;
	    }
	    }
	}
    }
  /*
  if(mflag)
    std::cout << '\n';
  */
#ifndef NDEBUG
  sum=0;
  for(ibeg=gametes->begin();ibeg!=gametes->end();++ibeg)
    {
      sum+=ibeg->n;
    }
  if(sum != twoN) std::cerr << NRECS << ' ' << sum << ' ' << twoN << '\n';
  assert(sum==twoN);
#endif

  return NRECS;
}
#endif /* _DIPLOID_FUNCTIONS_TCC_ */
