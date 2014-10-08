#ifndef _DIPLOID_FUNCTIONS_HPP_
#define _DIPLOID_FUNCTIONS_HPP_

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

template< typename gamete_type,
	  typename vector_type_allocator,
	  template<typename,typename> class vector_type,
	  typename diploid_fitness_function>
double sample_diploid(gsl_rng * r, vector_type<gamete_type,vector_type_allocator > * gametes, 
		      const unsigned & twoN, const diploid_fitness_function & ff,
		      const double & f = 0.);

template< typename gamete_type,
	  typename vector_type_allocator,
	  template<typename,typename> class vector_type,
	  typename diploid_fitness_function>
double sample_diploid_mc(gsl_rng * r, vector_type<gamete_type,vector_type_allocator > * gametes, 
		      const unsigned & twoN, const diploid_fitness_function & ff,
		      const double & f = 0.);


template< typename gamete_type,
	  typename vector_type_allocator,
	  template<typename,typename> class vector_type,
	  typename diploid_fitness_function >
	  double sample_diploid(gsl_rng * r, vector_type<gamete_type,vector_type_allocator > * gametes, 
				const unsigned & twoN_curr, const unsigned & twoN_next,
				const diploid_fitness_function & ff,
				const double & f = 0.);

template< typename gamete_type,
	  typename vector_type_allocator,
	  template<typename,typename> class vector_type,
	  typename diploid_fitness_function >
	  double sample_diploid_mc(gsl_rng * r, vector_type<gamete_type,vector_type_allocator > * gametes, 
				const unsigned & twoN_curr, const unsigned & twoN_next,
				const diploid_fitness_function & ff,
				const double & f = 0.);


template< typename gamete_type,
	  typename vector_type_allocator,
 	  template<typename,typename> class vector_type,
	  typename map_function>
	  unsigned recombine(gsl_rng * r, vector_type<gamete_type,vector_type_allocator > * gametes,
			     const unsigned & twoN, const double & littler,const map_function & mf,
			     const bool single,
			     const double & f = 0);

#endif /* _DIPLOID_FUNCTIONS_HPP_ */
#include <diploid_functions.tcc>
