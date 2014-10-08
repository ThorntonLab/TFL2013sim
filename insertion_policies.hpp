#ifndef _INSERTION_POLICIES_HPP_
#define _INSERTION_POLICIES_HPP_

#include <algorithm>

struct push_back_gamete
//! An insertion policy for the mutate function.
{
  template< typename gamete_type,
	    typename vector_type_allocator,
	    typename list_type_allocator,
	    template<typename,typename> class vector_type,
	    template<typename,typename> class list_type>
  inline void operator()( const gamete_type & ng, 
			  vector_type<gamete_type,vector_type_allocator > * gametes, 
			  const size_t current_num_gametes,
			  const list_type<typename gamete_type::mutation_type,list_type_allocator> * mutations) const
  {
    gametes->push_back(ng);
  }
};

struct insert_unique
//! An insertion policy
{
  template<typename gamete_type,
	   typename vector_type_allocator,
	   typename list_type_allocator,
	   template<typename,typename> class vector_type,
	   template<typename,typename> class list_type>
  inline void operator()(const gamete_type & ng, 
			 vector_type<gamete_type,vector_type_allocator > * gametes, 
			 const size_t current_num_gametes,
			 const list_type<typename gamete_type::mutation_type,list_type_allocator> * mutations) const
  {
    typedef typename  vector_type<gamete_type,vector_type_allocator >::iterator vtype_iterator;
    vtype_iterator itr=std::find(gametes->begin(),gametes->end(),ng);
    if(itr == gametes->end())
      {
	gametes->push_back(ng);
      }
    else
      {
	itr->n++;
      }
  }
  
  template< typename gamete_type, 
	    typename vector_type_allocator,	    
	    template<typename,typename> class vector_type >	    
  inline void operator()(const gamete_type & ng,
			 vector_type<gamete_type, vector_type_allocator > * gametes) const
  {
    typedef typename vector_type<gamete_type, vector_type_allocator >::iterator vtype_iterator;
    vtype_iterator itr=std::find(gametes->begin(),gametes->end(),ng);
    if(itr == gametes->end())
      {
	gametes->push_back(ng);
      }
    else
      {
	itr->n++;
      }
  }
};

template<typename mutation_type,typename list_type>
inline typename list_type::iterator insert_mutation_at_end(list_type * mutations, const mutation_type & m)
{
  return (mutations->insert(mutations->end(),m));
}

template<typename mutation_type,
	 typename list_type>
inline typename list_type::iterator
insert_unique_mutation(list_type * mutations, const mutation_type & m)
{
  typename list_type::iterator itr = std::find(mutations->begin(),mutations->end(),m);
  if( itr != mutations->end() ) return (itr);
  return (mutations->insert(mutations->end(),m));
}

#endif /* _INSERTION_POLICIES_HPP_ */
