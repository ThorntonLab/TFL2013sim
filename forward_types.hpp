#ifndef _FORWARD_TYPES_HPP_
#define _FORWARD_TYPES_HPP_

#include <limits>
#include <vector>
#include <list>
#include <cmath>
#include <boost/static_assert.hpp>
#include <boost/type_traits/is_base_and_derived.hpp>


struct mutation_base
/*!
  At minimum, a mutation must contain a position and a count in the population.	
  You can derive from this class, for instance to add selection coefficients,
  counts in different sexes, etc.
*/
{
  mutable double pos;
  unsigned n,g;
  bool neutral,checked;
  mutation_base(const double & position, const unsigned & count, const bool & isneutral = true)
    : pos(position),n(count),neutral(isneutral),checked(false)
  {	
  }
  virtual ~mutation_base(){}
};

struct mutation : public mutation_base
		  //!The simplest mutation type, adding just a selection coefficient and dominance to the interface
{
  mutable double s,h;
  mutation( const double & position, const double & sel_coeff,const unsigned & count,
	    const double dominance = 0.5) 
    : mutation_base(position,count,(sel_coeff==0)),s(sel_coeff),h(dominance)
  {
  }
  bool operator==(const mutation & rhs) const
  {
    return( fabs(this->pos-rhs.pos) <= std::numeric_limits<double>::epsilon() &&
	    this->s == rhs.s );
  }
};

template<typename mut_type,
	 typename list_type = std::list<mut_type> >
struct gamete_base
//! A gamete is a container of pointers (iterators) to mutations + a count in the population
{
  BOOST_STATIC_ASSERT( (boost::is_base_and_derived<mutation_base,mut_type>::value) );
  unsigned n;
  typedef mut_type mutation_type;
  typedef list_type mutation_list_type;
  typedef typename list_type::iterator mutation_list_type_iterator;
  typedef std::vector< mutation_list_type_iterator > mutation_container;
  //mutations is for neutral mutations, smutations for selected...
  mutation_container mutations,smutations;
  gamete_base(const unsigned & icount) : n(icount),mutations( mutation_container() )
  {
  }
  gamete_base(const unsigned & icount, const mutation_container & n,
	      const mutation_container & s) : n(icount),mutations(n),smutations(s)
  {
  }
  virtual ~gamete_base(){}
  
  inline bool operator==(const gamete_base<mut_type,list_type> & rhs) const
  {
    return(this->mutations == rhs.mutations && this->smutations == rhs.smutations);
  }
};

//! The simplest gamete adds nothing to the interface of the base class.
typedef gamete_base<mutation> gamete;

struct fake_less
{
  template<typename iterator_type>
  inline bool operator()( iterator_type  i,iterator_type  j) const
  {
    return (i->pos < j->pos);
  }
};

struct greater_pos
{
  typedef bool result_type;
  template<typename iterator_type>
  inline bool operator()( const iterator_type & i , const iterator_type & j ) const
  {
    return (i->pos > j->pos);
  }
};

#endif /* _FORWARD_TYPES_HPP_ */
