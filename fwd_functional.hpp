#ifndef _FWD_FUNCTIONAL_HPP_
#define _FWD_FUNCTIONAL_HPP_

#include <cmath>
#include <limits>

struct mutation_at_pos
{
  typedef bool result_type;
  template<typename mutation_type>
  inline bool operator()(const mutation_type & m, const double & d) const
  {
    return( std::fabs(m.pos-d) <= std::numeric_limits<double>::epsilon() );
  }
};

struct same_pos
{
  typedef bool result_type;
  template<typename mutation_type>
  inline bool operator()(const mutation_type & m1, const mutation_type & m2) const
  {
    return( std::fabs(m1.pos-m2.pos) <= std::numeric_limits<double>::epsilon() );
  }
};
#endif /* _FWD_FUNCTIONAL_HPP_ */
