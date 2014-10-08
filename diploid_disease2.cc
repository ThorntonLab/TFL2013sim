//same as diploid_disease1.cc, but sigma_s is a command-line parameter
#include <diploid.hh>
//#include <boost/tuple/tuple.hpp>
#include <Sequence/SimData.hpp>
#include <utility>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>

using namespace std;
using namespace boost::iostreams;
using Sequence::SimData;

struct mutation_with_age : public mutation_base
{
  mutable unsigned o;
  double s;
  char label;
  //mutable bool neutral;
  mutation_with_age( const double & position, const double & sel_coeff,
		     const unsigned & count,const unsigned & origin, const char & ch,
		     const bool & n=true) 
    : mutation_base(position,count,n),o(origin),s(sel_coeff),label(ch)
  {
  }
  bool operator==(const mutation_with_age & rhs) const
  {
    return( fabs(this->pos-rhs.pos) <= std::numeric_limits<double>::epsilon() &&
	    this->s == rhs.s );
  }	
};

typedef std::list<mutation_with_age > mlist;
typedef gamete_base<mutation_with_age, mlist> gtype;
typedef vector<gtype> gvector;
typedef vector<mutation_with_age> mvector;


struct disease_effect
{
  typedef double result_type;
  template< typename iterator_type >
  inline std::pair<double,double> operator()(const iterator_type & g1, const iterator_type & g2,
			   const double & sd, gsl_rng * r) const
  {
    //The effect of each allele is additive across mutations
    double e1 = 0.,e2=0.;
    //const double sd2=1.;
    //const double max_fitness = 1/sqrt(2.*M_PI*sd2);
    //cerr << g1->n << ' ' << g2->n << '\n';
    typename  iterator_type::value_type::mutation_container::const_iterator itr;
    for(itr = g1->smutations.begin() ; itr != g1->smutations.end() ; ++itr)
      {
	e1 += (*itr)->s;
      }
    for(itr = g2->smutations.begin() ; itr != g2->smutations.end() ; ++itr)
      {
	e2 += (*itr)->s;
      }
    //the effect size is 1/(1/e1 + 1/e2) = max(0.,e1e2/(e1+e2) + N(0,sd))
    //double effect = (e1>0.||e2>0.) ? 2.*(e1*e2)/(e1 + e2) : 0.;
    double effect = pow( e1*e2, 0.5 );
    //    cerr.precision(6);
    //     if( fabs(effect -  2.*(e1*e2)/(e1 + e2)) >= std::numeric_limits<double>::epsilon() )
    //       std::cerr << "muts: " << e1 << ' ' << e2 << ' ' << effect << ' '<< 2.*(e1*e2)/(e1 + e2) << ' ' <<  pow( e1*e2, 0.5 ) << '\n';
    //effect +=   gsl_ran_gaussian(r,sd);
    return make_pair(effect, gsl_ran_gaussian(r,sd));
  }
};
//calculates the fitess of a diploid
struct disease_effect_to_fitness
{
  typedef double result_type;
  template< typename iterator_type >
  inline double operator()(const iterator_type & g1, const iterator_type & g2,
			   const double & sd, const double & sd_s,gsl_rng * r) const
  {
    //const double sd2=1./3.;
    pair<double,double> effect = disease_effect()(g1,g2,sd,r);
    double fitness = exp( (-1. * pow(effect.first+effect.second,2.))/(2.*pow(sd_s,2)) );
    return ( fitness );
  }
};

struct mutation_model
{
  typedef mutation_with_age result_type;
  inline result_type operator()( gsl_rng * r, const unsigned int & ttl_generations,
				 const double & s, const double & ud, const double & un, const  mlist * mutations,
				 const bool dist_effects = false) const
  {
    double pos = gsl_rng_uniform(r);
    while( std::find_if(mutations->begin(),mutations->end(),boost::bind(mutation_at_pos(),_1,pos)) 
	   != mutations->end() )
      {
	pos = gsl_rng_uniform(r);
      }

    if( gsl_rng_uniform(r) <= ud/(ud+un) )
      {
	if( ! dist_effects )
	  {
	    return mutation_with_age(pos,s,1,ttl_generations,'A',false);
	  }
	else
	  {
	    return mutation_with_age(pos,gsl_ran_exponential(r,s),1,ttl_generations,'A',false);
	  }
      }
    return mutation_with_age(pos,0.,1,ttl_generations,'S',true);
  }
};

//typedef boost::tuples::tuple<unsigned,double,unsigned,unsigned,unsigned> logtype;
/*std::ostream & operator<<(std::ostream & s, const logtype & l)
{
  s << boost::get<0>(l) << ' '
    << boost::get<1>(l) << ' '
    << boost::get<2>(l) << ' '
    << boost::get<3>(l) << ' '
    << boost::get<4>(l);
  return s;
  }*/

#ifndef NDEBUG
bool check_sum(const std::vector<gtype> & gametes, const unsigned & twoN)
{
  unsigned check=0;
  for(std::vector<gtype>::const_iterator i=gametes.begin();i<gametes.end();++i)
    {
      check+=i->n;
    }
  return (check == twoN);
}
#endif

int main(int argc, char ** argv)
{
  int argument=1;
  const unsigned N = atoi(argv[argument++]);
  const double mu_disease = atof(argv[argument++]);
  const double mu_neutral = atof(argv[argument++]);
  const double littler = atof(argv[argument++])/2.; //convert per-diploid into per-gamete
  const double s = atof(argv[argument++]);
  const bool dist_effects = atoi(argv[argument++]);
  const double sd = atof(argv[argument++]);//sigma for "E"
  const double sd_s = atof(argv[argument++]);//sigma for selection
  const unsigned ngens_burnin = atoi(argv[argument++]);
  const unsigned ngens_evolve = atoi(argv[argument++]);
  const unsigned samplesize1 = atoi(argv[argument++]);
  const char * ofn = argv[argument++];
  const char * ofn2 = argv[argument++];
  const char * ofn3 = NULL;
  if( dist_effects )
    {
      ofn3 = argv[argument++];
    }
  int nreps=atoi(argv[argument++]);
  const unsigned seed = atoi(argv[argument++]);

  cerr << '#';
  for(unsigned i=0;i<argc;++i)
    {
      cerr << argv[i] << ' ';
    }
  cerr << endl;
  cerr << "#" << mu_disease << endl;
  /*
    const char * ofn = argv[argument++];
    const char * ofn2 = argv[argument++];
    const char * logfile = argv[argument++];
  */	
  //ofstream logstream(logfile);

  /*
    logstream <<"#theta = "<< 4.*double(N)*mu << '\n'
    << "#rho = " << 8.*double(N)*littler << '\n'
    << "#seed = " << seed << "\n#";
    copy(argv,argv+argc,ostream_iterator<char*>(logstream," "));
    logstream << '\n';
  */
  gsl_rng * r =  gsl_rng_alloc(gsl_rng_taus2);
  gsl_rng_set(r,seed);

//   ofstream msfile(ofn);
//   msfile.close();
//   ofstream effectfile(ofn2);
//   effectfile.close();
//   ofstream mutationfile(ofn3);
//   mutationfile.close();
  while(nreps--)
    {
      //the population begins with 1 gamete with no mutations
      gvector gametes(1,gtype(2*N));
      mlist mutations;
      mvector fixations;
      
      vector<unsigned> fixation_times;
      //vector<logtype> logdata;
      time_t timeval;
      (void)time(&timeval);

      unsigned generation;

      fixations.clear();
      fixation_times.clear();
      unsigned ttl_gen = 0;
      double wbar=1;
      for( generation = 0; generation < ngens_burnin; ++generation,++ttl_gen )
	{
	  if(generation%1000==0.)
	    {
	      cerr << generation << ' '
		   << wbar << ' ' << mutations.size() << ' ' << gametes.size()<< '\n';
	    }
	  //a generation is drift, then mutation, recombination
	  wbar = sample_diploid(r,&gametes,2*N,boost::bind(disease_effect_to_fitness(),_1,_2,sd,sd_s,r));
	  //cerr << wbar << ' ';
	  remove_fixed_lost(&mutations,&fixations,&fixation_times,generation,2*N);
	  assert(check_sum(gametes,2*N));
	  unsigned nmuts= mutate(r,&gametes,&mutations,mu_disease+mu_neutral,boost::bind(mutation_model(),r,ttl_gen,s,mu_disease,mu_neutral,&mutations,dist_effects),
				 push_back_gamete(),boost::bind(insert_mutation_at_end<mutation_with_age,list<mutation_with_age> >,_1,_2));
	  //cerr << nmuts << ' ' << mutations.size() << '\n';
	  assert(check_sum(gametes,2*N));
	  unsigned nrec = recombine(r, &gametes, 2*N, littler, boost::bind(gsl_rng_uniform,r));
	  assert(check_sum(gametes,2*N));


	}
      unsigned nfix = fixations.size();
      for( generation = 0; generation < ngens_evolve; ++generation,++ttl_gen )
	{
	  //a generation is drift, then mutation, recombination
	  double wbar = sample_diploid(r,&gametes,2*N,boost::bind(disease_effect_to_fitness(),_1,_2,sd,sd_s,r));
	  remove_fixed_lost(&mutations,&fixations,&fixation_times,generation,2*N);
	  assert(check_sum(gametes,2*N));
	  unsigned nmuts= mutate(r,&gametes,&mutations,mu_disease+mu_neutral,boost::bind(mutation_model(),r,ttl_gen,s,mu_disease,mu_neutral,&mutations,dist_effects),
				 push_back_gamete(),boost::bind(insert_mutation_at_end<mutation_with_age,list<mutation_with_age> >,_1,_2));
	  assert(check_sum(gametes,2*N));
	  unsigned nrec = recombine(r, &gametes, 2*N, littler, boost::bind(gsl_rng_uniform,r));
	  assert(check_sum(gametes,2*N));
	}
      wbar = sample_diploid(r,&gametes,2*N,boost::bind(disease_effect_to_fitness(),_1,_2,sd,sd_s,r)); //These are now gamete frequencies in the adult generation.
      remove_fixed_lost(&mutations,&fixations,&fixation_times,generation,2*N);


      //so now, make 2N diploids from the population as it stands, and assign a phenotype value to each individual
      vector< pair< vector<gtype>::iterator,vector<gtype>::iterator > > diploids;
      vector< unsigned > gamcounts;
      vector< double > phenotypes;
      unsigned gamsum=0;
      for(unsigned i=0;i<gametes.size();++i)
	{
	  gamcounts.push_back(gametes[i].n);
	  gamsum += gametes[i].n;
	}

      //generate the 2N diploids
      unsigned gam1,gam2,j,k,l;
      ofstream effectfile(ofn2);
      for(unsigned i=0;i<N;)
	{
	  gam1 = unsigned(gsl_ran_flat(r,0,gamsum))+1;
	  gam2 = unsigned(gsl_ran_flat(r,0,gamsum))+1;
	  j=0;
	  l=0;
	  for( ;l < gamcounts.size() ; ++l)
	    {
	      j += gamcounts[l];
	      if( j >= gam1 ) { j=l;break; }
	    }
	  k=0;
	  l=0;
	  for( ;l < gamcounts.size() ; ++l)
	    {
	      k += gamcounts[l];
	      if( k >= gam2 ) { k=l;break; }
	    }
 	  if( ( k==j && gamcounts[k]>1 ) ||
 	      ( k!=j && gamcounts[k]>0 && gamcounts[j]>0 ) )
 	    {
	      diploids.push_back( make_pair( gametes.begin()+j, gametes.begin()+k ) );
	      // 	  phenotypes.push_back( disease_effect()(gametes.begin()+j,gametes.begin()+k,sd,r) );
	      pair<double,double> effect = disease_effect()(gametes.begin()+j,gametes.begin()+k,sd,r);
	      effectfile << effect.first << '\t' << effect.second << '\n';
	      // 	  effectfile << phenotypes[phenotypes.size()-1] << '\n';
	      gamcounts[j]--;
	      gamcounts[k]--;
	      gamsum -= 2;
	      ++i;
	    }
	}
      effectfile.close();
      assert( gamsum == 0 );

      //print out num deleterious mutations vs effect
      vector< pair<double,string> > neutral,selected;
      SimData msneut,mssel;
      for(unsigned i=0;i<diploids.size();++i)
	{
	  //neutral mutations
	  for( unsigned mut = 0 ; mut < diploids[i].first->mutations.size() ; ++mut )
	    {
	      //assert(diploids[i].first->mutations[mut]->n);
	      vector< pair<double,string> >::iterator itr = find_if(neutral.begin(),neutral.end(),std::bind2nd(find_mut_pos(), diploids[i].first->mutations[mut]->pos));
	      if( itr == neutral.end() )
		{
		  neutral.push_back( make_pair(diploids[i].first->mutations[mut]->pos,std::string(2*N,'0')) );
		  neutral[neutral.size()-1].second[2*i] = '1';
		}
	      else
		{
		  itr->second[2*i] = '1';
		}
	    }
	  for( unsigned mut = 0 ; mut < diploids[i].second->mutations.size() ; ++mut )
	    {
	      //assert(diploids[i].second->mutations[mut]->n);
	      vector< pair<double,string> >::iterator itr = find_if(neutral.begin(),neutral.end(),std::bind2nd(find_mut_pos(), diploids[i].second->mutations[mut]->pos));
	      if( itr == neutral.end() )
		{
		  neutral.push_back( make_pair(diploids[i].second->mutations[mut]->pos,std::string(2*N,'0')) );
		  neutral[neutral.size()-1].second[2*i+1] = '1';
		}
	      else
		{
		  itr->second[2*i+1] = '1';
		}
	    }

	  //selected
	  for( unsigned mut = 0 ; mut < diploids[i].first->smutations.size() ; ++mut )
	    {
	      vector< pair<double,string> >::iterator itr = find_if(selected.begin(),selected.end(),std::bind2nd(find_mut_pos(), diploids[i].first->smutations[mut]->pos));
	      if( itr == selected.end() )
		{
		  selected.push_back( make_pair(diploids[i].first->smutations[mut]->pos,std::string(2*N,'0')) );
		  selected[selected.size()-1].second[2*i] = '1';
		}
	      else
		{
		  itr->second[2*i] = '1';
		}
	    }
	  for( unsigned mut = 0 ; mut < diploids[i].second->smutations.size() ; ++mut )
	    {
	      vector< pair<double,string> >::iterator itr = find_if(selected.begin(),selected.end(),std::bind2nd(find_mut_pos(), diploids[i].second->smutations[mut]->pos));
	      if( itr == selected.end() )
		{
		  selected.push_back( make_pair(diploids[i].second->smutations[mut]->pos,std::string(2*N,'0')) );
		  selected[selected.size()-1].second[2*i+1] = '1';
		}
	      else
		{
		  itr->second[2*i+1] = '1';
		}
	    }
	}
      std::sort(neutral.begin(),neutral.end(),sortpos());
      std::sort(selected.begin(),selected.end(),sortpos());
      msneut.assign(neutral.begin(),neutral.end());
      mssel.assign(selected.begin(),selected.end());
      //msfile.open(ofn,ios::app);
      filtering_ostream msfile;
      msfile.push(gzip_compressor());
      msfile.push(file_sink(ofn,ios_base::out|ios_base::binary));
      msfile << msneut << '\n'<<mssel<<'\n';
      //msfile.close();

      //mutationfile.open(ofn3,ios::app);
      filtering_ostream mutationfile;
      mutationfile.push(gzip_compressor());
      mutationfile.push(file_sink(ofn3,ios_base::out|ios_base::binary));
      mutationfile << "//\n";
      //mutationfile.precision(10);
      for(mlist::iterator i = mutations.begin() ; i != mutations.end() ; ++i )
	{
	  if (!i->neutral)
	    {
	      mutationfile << i->pos << '\t' << i->s << '\n';
	    }
	}
      //˛mutationfile.close();
	  
      //cerr << fixations.size()-nfix << ' ' << fixations.size() << ' ' <<nfix<< '\n';
      /*
      bool flag = false;
      unsigned nfixed_old = fixations.size();
      SimData sdata1,sdata2;
      msfile.open(ofn,ios::app);
      pair < vector<pair<double,string> >,
	vector<pair<double,string> > > mslike = ms_sample_separate(r,gametes,samplesize1,N,true);
      if(!mslike.first.empty())
	{
	  sdata1.assign(mslike.first.begin(),mslike.first.end());

	  msfile << sdata1 << '\n';

	}
      else
	{

	  msfile << "//\nsegsites: 0\n";
	}
      if(!mslike.second.empty())
	{
	  sdata2.assign(mslike.second.begin(),mslike.second.end());

	  msfile << sdata2 << '\n';

	}
      else
	{

	  msfile << "//\nsegsites: 0\n";
	}
      msfile.close();
      */
    }
}
