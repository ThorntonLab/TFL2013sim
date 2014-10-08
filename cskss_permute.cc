//tony's ad-hoc test
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <functional>
#include <cmath>
#include <cstdlib>
#include <numeric>
#include <sstream>
#include <utility>
#include <cassert>
#include <ctest.h>
#include <isbinary.hpp>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <boost/bind.hpp>
#include <Sequence/Portability/random_shuffle.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>

using namespace std;
//using namespace boost;
//using namespace boost::iostreams;

struct FETconfig
{
  unsigned minor_cases,minor_controls;
  double lnpv;
  FETconfig(const unsigned &mica,const unsigned&mico,const double & p):
    minor_cases(mica),minor_controls(mico),lnpv(p)
  {
  }
};

bool operator==(const FETconfig & lhs,
		const FETconfig & rhs)
{
  return (lhs.minor_cases == rhs.minor_cases &&
	  lhs.minor_controls == rhs.minor_controls);
}
bool operator<(const FETconfig & lhs,
		const FETconfig & rhs)
{
  return (lhs.minor_cases < rhs.minor_cases &&
	  lhs.minor_controls < rhs.minor_controls);
}

bool file_exists( const char * infile )
{
  ifstream in(infile);
  if(!in) return false;
  return true;
}

int main( int argc, char ** argv )
{
  int argn=1;
  const char * infile = argv[argn++]; //the case-control genotype file
  if(!file_exists(infile))
    {
      cerr << infile << " does not exist\n";
      exit(10);
    }
  const char * outfile = argv[argn++];//name of output file
  const unsigned ncases=atoi(argv[argn++]);//number of cases (assumes = no. controls)
  const unsigned mincount = atoi(argv[argn++]);//min. count of minor allele (set to 4 for paper)
  const double maxfreq = atof(argv[argn++]);//max MAF (0.05 for paper)
  const unsigned mmarker = atoi(argv[argn++]); //max no. markers to use
  const unsigned seed = atoi(argv[argn++]);
  vector< vector<int> > data;
  if(isbinary(infile))
    {
      boost::iostreams::filtering_istream in;
      in.push(boost::iostreams::gzip_decompressor());
      in.push(boost::iostreams::file_source(infile,ios_base::in|ios_base::binary));  
      string temp,temp2;
      getline(in,temp);
      
      //figure out how many columns there are.
      istringstream figure(temp);
      unsigned ncol=0;
      while( ! figure.eof() )
	{
	  figure >> temp2 >> ws;
	  ++ncol;
	}
     
      //now, read in the data, and store by column
      data = vector<vector<int> >(ncol-4,vector<int>());
      //deal with that first line of input again
      istringstream firstline(temp);
      int foo;

      for(unsigned i=0;i<ncol-4;++i)
	{
	  firstline >> foo >> ws;
	  data[i].push_back(foo);//atoi(temp.c_str()));
	}
      int genotype;
      while(! in.eof() )
	{
	  for(unsigned i=0;i<ncol-4;++i)
	    {
	      in >> genotype >> ws;
	      data[i].push_back(genotype);
	    }
	  if(!in.eof())
	    getline(in,temp); //skip rest of line
	  in >> ws;
	}
    }
  else
    {
      ifstream in;
      string temp,temp2;
      getline(in,temp);
      
      //figure out how many columns there are.
      istringstream figure(temp);
      unsigned ncol=0;
      while( ! figure.eof() )
	{
	  figure >> temp2 >> ws;
	  ++ncol;
	}
     
      //now, read in the data, and store by column
      data = vector<vector<int> >(ncol-4,vector<int>());
      //deal with that first line of input again
      istringstream firstline(temp);
      int(foo);
      for(unsigned i=0;i<ncol-4;++i)
	{
	  firstline >> foo >> ws;
	  data[i].push_back(foo);//atoi(temp.c_str()));
	}
      int genotype;
      while(! in.eof() )
	{
	  for(unsigned i=0;i<ncol-4;++i)
	    {
	      in >> genotype >> ws;
	      data[i].push_back(genotype);
	    }
	  if(!in.eof())
	    getline(in,temp); //skip rest of line
	  in >> ws;
	}
    }

  //stuff for fexact
  int a=2,y=2;
  int work=100000;
  int leading = 2;
  double expect,percnt,emin,prt,pre;

  //get the unique data columns, and do allele freq and mincount checks.
  vector<int> keep(data.size(),1);
  unsigned failcount = 0;
  for(unsigned i=0; i < data.size() ; ++i )
    {
      bool done = false;
      vector<int> tcol(data[i].size());
      transform(data[i].begin(),data[i].end(),
		tcol.begin(),
		bind2nd(plus<int>(),1));
      unsigned sum = accumulate(tcol.begin(),tcol.end(),0);

      if ( sum <= mincount ||
	   double(sum) >= maxfreq * double(4*ncases) )
	//it is 4*ncases because we assume ncases = ncontrols, and diploid!!!
	{
	  keep[i]=0;
	  ++failcount;
	  done = true;
	}
      for( unsigned j=i+1 ; !done&&j < data.size() ; ++j )
	{
	  if( data[i] == data[j] )
	    {
	      keep[i]=0;
	      done = true;
	    }
	}
    }

  vector< vector<int> > reduced_data;
  for( unsigned i=0;i<keep.size();++i )
    {
      if( keep[i] )
	{
	  reduced_data.push_back(data[i]);
	}
    }
  data.clear();


  //now, get the FET for each of these markers in the columns we kept
  vector<double> origFET;
  vector<FETconfig> vFETconfig;
  for( unsigned site=0;site<reduced_data.size() ;++site )
    {
      //vector<int> tcol;
      //add 1 to each genotype to make MAF counting easy
      transform(reduced_data[site].begin(),reduced_data[site].end(),
		reduced_data[site].begin(),
		//std::back_inserter(tcol),
		bind2nd(plus<int>(),1));
      assert(reduced_data[site].size() == 2*ncases);
      //calculate FET
      unsigned minor_control = accumulate(reduced_data[site].begin(),reduced_data[site].begin()+ncases,0);
      unsigned major_control = 6000-minor_control;
      unsigned minor_cases = accumulate(reduced_data[site].begin()+ncases,reduced_data[site].end(),0);
      unsigned major_cases = 6000 - minor_cases;
      vector<FETconfig>::iterator itr = find(vFETconfig.begin(),vFETconfig.end(),FETconfig(minor_cases,minor_control,0.));
      if( itr ==  vFETconfig.end() )
	{
	  double ctable[4];
	  ctable[0] = minor_control;
	  ctable[1] = minor_cases;
	  ctable[2] = major_control;
	  ctable[3] = major_cases;
	  
	  expect=percnt=emin=prt=pre=-1.;
	  
	  fexact(&a,&y,ctable,
		 &leading,&expect,&percnt,&emin,&prt,&pre,&work);
	  origFET.push_back( -log(pre) );//add -log(pval) to array
	  vFETconfig.push_back(FETconfig(minor_cases,minor_control,-log(pre)));
	}
      else
	{
	  origFET.push_back( itr->lnpv );
	}
    }
  //sort in descending order
  sort(origFET.begin(),origFET.end(),
       greater<double>());

  //log of expected pvals
  vector<double> logEpv;
  for(unsigned i=0; i < origFET.size() ; ++i )
    {
      logEpv.push_back( -log( double(i+1)/double(origFET.size()) ) );
    }

  vector<double> kss;

  for( unsigned i = 0 ; i < min(size_t(250),origFET.size()) ; ++i )//keep max of 1st 250
    {
      kss.push_back( (origFET[i]-logEpv[i])/log(10.) );
    }
  //keep track of running sum, which is the output˛
  // ofstream out(outfile);
  // for(unsigned j=0;j<min(size_t(250),kss.size());++j)
  //   {
  //     out << ( accumulate( kss.begin(), kss.begin()+j+1,0. ) ) << '\n';
  //   }
  const double observed = accumulate(kss.begin(),kss.begin()+mmarker+1,0.);

  if( ! isfinite(observed) )
    {
      ofstream out(outfile);
      for(unsigned i=0;i<1001;++i)
	{
	  out << "NA\n";
	}
      out.close();
      exit(1);
    }
  //now, permute
  vector<unsigned> indexes;
  for(unsigned i = 0 ; i < reduced_data[0].size() ; ++i )
    {
      indexes.push_back(i);
    }
  bool permuting = true;
  unsigned permutation = 0;
  unsigned winners = 0;
  //unsigned npermutes[] = {1000};//,10000,100000,250000, 500000, 750000, 1000000};
  //unsigned npermutes[] = {100,1000};
  unsigned PI = 0;
  //const unsigned MAXP = 6;
  vector<double> pdist;
  const unsigned npermutes = 1000;
  gsl_rng * r =  gsl_rng_alloc(gsl_rng_taus2);
  gsl_rng_set(r,seed);
  for( ; permutation < npermutes ; ++permutation )
    {
      vector<double> permuteFET;
      Sequence::random_shuffle(indexes.begin(),indexes.end(),
			       boost::bind(gsl_ran_flat,r,0.,_1));
      for( unsigned site=0;site<reduced_data.size() ;++site )
	{
	      //vector<int> tcol;
	      //add 1 to each genotype to make MAF counting easy
	      assert(reduced_data[site].size() == 2*ncases);
	      //calculate FET
	      unsigned minor_control = 0;//accumulate(reduced_data[site].begin(),reduced_data[site].begin()+ncases,0);
	      unsigned minor_cases = 0;//accumulate(reduced_data[site].begin()+ncases,reduced_data[site].end(),0);
	      for(unsigned j=0;j<reduced_data[site].size();++j)
		{
		  if(j < ncases)
		    {
		      minor_cases += reduced_data[site][indexes[j]];
		    }
		  else
		    {
		      minor_control += reduced_data[site][indexes[j]];
		    }
		}
	      unsigned major_control = 6000-minor_control;
	      unsigned major_cases = 6000 - minor_cases;
	      vector<FETconfig>::iterator itr = find(vFETconfig.begin(),vFETconfig.end(),
						     FETconfig(minor_cases,minor_control,0.));
	      if( itr ==  vFETconfig.end() )
		{
		  double ctable[4];
		  ctable[0] = minor_control;
		  ctable[1] = minor_cases;
		  ctable[2] = major_control;
		  ctable[3] = major_cases;
		  
		  expect=percnt=emin=prt=pre=-1.;
		  fexact(&a,&y,ctable,
			 &leading,&expect,&percnt,&emin,&prt,&pre,&work);
		  permuteFET.push_back( -log(pre) );//add -log(pval) to array
		  vFETconfig.push_back(FETconfig(minor_cases,minor_control,-log(pre)));
		}
	      else
		{
		  permuteFET.push_back(itr->lnpv);
		}
	    }
	  //sort in descending order
	  sort(permuteFET.begin(),permuteFET.end(),
	       greater<double>());
	  vector<double> permlogEpv;
	  for(unsigned i=0; i < permuteFET.size() ; ++i )
	    {
	      permlogEpv.push_back( -log( double(i+1)/double(permuteFET.size()) ) );
	    }
      
	  vector<double> pkss;
      
	  for( unsigned i = 0 ; i < min(size_t(250),permuteFET.size()) ; ++i )//keep max of 1st 250
	    {
	      pkss.push_back( (permuteFET[i]-permlogEpv[i])/log(10.) );
	    }
	  double permstat = accumulate(pkss.begin(),pkss.begin()+mmarker+1,0.);
	  pdist.push_back(permstat);
	  /*
	  if(permstat >= observed)
	    {
	      ++winners;
	    }
	  if( (npermutes[PI] > 10000 && winners > 4) || (npermutes[PI] == npermutes[MAXP] && winners > 1))
	    {
	      //terminate early, as we'll never make it to p <= 1e-6
	      permuting = false;
	    }
	  */
	}
      /*
      if (winners > 4) //out of 1 million reps, if the true p = 1e-6, P(observing > 4) is ~ 3e-3
	{
	  permuting = false;
	}
      else
	{
	  if(PI==MAXP)
	    {
	      permuting = false;
	    }
	  PI++;
	}
      */
      //out << ' ' << double(winners)/double(permutation) << endl;

  ofstream out(outfile);
  out << observed << '\n';
  for(unsigned i=0;i<pdist.size();++i)
    {
      out << pdist[i] << '\n';
    }
  out.close();
}
