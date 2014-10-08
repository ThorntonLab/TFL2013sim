#include <Sequence/SimData.hpp>
#include <Sequence/PolyTableFunctions.hpp>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <vector>
#include <iterator>
#include <algorithm>
#include <boost/bind.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <isbinary.hpp>
#include <gsl/gsl_statistics.h>

using namespace std;
using namespace boost::iostreams;
using namespace Sequence;

double within(const double & x ,
	      const double & minimum,
	      const double & maximum)	      
  {
    return( x >= minimum && x <= maximum );
  }

int main(int argc, char ** argv)
{
  int argn=1;
  const char * msfile = argv[argn++];
  const char * phenofile = argv[argn++];
  const char * nmfile = argv[argn++];
  const char * anovafile = argv[argn++];
  const char * msoutfile = argv[argn++];
  const double case_proportion = atof(argv[argn++]); //the proportion of extreme phenotypes called cases
  const unsigned ncases = atoi(argv[argn++]);
  const unsigned ncontrols = atoi(argv[argn++]);

  //read in phenotypes
  vector<double > phenotypes;
  vector<pair<double,double> > raw_phenotypes;
  if(!isbinary(phenofile))
    {
      filtering_istream phenostream;
      phenostream.push(gzip_decompressor());
      phenostream.push(file_source(msfile,ios_base::in | ios_base::binary));
      double g,e;
      while(! phenostream.eof() )
	{
	  phenostream >> g>>e >> ws;
	  phenotypes.push_back(g+e);
	  raw_phenotypes.push_back(make_pair(g,e));
	}
    }
  else
    {
      ifstream phenostream(phenofile);
      if(! phenostream)
	{
	  cerr << "cannot open " << phenofile << endl;
	  exit(10);
	}
      double g,e;
      while(! phenostream.eof() )
	{
	  phenostream >> g>>e >> ws;
	  phenotypes.push_back(g+e);
	  raw_phenotypes.push_back(make_pair(g,e));
	}
    }
  vector<double> sorted_phenotypes(phenotypes.begin(),phenotypes.end());
  sort(sorted_phenotypes.begin(),sorted_phenotypes.end());
  
  const double case_cutoff =  gsl_stats_quantile_from_sorted_data(&sorted_phenotypes[0],
								  1,sorted_phenotypes.size(),1.-case_proportion);
  //exit(10);
  //assemble a vector containing pointers to potential cases
  vector< vector<double>::const_iterator > putative_cases;
  vector<double>::iterator itr2,itr = 
    find_if(phenotypes.begin(),phenotypes.end(),
	    bind2nd(greater_equal<double>(),case_cutoff));
  while( itr != phenotypes.end() )
    {
      putative_cases.push_back(itr);
      itr2=itr+1;
      itr = 
	find_if(itr2,phenotypes.end(),
		bind2nd(greater_equal<double>(),case_cutoff));
    }
  //sometimes, due to ties, etc., the
  //quantile method won't give exactly the desired number of 
  //cases.  So ,we go back to the population, and find 
  //the next most-affected individuals, and make up the difference
  double min_case = case_cutoff;
  while( putative_cases.size() < ncases )
    {
      //find the least-extreme included case in the sorted list
      itr = find_if(sorted_phenotypes.begin(),
		    sorted_phenotypes.end(),
		    bind2nd(greater_equal<double>(),case_cutoff));
      assert(itr != sorted_phenotypes.end());
      //go one step back, and find that individual
      //in the population
      itr--;
      itr2 = find(phenotypes.begin(),
		  phenotypes.end(),
		  *itr);
      if( find(putative_cases.begin(),putative_cases.end(),
	       itr2) == putative_cases.end() )
	{
	  putative_cases.push_back(itr2);
	  min_case = *itr;
	}
    }


  //now, assemble a list of controls, +/- 0.5*sd of phenotype
  const double mean_pheno = gsl_stats_mean(&phenotypes[0],
						1,phenotypes.size());
  const double sd_pheno = gsl_stats_sd(&phenotypes[0],
				       1,phenotypes.size());

//   cerr <<'#'<< (mean_pheno - 0.5*sd_pheno) << ' '
//        << (mean_pheno + 0.5*sd_pheno) << '\n';
  vector< vector<double>::const_iterator > putative_controls;
  itr = find_if(phenotypes.begin(),
		phenotypes.end(),
		boost::bind(within,_1,
			    mean_pheno-0.5*sd_pheno,mean_pheno+0.5*sd_pheno));
  while(itr != phenotypes.end())
    {
      //      cerr << *itr << '\n';
      putative_controls.push_back(itr);
      itr2=itr+1;
      itr = find_if(itr2,
		    phenotypes.end(),
		    boost::bind(within,_1,
				mean_pheno-0.5*sd_pheno,mean_pheno+0.5*sd_pheno));
    }
  SimData neutral,causative;
  if ( isbinary(msfile) )
    {
      filtering_istream msstream;
      msstream.push(gzip_decompressor());
      msstream.push(file_source(msfile,ios_base::in | ios_base::binary));
      msstream >> neutral >> causative >> ws;
      msstream.reset();
    }
  else
    {
      ifstream msstream(msfile);
      if(!msstream)
	{
	  cerr << "cannot open " << msfile << endl;
	  exit(10);
	}
      msstream >> neutral >> causative >> ws;
      msstream.close();
    }




  //now, we can make an output file of cases and controls.
  //because the "msfile" is output in random order with respect to phenotype,
  //we just have to iterate over our putatives lists, and generate the output

  //first, let's make new SimData objects of cases and controls,
  //so that we can easily remove invariant sites

  vector<string> data_neut,data_caus;
  for( unsigned i = 0 ; i < ncontrols ; ++i )
    {
      unsigned a = 2*(putative_controls[i]-phenotypes.begin()), 
	b = a+1; //indexes of i-th control in data matrix
      data_neut.push_back(neutral[a]);
      data_neut.push_back(neutral[b]);
      if (! causative.empty())
	{
	  data_caus.push_back(causative[a]);
	  data_caus.push_back(causative[b]);
	}
    }
  for( unsigned i = 0 ; i < ncases ; ++i )
    {
      unsigned a = 2*(putative_cases[i]-phenotypes.begin()), 
	b = a+1; //indexes of i-th case in data matrix
      data_neut.push_back(neutral[a]);
      data_neut.push_back(neutral[b]);
      if (! causative.empty())
	{
	  data_caus.push_back(causative[a]);
	  data_caus.push_back(causative[b]);
	}
    }
  SimData ccblock_neutral,ccblock_causative;
  ccblock_neutral.assign(&*neutral.pbegin(),
			 neutral.numsites(),
			 &data_neut[0],data_neut.size());
  ccblock_causative.assign(&*causative.pbegin(),
			   causative.numsites(),
			   &data_caus[0],data_caus.size());
  data_neut.clear();
  data_caus.clear();
  RemoveInvariantColumns(&ccblock_neutral);
  RemoveInvariantColumns(&ccblock_causative);
  //ofstream msoutstream(msoutfile);
  filtering_ostream msoutstream;
  msoutstream.push(gzip_compressor());
  msoutstream.push(file_sink(msoutfile,ios_base::out|ios_base::binary));
  msoutstream << ccblock_neutral << '\n' << ccblock_causative << endl;
  //msoutstream.close();
  vector<char> neutral_minors,causative_minors;
  for(SimData::const_site_iterator itr = ccblock_neutral.sbegin();
      itr != ccblock_neutral.send() ; ++itr)
    {
      size_t c = count(itr->second.begin(),itr->second.end(),'1');
      if( c <= ccblock_neutral.size()/2 )
	{
	  neutral_minors.push_back('1');
	}
      else
	{
	  neutral_minors.push_back('0');
	}
    }

  for(SimData::const_site_iterator itr = ccblock_causative.sbegin();
      itr != ccblock_causative.send() ; ++itr)
    {
      size_t c = count(itr->second.begin(),itr->second.end(),'1');
      if( c <= ccblock_causative.size()/2 )
	{
	  causative_minors.push_back('1');
	}
      else
	{
	  cerr << "ancestral allele minor @ causal site!\n";
	  causative_minors.push_back('0');
	}
    }
  //ofstream anovastream(anovafile);
  filtering_ostream anovastream;
  anovastream.push( gzip_compressor() );
  anovastream.push( file_sink(anovafile,ios_base::out|ios_base::binary) );
  if(!anovastream)
    {
      cerr << "error: cannot open " << anovafile << " for writing" << endl;
      exit(10);
    }
  //controls first
  for( unsigned i = 0 ; i < ncontrols ; ++i )
    {
      unsigned a = 2*(putative_controls[i]-phenotypes.begin()), 
	b = a+1; //indexes of i-th control in data matrix

      for(unsigned site=0;site<ccblock_neutral.numsites();++site)
	{
	  char ch1 = ccblock_neutral[2*i][site],
	    ch2=ccblock_neutral[2*i+1][site];
	  if( ch1 != neutral_minors[site] &&
	      ch2 != neutral_minors[site] ) //major/major
	    {
	      anovastream << "-1\t";
	    }
	  else if ( ch1 == neutral_minors[site] &&
		    ch2 == neutral_minors[site] ) //minor/minoc
	    {
	      anovastream << "1\t";
	    }
	  else
	    {
	      anovastream << "0\t"; //het
	    }
	}
      for(unsigned site=0;site<ccblock_causative.numsites();++site)
	{
	  char ch1 = ccblock_causative[2*i][site],
	    ch2=ccblock_causative[2*i+1][site];
	  if( ch1 != causative_minors[site] &&
	      ch2 != causative_minors[site] ) //major/major
	    {
	      anovastream << "-1\t";
	    }
	  else if ( ch1 == causative_minors[site] &&
		    ch2 == causative_minors[site] ) //minor/minoc
	    {
	      anovastream << "1\t";
	    }
	  else
	    {
	      anovastream << "0\t";
	    }
	}
      size_t X1 = (ccblock_causative.empty()) ? 0 : count(ccblock_causative[2*i].begin(),
							  ccblock_causative[2*i].end(),'1');
      size_t X2 = (ccblock_causative.empty()) ? 0 : count(ccblock_causative[2*i+1].begin(),
							  ccblock_causative[2*i+1].end(),'1');
      anovastream << X1 << '\t'
		  << X2 << '\t'
		  << raw_phenotypes[putative_controls[i]-phenotypes.begin()].first << '\t'
		  << raw_phenotypes[putative_controls[i]-phenotypes.begin()].second << '\n';
    }

  //now, cases
  for( unsigned j = 0 ; j < ncases ; ++j )
    {
      unsigned i = j+ncontrols;
      unsigned a = 2*(putative_cases[j]-phenotypes.begin());

      for(unsigned site=0;site<ccblock_neutral.numsites();++site)
	{
	  char ch1 = ccblock_neutral[2*i][site],
	    ch2=ccblock_neutral[2*i+1][site];
	  if( ch1 != neutral_minors[site] &&
	      ch2 != neutral_minors[site] ) //major/major
	    {
	      anovastream << "-1\t";
	    }
	  else if ( ch1 == neutral_minors[site] &&
		    ch2 == neutral_minors[site] ) //minor/minor
	    {
	      anovastream << "1\t";
	    }
	  else
	    {
	      anovastream << "0\t";
	    }
	}
      for(unsigned site=0;!ccblock_causative.empty()&&
	    site<ccblock_causative.numsites();++site)
	{
	  char ch1 = ccblock_causative[2*i][site],
	    ch2=ccblock_causative[2*i+1][site];
	  if( ch1 != causative_minors[site] &&
	      ch2 != causative_minors[site] ) //major/major
	    {
	      anovastream << "-1\t";
	    }
	  else if ( ch1 == causative_minors[site] &&
		    ch2 == causative_minors[site] ) //minor/minoc
	    {
	      anovastream << "1\t";
	    }
	  else
	    {
	      anovastream << "0\t";
	    }
	}
      size_t X1 = (ccblock_causative.empty())?0:count(ccblock_causative[2*i].begin(),
						      ccblock_causative[2*i].end(),'1');
      size_t X2 = (ccblock_causative.empty())?0:count(ccblock_causative[2*i+1].begin(),
						      ccblock_causative[2*i+1].end(),'1');
      anovastream << X1 << '\t'
		  << X2 << '\t'
		  << raw_phenotypes[putative_cases[j]-phenotypes.begin()].first << '\t'
		  << raw_phenotypes[putative_cases[j]-phenotypes.begin()].second << '\n';
    }

  //anovastream.close();
  ofstream  nmstream(nmfile);
  nmstream << ccblock_neutral.numsites() << ' ' << ccblock_causative.numsites() << endl;
  nmstream.close();
}


