/*
 * 
 */
#include <vector>
#include <string>
#include <cstdlib>
#include "sitemarker.hpp"

//using namespace std;


SiteMarker::SiteMarker(int id, int num_alleles):site_id(id)
{
  alleles = boost::dynamic_bitset<>( num_alleles ) ;
}

void SiteMarker::zero_all_alleles()
{
  alleles.reset();
}
void SiteMarker::set_allele(int i)
{
  alleles[i]=1;
}

void SiteMarker::set_these_alleles(std::vector<int> v)
{
  for (int i : v)
    {
      alleles[i]=1;
    }
}

int SiteMarker::get_num_alleles()
{
  return alleles.size();
}


SiteMarkerArray::SiteMarkerArray(std::string filename_lin_prg)
{
  std::vector<int> num_alleles_in_site;
  num_alleles_in_site.reserve(1000);
  //open lin PRG file, parse it and count up number of sites and how many alleles each
  //as you go through, for each site do:      num_alleles_in_site[i]=whatever

  //then 
  
}



int main()
{
  return 1;
}
