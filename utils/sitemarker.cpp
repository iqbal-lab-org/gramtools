/*
 * 
 */
#include <vector>
#include <string>
#include <cstdlib>
#include "sitemarker.hpp"
#include <algorithm>
#include <fstream>

//using namespace std;



void SiteMarker::zero_all_alleles()
{
  alleles.reset();
}

void SiteMarker::set_allele(uint32_t i)
{
  if (i<alleles.size())
    {
      alleles[i]=1;
    }
}

void SiteMarker::set_these_alleles(std::vector<int> v)
{
  for (uint32_t i : v)
    {
      set_allele(i);
    }
}

int SiteMarker::get_num_alleles()
{
  return alleles.size();
}

int SiteMarker::get_allele_bit(uint32_t i)
{
  return alleles[i];
}

//pass in a file which has one column
//the i-th row = number of alleles in i-th site.
SiteMarkerArray::SiteMarkerArray(std::string sitefile)
{
  //count lines in file so know how much to alloc
  std::ifstream inFile(sitefile); 
  num_sites = std::count(std::istreambuf_iterator<char>(inFile), 
			 std::istreambuf_iterator<char>(), '\n');
  inFile.close();

  //collect info on how many alleles in each site from input file
  std::vector<int> num_alleles_in_each_site;
  num_alleles_in_each_site.reserve(num_sites);

  std::ifstream fs(sitefile);
  int i=0;
  for(std::string line; std::getline(fs, line); )
    {
      num_alleles_in_each_site.push_back(stoi(line));
      printf("Got %d alleles for site %d\n", num_alleles_in_each_site[i], i);
      i++;
    }
  
  //we now know how many sites there are and how many alleles for each


  //alloc enoough memory
  sites.reserve(num_sites);
  i=0;
  for (int num : num_alleles_in_each_site)
    {
      int odd_id = 2*i+5;
      sites[i]=new SiteMarker(odd_id, num);
      i++;
    }
  
}

SiteMarkerArray::~SiteMarkerArray()
{
  for (int i=0; i< num_sites; i++)
    {
      delete sites[i];
    }  
}

SiteMarker* SiteMarkerArray::get_site_and_set_allele(int site_id, int allele)
{
  if (site_id>num_sites-1)
    {
      printf("Calling site id %d which is too big. Expected to be < %d. Also - need to put proper exception handling in\n",
	     site_id, num_sites-1);
      exit(1);
    }
  SiteMarker* s = sites[site_id];
  s->set_allele(allele);
  return s;
}

SiteOverlapTracker::SiteOverlapTracker(SiteMarkerArray* SMA): sma(SMA)
{
  vec.reserve(100);
}

void SiteOverlapTracker::push(int site_id, int allele)
{
  vec.push_back(sma->get_site_and_set_allele(site_id, allele));
}

void SiteOverlapTracker::clear()
{
  vec.clear();
}

int main()
{
  //set up the structure that mirrors the PRG sites:
  SiteMarkerArray* sma = new SiteMarkerArray(std::string("testmarkers"));
  
  SiteOverlapTracker* sot = new SiteOverlapTracker(sma);

  //fake some data
  sot->push(0,1);
  sot->push(1,2);
  sot->push(2,1);
  sot->push(3,4);
  sot->push(4,21);
  

  //now let's take a look.
  printf("Start with site 0 - allele bits are\n");
  int j;
  for (j=0; j<sot->vec[0]->get_num_alleles(); j++)
    {
      printf("%d ", sot->vec[0]->get_allele_bit(j));
    }
  printf("\n");

  printf("Site 1 - allele bits are\n");

  for (j=0; j<sot->vec[1]->get_num_alleles(); j++)
    {
      printf("%d ", sot->vec[1]->get_allele_bit(j));
    }
  printf("\n");


  printf("Site 2 - allele bits are\n");
  for (j=0; j<sot->vec[2]->get_num_alleles(); j++)
    {
      printf("%d ", sot->vec[2]->get_allele_bit(j));
    }
  printf("\n");


  printf("Site 3 - allele bits are\n");
  for (j=0; j<sot->vec[3]->get_num_alleles(); j++)
    {
      printf("%d ", sot->vec[3]->get_allele_bit(j));
    }
  printf("\n");



  printf("Site 4 - allele bits are\n");
  for (j=0; j<sot->vec[4]->get_num_alleles(); j++)
    {
      printf("%d ", sot->vec[4]->get_allele_bit(j));
    }
  printf("\n");






  return 1;
}
