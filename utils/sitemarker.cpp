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

void SiteMarker::print_all_info()
{
  printf("Marker, site id %d, odd-id %d:\n", site_index, int_alphabet_odd_id);
  printf("Bits are: \n");
  uint32_t i;
  for (i=0; i<alleles.size(); i++)
    {
      printf("%d", int(alleles[i]));
    }
  printf("\n");
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
  printf("Found %d sites\n", num_sites);
  //collect info on how many alleles in each site from input file
  std::vector<int> num_alleles_in_each_site;
  num_alleles_in_each_site.reserve(num_sites);

  std::ifstream fs(sitefile);
  int i=0;
  for(std::string line; std::getline(fs, line); )
    {
      num_alleles_in_each_site.push_back(stoi(line));
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

int SiteMarkerArray::get_num_sites()
{
  return num_sites;
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
  // compile thus:
  // g++ -g -O0 -std=c++11 -I /data2/apps/boost_1_60_0/ -Wall -o zam sitemarker.cpp



  //set up the structure that mirrors the PRG sites:
  SiteMarkerArray* sma = new SiteMarkerArray(std::string("test/testmarkers"));

  SiteOverlapTracker* tracker = new SiteOverlapTracker(sma);

  //fake some data
  tracker->push(0,1);
  tracker->push(1,2);
  tracker->push(2,1);
  tracker->push(3,4);
  tracker->push(4,21);
  

  //now let's take a look.
  
  for (SiteMarker* v : tracker->vec)
    {
      v->print_all_info();
    }
  
  /* this is what you see 

  Marker, site id 0, odd-id 5:
  Bits are: 
  01
  Marker, site id 1, odd-id 7:
  Bits are: 
  001
  Marker, site id 2, odd-id 9:
  Bits are: 
  0100
  Marker, site id 3, odd-id 11:
  Bits are: 
  00001
  Marker, site id 4, odd-id 13:
  Bits are: 
  000000000000000000000100000000

  */




  return 1;
}
