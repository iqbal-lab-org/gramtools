
#include <boost/dynamic_bitset.hpp>



class SiteMarker
{
public:
  int site_id;//number from 5 to whatever
  //no need to store number of alleles - see below

  SiteMarker(int id, int num_alleles);//constructor
  ~SiteMarker();
  zero_all_alleles();
  set_allele(int i);//set it to 1 if we cross this allele
  set_these_alleles(std::vector<int>); //to set 5th allele, do    alleles[5]=1;
  get_num_alleles(); // this returns alleles.size();

private:
  boost::dynamic_bitset<> alleles;  


};


//this object gets allocated once - parse the linPRG file
// the sites vector has alphabet+1 entries. at index 0,1,2,3,4 - just have NULL pointers or something
// then you index each site with the number it has in the linPRG. 
class ArrayOfSiteMarkers
{
public:
  ArrayOfSiteMarkers(std::string filename_lin_prg);//constructor


  //suppose you want to get a pointer to the SiteMarker for site 56
  //and at the same time set alleles 1,5,6 to 1 (meaning this read overlaps alleles 1,5,6:
  SiteMarker* get_site_and_set_alleles(int site_id, std::vector<int> alleles);

  
private:
  std::vector<SiteMarker*> sites;
}


//what is the thing we use when tracking which sites/alleles a read crosses?
class SiteOverlapTracker
{
public:
  SiteOverlapTracker(int num); //constructor - allocs the vector so it has 100 ptrs
  vector<SiteMarker*> vec;
}
//  vector<SiteMarker*> would do fine - an array of pointers, pointing to the fixed
//  SiteMarkers sitting inside ArrayOfSiteMarkers.
// PREALLOC this so it has 100 pointers.
