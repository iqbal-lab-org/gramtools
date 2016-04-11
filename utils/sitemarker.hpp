#include <vector>
#include <boost/dynamic_bitset.hpp>


//note the number of alleles is not stored explicitly, can get via the alleles member data
class SiteMarker
{
public:
  int int_alphabet_odd_id; //odd number corresponding to this site
  int site_index;//which-th site is this. Runs from 0 to num_sites-1
  

  //site index is like 0,1,2,3,....
  //odd_id is like 5,7,9,11....
  //odd_id = 2*site_index+5
  SiteMarker(int odd_id, int num_alleles): int_alphabet_odd_id(odd_id), 
					   site_index((odd_id-5)/2), 
					   alleles(num_alleles) {};

  void zero_all_alleles();
  void set_allele(uint32_t i);//set it to 1 if we cross this allele
  void set_these_alleles(std::vector<int> v); //to set 5th allele, do    alleles[5]=1;
  int get_num_alleles(); // this returns alleles.size();
  int get_allele_bit(uint32_t i);

private:
  boost::dynamic_bitset<> alleles;  


};


//this object gets allocated once - parse a file containing one line per alpha
// the sites vector has alphabet+1 entries. at index 0,1,2,3,4 - just have NULL pointers or something
// then you index each site with the number it has in the linPRG. 
class SiteMarkerArray
{
public:
  SiteMarkerArray(std::string sitefile);
  ~SiteMarkerArray();
  
  //suppose you want to get a pointer to the SiteMarker for site 56
  //and at the same time set alleles 1,5,6 to 1 (meaning this read overlaps alleles 1,5,6:
  //note this will zero all other bits.
  SiteMarker* get_site_and_set_allele(int site_id, int allele);
  
  
private:
  std::vector<SiteMarker*> sites;
  int num_sites;//literally the number of sites, nothing to do with int-alphabet/even/odd stuff
};
  

//what is the thing we use when tracking which sites/alleles a read crosses?
class SiteOverlapTracker
{
public:
  SiteOverlapTracker(SiteMarkerArray*);
  void push(int site_id, int allele);
  void clear();
  std::vector<SiteMarker*> vec;
private:

  SiteMarkerArray* sma; //pointer to "PRG" object but not ownership
};
