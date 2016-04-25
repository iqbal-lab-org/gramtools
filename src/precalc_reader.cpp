#include "bwt_search.h"
#include <iostream> 
#include <fstream> 
#include <vector> 
#include <string.h>
#include <algorithm> 
#include <functional> 
#include <cctype>
#include <locale>

//trim from start
static inline std::string &ltrim(std::string &s) {
	s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
	return s;
}

// trim from end
static inline std::string &rtrim(std::string &s) {
	s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
	return s;
}

// trim from both ends
static inline std::string &trim(std::string &s) {
	return ltrim(rtrim(s));
}

std::vector<std::string> split(std::string cad,char *delim)
{
	std::vector<std::string> v;
	int p,q,d;

	q=cad.size();
	p=0;
	d=strlen(delim);
	while (p<q){
		int posfound=cad.find(delim,p);
		std::string token;

		if (posfound>=0){ token = cad.substr(p, posfound-p);}
		else{ token = cad.substr(p, cad.size()+1);}
		p+=token.size()+d;
		trim(token);
		v.push_back(token);
	}

	
	return v;
}

void readInput(char *fil)
{
	std::ifstream kfile;
	std::string line;
	kfile.open(fil);

	sequence_map<std::vector<uint8_t>, std::list<std::pair<uint64_t,uint64_t>>> kmer_idx;
	sequence_map<std::vector<uint8_t>, std::list<std::pair<uint64_t,uint64_t>>> kmer_idx_rev;
	sequence_map<std::vector<uint8_t>, std::list<std::vector<std::pair<uint32_t, std::vector<int>>>>> kmer_sites;

	while (std::getline(kfile,line))
	{
		std::vector<std::string> parts=split(line,"|");
		std::list<std::pair<uint64_t,uint64_t>> idx,idx_rev;
		std::list<std::vector<std::pair<uint32_t, std::vector<int>>>> sites;

		std::vector<uint8_t> kmer;
		for (auto c: split(parts[0]," ")) kmer.push_back(std::stoi(c));

		std::vector<std::string> idx_str=split(parts[1]," ");
		std::vector<std::string> idx_rev_str=split(parts[2]," ");
		for (int i=0;i<idx_str.size()/2;i++) idx.push_back(std::pair<uint64_t,uint64_t>(std::stoi(idx_str[i*2]),std::stoi(idx_str[i*2+1])));
		for (int i=0;i<idx_rev_str.size()/2;i++) idx_rev.push_back(std::pair<uint64_t,uint64_t>(std::stoi(idx_rev_str[i*2]),std::stoi(idx_rev_str[i*2+1])));

		for (int i=3;i<parts.size();i++)
		{
			std::vector<std::pair<uint32_t, std::vector<int>>> v;
			for (auto pair_i_v: split(parts[i],"@")){
				int first=std::stoi(pair_i_v[0]);
				std::vector<int> rest;
				for (int i=1;i<pair_i_v.size();i++) rest.push_back(std::stoi(pair_i_v[i]));
				v.push_back(std::pair<uint32_t, std::vector<int>>(first,rest));
			}
			sites.push_back(v);
		}

//		for (auto c: parts)
//		{
//			std::cout << c << '~';
//		}
//		std::cout << std::endl;
//		break;
	}
}

int main (int argc, char **argv)
{
	readInput("TESTRESULT");
}
