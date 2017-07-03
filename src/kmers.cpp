#include <algorithm>

#include "bwt_search.h"
#include "kmers.hpp"

#define THREADS 25

inline bool fexists (const std::string& name) {
    ifstream f(name.c_str());
    return f.good();
}

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

std::vector<std::string> split(std::string cad,std::string delim)
{
    std::vector<std::string> v;
    int p,q,d;

    q=cad.size();
    p=0;
    d=strlen(delim.c_str());
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

void * worker (void *st)
{
    thread_data* th=(thread_data *) st;
    precalc_kmer_matches(*(th->csa),th->k,*(th->kmer_idx),*(th->kmer_idx_rev),*(th->kmer_sites),
                         *(th->mask_a),th->maxx,*(th->kmers_in_ref),*(th->kmers), th->thread_id);
    return NULL;
}


void gen_precalc_kmers(
        csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,16777216> &csa,
//		sequence_map<std::vector<uint8_t>, std::list<std::pair<uint64_t,uint64_t>>>& kmer_idx,
//		sequence_map<std::vector<uint8_t>, std::list<std::pair<uint64_t,uint64_t>>>& kmer_idx_rev,
//		sequence_map<std::vector<uint8_t>, std::list<std::vector<std::pair<uint32_t, std::vector<int>>>>>& kmer_sites,
//		sequence_set<std::vector<uint8_t>> &kmers_in_ref,
//		std::vector<std::vector<uint8_t>> &kmers,
        std::vector<int> &mask_a,
        std::string kmer_fname,
        uint64_t maxx,
        int k
)
{


    pthread_t threads[THREADS];
    struct thread_data td[THREADS];
    std::vector<std::vector<uint8_t>> kmers[THREADS];

    std::ifstream kfile;
    std::string line;
    kfile.open(kmer_fname);

    int i=0;
    while (std::getline(kfile,line))
    {
        std::vector<uint8_t> kmer;
        for (auto c: line)
            switch (c)
            {
                case 'A': case 'a': kmer.push_back(1);break;
                case 'C': case 'c': kmer.push_back(2);break;
                case 'G': case 'g': kmer.push_back(3);break;
                case 'T': case 't': kmer.push_back(4);break;
            }
        kmers[i++].push_back(kmer);
        i%=THREADS;
    }

    sequence_map<std::vector<uint8_t>, std::list<std::pair<uint64_t,uint64_t>>> kmer_idx[THREADS],kmer_idx_rev[THREADS];
    sequence_map<std::vector<uint8_t>, std::list<std::vector<std::pair<uint32_t, std::vector<int>>>>> kmer_sites[THREADS];
    sequence_set<std::vector<uint8_t>> kmers_in_ref[THREADS];

    for (int i=0;i<THREADS;i++)
    {
        td[i].csa=&csa;
        td[i].k=k;
        td[i].kmer_idx=&kmer_idx[i];
        td[i].kmer_idx_rev=&kmer_idx_rev[i];
        td[i].kmer_sites=&kmer_sites[i];
        td[i].mask_a=&mask_a;
        td[i].maxx=maxx;
        td[i].kmers_in_ref=&kmers_in_ref[i];
        td[i].kmers=&kmers[i];
        td[i].thread_id = i;

        std::cout << "Starting thread: " << i << std::endl;
        pthread_create(&threads[i], NULL, worker, &td[i]);
    }

    std::ofstream precalc_file;
    precalc_file.open (std::string(kmer_fname)+".precalc");

    for (int i=0;i<THREADS;i++){
        void * status;
        pthread_join(threads[i],&status);
        for (auto obj: kmer_idx[i])
        {
            auto k=obj.first;
            for (auto n: k) precalc_file<<(int) n << ' ';
            precalc_file << '|';

            if (kmers_in_ref[i].count(k))
            {
                precalc_file << 1;
            }
            else
            {
                precalc_file << 0;
            }
            precalc_file << '|';
            for (auto o: kmer_idx[i][k]) precalc_file << o.first << ' ' << o.second << ' ';
            precalc_file << '|';
            for (auto o: kmer_idx_rev[i][k]) precalc_file << o.first << ' ' << o.second << ' ';
            precalc_file << '|';
            for (auto o: kmer_sites[i][k])
            {
                for (auto v: o)
                {
                    precalc_file << v.first << ' ';
                    for (auto n: v.second) precalc_file<<(int) n << ' ';
                    precalc_file << '@';
                }
                precalc_file << '|';
            }

            precalc_file << std::endl;
        }
    }

}

void read_precalc_kmers(std::string fil, sequence_map<std::vector<uint8_t>,
        std::list<std::pair<uint64_t,uint64_t>>> &kmer_idx,
                        sequence_map<std::vector<uint8_t>, std::list<std::pair<uint64_t,uint64_t>>> &kmer_idx_rev,
                        sequence_map<std::vector<uint8_t>, std::list<std::vector<std::pair<uint32_t,
                                std::vector<int>>>>> &kmer_sites,
                        sequence_set<std::vector<uint8_t>>&kmers_in_ref
)
{
    std::ifstream kfile;
    std::string line;
    kfile.open(fil);

    while (std::getline(kfile,line))
    {
        std::vector<std::string> parts=split(line,"|");
        std::list<std::pair<uint64_t,uint64_t>> idx,idx_rev;
        std::list<std::vector<std::pair<uint32_t, std::vector<int>>>> sites;

        std::vector<uint8_t> kmer;
        for (auto c: split(parts[0]," ")) kmer.push_back(std::stoi(c));

        if (!parts[1].compare("1"))
            kmers_in_ref.insert(kmer);

        std::vector<std::string> idx_str=split(parts[2]," ");
        std::vector<std::string> idx_rev_str=split(parts[3]," ");
        for (unsigned int i=0;i<idx_str.size()/2;i++) idx.push_back(std::pair<uint64_t,uint64_t>(std::stoi(idx_str[i*2]),std::stoi(idx_str[i*2+1])));
        for (unsigned int i=0;i<idx_rev_str.size()/2;i++) idx_rev.push_back(std::pair<uint64_t,uint64_t>(std::stoi(idx_rev_str[i*2]),std::stoi(idx_rev_str[i*2+1])));

        int flag=0;
        if (! idx.empty())
        {
            kmer_idx[kmer]=idx;
            flag=1;
        }
        if (! idx_rev.empty() )
        {
            kmer_idx_rev[kmer]=idx_rev;
        }
        // 3 1 1 3 1 1 1 1 4 | 1809810 1809950 2244456 2244457 2244471 2244472 |2278258 2278407 3409934 3409934 3410007 3410009 ||9 @9 1 @7 @7 1 @5 @5 1 @|53 4 @51 @|

        if (flag==1)
        {
            for (unsigned int i=4;i<parts.size();i++)
            {
                std::vector<std::pair<uint32_t, std::vector<int>>> v;
                for (auto pair_i_v: split(parts[i],"@")){
                    std::vector<std::string> strvec=split(pair_i_v," ");
                    if (strvec.size())
                    {
                        int first=std::stoi(strvec[0]);

                        std::vector<int> rest;
                        for (unsigned int i=1;i<strvec.size();i++)
                            if (strvec[i].size())
                                rest.push_back(std::stoi(strvec[i]));

                        v.push_back(std::pair<uint32_t, std::vector<int>>(first,rest));
                    }
                }
                sites.push_back(v);
            }
            kmer_sites[kmer]=sites;
        }
    }


//		for (auto items : kmer_idx)
//		{
//			std::vector<uint8_t> k=items.first;
//
//			for (auto n: k) std::cout<<(int) n << ' ';
//			std::cout << '|';
//
//			if (kmers_in_ref.count(k))
//			{
//				std::cout << 1;
//			}
//			else
//			{
//				std::cout << 0;
//			}
//			std::cout << '|';
//			for (auto o: kmer_idx[k]) std::cout << o.first << ' ' << o.second << ' ';
//			std::cout << '|';
//			for (auto o: kmer_idx_rev[k]) std::cout << o.first << ' ' << o.second << ' ';
//			std::cout << '|';
//			for (auto o: kmer_sites[k])
//			{
//				for (auto v: o)
//				{
//					std::cout << v.first << ' ';
//					for (auto n: v.second) std::cout<<(int) n << ' ';
//					std::cout << '@';
//				}
//				std::cout << '|';
//			}
//
//			std::cout << std::endl;
//		}
}


KmersData get_kmers(csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,16777216> &csa,
                   std::vector<int> &mask_a, std::string kmer_fname,
                   uint64_t maxx, int k){

    if (!fexists(std::string(kmer_fname)+".precalc"))
    {
        std::cout << "Precalculated kmers not found, calculating them using "
                  << THREADS << " threads" << std::endl;
        gen_precalc_kmers(csa, mask_a, kmer_fname, maxx, k);
        std::cout << "Finished precalculating kmers" << std::endl;
    }

    std::cout << "Reading K-mers" << std::endl;
    KmersData kmers;
    read_precalc_kmers(std::string(kmer_fname)+".precalc" , kmers.index,
                       kmers.index_reverse, kmers.sites,
                       kmers.in_reference);
    return kmers;
}
