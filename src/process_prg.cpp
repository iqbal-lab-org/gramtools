#include <inttypes.h>
#include <sdsl/suffix_arrays.hpp>

#include "process_prg.hpp"


std::string read_prg_file(std::string prg_fpath){
    std::string prg;
    std::ifstream fhandle(prg_fpath, std::ios::in | std::ios::binary);

    if (!fhandle){
        std::cout << "Problem reading PRG input file" << std::endl;
        exit(1);
    }

    fhandle.seekg(0, std::ios::end);
    prg.resize(fhandle.tellg());
    fhandle.seekg(0, std::ios::beg);
    fhandle.read(&prg[0], prg.size());
    fhandle.close();

    return prg;
}


uint64_t *encode_prg(std::string prg, uint32_t &i, uint64_t &base_index){
    // TODO: needs to be a smart pointer to an array
    uint64_t *prg_int = (uint64_t *) malloc(prg.length() * sizeof(uint64_t));
    if (prg_int == nullptr) {
        exit(1);
    }

    while (i < prg.length()) {
        if (isdigit(prg[i])) {
            int j = 1;
            while (isdigit(prg[i + 1])) {
                j++;
                i++;
            }
            auto al_ind = prg.substr(i - j + 1, j);
            //uint64_t l=(uint64_t) stoull(al_ind,NULL,0);
            auto l = stoull(al_ind, NULL, 0);
            //uint64_t l=boost::lexical_cast<uint64_t>(al_ind);
            prg_int[base_index] = l;
        } else {
            if (prg[i] == 'A' or prg[i] == 'a') prg_int[base_index] = 1;
            if (prg[i] == 'C' or prg[i] == 'c') prg_int[base_index] = 2;
            if (prg[i] == 'G' or prg[i] == 'g') prg_int[base_index] = 3;
            if (prg[i] == 'T' or prg[i] == 't') prg_int[base_index] = 4;
        }
        i++;
        base_index++;// base_index keeps track of actual base position - it's aware of numbers with more than one digit
    }
    return prg_int;
}


//make SA sampling density and ISA sampling density customizable
//make void fcn and pass csa by reference? return ii?
FM_Index construct_fm_index(std::string prg_fpath,
                            std::string prg_int_fpath,
                            std::string memory_log_fname,
                            std::string fm_index_fpath,
                            bool fwd) {

    std::string prg = read_prg_file(prg_fpath);

    uint32_t i = 0;
    uint64_t base_index = 0;
    uint64_t *prg_int = encode_prg(prg, i, base_index);

    uint64_t *prg_int_to_write = prg_int;
    if (!fwd){
        prg_int_fpath = prg_int_fpath + "_rev";
        uint64_t prg_int_rev[base_index];
        std::reverse_copy(prg_int, prg_int + base_index, prg_int_rev);
        prg_int_to_write = prg_int_rev;
    }

    // write prg_int to file
    FILE *fp = fopen(prg_int_fpath.c_str(), "wb");
    fwrite(prg_int_to_write, sizeof(uint64_t), base_index, fp);
    fclose(fp);
    free(prg_int);

    // construct fm_index
    FM_Index fm_index;

    std::ofstream out(memory_log_fname);
    std::cout.rdbuf(out.rdbuf());

    sdsl::memory_monitor::start();
    sdsl::construct(fm_index, prg_int_fpath.c_str(), 8);
    sdsl::memory_monitor::stop();
    sdsl::memory_monitor::write_memory_log<sdsl::HTML_FORMAT>(std::cout);

    std::streambuf *coutbuf = std::cout.rdbuf();
    std::cout.rdbuf(coutbuf);
    sdsl::store_to_file(fm_index, fm_index_fpath);

    return fm_index;
}
