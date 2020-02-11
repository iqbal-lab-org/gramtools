/**
 * @file Print out fm index information: index, BWT, SA, full text suffixes.
 * Used in testing code in gramtools,
 * usable more generally for illustrative purposes.
 */
#include <iostream>
#include "src_common/generate_prg.hpp"
#include "prg/linearised_prg.hpp"

std::string decode(const uint64_t base){
    switch(base){
        case 1: return "A";
        case 2: return "C";
        case 3: return "G";
        case 4: return "T";
        default: return std::to_string(base);
    }
}

void usage(const char* argv[]){
    std::cout << "Usage: " << argv[0] << " prg_string"
        << " --make_prg PATH" << std::endl;
    std::cout << "The prg string variant markers should be encoded using '[', ']' and ','" << std::endl;
    std::cout << "Use --make_prg to write the produced binary encoded linearised prg." << std::endl;
    exit(1);
}

int main(int argc, const char* argv[]){
    if (argc != 2 && argc != 4){
        usage(argv);
    }

    std::string prg_string = argv[1];
    std::transform(prg_string.begin(), prg_string.end(), prg_string.begin(), ::toupper);

    auto as_marker_vec = prg_string_to_ints(prg_string);

    if (argc == 4){
        std::ofstream fout;
        std::string flag{argv[2]};
        if (flag != "--make_prg") usage(argv);
        std::string fout_path = argv[3];
        PRG_String handler(as_marker_vec);
        handler.write(fout_path);

        std::cout << "Made binary linearised PRG." << std::endl;
        std::cout << "Use:  hexdump -v -e '1/4 \"%d \"' " << fout_path <<
            " to get textual representation." << std::endl;
    }

    gram::PRG_Info prg_info = generate_prg_info(as_marker_vec);

    const auto &fm_index = prg_info.fm_index;

    std::cout << std::endl << "PRG: " << prg_string << std::endl;
    std::cout << "i\tBWT\tSA\ttext_suffix" << std::endl;
    for (int i = 0; i < fm_index.size(); ++i){
        std::cout << i << "\t" << decode(fm_index.bwt[i]) << "\t" << fm_index[i] << "\t";
        for (auto j = fm_index[i]; j < fm_index.size(); ++j)
            // Note: we do not use the prg_string here, because it does not encode each variant marker as its own entity.
            // TODO: find how to use prg_info.encoded_prg
            std::cout << decode(fm_index.text[j]) << " ";
        std::cout << std::endl;
    }
}

