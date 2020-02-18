#ifndef SIMU_PARAMETERS_HPP
#define SIMU_PARAMETERS_HPP

#include "common/parameters.hpp"

namespace gram{

    using Seed = uint32_t;

    class SimulateParams : public CommonParameters {
    public:
        std::string json_out_fpath;
        std::string fasta_out_fpath;
        std::string sample_id;
        uint64_t max_num_paths;
        Seed seed;
    };

    namespace commands::simulate {

        SimulateParams parse_parameters(po::variables_map &vm,
                                        const po::parsed_options &parsed);
    }
}

#endif //SIMU_PARAMETERS_HPP
