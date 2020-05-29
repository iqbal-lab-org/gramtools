#ifndef GTYPE_FIELDS_HPP
#define GTYPE_FIELDS_HPP

#include <string>
#include <sstream>
#include <map>
#include <vector>

namespace gram::genotype::output_spec {

    template <typename T>
        struct site_entry {
        std::string meta_type, ID;
        std::vector<T> vals;
        bool single_val;
        };

    /*
     * I cannot figure out a nice way to store a set of arbitrarily typed site_entry objects,
     * so defaulted to explicitly using those that have come up
     */
    struct site_entries {
       std::vector<site_entry<double>> doubles = {};
    };

    struct vcf_meta_info_line {
        std::string meta_type, ID = "", desc = "",
                flat_value = "", num = "", type = "";
        std::size_t length{0};

        /**
         * Constructor for generic header, simple key-val. Eg: ##source=my_source
         */
        vcf_meta_info_line(std::string m_t, std::string flat_val) :
                meta_type(m_t), flat_value(flat_val) {}

        /**
         * Constructor for structured header of type: FORMAT, INFO
         */
        vcf_meta_info_line(std::string m_t, std::string id,
                           std::string desc, std::string num, std::string type) :
                meta_type(m_t), ID(id), desc(desc), num(num), type(type) {}

        /**
         * Constructor for structured header of type: FILTER, ALT
         */
        vcf_meta_info_line(std::string m_t, std::string id, std::string desc) :
                meta_type(m_t), ID(id), desc(desc) {}

        /**
         * Constructor for structured header of type: contig
         */
        vcf_meta_info_line(std::string m_t, std::string id, std::size_t len) :
                meta_type(m_t), ID(id), length(len) {}

        std::string const to_string() const{
            std::stringstream out;
            out << "##" << meta_type << "=";
            if (flat_value != "") {
                out << flat_value;
                return out.str();
            }
            out << "<ID=" << ID;
            if (num != "") out << ",Number=" << num;
            if (type != "") out << ",Type=" << type;
            if (desc != "") out << ",Description=\"" << desc << "\"";
            if (length != 0) out << ",length=" << std::to_string(length);
            out << ",Source=\"gramtools\"";
            out << ">";
            return out.str();
        }
    };

    using header_vec = std::vector<vcf_meta_info_line>;

    static header_vec const vcf_format_headers() {
        header_vec result{
                vcf_meta_info_line{
                        "FORMAT",
                        "GT",
                        "Genotype",
                        "1",
                        "String"},
                vcf_meta_info_line{
                        "FORMAT",
                        "DP",
                        "Total read depth on variant site",
                        "1",
                        "Integer"},
                vcf_meta_info_line{
                        "FORMAT",
                        "COV",
                        "Read coverage on each allele",
                        "R",
                        "Integer"}
        };
    return result;
    }
}

#endif //GTYPE_FIELDS_HPP
