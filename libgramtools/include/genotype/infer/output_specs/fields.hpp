#ifndef GTYPE_FIELDS_HPP
#define GTYPE_FIELDS_HPP

#include <string>
#include <sstream>
#include <map>

namespace gram::genotype::output_spec {

    struct vcf_meta_info_line {
        std::string meta_type, ID, desc,
                flat_value = "", num = "", type = "";

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

        std::string const to_string() {
            std::stringstream out;
            out << "##" << meta_type << "=";
            if (flat_value != "") {
                out << flat_value;
                return out.str();
            }
            out << "<ID=" << ID;
            if (desc != "") out << ",Description=\"" << desc << "\"";
            if (num != "") out << ",Number=" << num;
            if (type != "") out << ",Type=" << type;
            out << ",Source=\"gramtools\"";
            out << ">";
            return out.str();
        }

        char const *const c_str() {
            return this->to_string().c_str();
        }
    };

    using headers = std::map<std::string, vcf_meta_info_line>;

    static headers const vcf_format_headers() {
        headers result;
        result.insert({
              {"FORMAT_GT",
                      vcf_meta_info_line{
                          "FORMAT",
                          "GT",
                          "Genotype",
                          "1",
                          "String"}
              },
              {"FORMAT_DP",
                      vcf_meta_info_line{
                          "FORMAT",
                          "DP",
                          "Total read depth on variant site",
                           "1",
                          "Integer"}
              },
              {"FORMAT_COV",
                      vcf_meta_info_line{
                          "FORMAT",
                          "COV",
                          "Read coverage on each allele",
                          "R",
                          "Integer"}
                          }
        });
        return result;
    }

}

#endif //GTYPE_FIELDS_HPP
