#ifndef MAKE_VCF_HPP
#define MAKE_VCF_HPP

#include <htslib/vcf.h>
#include "genotype/infer/interfaces.hpp"
#include "json_prg_spec.hpp"

struct vcf_meta_info_line{
   std::string meta_type, ID,
   desc = "", num = "", type = "";

   vcf_meta_info_line(std::string m_t, std::string id,
           std::string desc, std::string num, std::string type) :
   meta_type(m_t), ID(id), desc(desc), num(num), type(type){}

   vcf_meta_info_line(std::string m_t, std::string id, std::string desc) :
        meta_type(m_t), ID(id), desc(desc){}

   std::string const to_string(){
       std::stringstream out;
       out << "##" << meta_type << "=<ID=" << ID;
       if (desc != "") out << ",Description=\"" << desc << "\"";
       if (num != "") out << ",Number=" << num;
       if (type != "") out << ",Type=" << type;
       out << "Source=\"gramtools\"";
       out << ">";
       return out.str();
   }
};

using headers = std::map<std::string, vcf_meta_info_line>;

static headers const get_format_header_spec(){
   headers result;
   result.insert({
                         {"FORMAT_GT",
                                 vcf_meta_info_line{"FORMAT", "GT", "Genotype", "1", "String"}
                         },
                         {"FORMAT_DP",
                          vcf_meta_info_line{"FORMAT", "DP", "Total read depth on variant site",
                                             "1", "Integer"}
                          },
                         {"FORMAT_COV",
                          vcf_meta_info_line{"FORMAT", "COV", "Read coverage on each allele",
                                             "R", "Integer"}},
                 });
    return result;
}

void populate_vcf_hdr(bcf_hdr_t * hdr, gram::genotype::infer::Genotyper gtyper);

#endif //MAKE_VCF_HPP
