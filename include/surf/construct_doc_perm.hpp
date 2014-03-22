#ifndef SURF_CONSTRUCT_DOC_PERM_HPP
#define SURF_CONSTRUCT_DOC_PERM_HPP

#include "doc_perm.hpp"
#include <sdsl/int_vector.hpp>
#include <algorithm>
#include <utility>

namespace surf{

template<uint8_t t_width>
void construct_doc_perm(sdsl::cache_config& cc)
{
    using namespace sdsl;
    using namespace std;
    static_assert(t_width == 0 or t_width == 8 , 
            "construct_doc_perm: width must be `0` for integer alphabet and `8` for byte alphabet");

    if ( !cache_file_exists(KEY_DOCPERM, cc) ) {
        const char* KEY_TEXT  = key_text_trait<t_width>::KEY_TEXT;
        std::string text_file = cache_file_name(KEY_TEXT, cc);
        if (!cache_file_exists(KEY_TEXT, cc)) {
            std::cerr << "ERROR: construct_doc_perm: " << text_file
                      << " does not exist. Abort." << std::endl;
            return;
        }
        int_vector_buffer<t_width> text(text_file);

        std::cout<<"constructing doc_perm start"<<std::endl;
        typedef std::pair<uint64_t, uint64_t> tPII;
        std::vector<tPII> len_id; 
        for (uint64_t i=0, doc_len=0,id=0; i < text.size(); ++i){
            ++doc_len;
            if ( 1 == text[i] ){
                len_id.emplace_back(doc_len, id);
                ++id;
                doc_len = 0;
            }
        }
        std::cout<<"now sorting..."<<std::endl;
        std::sort(len_id.begin(),len_id.end());
        std::cout<<"end sorting"<<std::endl;
        doc_perm dp;
        dp.id2len = int_vector<>(len_id.size(), 0, sdsl::bits::hi(len_id.size()-1)+1);
        dp.len2id = dp.id2len;
        for (size_t i=0; i<len_id.size(); ++i){
            dp.id2len[len_id[i].second] = i;
        }
        std::cout << "inv perm..." << std::endl;
        for (size_t i=0; i<len_id.size(); ++i){
            dp.len2id[dp.id2len[i]] = i;
        }
        std::cout<<"constructing doc_perm end"<<std::endl;
        store_to_cache(dp, KEY_DOCPERM, cc);
    }
}

}// end namespace

#endif
