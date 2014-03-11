#ifndef SURF_CONSTRUCT_DOC_LENGTHS_HPP
#define SURF_CONSTRUCT_DOC_LENGTHS_HPP

#include <sdsl/int_vector.hpp>
#include <algorithm>

namespace surf{

template<uint8_t t_width>
void construct_doc_lengths(sdsl::cache_config& cconfig)
{
    using namespace sdsl;
    using namespace std;
    static_assert(t_width == 0 or t_width == 8 , 
            "construct_doc_border: width must be `0` for integer alphabet and `8` for byte alphabet");
    const char* KEY_TEXT  = key_text_trait<t_width>::KEY_TEXT;
    std::string text_file = cache_file_name(KEY_TEXT, cconfig);
    if (!cache_file_exists(KEY_TEXT, cconfig)) {
        std::cerr << "ERROR: construct_doc_cnt: " << text_file
                  << " does not exist. Abort." << std::endl;
        return;
    }
    int_vector_buffer<t_width> text(text_file);
    std::vector<uint64_t> doc_lengths;
    size_t len = 0;
    for (uint64_t i=0; i < text.size(); ++i){
        if ( 1 == text[i] ){
            doc_lengths.push_back(len);
            len = 0;
        } else {
            len++;
        }
    }
    sdsl::int_vector<> sdsl_doc_len(doc_lengths.size());
    for(size_t i=0;i<doc_lengths.size();i++) {
        sdsl_doc_len[i] = doc_lengths[i];
    }
    sdsl::util::bit_compress(sdsl_doc_len);
    store_to_cache(sdsl_doc_len, KEY_DOC_LENGTHS, cconfig);
}

}// end namespace

#endif
