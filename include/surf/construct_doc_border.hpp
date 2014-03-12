#ifndef SURF_CONSTRUCT_DOC_BORDER_HPP
#define SURF_CONSTRUCT_DOC_BORDER_HPP

#include <sdsl/int_vector.hpp>
#include <algorithm>

namespace surf{

template<uint8_t t_width>
void construct_doc_border(sdsl::cache_config& cc)
{
    using namespace sdsl;
    using namespace std;
    static_assert(t_width == 0 or t_width == 8 , 
            "construct_doc_border: width must be `0` for integer alphabet and `8` for byte alphabet");

    if ( !cache_file_exists(KEY_DOCBORDER, cc) ) {
        const char* KEY_TEXT  = key_text_trait<t_width>::KEY_TEXT;
        std::string text_file = cache_file_name(KEY_TEXT, cc);
        if (!cache_file_exists(KEY_TEXT, cc)) {
            std::cerr << "ERROR: construct_doc_cnt: " << text_file
                      << " does not exist. Abort." << std::endl;
            return;
        }
        int_vector_buffer<t_width> text(text_file);
        bit_vector doc_border(text.size(), 0);
        for (uint64_t i=0; i < text.size(); ++i){
            if ( 1 == text[i] ){
                doc_border[i] = 1;
            }
        }
        store_to_cache(doc_border, KEY_DOCBORDER, cc);
    }
}

}// end namespace

#endif
