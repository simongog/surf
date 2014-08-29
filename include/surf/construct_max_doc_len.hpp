#ifndef SURF_CONSTRUCT_MAX_DOC_LEN_HPP
#define SURF_CONSTRUCT_MAX_DOC_LEN_HPP

#include <sdsl/int_vector.hpp>

namespace surf{

template<uint8_t t_width>
void construct_max_doc_len(sdsl::cache_config& cc)
{
    using namespace sdsl;
    using namespace std;
    static_assert(t_width == 0 or t_width == 8 , 
        "construct_max_doc_len: width must be `0` for integer alphabet and `8` for byte alphabet");

    if ( !cache_file_exists(KEY_MAXDOCLEN, cc) ){
        const char* KEY_TEXT  = key_text_trait<t_width>::KEY_TEXT;
        std::string text_file = cache_file_name(KEY_TEXT, cc);
        if (!cache_file_exists(KEY_TEXT, cc)) {
            std::cerr << "ERROR: construct_max_doc_len: " << text_file
                      << " does not exist. Abort. t_width=" << (int)t_width << std::endl;
            return;
        }
        uint64_t n = 0;
        std::cerr<<"t_width="<<(int)t_width<<std::endl;
        int_vector_buffer<t_width> text(text_file);
        n = text.size();
        uint64_t len = 0, max_len=0;
        for(size_t i=0; i<n;++i){
            ++len;
            if ( text[i]==1 ){
                if ( len > max_len )
                    max_len = len;
                len = 0;
            }
        }
        std::cerr<<"max_doc_len = "<<max_len<<std::endl;
        store_to_cache(max_len, surf::KEY_MAXDOCLEN, cc);
    }
}

}// end namespace

#endif
