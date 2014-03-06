#ifndef SURF_CONSTRUCT_DOC_CNT_HPP
#define SURF_CONSTRUCT_DOC_CNT_HPP

#include <sdsl/int_vector.hpp>
#include <algorithm>

namespace surf{

template<uint8_t t_width>
void construct_doc_cnt(sdsl::cache_config& cconfig)
{
    using namespace sdsl;
    using namespace std;
    static_assert(t_width == 0 or t_width == 8 , 
            "construct_doc_cnt: width must be `0` for integer alphabet and `8` for byte alphabet");
    const char* KEY_TEXT  = key_text_trait<t_width>::KEY_TEXT;
    std::string text_file = cache_file_name(KEY_TEXT, cconfig);
    if (!cache_file_exists(KEY_TEXT, cconfig)) {
        std::cerr << "ERROR: construct_doc_cnt: " << text_file
                  << " does not exist. Abort." << std::endl;
        return;
    }
    uint64_t doc_cnt = 0;
    int_vector_buffer<t_width> text(text_file);
    doc_cnt = count_if(text.begin(), text.end(),
                 [](decltype(*(text.begin())) y){
                    return y==1;       
                 });
    store_to_cache(doc_cnt, KEY_DOCCNT, cconfig);
}

}// end namespace

#endif
