#ifndef SURF_CONSTRUCT_COL_LEN_HPP
#define SURF_CONSTRUCT_COL_LEN_HPP

#include <sdsl/int_vector.hpp>

namespace surf{

template<uint8_t t_width>
void construct_col_len(sdsl::cache_config& cc)
{
    using namespace sdsl;
    using namespace std;
    static_assert(t_width == 0 or t_width == 8 , 
        "construct_col_len: width must be `0` for integer alphabet and `8` for byte alphabet");

    if ( !cache_file_exists(KEY_COLLEN, cc) ){
        const char* KEY_TEXT  = key_text_trait<t_width>::KEY_TEXT;
        std::string text_file = cache_file_name(KEY_TEXT, cc);
        if (!cache_file_exists(KEY_TEXT, cc)) {
            std::cerr << "ERROR: construct_col_len: " << text_file
                      << " does not exist. Abort. t_width=" << (int)t_width << std::endl;
            return;
        }
        uint64_t n = 0;
        std::cerr<<"t_width="<<(int)t_width<<std::endl;
        int_vector_buffer<t_width> text(text_file);
        n = text.size();
        std::cerr<<"construct_col_len = "<<n<<std::endl;
        store_to_cache(n, KEY_COLLEN, cc);
    }
}

}// end namespace

#endif
