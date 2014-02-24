#ifndef SURF_DARRAY_HPP
#define SURF_DARRAY_HPP

#include "config.hpp"
#include "construct_doc_border.hpp"
#include <sdsl/suffix_arrays.hpp>
#include <algorithm>

namespace surf{

template<uint8_t t_width>
void construct_darray(sdsl::cache_config& cconfig)
{
    using namespace sdsl;
    using namespace std;

    bit_vector doc_border;
    load_from_cache(doc_border, KEY_DOCBORDER, cconfig);
    
    int_vector_buffer<> sa(cache_file_name(conf::KEY_SA, cconfig));

    rank_support_v<> doc_border_rank(&doc_border);
    uint64_t doc_cnt = doc_border_rank(doc_border.size());

    int_vector<> darray(sa.size(), 0, bits::hi(doc_cnt)+1);
    for (uint64_t i=0; i<sa.size(); ++i){
        darray[i] = doc_border_rank(sa[i]);
    }
    store_to_cache(darray, KEY_DARRAY, cconfig);
}

}// end namespace

#endif
