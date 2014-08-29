#ifndef SURF_DARRAY_HPP
#define SURF_DARRAY_HPP

#include "config.hpp"
#include "construct_doc_perm.hpp"
#include "construct_doc_border.hpp"
#include <sdsl/suffix_arrays.hpp>
#include <algorithm>

namespace surf{

template<uint8_t t_width>
void construct_darray(sdsl::cache_config& cc, bool permute=true)
{
    using namespace sdsl;
    using namespace std;
    if ( !cache_file_exists(KEY_DARRAY, cc) ) {
        bit_vector doc_border;
        construct_doc_border<t_width>(cc);
        load_from_cache(doc_border, KEY_DOCBORDER, cc);
        
        int_vector_buffer<> sa(cache_file_name(conf::KEY_SA, cc));

        rank_support_v<> doc_border_rank(&doc_border);
        uint64_t doc_cnt = doc_border_rank(doc_border.size());

        int_vector<> darray(sa.size(), 0, bits::hi(doc_cnt)+1);
        if ( permute ){
            construct_doc_perm<t_width>(cc);
            doc_perm dp;
            load_from_cache(dp, KEY_DOCPERM,cc);
            for (uint64_t i=0; i<sa.size(); ++i){
                darray[i] = dp.id2len[doc_border_rank(sa[i])];
            }
        } else {
            for (uint64_t i=0; i<sa.size(); ++i){
                darray[i] = doc_border_rank(sa[i]);
            }       
        }
        store_to_cache(darray, KEY_DARRAY, cc);
    }
}

}// end namespace

#endif
