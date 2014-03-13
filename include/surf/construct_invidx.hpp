#ifndef SURF_CONSTRUCT_INVIDX_HPP
#define SURF_CONSTRUCT_INVIDX_HPP

#include "surf/config.hpp"
#include "sdsl/config.hpp"
#include "sdsl/int_vector.hpp"

namespace surf{


void construct_term_ranges(sdsl::int_vector<>& ids, sdsl::int_vector<>& sp, 
                            sdsl::int_vector<>& ep,sdsl::cache_config& cconfig)
{
    sdsl::int_vector_buffer<> sa(cache_file_name(sdsl::conf::KEY_SA,cconfig));
    sdsl::int_vector<> T;
    load_from_cache(T,sdsl::conf::KEY_TEXT_INT,cconfig);
    size_t range_start = 0;
    std::vector<std::tuple<size_t,size_t,size_t>> ranges;
    std::cout << "determine term ranges"<< std::endl;
    for(size_t i=1;i<T.size();i++) {
        if(T[sa[i]] != T[sa[i-1]]) {
            ranges.emplace_back(T[sa[i-1]],range_start,i-1);
            range_start = i;
        }
    }
    ranges.emplace_back(T[sa[T.size()-1]],range_start,T.size()-1);
    sp.resize(ranges.size());
    ep.resize(ranges.size());
    ids.resize(ranges.size());
    size_t num_sym = 0;
    for(const auto& range : ranges) {
        ids[num_sym] = std::get<0>(range);
        sp[num_sym] = std::get<1>(range);
        ep[num_sym++] = std::get<2>(range);
    }
}


template<class t_pl,class t_rank>
void construct_postings_lists(std::vector<t_pl>& postings_lists,sdsl::cache_config& cconfig)
{
    using namespace sdsl;
    using namespace std;

    // load term ranges 
    sdsl::int_vector<> ids; sdsl::int_vector<> sp; sdsl::int_vector<> ep;
    if( cache_file_exists(surf::KEY_INVFILE_TERM_RANGES,cconfig) ) {
        std::ifstream ifs(cache_file_name(surf::KEY_INVFILE_TERM_RANGES,cconfig));
        ids.load(ifs);
        sp.load(ifs);
        ep.load(ifs);
    } else {
        construct_term_ranges(ids,sp,ep,cconfig);
        std::ofstream ofs(cache_file_name(surf::KEY_INVFILE_TERM_RANGES,cconfig));
        serialize(ids,ofs);
        serialize(sp,ofs);
        serialize(ep,ofs);
    }

    // load or construct D array
    std::cout << "stream D"<< std::endl;
    int_vector_buffer<> D(cache_file_name(surf::KEY_DARRAY,cconfig));

    // load or construct rank function
    std::cout << "load rank"<< std::endl;
    t_rank ranker(cconfig);

    // construct plist for each range
    std::cout << "create postings lists"<< endl;
    size_t max_id = ids[ids.size()-1];
    postings_lists.resize(max_id+1);
    for(size_t i=2;i<ids.size();i++) { // skip \0 and \1
        size_t range_size = ep[i] - sp[i] + 1;
        int_vector<> tmpD(range_size);
        std::copy(D.begin()+sp[i],D.begin()+ep[i],tmpD.begin());
        std::cout << "(" << i << ") |<" << sp[i] << "," << ep[i] << ">| = " << range_size << std::endl;
        postings_lists[ids[i]] = t_pl(ranker,T,0,range_size);
    }
}

}// end namespace

#endif
