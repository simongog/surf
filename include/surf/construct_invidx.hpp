#ifndef SURF_CONSTRUCT_INVIDX_HPP
#define SURF_CONSTRUCT_INVIDX_HPP

#include "surf/config.hpp"
#include "sdsl/config.hpp"
#include "construct_doc_cnt.hpp"
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

void construct_invidx_doc_permuations(sdsl::int_vector<>& id_mapping,sdsl::cache_config& cconfig)
{
    surf::construct_doc_cnt<sdsl::int_alphabet_tag::WIDTH>(cconfig);
    uint64_t doc_cnt = 0;
    load_from_cache(doc_cnt, surf::KEY_DOCCNT, cconfig);
    sdsl::int_vector<> doc_mapping(doc_cnt);
    {
        auto url_file = cconfig.dir + "/../" + surf::URL2ID_FILENAME;
        std::ifstream ufs(url_file);
        if(ufs.is_open()) {
            /* load current/indri order */
            std::unordered_map<std::string,uint64_t> id_mapping;
            auto docnames_file = cconfig.dir + "/../" + surf::DOCNAMES_FILENAME;
            std::ifstream dfs(docnames_file);
            std::string name_mapping;
            size_t j=0;
            while( std::getline(dfs,name_mapping) ) {
                id_mapping[name_mapping] = j;
                j++;
            }
            /* load url sorted order */
            std::string url_mapping;
            j=0;
            while( std::getline(ufs,url_mapping) ) {
                auto doc_name = url_mapping.substr(url_mapping.find(' ')+1);
                auto itr = id_mapping.find(doc_name);
                if(itr != id_mapping.end()) {
                    doc_mapping[itr->second] = j;
                } else {
                    std::cerr << "could not find mapping for '" << doc_name << "'" << std::endl;
                }
                j++;
            }
        } else {
            // identity permutation
            for(size_t i=0;i<doc_mapping.size();i++) doc_mapping[i] = i;
        }
    }
    // create the inverse permutation
    id_mapping.resize(doc_mapping.size());
    for(size_t i=0;i<doc_mapping.size();i++) {
        id_mapping[doc_mapping[i]] = i;
    }

    // store the forward to disk
    store_to_cache(doc_mapping, KEY_INVFILE_DOCPERM, cconfig);
}

void construct_F_t(sdsl::int_vector<>& F_t,sdsl::cache_config& cconfig)
{
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

    F_t.resize(ids.size());
    for(size_t i=0;i<ids.size();i++) {
        F_t[i] = sp[i] - ep[i] + 1;
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

    // load mapping if it exists
    std::cout << "load docid mapping" << std::endl;
    sdsl::int_vector<> doc_mapping;
    load_from_cache(doc_mapping, KEY_INVFILE_DOCPERM, cconfig);

    // construct plist for each range
    std::cout << "create postings lists"<< endl;
    size_t max_id = ids[ids.size()-1];
    postings_lists.resize(max_id+1);
    for(size_t i=2;i<ids.size();i++) { // skip \0 and \1
        size_t range_size = ep[i] - sp[i] + 1;
        int_vector<> tmpD(range_size);
        for(size_t j=sp[i];j<=ep[i];j++) tmpD[j-sp[i]] = doc_mapping[D[j]];
        if(range_size>1000) std::cout << "(" << i << ") |<" << sp[i] << "," << ep[i] << ">| = " << range_size << std::endl;
        postings_lists[ids[i]] = t_pl(ranker,tmpD,0,range_size-1);
    }
}

}// end namespace

#endif
