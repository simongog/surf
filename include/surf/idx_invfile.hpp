
#ifndef SURF_IDX_INVFILE_HPP
#define SURF_IDX_INVFILE_HPP

#include "sdsl/config.hpp"
#include "sdsl/int_vector.hpp"
#include "surf/construct_invidx.hpp"
#include "surf/construct_darray.hpp"
#include "surf/construct_doc_border.hpp"
#include "surf/invfile_postings_list.hpp"
#include "surf/util.hpp"
#include "surf/rank_functions.hpp"

using namespace sdsl;

namespace surf {

template<class t_pl = postings_list<compression_codec::optpfor,128>,class t_rank = rank_bm25<120,75>>
class idx_invfile {
public:
    using size_type = sdsl::int_vector<>::size_type;
    using plist_type = t_pl;
    using rank_type = t_rank;
private:
    std::vector<plist_type> m_postings_lists;
public:
	idx_invfile() = default;
    idx_invfile(cache_config& config)
    {
    	if( cache_file_exists<plist_type>(KEY_INVFILE_PLISTS,config) ) {
    		std::ifstream ifs(cache_file_name<plist_type>(KEY_INVFILE_PLISTS,config));
    		load(ifs);
    	} else {
    		construct_postings_lists<plist_type,rank_type>(m_postings_lists,config);
    		std::ofstream ofs(cache_file_name<plist_type>(KEY_INVFILE_PLISTS,config));
    		serialize(ofs);
    	}
    }
    auto serialize(std::ostream& out, sdsl::structure_tree_node* v=NULL, std::string name="") const -> size_type {
    	structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
        size_type written_bytes = 0;
        size_t num_lists = m_postings_lists.size();
        written_bytes += sdsl::serialize(num_lists,out,child,"num postings lists");
        for(const auto& pl : m_postings_lists) {
        	written_bytes += sdsl::serialize(pl,out,child,"postings list");
        }
        structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    void load(std::istream & in){
    	size_t num_lists;
    	read_member(num_lists,in);
    	m_postings_lists.resize(num_lists);
    	for(size_t i=0;i<num_lists;i++) {
    		m_postings_lists[i].load(in);
    	}
    }
};

template<class t_pl,class t_rank>
void construct(idx_invfile<t_pl,t_rank> &idx, const std::string& file,
               sdsl::cache_config& cconfig, uint8_t num_bytes)
{
    using namespace sdsl;

    cout << "construct(idx_invfile)"<< endl;
    register_cache_file(sdsl::conf::KEY_TEXT_INT, cconfig);

    if (!cache_file_exists(sdsl::conf::KEY_SA, cconfig)) {
        construct_sa<sdsl::int_alphabet_tag::WIDTH>(cconfig);
    }
    register_cache_file(sdsl::conf::KEY_SA, cconfig);
    cout << "sa constructed"<< endl;

    if (!cache_file_exists(surf::KEY_DOCBORDER, cconfig)){
        construct_doc_border<sdsl::int_alphabet_tag::WIDTH>(cconfig);
    }
    cout << "docborders constructed"<< endl;
    if (!cache_file_exists(surf::KEY_DARRAY, cconfig)){
        construct_darray<sdsl::int_alphabet_tag::WIDTH>(cconfig);
    }
    cout << "darray constructed"<< endl;

    cout << "call idx_invfile construct" << endl;
    idx = idx_invfile<t_pl,t_rank>(cconfig);
}

}

#endif
