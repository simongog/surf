
#ifndef SURF_IDX_INVFILE_HPP
#define SURF_IDX_INVFILE_HPP

#include "surf/query.hpp"
#include "sdsl/config.hpp"
#include "sdsl/int_vector.hpp"
#include "surf/construct_invidx.hpp"
#include "surf/construct_darray.hpp"
#include "surf/construct_doc_border.hpp"
#include "construct_doc_cnt.hpp"
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
    // determine lists
    struct plist_wrapper {
        typename plist_type::const_iterator cur;
        typename plist_type::const_iterator end;
        double f_qt;
        double f_t;
        double F_t;
        double list_max_score;
        double max_doc_weight;
        plist_wrapper() = default;
        plist_wrapper(plist_type& pl,double _F_t,double _f_qt) {
            cur = pl.begin();
            end = pl.end();
            list_max_score = pl.list_max_score();
            max_doc_weight = pl.max_doc_weight();
            f_t = pl.size();
            F_t = _F_t;
            f_qt = _f_qt;
        }
    };
private:
    std::vector<plist_type> m_postings_lists;
    sdsl::int_vector<> m_F_t;
    sdsl::int_vector<> m_id_mapping;
    rank_type ranker;
public:
	idx_invfile() = default;
    idx_invfile(cache_config& config)
    {
        if( cache_file_exists(KEY_INVFILE_IDOCPERM,config) ) {
            std::ifstream ifs(cache_file_name(KEY_INVFILE_IDOCPERM,config));
            m_id_mapping.load(ifs);
        } else {
            construct_invidx_doc_permuations(m_id_mapping,config);
            std::ofstream ofs(cache_file_name(KEY_INVFILE_IDOCPERM,config));
            m_id_mapping.serialize(ofs);
        }

        if( cache_file_exists(KEY_F_T,config) ) {
            std::ifstream ifs(cache_file_name(KEY_F_T,config));
            m_F_t.load(ifs);
        } else {
            construct_F_t(m_F_t,config);
            std::ofstream ofs(cache_file_name(KEY_F_T,config));
            m_F_t.serialize(ofs);
        }
    	if( cache_file_exists<plist_type>(KEY_INVFILE_PLISTS,config) ) {
    		std::ifstream ifs(cache_file_name<plist_type>(KEY_INVFILE_PLISTS,config));
            size_t num_lists;
            read_member(num_lists,ifs);
            m_postings_lists.resize(num_lists);
            for(size_t i=0;i<num_lists;i++) {
                m_postings_lists[i].load(ifs);
            }
    	} else {
    		construct_postings_lists<plist_type,rank_type>(m_postings_lists,config);
    		std::ofstream ofs(cache_file_name<plist_type>(KEY_INVFILE_PLISTS,config));
            size_t num_lists = m_postings_lists.size();
            sdsl::serialize(num_lists,ofs);
            for(const auto& pl : m_postings_lists) {
                sdsl::serialize(pl,ofs);
            }
    	}
    }
    auto serialize(std::ostream& out, sdsl::structure_tree_node* v=NULL, std::string name="") const -> size_type {
    	structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
        size_type written_bytes = 0;
        written_bytes += m_F_t.serialize(out,child,"F_t");
        written_bytes += m_id_mapping.serialize(out,child,"id mapping");
        size_t num_lists = m_postings_lists.size();
        written_bytes += sdsl::serialize(num_lists,out,child,"num postings lists");
        for(const auto& pl : m_postings_lists) {
        	written_bytes += sdsl::serialize(pl,out,child,"postings list");
        }
        structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    void load(sdsl::cache_config& cc){
        ranker = t_rank(cc);
    }

    typename std::vector<plist_wrapper*>::iterator
    find_shortest_list(std::vector<plist_wrapper*>& postings_lists,
                       const typename std::vector<plist_wrapper*>::iterator& end,
                       uint64_t id) 
    {
        auto itr = postings_lists.begin();
        if (itr != end) {
            size_t smallest = std::numeric_limits<size_t>::max();
            auto smallest_itr = itr;
            while (itr != end) {
                if ((*itr)->cur.remaining() < smallest && (*itr)->cur.docid() != id) {
                    smallest = (*itr)->cur.remaining();
                    smallest_itr = itr;
                }
                ++itr;
            }
            return smallest_itr;
        }
        return end;
    }

    void sort_list_by_id(std::vector<plist_wrapper*>& plists) {
        // delete if necessary
        auto del_itr = plists.begin();
        while(del_itr != plists.end()) {
            if((*del_itr)->cur == (*del_itr)->end) {
                del_itr = plists.erase(del_itr);
            } else {
                del_itr++;
            }
        }
        // sort
        auto id_sort = [](const plist_wrapper* a,const plist_wrapper* b) {
            return a->cur.docid() < b->cur.docid();
        };
        std::sort(plists.begin(),plists.end(),id_sort);
    }

    void forward_lists(std::vector<plist_wrapper*>& postings_lists,
                       const typename std::vector<plist_wrapper*>::iterator& pivot_list,
                       uint64_t id)
    {
        auto smallest_itr = find_shortest_list(postings_lists,pivot_list+1,id);

        // advance the smallest list to the new id
        (*smallest_itr)->cur.skip_to_id(id);

        if ((*smallest_itr)->cur == (*smallest_itr)->end) {
            // list is finished! reorder list by id
            sort_list_by_id(postings_lists);
            return;
        }

        // bubble it down!
        auto next = smallest_itr + 1;
        auto list_end = postings_lists.end();
        while (next != list_end && (*smallest_itr)->cur.docid() > (*next)->cur.docid()) {
            std::swap(*smallest_itr,*next);
            smallest_itr = next;
            next++;
        }
    }

    std::pair<typename std::vector<plist_wrapper*>::iterator,double>
    determine_candidate(std::vector<plist_wrapper*>& postings_lists,double threshold,size_t initial_lists) {
        double score = 0.0;
        double max_doc_weight = std::numeric_limits<double>::lowest();
        double total_score = 0.0;
        auto itr = postings_lists.begin();
        auto end = postings_lists.end();
        while(itr != end) {
            score += (*itr)->list_max_score;
            max_doc_weight = std::max(max_doc_weight,(*itr)->max_doc_weight);
            total_score = score + (max_doc_weight*initial_lists);
            if(total_score > threshold) {
                // forward to last list equal to pivot
                auto pivot_id = (*itr)->cur.docid();
                auto next = itr+1;
                while(next != end && (*next)->cur.docid() == pivot_id) {
                    itr = next;
                    score += (*itr)->list_max_score;
                    max_doc_weight = std::max(max_doc_weight,(*itr)->max_doc_weight);
                    total_score = score + (max_doc_weight*initial_lists);
                    next++;
                }
                return {itr,score};
            }
            itr++;
        }
        return {end,score};
    }

    double evaluate_pivot(std::vector<plist_wrapper*>& postings_lists,
                        std::priority_queue<doc_score,std::vector<doc_score>,std::greater<doc_score>>& heap,
                        double potential_score,
                        double threshold,
                        size_t initial_lists,
                        size_t k)
    {
        auto doc_id = postings_lists[0]->cur.docid();
        double W_d = ranker.doc_length(doc_id);
        double doc_score = initial_lists * ranker.calc_doc_weight(W_d);
        potential_score -= doc_score;

        auto itr = postings_lists.begin();
        auto end = postings_lists.end();
        while(itr != end) {
            if((*itr)->cur.docid() == doc_id) {
                double contrib = ranker.calculate_docscore((*itr)->f_qt,
                                                         (*itr)->cur.freq(),
                                                         (*itr)->f_t,
                                                         (*itr)->F_t,
                                                         W_d);
                doc_score += contrib;
                potential_score += contrib;
                potential_score -= (*itr)->list_max_score;
                ++((*itr)->cur); // move to next larger doc_id
                if(potential_score < threshold) {
                    break;
                }
            } else {
                break;
            }
            itr++;
        }

        // add if it is in the top-k
        if(heap.size() < k) {
            heap.push({doc_id,doc_score});
        } else {
            if( heap.top().score < doc_score ) {
                heap.pop();
                heap.push({doc_id,doc_score});
            }
        }

        // resort
        sort_list_by_id(postings_lists);

        if(heap.size()) {
            return heap.top().score;
        }
        return 0.0f;
    }

    void
    print_lists(std::vector<plist_wrapper*>& postings_lists,double thres) {
        double lm_sum = 0.0;
        std::cout << thres << " ==================================================================\n";
        for(size_t i=0;i<postings_lists.size();i++) {
            const auto& pl = postings_lists[i];
            std::cout << " (" << i << ") "
                      << " n=" << pl->cur.size() 
                      << " did=" << pl->cur.offset() 
                      << " fdt=" << pl->cur.freq()
                      << " rem=" << pl->cur.remaining()
                      << " lm=" << pl->list_max_score
                      << " sum=" << lm_sum + pl->list_max_score << std::endl;
            lm_sum += pl->list_max_score; 
        }
    }

    result_t process_wand(std::vector<plist_wrapper*>& postings_lists,size_t k) {
        // heap containing the top-k docs
        std::priority_queue<doc_score,std::vector<doc_score>,std::greater<doc_score>> score_heap;

        // init list processing 
        auto threshold = 0.0f;
        size_t initial_lists = postings_lists.size();
        sort_list_by_id(postings_lists);
        auto pivot_and_score = determine_candidate(postings_lists,threshold,initial_lists);
        auto pivot_list = std::get<0>(pivot_and_score);
        auto potential_score = std::get<1>(pivot_and_score);

        while(pivot_list != postings_lists.end()) {
            //print_lists(postings_lists,threshold);
            if (postings_lists[0]->cur.docid() == (*pivot_list)->cur.docid()) {
                threshold = evaluate_pivot(postings_lists,score_heap,potential_score,threshold,initial_lists,k);
            } else {
                forward_lists(postings_lists,pivot_list-1,(*pivot_list)->cur.docid());
            }
            pivot_and_score = determine_candidate(postings_lists,threshold,initial_lists);
            pivot_list = std::get<0>(pivot_and_score);
            potential_score = std::get<1>(pivot_and_score);
        }

        // return the top-k results
        result_t res(score_heap.size());
        for(size_t i=0;i<res.size();i++) {
            auto min = score_heap.top(); score_heap.pop();
            min.doc_id = m_id_mapping[min.doc_id];
            res[res.size()-1-i] = min;
        }

        return res;
    }

    result_t search(const std::vector<query_token>& qry,size_t k) {
        std::vector<plist_wrapper> pl_data(qry.size());
        std::vector<plist_wrapper*> postings_lists;
        size_t j=0;
        for(const auto& qry_token : qry) {
            pl_data[j++] = plist_wrapper(m_postings_lists[qry_token.token_id],(double)m_F_t[qry_token.token_id],(double)qry_token.f_qt);
            if(pl_data[j-1].list_max_score > 0) {
                postings_lists.emplace_back(&(pl_data[j-1]));
            }
        }
        return process_wand(postings_lists,k);
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

    if (!cache_file_exists(surf::KEY_DOCBORDER, cconfig)){
        construct_doc_border<sdsl::int_alphabet_tag::WIDTH>(cconfig);
    }
    if (!cache_file_exists(surf::KEY_DARRAY, cconfig)){
        construct_darray<sdsl::int_alphabet_tag::WIDTH>(cconfig);
    }

    idx = idx_invfile<t_pl,t_rank>(cconfig);
}

}

#endif
