#ifndef SURF_RANK_FUNCTIONS_HPP
#define SURF_RANK_FUNCTIONS_HPP

#include "construct_doc_lengths.hpp"
#include "surf/config.hpp"
#include <sdsl/suffix_trees.hpp>
#include "sdsl/int_vector.hpp"
#include "surf/util.hpp"

using namespace sdsl;

namespace surf {

template<uint32_t t_k1=120,uint32_t t_b=75>
struct rank_bm25 {
    static const double k1;
    static const double b;
    static const double epsilon_score;
	size_t num_docs;
	size_t num_terms;
	double avg_doc_len;
	double min_doc_len;
	sdsl::int_vector<> doc_lengths;

    rank_bm25(){}

    rank_bm25& operator=(const rank_bm25&) = default;

	rank_bm25(cache_config& cconfig) {
		uint64_t num_terms;
        load_from_cache(num_terms, surf::KEY_COLLEN, cconfig);
        if (!cache_file_exists(surf::KEY_DOC_LENGTHS, cconfig)){
            surf::construct_doc_lengths<sdsl::int_alphabet_tag::WIDTH>(cconfig);
        }
        load_from_cache(doc_lengths, surf::KEY_DOC_LENGTHS, cconfig);
		num_docs = doc_lengths.size();
        std::cerr<<"num_docs = "<<num_docs<<std::endl;
	    avg_doc_len = (double)num_terms / (double)num_docs;
        std::cerr<<"avg_doc_len = "<<avg_doc_len<<std::endl;
	}
	double doc_length(size_t doc_id) const {
		return (double) doc_lengths[doc_id];
	}
	double calc_doc_weight(size_t doc_id) const {
		return 0;
	}
	double calculate_docscore(const double f_qt,const double f_dt,const double f_t,
							  const double F_t,const double W_d,bool) const 
	{
        double w_qt = std::max(epsilon_score, log((num_docs - f_t + 0.5) / (f_t+0.5)) * f_qt);
        double K_d = k1*((1-b) + (b*(W_d/avg_doc_len)));
        double w_dt = ((k1+1)*f_dt) / (K_d + f_dt);
        return w_dt*w_qt;
    }
};


template<uint32_t t_k1=120,uint32_t t_b=75>
struct rank_bm25_simple_est {
    static const double k1;
    static const double b;
    static const double epsilon_score;
	size_t num_docs;
	size_t num_terms;
	double avg_doc_len;
	double min_doc_len;
	sdsl::int_vector<> doc_lengths;

    rank_bm25_simple_est(){}

    rank_bm25_simple_est& operator=(const rank_bm25_simple_est&) = default;

	rank_bm25_simple_est(cache_config& cconfig) {
		uint64_t num_terms;
        load_from_cache(num_terms, surf::KEY_COLLEN, cconfig);
        if (!cache_file_exists(surf::KEY_DOC_LENGTHS, cconfig)){
            surf::construct_doc_lengths<sdsl::int_alphabet_tag::WIDTH>(cconfig);
        }
        load_from_cache(doc_lengths, surf::KEY_DOC_LENGTHS, cconfig);
		num_docs = doc_lengths.size();
        std::cerr<<"num_docs = "<<num_docs<<std::endl;
	    avg_doc_len = (double)num_terms / (double)num_docs;
	    auto min_itr = std::min_element(doc_lengths.begin(),doc_lengths.end());
	    min_doc_len = *min_itr;
        std::cerr<<"avg_doc_len = "<<avg_doc_len<<std::endl;
        std::cerr<<"min_doc_len = "<<min_doc_len<<std::endl;
	}
	double doc_length(size_t doc_id) const {
		return (double) doc_lengths[doc_id];
	}
	double calc_doc_weight(size_t doc_id) const {
		return 0;
	}
	double calculate_docscore(const double f_qt,const double f_dt,const double f_t,
							  const double F_t,double W_d,bool use_W_d = true) const 
	{
		if(!use_W_d) {
			W_d = min_doc_len;
		}
        double w_qt = std::max(epsilon_score, log((num_docs - f_t + 0.5) / (f_t+0.5)) * f_qt);
        double K_d = k1*((1-b) + (b*(W_d/avg_doc_len)));
        double w_dt = ((k1+1)*f_dt) / (K_d + f_dt);
        return w_dt*w_qt;
    }
};


template<uint32_t t_k1,uint32_t t_b>
const double rank_bm25_simple_est<t_k1,t_b>::k1 = (double)t_k1/100.0;

template<uint32_t t_k1,uint32_t t_b>
const double rank_bm25_simple_est<t_k1,t_b>::b = (double)t_b/100.0;

template<uint32_t t_k1,uint32_t t_b>
const double rank_bm25_simple_est<t_k1,t_b>::epsilon_score = 1e-6;

template<uint32_t t_k1,uint32_t t_b>
const double rank_bm25<t_k1,t_b>::k1 = (double)t_k1/100.0;

template<uint32_t t_k1,uint32_t t_b>
const double rank_bm25<t_k1,t_b>::b = (double)t_b/100.0;

template<uint32_t t_k1,uint32_t t_b>
const double rank_bm25<t_k1,t_b>::epsilon_score = 1e-6;


} // end surf namespace

#endif
