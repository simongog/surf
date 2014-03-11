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
    double k1 = (double)t_k1/100.0;
    double b = (double)t_b/100.0;
	size_t num_docs;
	size_t num_terms;
	double avg_doc_len;
	sdsl::int_vector<> doc_lengths;
	rank_bm25(cache_config& cconfig) {
        if (!cache_file_exists(surf::KEY_DOC_LENGTHS, cconfig)){
            surf::construct_doc_lengths<sdsl::int_alphabet_tag::WIDTH>(cconfig);
        }
        load_from_cache(doc_lengths, surf::KEY_DOC_LENGTHS, cconfig);
		num_docs = doc_lengths.size();
	    std::string d_file = cache_file_name(surf::KEY_DARRAY, cconfig);
	    int_vector_buffer<> D(d_file);
	    num_terms = D.size() - num_docs;
	    avg_doc_len = (double)num_terms / (double)num_docs;
	}
	double doc_length(size_t doc_id) const {
		return (double) doc_lengths[doc_id];
	}
	double calc_doc_weight(size_t doc_id) const {
		return 0;
	}
	double calculate_docscore(const double f_qt,const double f_dt,const double f_t,
							  const double F_t,const double W_d) const 
	{
        double w_qt = log((num_docs - f_t + 0.5) / (f_t+0.5)) * f_qt;
        double K_d = k1*((1-b) + (b*(W_d/avg_doc_len)));
        double w_dt = ((k1+1)*f_dt) / (K_d + f_dt);
        return w_dt*w_qt;
    }
};

} // end surf namespace

#endif
