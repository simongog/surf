
#ifndef SURF_QUERY_HPP
#define SURF_QUERY_HPP

#include <vector>

namespace surf {

struct term_info{
    std::vector<uint64_t> t; // term_id
    uint64_t f_qt; // term_frequency
    uint64_t sp_Dt; // start of interval for term t in the suffix array
    uint64_t ep_Dt; // end of interval for term t in the suffix array
    uint64_t f_Dt;  // number of distinct document the term occurs in

    term_info() = default;
    term_info(const std::vector<uint64_t>& t, uint64_t f_qt, uint64_t sp_Dt, uint64_t ep_Dt, uint64_t f_Dt) :
            t(t), f_qt(f_qt), sp_Dt(sp_Dt), ep_Dt(ep_Dt), f_Dt(f_Dt) {

    }

    term_info(term_info&&) = default;
    term_info(const term_info&) = default;
    term_info& operator=(term_info&&) = default;
    term_info& operator=(const term_info&) = default;

    uint64_t F_Dt() const{
        return ep_Dt-sp_Dt+1;
    }
};

struct prox {
    uint64_t idx;
    term_info ti;

    prox() {};
    prox(uint64_t i, term_info t) : idx(i), ti(t) {};
};

struct doc_score {
	uint64_t doc_id;
	double score;
    std::vector<prox> query_proximities;
    bool operator>(const doc_score& rhs) const {
    	if(score == rhs.score)
    		return doc_id > rhs.doc_id;
        return score > rhs.score;
    }
    doc_score() {};
    doc_score(uint64_t did, double s, std::vector<prox> q_p) : doc_id(did), score(s), query_proximities(q_p) {};
};

struct result {
    std::vector<doc_score> list;
    std::vector<std::string> autocompletes;
    uint64_t wt_search_space = 0;
    uint64_t wt_nodes = 0;
    uint64_t postings_evaluated = 0;
    uint64_t postings_total = 0;
};

struct query_token{
    std::vector<uint64_t> token_ids;
    std::vector<std::string> token_strs;
	uint64_t f_qt;
	query_token(const std::vector<uint64_t>& ids,
                const std::vector<std::string>& strs,
                uint64_t f) : token_ids(ids), token_strs(strs), f_qt(f) 
    {
    }
    bool operator<(const query_token& qt) const {
        return std::lexicographical_compare(token_ids.begin(), token_ids.end(),
                                            qt.token_ids.begin(), qt.token_ids.end());
    }
};

using query_t = std::tuple<uint64_t,std::vector<query_token>>;


}

#endif