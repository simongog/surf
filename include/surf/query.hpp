
#ifndef SURF_QUERY_HPP
#define SURF_QUERY_HPP

#include <vector>

namespace surf {

struct doc_score {
	uint64_t doc_id;
	double score;
    bool operator>(const doc_score& rhs) const {
    	if(score == rhs.score)
    		return doc_id > rhs.doc_id;
        return score > rhs.score;
    }
    doc_score() {};
    doc_score(uint64_t did,double s) : doc_id(did) , score(s) {};
};

struct result {
    std::vector<doc_score> list;
    uint64_t wt_search_space = 0;
    uint64_t wt_nodes = 0;
    uint64_t postings_evaluated = 0;
    uint64_t postings_total = 0;
};

struct query_token{
    std::vector<uint64_t> token_ids;
    std::vector<std::string> token_strs;
	uint64_t f_qt = 1;
    query_token() {};
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