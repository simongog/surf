#ifndef SURF_IDX_D_HPP
#define SURF_IDX_D_HPP

#include "sdsl/suffix_trees.hpp"
#include "surf/df_sada.hpp"
#include "surf/rank_functions.hpp"
#include "surf/construct_col_len.hpp"
#include "surf/query.hpp"
#include <algorithm>
#include <limits>
#include <queue>
#include <sdsl/wt_algorithm.hpp>

namespace surf {

using range_type = sdsl::range_type;

template<typename t_wt_node>
struct s_state_t{
    double score;
    t_wt_node v;
    std::vector<term_info*> t_ptrs; // pointers to term_info array
    std::vector<range_type> r; // ranges

    s_state_t() = default;

    s_state_t(double score, const t_wt_node& v, 
              const std::vector<term_info*>& t_ptrs,
              const std::vector<range_type>& r):
        score(score), v(v), t_ptrs(t_ptrs),
        r(r){}

    s_state_t(s_state_t&&) = default;
    s_state_t(const s_state_t&) = default;

    s_state_t& operator=(s_state_t&&) = default;
    s_state_t& operator=(const s_state_t&) = default;

    bool operator<(const s_state_t& s)const{
        if ( score != s.score ){
            return score < s.score;
        }
        return v < s.v;
    }
};

/*! Class idx_d consists of a 
 *   - CSA over the collection concatenation
 *   - document frequency structure
 *   - a WT over the D array
 */
template<typename t_csa,
         typename t_wtd,
         typename t_df,
         typename t_ranker=rank_bm25<>>
class idx_d{
public:
    using size_type = sdsl::int_vector<>::size_type;
    typedef t_csa    csa_type;
    typedef t_wtd    wtd_type;
    typedef typename wtd_type::node_type node_type;
    typedef t_df     df_type;
    typedef t_ranker ranker_type;
public:
    csa_type    m_csa;
    wtd_type    m_wtd;
    df_type     m_df;
    doc_perm    m_docperm;
    ranker_type m_ranker;

    using state_type = s_state_t<typename t_wtd::node_type>;

private:

    void autocompletee(std::vector<std::vector<uint64_t>> &tokens, std::vector<uint64_t> token_ids, uint64_t& ws) {

        for (uint64_t t : token_ids) {
            std::cout << t << " ";
        }
        std::cout << " -> ";

        size_type occ_begin, occ_end;
        backward_search(m_csa, 0, m_csa.size()-1, token_ids.begin(), token_ids.end(), occ_begin, occ_end);

        size_type k = 10;
        std::vector<uint64_t> cs(m_csa.wavelet_tree.sigma);
        std::vector<size_type> rank_c_i(m_csa.wavelet_tree.sigma);
        std::vector<size_type> rank_c_j(m_csa.wavelet_tree.sigma);

        interval_symbols(m_csa.wavelet_tree, occ_begin, occ_end, k, cs, rank_c_i, rank_c_j);

        for (uint64_t t : cs) {
            std::cout << t << " ";
        }
        std::cout << std::endl;

        for(std::vector<uint64_t>::const_iterator i = cs.begin(); i != cs.end(); ++i)
        {
            const uint64_t c = *i;

            if (c == ws || token_ids.size() > 20)
            {
                if ((std::find(tokens.begin(), tokens.end(), token_ids) == tokens.end()))
                    tokens.push_back(token_ids);

            } else if (c > 1)
            {
                std::vector<uint64_t> token_ids_copy = token_ids;
                token_ids_copy.insert(token_ids_copy.begin(), c);
                autocompletee(tokens, token_ids_copy, ws);
            }
        }
    }

public:

    result search(const std::vector<query_token>& qry,size_t k,bool ranked_and = false,bool profile = false) const {
        typedef std::priority_queue<state_type> pq_type;
        typedef std::priority_queue<double, std::vector<double>, std::greater<double>> pq_min_type;
        std::vector<term_info> terms;
        std::vector<term_info*> term_ptrs;
        std::vector<range_type> ranges;
        result res;

        uint x = 0;
        /*std::cout << "csa: ";
        for (uint i : m_csa) {
            std::cout << x << ":" << i << " ";
            x++;
        }
        std::cout << std::endl;*/

        if (profile) {
            res.wt_nodes = 2*m_wtd.sigma-1;
        }

        for (size_t i=0; i<qry.size(); ++i){
            size_type sp=1, ep=0;
            if (backward_search(m_csa, 0, m_csa.size()-1,
                                 qry[i].token_ids.begin(),
                                 qry[i].token_ids.end(),
                                 sp, ep) > 0) {
                auto f_Dt = std::get<0>(m_df(sp,ep)); // document frequency
                terms.emplace_back(qry[i].token_ids, qry[i].f_qt, sp, ep,  f_Dt);
                ranges.emplace_back(sp, ep);
            }
        }
        term_ptrs.resize(terms.size()); 
        for (size_type i=0; i<terms.size(); ++i){
            term_ptrs[i] = &terms[i];
        }
        double initial_term_num = terms.size();

        auto push_node = [this,&initial_term_num, &res,&profile,&ranked_and]
                         (pq_type& pq, const std::vector<term_info*>& t_ptrs,node_type& v,
                          std::vector<range_type>& r,
                          pq_min_type& pq_min, const size_t& k){
            auto min_idx = m_wtd.sym(v) << (m_wtd.max_level - v.level);  
            auto min_doc_len = m_ranker.doc_length(m_docperm.len2id[min_idx]);
            state_type t; // new state
            t.v = v;

            t.score = initial_term_num * m_ranker.calc_doc_weight(min_doc_len);

            bool eval = false;
            bool is_leaf = m_wtd.is_leaf(v);
            for (size_t i = 0; i < r.size(); ++i){
                if ( !empty(r[i]) ){
                    eval = true;
                    t.r.push_back(r[i]);
                    t.t_ptrs.push_back(t_ptrs[i]);

                    auto score = m_ranker.calculate_docscore(
                                 t.t_ptrs.back()->f_qt,
                                 size(t.r.back()),
                                 t.t_ptrs.back()->f_Dt,
                                 t.t_ptrs.back()->F_Dt(),
                                 min_doc_len,
                                 is_leaf
                               );
                    t.score += score;
                } else if ( ranked_and ){
                    return;
                }
            }
            if (!eval){
                return;
            }
            if ( pq_min.size() < k ){ // not yet k leaves in score queue
                pq.emplace(t);
                if (profile) res.wt_search_space++;
                if ( m_wtd.is_leaf(t.v) )
                    pq_min.push(t.score);
            } else { // more than k leaves in score queue
                if ( t.score > pq_min.top() ){
                    pq.emplace(t);
                    if (profile) res.wt_search_space++;
                    if ( m_wtd.is_leaf(t.v) ){
                        pq_min.pop();
                        pq_min.push(t.score);
                    }
                } 
            }
        };

        constexpr double max_score = std::numeric_limits<double>::max();

        /*for(unsigned int i = 0; i<m_csa.wavelet_tree.size(); i++) {
            cout << i;
            cout << ":";
            cout << m_csa.wavelet_tree[i];
            cout << " ";
        }
        cout << endl;*/

        pq_min_type pq_min;
        pq_type pq;
        pq.emplace(max_score, m_wtd.root(), term_ptrs, ranges);
        if(profile) res.wt_search_space++;

        while ( !pq.empty() and res.list.size() < k ) {
            state_type s = pq.top();
            pq.pop();
            if ( m_wtd.is_leaf(s.v) ) {

                std::vector<prox> query_proximities;
                for (uint i=0; i<s.r.size(); i++) {
                    for (uint j=s.r[i].first; j<=s.r[i].second; j++) {

                        /*std::cout << s.t_ptrs.size() << ": ";
                        for (term_info* ti : s.t_ptrs) {
                            for (uint64_t t : ti->t) {
                                std::cout << t << " ";
                            }
                        }*/

                        uint64_t idx = m_csa[m_wtd.select(j+1, m_wtd.sym(s.v))];
                        term_info* ti = s.t_ptrs[i];
                        /*for (uint64_t t : ti->t) {
                            std::cout << t << " ";
                        }
                        std::cout << std::endl;*/

                        //std:cout << std::endl;
                        //std::cout << "range: " << s.r[i].first << " to: " << s.r[i].second << std::endl;
                        //uint64_t idx = m_csa[m_wtd.select(j+1, m_wtd.sym(s.v))];

                        prox p(idx, *ti);
                        query_proximities.push_back(p);
                    }
                }
                res.list.emplace_back(m_docperm.len2id[m_wtd.sym(s.v)], s.score, query_proximities);

            } else {
//fast_expand:               
                auto exp_v = m_wtd.expand(s.v);
                bool left_empty = m_wtd.empty(std::get<0>(exp_v));
                bool right_empty = m_wtd.empty(std::get<1>(exp_v));
                auto exp_r = m_wtd.expand(s.v, std::move(s.r));
                if ( std::get<1>(exp_r).size() == 0 and std::get<0>(exp_r).size() > 0 and !m_wtd.is_leaf(std::get<0>(exp_v) )){
                    std::cout<<"easy"<<std::endl;
                } 

                if ( !left_empty ) {
                    push_node(pq, s.t_ptrs, std::get<0>(exp_v), std::get<0>(exp_r), pq_min, k);
                } else{
                    //std::cout<<"left_empty"<<std::endl;
                }
                if ( !right_empty ) {
                    push_node(pq, s.t_ptrs, std::get<1>(exp_v), std::get<1>(exp_r), pq_min, k);
                } else{
                    //std::cout<<"right_empty"<<std::endl;
                }
            }
        }
        return res;
    }



    void load(sdsl::cache_config& cc){
        load_from_cache(m_csa, surf::KEY_CSA, cc, true);
        load_from_cache(m_wtd, surf::KEY_WTD, cc, true);
        load_from_cache(m_df, surf::KEY_SADADF, cc, true);
        load_from_cache(m_docperm, surf::KEY_DOCPERM, cc); 
        m_ranker = ranker_type(cc);
    }

    size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const {
        structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
        size_type written_bytes = 0;
        written_bytes += m_csa.serialize(out, child, "CSA");
        written_bytes += m_wtd.serialize(out, child, "WTD");
        written_bytes += m_df.serialize(out, child, "DF");
        written_bytes += m_docperm.serialize(out, child, "DOCPERM");
        structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    void mem_info(){
        std::cout << sdsl::size_in_bytes(m_csa) << ";"; // CSA
        std::cout << sdsl::size_in_bytes(m_wtd) << ";"; // WTD^\ell 
        std::cout << sdsl::size_in_bytes(m_df) << ";";  // DF
        std::cout << 0 << ";"; // WTR^\ell
        std::cout << sdsl::size_in_bytes(m_docperm) << std::endl;  // DOCPERM
    }

    // TODO this to actual search
    std::vector<std::vector<uint64_t>> autocomplete(query_token& qry, uint64_t ws) {
        std::vector<std::vector<uint64_t>> tokens;
        autocompletee(tokens, qry.token_ids, ws);
        return tokens;
    }

};

template<typename t_csa,
         typename t_wtd,
         typename t_df,
         typename t_ranker
        >
void construct(idx_d<t_csa,t_wtd,t_df,t_ranker>& idx,
               const std::string&,
               sdsl::cache_config& cc, uint8_t num_bytes)
{    
    using namespace sdsl;
    using namespace std;

    construct_col_len<t_df::alphabet_category::WIDTH>(cc);

    cout<<"...CSA"<<endl;
    if ( !cache_file_exists<t_csa>(surf::KEY_CSA, cc) )
    {
        t_csa csa;
        construct(csa, "", cc, 0);
        store_to_cache(csa, surf::KEY_CSA, cc, true);
    }
    cout<<"...WTD"<<endl;
    if (!cache_file_exists<t_wtd>(surf::KEY_WTD, cc) ){
        construct_doc_perm<t_csa::alphabet_type::int_width>(cc);
        construct_darray<t_csa::alphabet_type::int_width>(cc);
        t_wtd wtd;
        construct(wtd, cache_file_name(surf::KEY_DARRAY, cc), cc);
        cout << "wtd.size() = " << wtd.size() << endl;
        cout << "wtd.sigma = " << wtd.sigma << endl;
        store_to_cache(wtd, surf::KEY_WTD, cc, true);
    }
    cout<<"...DF"<<endl;
    if (!cache_file_exists<t_df>(surf::KEY_SADADF, cc))
    {
        t_df df;
        construct(df, "", cc, 0);
        store_to_cache(df, surf::KEY_SADADF, cc, true);
    }
}

} // end namespace surf

#endif
