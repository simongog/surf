#ifndef SURF_IDX_SAWIT_HPP
#define SURF_IDX_SAWIT_HPP

#include "sdsl/suffix_trees.hpp"
#include "surf/df_sada.hpp"
#include "surf/rank_functions.hpp"
#include "surf/construct_col_len.hpp"
#include <algorithm>
#include <limits>
#include <queue>

namespace surf{

using range_type = sdsl::range_type;

struct term_info{
    uint64_t t; // term_id
    uint64_t f_qt; // term_frequency
    uint64_t sp_Dt; // start of interval for term t in the suffix array
    uint64_t ep_Dt; // end of interval for term t in the suffix array
    uint64_t f_Dt;  // number of distinct document the term occurs in 

    term_info() = default;
    term_info(uint64_t t, uint64_t f_qt, uint64_t sp_Dt, uint64_t ep_Dt, uint64_t f_Dt) : 
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

/*! Class sawit (Suffix Array Wavelet tree Index Type) consists of a
 *  (compressed) suffix array, a wavelet tree over the document array,
 *  a (succinct) document frequency structure, and a wavelet tree
 *  over the duplication array.
 *  
 */
template<typename t_csa,
         typename t_wtd,
         typename t_df>
class idx_sawit{
public:
    using size_type = sdsl::int_vector<>::size_type;
    typedef t_csa   csa_type;
    typedef t_wtd   wtd_type;
    typedef typename wtd_type::node_type node_type;
    typedef t_df    df_type;
    typedef rank_bm25<> ranker_type;
private:
    csa_type    m_csa;
    wtd_type    m_wtd;
    df_type     m_df;
    doc_perm    m_docperm;
    ranker_type m_r;

    using state_type = s_state_t<typename t_wtd::node_type>;
public:

    result search(const std::vector<query_token>& qry,size_t k,bool ranked_and = false,bool profile = false) {
        typedef std::priority_queue<state_type> pq_type;
        typedef std::priority_queue<double, std::vector<double>, std::greater<double>> pq_min_type;
        std::vector<term_info> terms;
        std::vector<term_info*> term_ptrs;
        std::vector<range_type> ranges;
        result res;

        if(profile) {
            res.wt_nodes = 2*m_wtd.sigma-1;
        }

        for (size_t i=0; i<qry.size(); ++i){
            size_type sp=1, ep=0;
            if ( backward_search(m_csa, 0, m_csa.size()-1, qry[i].token_id, sp, ep) > 0 ) {
                auto f_Dt = std::get<0>(m_df(sp,ep)); // document frequency
                terms.emplace_back(qry[i].token_id, qry[i].f_qt, sp, ep,  f_Dt);
                ranges.emplace_back(sp, ep);
            }
        }
        term_ptrs.resize(terms.size()); 
        for (size_type i=0; i<terms.size(); ++i){
            term_ptrs[i] = &terms[i];
        }

        auto push_node = [this, &res,&profile,&ranked_and]
                         (pq_type& pq, state_type& s,node_type& v,
                          std::vector<range_type>& r,
                          pq_min_type& pq_min, const size_t& k){
            auto min_idx = m_wtd.sym(v) << (m_wtd.max_level - v.level);  
            auto min_doc_len = m_r.doc_length(m_docperm.len2id[min_idx]);
            state_type t; // new state
            t.v = v;
            t.score = 0;
            bool eval = false;
            for (size_t i = 0; i < r.size(); ++i){
                if ( !empty(r[i]) ){
                    eval = true;
                    t.r.push_back(r[i]);
                    t.t_ptrs.push_back(s.t_ptrs[i]);

                    auto score = m_r.calculate_docscore(
                                 t.t_ptrs.back()->f_qt,
                                 size(t.r.back()),
                                 t.t_ptrs.back()->f_Dt,
                                 t.t_ptrs.back()->F_Dt(),
                                 min_doc_len
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
        
        pq_min_type pq_min;
        pq_type pq;
        pq.emplace(max_score, m_wtd.root(), term_ptrs, ranges);
        if(profile) res.wt_search_space++;

        while ( !pq.empty() and res.list.size() < k ) {
            state_type s = pq.top();
            pq.pop();
            if ( m_wtd.is_leaf(s.v) ){
                res.list.emplace_back(m_docperm.len2id[m_wtd.sym(s.v)], s.score);
            } else {
                auto exp_v = m_wtd.expand(s.v);
                auto exp_r = m_wtd.expand(s.v, s.r);
                if ( !m_wtd.empty(std::get<0>(exp_v)) ) {
                    push_node(pq, s, std::get<0>(exp_v), std::get<0>(exp_r), pq_min, k);
                }
                if ( !m_wtd.empty(std::get<1>(exp_v)) ) {
                    push_node(pq, s, std::get<1>(exp_v), std::get<1>(exp_r), pq_min, k);
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
        m_r = ranker_type(cc);
    }

    size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const {
        structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
        size_type written_bytes = 0;
        written_bytes += m_csa.serialize(out, child, "csa");
        written_bytes += m_wtd.serialize(out, child, "wtd");
        written_bytes += m_df.serialize(out, child, "df");
        written_bytes += m_docperm.serialize(out, child, "docperm");
        structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

};

template<typename t_csa,
         typename t_wtd,
         typename t_df
        >
void construct(idx_sawit<t_csa,t_wtd,t_df>& idx,
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
