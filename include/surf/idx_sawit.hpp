#ifndef SURF_IDX_SAWIT_HPP
#define SURF_IDX_SAWIT_HPP

#include "sdsl/suffix_trees.hpp"
#include "surf/df_sada.hpp"
#include "surf/rank_functions.hpp"
#include <algorithm>
#include <limits>
#include <queue>

namespace surf{

using range_type = sdsl::range_type;

template<class ForwardIterator>
std::vector<std::pair<typename ForwardIterator::value_type, uint64_t>> 
unique_and_freq(ForwardIterator first, ForwardIterator last){
    std::sort(first, last);
    std::vector<std::pair<typename ForwardIterator::value_type, uint64_t>> res;
    if ( first == last ){
        return res;
    }
    ForwardIterator result = first;
    res.emplace_back(*first, 1);
    while (++first != last){
        if (!(*result == *first)){
            *(++result)=*first;
            res.emplace_back(*first, 1);
        } else {
            ++(std::get<1>(res.back()));
        }
    }
    return res;
}

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
            return score != s.score;
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
         typename t_df,
         typename t_wtdup>
class idx_sawit{
public:
    using size_type = sdsl::int_vector<>::size_type;
    typedef t_csa   csa_type;
    typedef t_wtd   wtd_type;
    typedef typename wtd_type::node_type node_type;
    typedef t_df    df_type;
    typedef t_wtdup wtdup_type;
    typedef rank_bm25<> ranker_type;
private:
    csa_type    m_csa;
    wtd_type    m_wtd;
    df_type     m_df;
    wtdup_type  m_wtdup; 
    ranker_type m_r;

    using state_type = s_state_t<typename t_wtd::node_type>;
public:
    result_t search(std::vector<uint64_t> qry,size_t k) {
        typedef std::priority_queue<state_type> pq_type;

        auto qry_frq = unique_and_freq(qry.begin(), qry.end());
        std::vector<term_info> terms;
        std::vector<term_info*> term_ptrs;
        std::vector<range_type> ranges;

        for (size_t i=0; i<qry_frq.size(); ++i){
            size_type sp=1, ep=0;
            if ( backward_search(m_csa, 0, m_csa.size()-1, qry_frq[i].first, sp, ep) > 0 ) {
                terms.emplace_back(qry_frq[i].first, qry_frq[i].second, sp, ep, std::get<0>(m_df(sp,ep)) );
                term_ptrs.emplace_back(&terms.back());
                ranges.emplace_back(sp, ep);
                cout << "interval of " << qry[i] << " ["
                     << sp << "," << ep << "]" << endl;
            }
        }

        auto push_node = [this/*,&m_wtd, &m_r*/](pq_type& pq, state_type& s,node_type& v,std::vector<range_type>& r){
            auto min_idx = m_wtd.sym(v) << (m_wtd.max_level - v.level);  
            auto min_doc_len = m_r.doc_length(min_idx);
            state_type t; // new state
            t.v = v;
            t.score = 0;
            for (size_t i = 0; i < r.size(); ++i){
                if ( !empty(r[i]) ){
                    t.r.push_back(r[i]);
                    t.t_ptrs.push_back(s.t_ptrs[i]);
                    t.score += m_r.calculate_docscore(
                                 t.t_ptrs.back()->f_qt,
                                 size(t.r.back()),
                                 t.t_ptrs.back()->f_Dt,
                                 t.t_ptrs.back()->F_Dt(),
                                 min_doc_len
                               );
                } 
            }
            pq.emplace(t);       
        };

        constexpr double max_score = std::numeric_limits<double>::max();
        
        pq_type pq;
        pq.emplace(max_score, m_wtd.root(), term_ptrs, ranges);

        result_t res;
        while ( !pq.empty() and res.size() < k ) {
            state_type s = pq.top();
            pq.pop();
            if ( m_wtd.is_leaf(s.v) ){
                res.emplace_back(s.score, m_wtd.sym(s.v));//TODO: symbol ok?
            } else {
                auto exp_v = m_wtd.expand(s.v);
                auto exp_r = m_wtd.expand(s.v, s.r);
                if ( !m_wtd.empty(std::get<0>(exp_v)) ) {
                    push_node(pq, s, std::get<0>(exp_v), std::get<0>(exp_r));
                }
                if ( !m_wtd.empty(std::get<1>(exp_v)) ) {
                    push_node(pq, s, std::get<1>(exp_v), std::get<1>(exp_r));
                }
            }
        }
        return res;
    }

    void load(sdsl::cache_config& cc){
        load_from_cache(m_csa, surf::KEY_CSA, cc, true);
        load_from_cache(m_wtd, surf::KEY_WTD, cc, true);
        load_from_cache(m_df, surf::KEY_SADADF, cc, true);
        load_from_cache(m_wtdup, surf::KEY_WTDUP, cc, true);
        m_r = ranker_type(cc);
    }

    size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const {
        structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
        size_type written_bytes = 0;
        written_bytes += m_csa.serialize(out, child, "csa");
        written_bytes += m_wtd.serialize(out, child, "wtd");
        written_bytes += m_df.serialize(out, child, "df");
        written_bytes += m_wtdup.serialize(out, child, "wtdup");
        structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

};

template<typename t_csa,
         typename t_wtd,
         typename t_df,
         typename t_wtdup>
void construct(idx_sawit<t_csa,t_wtd,t_df,t_wtdup>& idx,
               const std::string&,
               sdsl::cache_config& cc, uint8_t num_bytes)
{    
    using namespace sdsl;
    using namespace std;
    cout<<"...CSA"<<endl;
    if ( !cache_file_exists<t_csa>(surf::KEY_CSA, cc) )
    {
        t_csa csa;
        construct(csa, "", cc, 0);
        store_to_cache(csa, surf::KEY_CSA, cc, true);
    }
    cout<<"...WTD"<<endl;
    if (!cache_file_exists<t_wtd>(surf::KEY_WTD, cc) ){
        t_wtd wtd;
        construct_darray<t_csa::alphabet_type::int_width>(cc);
        construct(wtd, cache_file_name(surf::KEY_DARRAY, cc));
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
    cout<<"...WTDUP"<<endl;
    if (!cache_file_exists<t_wtdup>(surf::KEY_WTDUP,cc)){
        t_wtdup wtdup;
        construct(wtdup, cache_file_name(surf::KEY_DUP, cc));
        store_to_cache(wtdup, surf::KEY_WTDUP, cc, true);
        cout << "wtdup.size() = " << wtdup.size() << endl;
        cout << "wtdup.sigma = " << wtdup.sigma << endl;
    }
}

} // end namespace surf

#endif
