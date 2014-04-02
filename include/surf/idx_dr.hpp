#ifndef SURF_IDX_DR_HPP
#define SURF_IDX_DR_HPP

#include "sdsl/suffix_trees.hpp"
#include "surf/df_sada.hpp"
#include "surf/rank_functions.hpp"
#include "surf/idx_d.hpp"
#include "surf/construct_col_len.hpp"
#include "surf/construct_DUP2.hpp"
#include <algorithm>
#include <limits>
#include <queue>

namespace surf{

using range_type = sdsl::range_type;


template<typename t_wtd_node, typename t_wtr_node>
struct s_state2_t{
    double score;
    t_wtd_node v; // node in document array wavelet tree
    std::vector<term_info*> t_ptrs; // pointers to term_info array
    std::vector<range_type> r_v; // ranges in v
    t_wtr_node w; // node in repetition array wavelet tree
    std::vector<range_type> r_w; // ranges in w

    s_state2_t() = default;

    s_state2_t(double score, const t_wtd_node& v, 
              const std::vector<term_info*>& t_ptrs,
              const std::vector<range_type>& r_v,
              const t_wtr_node& w,
              const std::vector<range_type>& r_w):
        score(score), v(v), t_ptrs(t_ptrs),
        r_v(r_v),w(w),r_w(r_w)
    {}

    s_state2_t(s_state2_t&&) = default;
    s_state2_t(const s_state2_t&) = default;

    s_state2_t& operator=(s_state2_t&&) = default;
    s_state2_t& operator=(const s_state2_t&) = default;

    bool operator<(const s_state2_t& s)const{
        if ( score != s.score ){
            return score < s.score;
        }
        return v < s.v;
    }
};

/*! Class idx_dr consists of a 
 *   - CSA over the collection concatenation
 *   - document frequency structure
 *   - a WT over the D array
 *   - a WT over the repetition array
 */
template<typename t_csa,
         typename t_df,
         typename t_wtd,
         typename t_wtr,
         typename t_ranker=rank_bm25<>,
         typename t_rbv=sdsl::rrr_vector<63>,
         typename t_rrank=typename t_rbv::rank_1_type>
class idx_dr{
public:
    using size_type = sdsl::int_vector<>::size_type;
    typedef t_csa                        csa_type;
    typedef t_wtd                        wtd_type;
    typedef typename wtd_type::node_type node_type;
    typedef t_df                         df_type;
    typedef t_wtr                        wtr_type;
    typedef typename wtr_type::node_type node2_type;
    typedef t_rbv                        rbv_type;
    typedef t_rrank                      rrank_type;
    typedef t_ranker                  ranker_type;
public:
    csa_type    m_csa;
    df_type     m_df;
    wtr_type    m_wtr; 
    t_wtd       m_wtd;
    doc_perm    m_docperm;
    ranker_type m_ranker;

    using state_type = s_state2_t<node_type, node2_type>;
public:

    result search(const std::vector<query_token>& qry,size_t k,bool ranked_and = false,bool profile = false) {
        typedef std::priority_queue<state_type> pq_type;
        std::vector<term_info> terms;
        std::vector<term_info*> term_ptrs;
        std::vector<range_type> v_ranges; // ranges in wtd
        std::vector<range_type> w_ranges; // ranges in wtdup
        result res;

        for (size_t i=0; i<qry.size(); ++i){
            size_type sp=1, ep=0;
            if ( backward_search(m_csa, 0, m_csa.size()-1, 
                                qry[i].token_ids.begin(),
                                qry[i].token_ids.end(),
                                sp, ep) > 0 ) {
                auto df_info = m_df(sp,ep);
                auto f_Dt = std::get<0>(df_info); // document frequency
                terms.emplace_back(qry[i].token_ids, qry[i].f_qt, sp, ep,  f_Dt);
                v_ranges.emplace_back(sp, ep);
                w_ranges.emplace_back(std::get<1>(df_info),std::get<2>(df_info));
            }
        }
        term_ptrs.resize(terms.size()); 
        for (size_type i=0; i<terms.size(); ++i){
            term_ptrs[i] = &terms[i];
        }

        auto push_node = [this,&res,&profile,&ranked_and](pq_type& pq, state_type& s,node_type& v,std::vector<range_type>& r_v,
                                node2_type& w, std::vector<range_type>& r_w){
            auto min_idx = m_wtd.sym(v) << (m_wtd.max_level - v.level);  
            auto min_doc_len = m_ranker.doc_length(m_docperm.len2id[min_idx]);
            state_type t; // new state
            t.v = v;
            t.w = w;
            t.score = 0;
            bool eval = false;
            for (size_t i = 0; i < r_v.size(); ++i){
                if ( !empty(r_v[i]) ){
                    eval = true;
                    t.r_v.push_back(r_v[i]);
                    t.r_w.push_back(r_w[i]);
                    t.t_ptrs.push_back(s.t_ptrs[i]);
                    auto score = m_ranker.calculate_docscore(
                                 t.t_ptrs.back()->f_qt,
                                 size(t.r_w.back())+1,
                                 t.t_ptrs.back()->f_Dt,
                                 t.t_ptrs.back()->F_Dt(),
                                 min_doc_len
                               );
                    t.score += score;
                } else if ( ranked_and ) {
                    return;
                }
            }
            if (eval){ 
                if (profile) res.wt_search_space++;
                pq.emplace(t);       
            }
        };

        constexpr double max_score = std::numeric_limits<double>::max();
        
        pq_type pq;
        size_type search_space=0;
        pq.emplace(max_score, m_wtd.root(), term_ptrs, v_ranges, m_wtr.root(), w_ranges);
        if(profile) res.wt_search_space++;

        while ( !pq.empty() and res.list.size() < k ) {
            state_type s = pq.top();
            pq.pop();
            if ( m_wtd.is_leaf(s.v) ){
                res.list.emplace_back(m_docperm.len2id[m_wtd.sym(s.v)], s.score);
            } else {
                auto exp_v = m_wtd.expand(s.v);
                auto exp_r_v = m_wtd.expand(s.v, s.r_v);
                auto exp_w = m_wtr.expand(s.w);
                auto exp_r_w = m_wtr.expand(s.w, s.r_w);

                if ( !m_wtd.empty(std::get<0>(exp_v)) ) {
                    push_node(pq, s, std::get<0>(exp_v), std::get<0>(exp_r_v), 
                                     std::get<0>(exp_w), std::get<0>(exp_r_w));
                }
                if ( !m_wtd.empty(std::get<1>(exp_v)) ) {
                    push_node(pq, s, std::get<1>(exp_v), std::get<1>(exp_r_v),
                                     std::get<1>(exp_w), std::get<1>(exp_r_w));
                }
            }
        }
        return res;
    }

    void load(sdsl::cache_config& cc){
        load_from_cache(m_csa, surf::KEY_CSA, cc, true);
        load_from_cache(m_df, surf::KEY_SADADF, cc, true);
        load_from_cache(m_wtr, surf::KEY_WTDUP, cc, true);
        std::cerr<<"m_wtr.size()="<<m_wtr.size()<<std::endl;
        std::cerr<<"m_wtr.sigma()="<<m_wtr.sigma<<std::endl;
        load_from_cache(m_wtd, surf::KEY_WTD, cc, true);
        std::cerr<<"m_wtd.size()="<<m_wtd.size()<<std::endl;
        std::cerr<<"m_wtd.sigma()="<<m_wtd.sigma<<std::endl;
        load_from_cache(m_docperm, surf::KEY_DOCPERM, cc); 
        m_ranker = ranker_type(cc);
    }

    size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const {
        structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
        size_type written_bytes = 0;
        written_bytes += m_csa.serialize(out, child, "csa");
        written_bytes += m_df.serialize(out, child, "df");
        written_bytes += m_wtd.serialize(out, child, "wtd");
        written_bytes += m_wtr.serialize(out, child, "wtr");
        written_bytes += m_docperm.serialize(out, child, "docperm");
        structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

};

template<typename t_csa,
         typename t_df,
         typename t_wtd,
         typename t_wtr,
         typename t_ranker,
         typename t_rbv,
         typename t_rrank>
void construct(idx_dr<t_csa,t_df,t_wtd,t_wtr, t_ranker, t_rbv, t_rrank>& idx,
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
    cout<<"...WTR"<<endl;
    if (!cache_file_exists<t_wtr>(surf::KEY_WTDUP2,cc)){
        construct_dup2<t_df>(cc); // construct DUP2 and DUPMARK
        t_wtr wtr;
        construct(wtr, cache_file_name(surf::KEY_DUP2, cc), cc);
        store_to_cache(wtr, surf::KEY_WTDUP2, cc, true);
        cout << "wtr.size() = " << wtr.size() << endl;
        cout << "wtr.sigma = " << wtr.sigma << endl;
    }
    cout<<"...WTR_BV"<<endl;
    if (!cache_file_exists<t_rbv>(surf::KEY_UMARK, cc) ){
        bit_vector bv;
        load_from_cache(bv, surf::KEY_UMARK, cc);
        t_rbv rbv(bv);
        store_to_cache(rbv, surf::KEY_UMARK, cc, true);
        t_rrank rrank(&rbv);
        store_to_cache(rrank, surf::KEY_URANK, cc, true);
    }
}

} // end namespace surf

#endif
