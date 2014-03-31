#ifndef SURF_IDX_SAWIT2_HPP
#define SURF_IDX_SAWIT2_HPP

#include "sdsl/suffix_trees.hpp"
#include "surf/df_sada.hpp"
#include "surf/rank_functions.hpp"
#include "surf/idx_sawit.hpp"
#include "surf/construct_col_len.hpp"
#include <algorithm>
#include <limits>
#include <queue>

namespace surf{

using range_type = sdsl::range_type;


template<typename t_wtu_node, typename t_wtp_node>
struct s_state2_t{
    double score;
    t_wtu_node v; // node in document array wavelet tree
    std::vector<term_info*> t_ptrs; // pointers to term_info array
    std::vector<range_type> r_v; // ranges in v
    t_wtp_node w; // node in duplication array wavelet tree
    std::vector<range_type> r_w; // ranges in w

    s_state2_t() = default;

    s_state2_t(double score, const t_wtu_node& v, 
              const std::vector<term_info*>& t_ptrs,
              const std::vector<range_type>& r_v,
              const t_wtp_node& w,
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

/*! Class sawit2 (Suffix Array Wavelet tree Index Type) consists of a
 *  (compressed) suffix array, a wavelet tree over the duplicate array,
 *  a (succinct) document frequency structure, and a wavelet tree
 *  over the unique array.
 *  
 */
template<typename t_csa,
         typename t_df,
         typename t_wtp,
         typename t_wtu,
         typename t_ubv=sdsl::rrr_vector<63>,
         typename t_urank=typename t_ubv::rank_1_type>
class idx_sawit2{
public:
    using size_type = sdsl::int_vector<>::size_type;
    typedef t_csa                        csa_type;
    typedef t_wtu                        wtu_type;
    typedef typename wtu_type::node_type node_type;
    typedef t_df                         df_type;
    typedef t_wtp                        wtp_type;
    typedef typename wtp_type::node_type node2_type;
    typedef t_ubv                        ubv_type;
    typedef t_urank                      urank_type;
    typedef rank_bm25<>                  ranker_type;
private:
    csa_type    m_csa;
    df_type     m_df;
    wtp_type    m_wtp; 
    t_wtu       m_wtu;
    ubv_type    m_ubv;
    urank_type  m_urank;
    doc_perm    m_docperm;
    ranker_type m_r;

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
                sp = m_urank(sp);
                ep = sp+f_Dt-1;
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
            auto min_idx = m_wtu.sym(v) << (m_wtu.max_level - v.level);  
            auto min_doc_len = m_r.doc_length(m_docperm.len2id[min_idx]);
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
                    auto score = m_r.calculate_docscore(
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
        pq.emplace(max_score, m_wtu.root(), term_ptrs, v_ranges, m_wtp.root(), w_ranges);
        if(profile) res.wt_search_space++;

        while ( !pq.empty() and res.list.size() < k ) {
            state_type s = pq.top();
            pq.pop();
            if ( m_wtu.is_leaf(s.v) ){
                res.list.emplace_back(m_docperm.len2id[m_wtu.sym(s.v)], s.score);
            } else {
                auto exp_v = m_wtu.expand(s.v);
                auto exp_r_v = m_wtu.expand(s.v, s.r_v);
                auto exp_w = m_wtp.expand(s.w);
                auto exp_r_w = m_wtp.expand(s.w, s.r_w);

                if ( !m_wtu.empty(std::get<0>(exp_v)) ) {
                    push_node(pq, s, std::get<0>(exp_v), std::get<0>(exp_r_v), 
                                     std::get<0>(exp_w), std::get<0>(exp_r_w));
                    if(profile) res.wt_search_space++;
                }
                if ( !m_wtu.empty(std::get<1>(exp_v)) ) {
                    push_node(pq, s, std::get<1>(exp_v), std::get<1>(exp_r_v),
                                     std::get<1>(exp_w), std::get<1>(exp_r_w));
                    if(profile) res.wt_search_space++;
                }
            }
        }
        return res;
    }

    void load(sdsl::cache_config& cc){
        load_from_cache(m_csa, surf::KEY_CSA, cc, true);
        load_from_cache(m_df, surf::KEY_SADADF, cc, true);
        load_from_cache(m_wtp, surf::KEY_WTDUP, cc, true);
        std::cerr<<"m_wtp.size()="<<m_wtp.size()<<std::endl;
        std::cerr<<"m_wtp.sigma()="<<m_wtp.sigma<<std::endl;
        load_from_cache(m_wtu, surf::KEY_WTU, cc, true);
        std::cerr<<"m_wtu.size()="<<m_wtu.size()<<std::endl;
        std::cerr<<"m_wtu.sigma()="<<m_wtu.sigma<<std::endl;
        load_from_cache(m_ubv, surf::KEY_UMARK, cc, true);
        std::cerr<<"m_ubv.size()="<<m_ubv.size()<<std::endl;
        load_from_cache(m_urank, surf::KEY_URANK, cc, true);
        m_urank.set_vector(&m_ubv);
        load_from_cache(m_docperm, surf::KEY_DOCPERM, cc); 
        m_r = ranker_type(cc);
    }

    size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const {
        structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
        size_type written_bytes = 0;
        written_bytes += m_csa.serialize(out, child, "csa");
        written_bytes += m_df.serialize(out, child, "df");
        written_bytes += m_wtp.serialize(out, child, "wtp");
        written_bytes += m_wtu.serialize(out, child, "wtu");
        written_bytes += m_ubv.serialize(out, child, "ubv");
        written_bytes += m_urank.serialize(out, child, "urank");
        written_bytes += m_docperm.serialize(out, child, "docperm");
        structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

};

template<typename t_csa,
         typename t_df,
         typename t_wtp,
         typename t_wtu,
         typename t_ubv,
         typename t_urank>
void construct(idx_sawit2<t_csa,t_df,t_wtp,t_wtu, t_ubv, t_urank>& idx,
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
    cout<<"...DF"<<endl;
    if (!cache_file_exists<t_df>(surf::KEY_SADADF, cc))
    {
        t_df df;
        construct(df, "", cc, 0);
        store_to_cache(df, surf::KEY_SADADF, cc, true);
    }
    cout<<"...WTP"<<endl;
    if (!cache_file_exists<t_wtp>(surf::KEY_WTDUP,cc)){
        t_wtp wtp;
        construct(wtp, cache_file_name(surf::KEY_DUP, cc), cc);
        store_to_cache(wtp, surf::KEY_WTDUP, cc, true);
        cout << "wtp.size() = " << wtp.size() << endl;
        cout << "wtp.sigma = " << wtp.sigma << endl;
    }
    cout<<"...WTU"<<endl;
    if (!cache_file_exists<t_wtu>(surf::KEY_WTU, cc) ){
        t_wtu wtu;
        construct(wtu, cache_file_name(surf::KEY_U, cc), cc);
        cout << "wtu.size() = " << wtu.size() << endl;
        cout << "wtu.sigma = " << wtu.sigma << endl;
        store_to_cache(wtu, surf::KEY_WTU, cc, true);
    }
    if (!cache_file_exists<t_ubv>(surf::KEY_UMARK, cc) ){
        bit_vector bv;
        load_from_cache(bv, surf::KEY_UMARK, cc);
        t_ubv ubv(bv);
        store_to_cache(ubv, surf::KEY_UMARK, cc, true);
        t_urank urank(&ubv);
        store_to_cache(urank, surf::KEY_URANK, cc, true);
    }
}

} // end namespace surf

#endif
