#ifndef SURF_IDX_SAWIT2_HPP
#define SURF_IDX_SAWIT2_HPP

#include "sdsl/suffix_trees.hpp"
#include "surf/df_sada.hpp"
#include "surf/rank_functions.hpp"
#include "surf/idx_sawit.hpp"
#include <algorithm>
#include <limits>
#include <queue>

namespace surf{

using range_type = sdsl::range_type;

template<typename t_wt_node, typename t_wtdup_node>
struct s_state2_t{
    double score;
    t_wt_node v; // node in document array wavelet tree
    std::vector<term_info*> t_ptrs; // pointers to term_info array
    std::vector<range_type> r_v; // ranges in v
    t_wtdup_node w; // node in duplication array wavelet tree
    std::vector<range_type> r_w; // ranges in w

    s_state2_t() = default;

    s_state2_t(double score, const t_wt_node& v, 
              const std::vector<term_info*>& t_ptrs,
              const std::vector<range_type>& r_v,
              const t_wtdup_node& w,
              const std::vector<range_type>& r_w):
        score(score), v(v), t_ptrs(t_ptrs),
        r_v(r_v),w(w),r_w(r_w){}

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
class idx_sawit2{
public:
    using size_type = sdsl::int_vector<>::size_type;
    typedef t_csa   csa_type;
    typedef t_wtd   wtd_type;
    typedef typename wtd_type::node_type node_type;
    typedef t_df    df_type;
    typedef t_wtdup wtdup_type;
    typedef typename wtdup_type::node_type node2_type;
    typedef rank_bm25<> ranker_type;
private:
    csa_type    m_csa;
    wtd_type    m_wtd;
    df_type     m_df;
    wtdup_type  m_wtdup; 
    doc_perm    m_docperm;
    ranker_type m_r;

    using state_type = s_state2_t<node_type, node2_type>;
public:

    result_t search(const std::vector<query_token>& qry,size_t k) {
        typedef std::priority_queue<state_type> pq_type;
        std::vector<term_info> terms;
        std::vector<term_info*> term_ptrs;
        std::vector<range_type> v_ranges; // ranges in wtd
        std::vector<range_type> w_ranges; // ranges in wtdup

        for (size_t i=0; i<qry.size(); ++i){
            size_type sp=1, ep=0;
            if ( backward_search(m_csa, 0, m_csa.size()-1, qry[i].token_id, sp, ep) > 0 ) {
                auto df_info = m_df(sp,ep);
                auto f_Dt = std::get<0>(df_info); // document frequency
                terms.emplace_back(qry[i].token_id, qry[i].f_qt, sp, ep,  f_Dt);
                v_ranges.emplace_back(sp, ep);
                w_ranges.emplace_back(std::get<1>(df_info),std::get<2>(df_info));
/*
                std::cout<<"f_Dt="<<f_Dt<<std::endl;
                std::set<uint64_t> sss;
                for (size_t ii = sp; ii<=ep; ++ii){
                    sss.insert(m_wtd[ii]);
                }
                std::cout<<"sss.size()="<<sss.size()<<std::endl;
                std::cout<<"wtd=";
                for(size_t ii=sp; ii<=ep; ++ii){
                    std::cout<<" "<<m_wtd[ii];
                }
                std::cout<<std::endl;
                std::cout<<"wtdup=";
                for(size_t ii=std::get<1>(df_info); ii<=std::get<2>(df_info); ++ii){
                    std::cout<<" "<<m_wtdup[ii];
                }
                std::cout<<std::endl;
*/
/*                
                std::map<uint64_t, uint64_t> vm, wm;
                for(size_t ii=sp; ii<=ep; ++ii){
                    vm[m_wtd[ii]]++;
                }
                for(size_t ii=std::get<1>(df_info); ii<=std::get<2>(df_info); ++ii){
                    wm[m_wtdup[ii]]++;
                }
                for(auto x : vm){
                    if ( x.second > 1 ){
                        if ( wm.find(x.first) == wm.end() ){
                            std::cout<<"error: "<<x.first<<" is in D with count "<<x.second<<" but not in DUP"<<std::endl;
                        } else{
                            if ( wm[x.first]+1 != x.second ){
                                std::cout<<"x=("<<x.first<<","<<x.second<<") wm[x.first]="<<wm[x.first]<<std::endl;
                            }
                        }
                    }
                }
*/                
            }
        }
        term_ptrs.resize(terms.size()); 
        for (size_type i=0; i<terms.size(); ++i){
            term_ptrs[i] = &terms[i];
        }

        auto push_node = [this](pq_type& pq, state_type& s,node_type& v,std::vector<range_type>& r_v,
                                node2_type& w, std::vector<range_type>& r_w){
            auto min_idx = m_wtd.sym(v) << (m_wtd.max_level - v.level);  
            auto min_doc_len = m_r.doc_length(m_docperm.len2id[min_idx]);
            state_type t; // new state
            t.v = v;
            t.w = w;
            t.score = 0;
            for (size_t i = 0; i < r_v.size(); ++i){
                if ( !empty(r_v[i]) ){
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
                } 
            }
            pq.emplace(t);       
        };

        constexpr double max_score = std::numeric_limits<double>::max();
        
        pq_type pq;
        size_type search_space=0;
        pq.emplace(max_score, m_wtd.root(), term_ptrs, v_ranges, m_wtdup.root(), w_ranges);
        ++search_space;

        result_t res;
        while ( !pq.empty() and res.size() < k ) {
            state_type s = pq.top();
            pq.pop();
            if ( m_wtd.is_leaf(s.v) ){
                res.emplace_back(m_docperm.len2id[m_wtd.sym(s.v)], s.score);
            } else {
                auto exp_v = m_wtd.expand(s.v);
                auto exp_r_v = m_wtd.expand(s.v, s.r_v);
                auto exp_w = m_wtdup.expand(s.w);
                auto exp_r_w = m_wtdup.expand(s.w, s.r_w);
                if ( !m_wtd.empty(std::get<0>(exp_v)) ) {
                    push_node(pq, s, std::get<0>(exp_v), std::get<0>(exp_r_v), 
                                     std::get<0>(exp_w), std::get<0>(exp_r_w));
                    ++search_space;
                }
                if ( !m_wtd.empty(std::get<1>(exp_v)) ) {
                    push_node(pq, s, std::get<1>(exp_v), std::get<1>(exp_r_v),
                                     std::get<1>(exp_w), std::get<1>(exp_r_w));
                    ++search_space;
                }
            }
        }
        std::cerr << "search_space = " << search_space << std::endl;
        for(size_t i=0; i<res.size(); ++i){
            std::cerr<<res[i].score<<","<<res[i].doc_id<<std::endl;
        }
        return res;
    }

    void load(sdsl::cache_config& cc){
        load_from_cache(m_csa, surf::KEY_CSA, cc, true);
        load_from_cache(m_wtd, surf::KEY_WTD, cc, true);
        load_from_cache(m_df, surf::KEY_SADADF, cc, true);
        load_from_cache(m_wtdup, surf::KEY_WTDUP, cc, true);
        load_from_cache(m_docperm, surf::KEY_DOCPERM, cc); 
        m_r = ranker_type(cc);
    }

    size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const {
        structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
        size_type written_bytes = 0;
        written_bytes += m_csa.serialize(out, child, "csa");
        written_bytes += m_wtd.serialize(out, child, "wtd");
        written_bytes += m_df.serialize(out, child, "df");
        written_bytes += m_wtdup.serialize(out, child, "wtdup");
        written_bytes += m_docperm.serialize(out, child, "docperm");
        structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

};

template<typename t_csa,
         typename t_wtd,
         typename t_df,
         typename t_wtdup>
void construct(idx_sawit2<t_csa,t_wtd,t_df,t_wtdup>& idx,
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
        construct_doc_perm<t_csa::alphabet_type::int_width>(cc);
        construct_darray<t_csa::alphabet_type::int_width>(cc);
        t_wtd wtd;
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
