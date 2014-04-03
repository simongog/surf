#ifndef SURF_IDX_SAWIT3_HPP
#define SURF_IDX_SAWIT3_HPP

#include "sdsl/suffix_trees.hpp"
#include "surf/df_sada.hpp"
#include "surf/rank_functions.hpp"
#include "surf/idx_d.hpp"
#include "surf/idx_dr.hpp"
#include "surf/construct_col_len.hpp"
#include "surf/construct_U.hpp"
#include <algorithm>
#include <limits>
#include <queue>

namespace surf{

using range_type = sdsl::range_type;


/*! Class idx_d1r1 consists of a 
 *   - CSA over the collection concatenation
 *   - document frequency structure
 *   - a WT over the reduced D array (only 1-phrases)
 *   - a WT over the (1-phrases) sorted repetition array
 */
template<typename t_csa,
         typename t_df,
         typename t_wtr,
         typename t_wtd1,
         typename t_ranker=rank_bm25<>,
         typename t_d1bv=sdsl::rrr_vector<63>,
         typename t_d1rank=typename t_d1bv::rank_1_type,
         typename t_rbv=sdsl::rrr_vector<63>,
         typename t_rrank=typename t_rbv::rank_1_type
         >
class idx_d1r1{
public:
    using size_type = sdsl::int_vector<>::size_type;
    typedef t_csa                        csa_type;
    typedef t_wtd1                        wtd1_type;
    typedef typename wtd1_type::node_type node_type;
    typedef t_df                         df_type;
    typedef t_wtr                        wtr_type;
    typedef typename wtr_type::node_type node2_type;
    typedef t_d1bv                       d1bv_type;
    typedef t_d1rank                     d1rank_type;
    typedef t_rbv                        rbv_type;
    typedef t_rrank                      rrank_type;
    typedef t_ranker                     ranker_type;
private:
    csa_type    m_csa;
    df_type     m_df;
    wtr_type    m_wtr; 
    t_wtd1       m_wtd1;
    d1bv_type   m_d1bv;
    d1rank_type m_d1rank;
    rbv_type    m_rbv;
    rrank_type  m_rrank;
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
                sp = m_d1rank(sp);
                ep = sp+f_Dt-1;
                v_ranges.emplace_back(sp, ep);
                w_ranges.emplace_back(m_rrank(std::get<1>(df_info)),
                                      m_rrank(std::get<2>(df_info)+1)-1);
            }
        }
        term_ptrs.resize(terms.size()); 
        for (size_type i=0; i<terms.size(); ++i){
            term_ptrs[i] = &terms[i];
        }

        auto push_node = [this,&res,&profile,&ranked_and](pq_type& pq, state_type& s,node_type& v,std::vector<range_type>& r_v,
                                node2_type& w, std::vector<range_type>& r_w){
            auto min_idx = m_wtd1.sym(v) << (m_wtd1.max_level - v.level);  
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
                                 min_doc_len,
                                 m_wtd1.is_leaf(v)
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
        pq.emplace(max_score, m_wtd1.root(), term_ptrs, v_ranges, m_wtr.root(), w_ranges);
        if(profile) res.wt_search_space++;

        while ( !pq.empty() and res.list.size() < k ) {
            state_type s = pq.top();
            pq.pop();
            if ( m_wtd1.is_leaf(s.v) ){
                res.list.emplace_back(m_docperm.len2id[m_wtd1.sym(s.v)], s.score);
            } else {
                auto exp_v = m_wtd1.expand(s.v);
                auto exp_r_v = m_wtd1.expand(s.v, s.r_v);
                auto exp_w = m_wtr.expand(s.w);
                auto exp_r_w = m_wtr.expand(s.w, s.r_w);

                if ( !m_wtd1.empty(std::get<0>(exp_v)) ) {
                    push_node(pq, s, std::get<0>(exp_v), std::get<0>(exp_r_v), 
                                     std::get<0>(exp_w), std::get<0>(exp_r_w));
                }
                if ( !m_wtd1.empty(std::get<1>(exp_v)) ) {
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
    	std::string WTR_KEY = surf::KEY_WTR+"-"+std::to_string(1);
        load_from_cache(m_wtr, WTR_KEY, cc, true);
        std::cerr<<"m_wtr.size()="<<m_wtr.size()<<std::endl;
        std::cerr<<"m_wtr.sigma()="<<m_wtr.sigma<<std::endl;
        load_from_cache(m_wtd1, surf::KEY_WTU, cc, true);
        std::cerr<<"m_wtd1.size()="<<m_wtd1.size()<<std::endl;
        std::cerr<<"m_wtd1.sigma()="<<m_wtd1.sigma<<std::endl;
        load_from_cache(m_d1bv, surf::KEY_UMARK, cc, true);
        std::cerr<<"m_d1bv.size()="<<m_d1bv.size()<<std::endl;
        load_from_cache(m_d1rank, surf::KEY_URANK, cc, true);
        m_d1rank.set_vector(&m_d1bv);
        load_from_cache(m_rbv, surf::KEY_DUPMARK, cc, true);
        std::cerr<<"m_rbv.size()="<<m_rbv.size()<<std::endl;
        load_from_cache(m_rrank, surf::KEY_DUPRANK, cc, true);
        m_rrank.set_vector(&m_rbv);
        std::cerr<<"m_rrank(m_rbv.size())="<<m_rrank(m_rbv.size())<<std::endl;
        load_from_cache(m_docperm, surf::KEY_DOCPERM, cc); 
        m_ranker = ranker_type(cc);
    }

    size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const {
        structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
        size_type written_bytes = 0;
        written_bytes += m_csa.serialize(out, child, "CSA");
        written_bytes += m_df.serialize(out, child, "DF");
        written_bytes += m_wtr.serialize(out, child, "WTR");
        written_bytes += m_wtd1.serialize(out, child, "WTD1");
        written_bytes += m_d1bv.serialize(out, child, "D1_BV");
        written_bytes += m_d1rank.serialize(out, child, "D1_RANK");
        written_bytes += m_rbv.serialize(out, child, "R_BV");
        written_bytes += m_rrank.serialize(out, child, "R_RANK");
        written_bytes += m_docperm.serialize(out, child, "DOCPERM");
        structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

};

template<typename t_csa,
         typename t_df,
         typename t_wtr,
         typename t_wtd1,
         typename t_ranker,
         typename t_d1bv,
         typename t_d1rank,
         typename t_rbv,
         typename t_rrank>
void construct(idx_d1r1<t_csa,t_df,t_wtr,t_wtd1, t_ranker, t_d1bv, t_d1rank, t_rbv, t_rrank>& idx,
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
    cout<<"...WTR"<<endl;
    if (!cache_file_exists<t_wtr>(surf::KEY_WTDUP,cc)){
        t_wtr wtr;
        construct(wtr, cache_file_name(surf::KEY_DUP, cc), cc);
        store_to_cache(wtr, surf::KEY_WTDUP, cc, true);
        cout << "wtr.size() = " << wtr.size() << endl;
        cout << "wtr.sigma = " << wtr.sigma << endl;
    }
    cout<<"...U and Umark"<<endl;
    if (!cache_file_exists(surf::KEY_U,cc) or !cache_file_exists(surf::KEY_UMARK,cc)){
        construct_u<t_df>(cc);
    }
    cout<<"...WTU"<<endl;
    if (!cache_file_exists<t_wtd1>(surf::KEY_WTU, cc) ){
        t_wtd1 wtd1;
        construct(wtd1, cache_file_name(surf::KEY_U, cc), cc);
        cout << "wtd1.size() = " << wtd1.size() << endl;
        cout << "wtd1.sigma = " << wtd1.sigma << endl;
        store_to_cache(wtd1, surf::KEY_WTU, cc, true);
    }
    cout<<"...D1_BV"<<endl;
    if (!cache_file_exists<t_d1bv>(surf::KEY_UMARK, cc) ){
        bit_vector bv;
        load_from_cache(bv, surf::KEY_UMARK, cc);
        t_d1bv d1bv(bv);
        store_to_cache(d1bv, surf::KEY_UMARK, cc, true);
        t_d1rank d1rank(&d1bv);
        store_to_cache(d1rank, surf::KEY_URANK, cc, true);
    }
    cout<<"...WTR2"<<endl;
    const uint64_t depth=1; // depth of sorting in the repetition structure
    std::string R_KEY = surf::KEY_R+"-"+to_string(depth);
    std::string WTR_KEY = surf::KEY_WTR+"-"+to_string(depth);
    if (!cache_file_exists<t_wtr>(WTR_KEY,cc)){
        string dup2_file = cache_file_name(surf::KEY_DUP2,cc);
        if (!cache_file_exists(surf::KEY_DUP2,cc)){
            construct_dup2<t_df>(cc); // construct DUP2 and DUPMARK
        }
        if (!cache_file_exists<t_rbv>(surf::KEY_DUPMARK, cc) ){
            bit_vector bv;
            load_from_cache(bv, surf::KEY_DUPMARK, cc);
            t_rbv rbv(bv);
            store_to_cache(rbv, surf::KEY_DUPMARK, cc, true);
            t_rrank rrank(&rbv);
            store_to_cache(rrank, surf::KEY_DUPRANK, cc, true);
        }
        if ( !cache_file_exists(R_KEY, cc) ) {
            std::cout<<"generate "<<R_KEY<<" file"<<std::endl;
            cout<<".........load df"<<endl;
            t_df df;
            load_from_cache(df, surf::KEY_SADADF, cc, true);
            cout<<".........load cst"<<endl;
            using cst_type =  typename t_df::cst_type;
            cst_type cst;
            load_from_file(cst, cache_file_name<cst_type>(surf::KEY_TMPCST, cc));
            cout<<".........cst.size()="<<cst.size()<<endl;
            int_vector_buffer<> dup(cache_file_name(surf::KEY_DUP, cc));
            cout<<"dup.size()="<<dup.size()<<" dup.width()="<<(int)dup.width()<<endl;
            int_vector_buffer<> R(cache_file_name(R_KEY,cc), 
                                  std::ios::out, 1024*1024, dup.width());
            cout<<"R intialized"<<endl;
            auto root = cst.root();
            auto left_most_leaf = cst.select_leaf(1);
            for (const auto& v : cst.children(root) ){
                if ( v == left_most_leaf )
                    continue;
                auto lb = cst.lb(v);
                auto rb = cst.rb(v);
                auto df_info = df(lb, rb);
                std::vector<uint64_t> buf;
                for (uint64_t i = std::get<1>(df_info); i <= std::get<2>(df_info); ++i) {
                    buf.push_back(dup[i]); 
                }
                std::sort(buf.begin(), buf.end());
                for (size_t i=0; i < buf.size(); ++i){
                    R.push_back(buf[i]);
                }
            }
        }
        {
            cout<<"......generate WT"<<endl;
            t_wtr wtr2;
            construct(wtr2, cache_file_name(R_KEY, cc), cc);
            store_to_cache(wtr2, WTR_KEY,cc,true);
        }
    }
}

} // end namespace surf

#endif
