#ifndef SURF_IDX_NN_HPP
#define SURF_IDX_NN_HPP

#include "sdsl/suffix_trees.hpp"
#include "sdsl/k2_treap.hpp"
#include "surf/df_sada.hpp"
#include "surf/rank_functions.hpp"
#include "surf/idx_d.hpp"
#include "surf/construct_col_len.hpp"
#include "surf/construct_max_doc_len.hpp"
#include "surf/construct_DUP2.hpp"
#include <algorithm>
#include <limits>
#include <queue>
#include <set>

namespace surf{

using range_type = sdsl::range_type;


template<typename t_select>
struct map_to_dup3_type{
    const t_select* m_sel;

    map_to_dup3_type(const t_select* select=nullptr) : 
        m_sel(select)
    {}

    range_type
    operator()(size_t sp, size_t ep)const{
        uint64_t y = (*m_sel)(ep);
        uint64_t ep_ = (y + 1) - ep - 1; // # of 0 left to y - 1
        uint64_t sp_ = 0;          // # of 0 left to x
        if (0 == sp) {
        } else {
            uint64_t x = (*m_sel)(sp);
            sp_ = (x+1) - sp;
        }
        return std::make_pair(sp_, ep_);       
    }
};



/*! Class idx_nn consists of a 
 *   - CSA over the collection concatenation
 *   - H 
 */
// TODO:
//  - sort node-ranges in DOC
template<typename t_csa,
         typename t_df,
         typename t_wtd,
         typename t_wtr,
         typename t_rbv=sdsl::rrr_vector<63>,
         typename t_rrank=typename t_rbv::rank_1_type>
class idx_nn{
public:
    using size_type = sdsl::int_vector<>::size_type;
    typedef t_csa                                      csa_type;
    typedef t_df                                       df_type;
    typedef t_wtd                                      wtd_type;
    typedef t_wtr                                      wtr_type;
    typedef t_rbv                                      rbv_type;
    typedef t_rrank                                    rrank_type;
    typedef sd_vector<>                                border_type;
    typedef sd_vector<>::rank_1_type                   border_rank_type;
    typedef sd_vector<>::select_1_type                 border_select_type;
    typedef rrr_vector<63>                             h3_type;
    typedef rrr_vector<63>::select_1_type              h3_select_type;
    typedef rmq_succinct_sct<>                         rmqc_type;
    typedef k2_treap<2,rrr_vector<63>>                 k2treap_type;
    typedef k2_treap_ns::top_k_iterator<k2treap_type>  k2treap_iterator;
    typedef typename t_csa::alphabet_category          alphabet_category;

    typedef map_to_dup3_type<h3_select_type> map_to_h3_type;
public:
    csa_type           m_csa;
    border_type        m_border;
    border_rank_type   m_border_rank;
    border_select_type m_border_select;
    h3_type            m_h3;
    h3_select_type     m_h3_select;
    int_vector<>       m_doc; // documents in node lists
    rmqc_type          m_rmqc;
    k2treap_type       m_k2treap;
    map_to_h3_type     m_map_to_h3;

public:

    class top_k_iterator{
        public:
            typedef void(*t_mfptr)();
            typedef std::pair<uint64_t, double> t_doc_val;
            typedef std::stack<std::array<uint64_t,2>> t_stack_array;
        private:
            const idx_nn*      m_idx;
            uint64_t           m_sp;  // start point of lex interval
            uint64_t           m_ep;  // end point of lex interval
            t_doc_val          m_doc_val;  // stores the current result
            bool               m_valid = false;
            k2treap_iterator   m_k2_iter;
            std::set<uint64_t> m_reported;
            std::set<uint64_t> m_singletons;
            t_stack_array      m_states;
        public:
            top_k_iterator() = delete;
            top_k_iterator(const top_k_iterator&) = default;
            top_k_iterator(top_k_iterator&&) = default;
            top_k_iterator& operator=(const top_k_iterator&) = default;
            top_k_iterator& operator=(top_k_iterator&&) = default;
            template<typename t_pat_iter>
            top_k_iterator(const idx_nn* idx, t_pat_iter begin, t_pat_iter end) : m_idx(idx) {
                m_valid = backward_search(m_idx->m_csa, 0, m_idx->m_csa.size()-1, begin, end, m_sp, m_ep) > 0;
                if ( m_valid ){
                    auto h3_range = m_idx->m_map_to_h3(m_sp, m_ep);
                    if ( !empty(h3_range) ) {
                        uint64_t depth = 0;
                        // determine depth
                        size_type _sp=0, _ep = m_idx->m_csa.size()-1;
                        for (auto it = begin+1; it <= end; ++it) {
                            size_type __sp, __ep;
                            backward_search(m_idx->m_csa, 0, m_idx->m_csa.size()-1, begin, it, __sp, __ep);
                            if ( __ep-__sp+1 < _ep-_sp+1 ){
                                ++depth;
                            }
                            _sp = __sp;
                            _ep = __ep;
                        } 
                        m_k2_iter = top_k(m_idx->m_k2treap, {std::get<0>(h3_range) ,0}, {std::get<1>(h3_range), depth-1});
                    }
                    m_states.push({m_sp, m_ep});
                    ++(*this);
                }
            }

            top_k_iterator& operator++(){
                if ( m_valid ){
                    m_valid = false;
                    if ( m_k2_iter ) { // multiple occurrence result exists
                        auto xy_w       = *m_k2_iter;
                        uint64_t doc_id = m_idx->m_doc[real(xy_w.first)]; 
                        m_doc_val = t_doc_val(doc_id, xy_w.second+1);
                        m_reported.insert(doc_id);
                        m_valid = true;
                        ++m_k2_iter;
                    } else { // search for singleton results
                        while ( !m_states.empty() ) {
                            auto state = m_states.top();
                            m_states.pop();
                            uint64_t min_idx = m_idx->m_rmqc(state[0], state[1]);
                            uint64_t doc_id  = m_idx->m_border_rank(m_idx->m_csa[min_idx]);
                            if ( m_singletons.find(doc_id) == m_singletons.end() ){
                                m_singletons.insert(doc_id);
                                if ( min_idx + 1 <= state[1] )
                                    m_states.push({min_idx+1, state[1]});
                                if ( state[0] + 1 <= min_idx )
                                    m_states.push({state[0], min_idx-1});
                                if ( m_reported.find(doc_id) == m_reported.end() ){
                                    m_doc_val = t_doc_val(doc_id, 1);
                                    m_reported.insert(doc_id);
                                    m_valid = true;
                                    break;
                                }
                            }
                        }
                    }
                }
                return *this;
            }

            t_doc_val operator*() const{
                return m_doc_val; 
            }

            operator t_mfptr() const{
                return (t_mfptr)(m_valid);
            }
    };

    template<typename t_pat_iter>
    top_k_iterator topk(t_pat_iter begin, t_pat_iter end) const{
        return top_k_iterator(this, begin, end);
    }

    result search(const std::vector<query_token>& qry,size_t k,bool ranked_and = false,bool profile = false) const {
        result res;
        if ( qry.size() > 0 ){
            auto res_iter = topk(qry[0].token_ids.begin(), qry[0].token_ids.end());
            size_t i = 0;
            while ( i < k and res_iter ){
                ++i;
                auto docid_weight = *res_iter;
                res.list.emplace_back(docid_weight.first, docid_weight.second);
                ++res_iter;
            }
        }
        return res;
    }

/*                    
                    size_type doc_begin=0;
                    if ( doc_id ) {
                        doc_begin = m_border_select(doc_id)+1;
                    }
                    size_type doc_end=m_border_select(doc_id+1)-1;
                    cout<<"doc text range=["<<doc_begin<<","<<doc_end<<"]"<<endl;
                    cout << extract(m_csa, doc_begin, doc_end) << endl;
*/                  


    void load(sdsl::cache_config& cc){
        load_from_cache(m_csa, surf::KEY_CSA, cc, true);
        load_from_cache(m_doc, surf::KEY_DUP3, cc);
        load_from_cache(m_border, surf::KEY_DOCBORDER, cc, true); 
        load_from_cache(m_border_rank, surf::KEY_DOCBORDER_RANK, cc, true); 
        m_border_rank.set_vector(&m_border);
        load_from_cache(m_border_select, surf::KEY_DOCBORDER_SELECT, cc, true); 
        m_border_select.set_vector(&m_border);
        load_from_cache(m_h3, surf::KEY_H3, cc, true);
        load_from_cache(m_h3_select, surf::KEY_H3_SELECT, cc, true);
        m_h3_select.set_vector(&m_h3);
        m_map_to_h3 = map_to_h3_type(&m_h3_select);
        load_from_cache(m_rmqc, surf::KEY_RMQC, cc, true); 
        load_from_cache(m_k2treap, surf::KEY_W_AND_P, cc, true); 
//        int_vector<> darray;
//        load_from_cache(darray, surf::KEY_DARRAY, cc);
//        std::cout<<"DARRAY="<<darray <<std::endl;
    }

    size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const {
        structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
        size_type written_bytes = 0;
        written_bytes += m_csa.serialize(out, child, "CSA");
        written_bytes += m_doc.serialize(out, child, "DOC");
        written_bytes += m_border.serialize(out, child, "BORDER");
        written_bytes += m_border_rank.serialize(out, child, "BORDER_RANK");
        written_bytes += m_h3.serialize(out, child, "H3");
        written_bytes += m_h3_select.serialize(out, child, "H3_SELECT");
        written_bytes += m_rmqc.serialize(out, child, "RMQ_C");
        written_bytes += m_k2treap.serialize(out, child, "W_AND_P");
        structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    void mem_info()const{}
};

template<typename t_df,
         typename t_rank>
struct map_to_dup2_type{
    const t_df* m_df;
    const t_rank* m_rank;

    map_to_dup2_type(const t_df* df, const t_rank* rank) : 
        m_df(df), m_rank(rank)
    {}

    range_type
    operator()(size_t sp, size_t ep)const{
        auto df_info = (*m_df)(sp,ep); // select on H
        return std::make_pair((*m_rank)(std::get<1>(df_info)),
                              (*m_rank)(std::get<2>(df_info)+1)-1);
    }
};

template<typename t_df,
         typename t_rank>
struct map_node_to_dup2_type{
    typedef typename t_df::cst_type t_cst;
    typedef typename t_cst::node_type t_node;
    const map_to_dup2_type<t_df,t_rank> m_map;
    const t_cst* m_cst;

    map_node_to_dup2_type(const t_df* df, const t_rank* rank, const t_cst* cst):
        m_map(df, rank), m_cst(cst)
    { }

    range_type
    operator()(const t_node& v, size_t l_child, size_t r_child)const{
        auto left    = l_child + (r_child - l_child) / 2;
        auto left_rb = m_cst->rb(m_cst->select_child(v, left));
        return m_map(left_rb, left_rb+1);
    }
    // id \in [1..n-1]
    uint64_t id(const t_node& v, size_t l_child, size_t r_child)const{
        auto left    = l_child + (r_child - l_child) / 2;
        return m_cst->rb(m_cst->select_child(v, left))+1;
    }
};
template<typename t_cst,
         typename t_select>
struct map_node_to_dup3_type{
    typedef typename t_cst::node_type t_node;
    const map_to_dup3_type<t_select> m_map;
    const t_cst* m_cst;

    map_node_to_dup3_type(const t_select* select, const t_cst* cst):
        m_map(select), m_cst(cst)
    { }

    range_type
    operator()(const t_node& v)const{
        auto left    = 1 + (m_cst->degree(v) - 1) / 2;
        auto left_rb = m_cst->rb(m_cst->select_child(v, left));
        return m_map(left_rb, left_rb+1);
    }
    // id \in [1..n-1]
    uint64_t id(const t_node& v)const{
        auto left    = 1 + (m_cst->degree(v) - 1) / 2;
        return m_cst->rb(m_cst->select_child(v, left))+1;
    }
};



template<typename t_csa,
         typename t_df,
         typename t_wtd,
         typename t_wtr,
         typename t_rbv,
         typename t_rrank>
void construct(idx_nn<t_csa,t_df,t_wtd,t_wtr, t_rbv, t_rrank>& idx,
               const std::string&,
               sdsl::cache_config& cc, uint8_t num_bytes)
{    
    using namespace sdsl;
    using namespace std;
    using cst_type = typename t_df::cst_type;
    using k2treap_type = typename idx_nn<t_csa,t_df,t_wtd,t_wtr,t_rbv,t_rrank>::k2treap_type;

    construct_col_len<t_df::alphabet_category::WIDTH>(cc);

    cout<<"...CSA"<<endl; // CSA to get the lex. range
    if ( !cache_file_exists<t_csa>(surf::KEY_CSA, cc) )
    {
        t_csa csa;
        construct(csa, "", cc, 0);
        store_to_cache(csa, surf::KEY_CSA, cc, true);
    }
    cout<<"...WTD"<<endl; // Document array and wavelet tree of it
    if (!cache_file_exists<t_wtd>(surf::KEY_WTD, cc) ){
        construct_darray<t_csa::alphabet_type::int_width>(cc, false);
        t_wtd wtd;
        construct(wtd, cache_file_name(surf::KEY_DARRAY, cc), cc);
        cout << "wtd.size() = " << wtd.size() << endl;
        cout << "wtd.sigma = " << wtd.sigma << endl;
        store_to_cache(wtd, surf::KEY_WTD, cc, true);
    }
    cout<<"...DF"<<endl; // 
    if (!cache_file_exists<t_df>(surf::KEY_SADADF, cc))
    {
        t_df df;
        construct(df, "", cc, 0);
        store_to_cache(df, surf::KEY_SADADF, cc, true);
    }
    cout<<"...REP2"<<endl;
    if (!cache_file_exists<t_wtr>(surf::KEY_WTDUP2,cc)){
        construct_dup2<t_df>(cc); // construct DUP2 and DUPMARK
    }
    cout<<"...REP_BV"<<endl;
    if (!cache_file_exists<t_rbv>(surf::KEY_DUPMARK, cc) ){
        bit_vector bv;
        load_from_cache(bv, surf::KEY_DUPMARK, cc);
        t_rbv rbv(bv);
        store_to_cache(rbv, surf::KEY_DUPMARK, cc, true);
        t_rrank rrank(&rbv);
        store_to_cache(rrank, surf::KEY_DUPRANK, cc, true);
    }
    cout<<"...DOC_BORDER"<<endl;
    {
        bit_vector doc_border;
        load_from_cache(doc_border, surf::KEY_DOCBORDER,cc);
        sd_vector<> sd_doc_border(doc_border);
        store_to_cache(sd_doc_border, surf::KEY_DOCBORDER, cc, true);
        sd_vector<>::rank_1_type doc_border_rank(&sd_doc_border);
        store_to_cache(doc_border_rank, surf::KEY_DOCBORDER_RANK, cc, true);
        sd_vector<>::select_1_type doc_border_select(&sd_doc_border);
        store_to_cache(doc_border_select, surf::KEY_DOCBORDER_SELECT, cc, true);
    }
    cout<<"...DOC"<<endl; 
    size_t max_depth = 0;
    {
        t_rbv rbv;
        load_from_cache(rbv, surf::KEY_DUPMARK, cc, true);
        std::cerr<<"m_rbv.size()="<<rbv.size()<<std::endl;
        t_rrank rrank(&rbv);
        load_from_cache(rrank, surf::KEY_DUPRANK, cc, true);
        rrank.set_vector(&rbv);
        std::cerr<<"rrank(rbv.size())="<<rrank(rbv.size())<<std::endl;
        int_vector<> dup2;
        load_from_cache(dup2, surf::KEY_DUP2, cc);
        std::cerr<<"dup2.size()="<<dup2.size()<<std::endl;
        t_df df;
        load_from_cache(df, surf::KEY_SADADF, cc, true);
        std::cerr<<"df loaded"<<std::endl;
        map_to_dup2_type<t_df,t_rrank> map_to_dup2(&df,&rrank);
        cst_type cst;
        load_from_file(cst, cache_file_name<cst_type>(surf::KEY_TMPCST, cc));
        map_node_to_dup2_type<t_df,t_rrank> map_node_to_dup2(&df,&rrank,&cst);
        // tuple < node, leftmost child, rightmost child, first visit, [sp_h, ep_h],depth>
        using node_type = typename cst_type::node_type;
        using n_type = std::tuple<node_type, size_t, size_t, bool, range_type, size_t>;
        std::stack<n_type> s;
        std::vector<std::set<uint64_t>> data(1);
        cout<<"interval   ["<<cst.lb(cst.root())<<","<<cst.rb(cst.root())<<"]"<<endl;
        auto x = map_to_dup2(cst.lb(cst.root()),cst.rb(cst.root()));
        cout<<"h-interval ["<<get<0>(x)<<","<<get<1>(x)<<"]"<<endl;
        auto y = map_node_to_dup2(cst.root(), 1, cst.degree(cst.root()));
        cout<<"h-node     ["<<get<0>(y)<<","<<get<1>(y)<<"]"<<endl;
        uint64_t real_nodes=0;
        uint64_t virtual_nodes=0;

        auto add_child = [&s,&cst,&map_node_to_dup2,&real_nodes,&virtual_nodes]
                         (const node_type& v, size_t l_child, size_t r_child, size_t depth){
            if ( r_child == l_child ){ // real node
                auto w = cst.select_child(v, l_child);
                if (!cst.is_leaf(w)){
                    auto degree = cst.degree(w);
                    auto r = map_node_to_dup2(w, 1, degree);
                    s.emplace(w, 1, degree, true, r,depth+1);
                    ++real_nodes;
                } 
            } else { // virtual node
                auto r = map_node_to_dup2(v, l_child, r_child);
                s.emplace(v, l_child, r_child, true, r,depth);
                ++virtual_nodes;
            }
        };

        auto is_real_node = [&cst](const node_type&v, size_t l_child, size_t r_child){
            return l_child==1 and cst.degree(v)==r_child;
        };

        add_child(cst.root(), 1, cst.degree(cst.root()),1);
        size_t count=0;
        uint64_t R_size = 0;
        uint64_t R_red_size = 0;
        vector<vector<uint64_t>> reps(cst.size());

// FROM binary suffix tree back to suffix tree
// DUP3 will contain the duplicated documents
// from real ST nodes
//        
// Select on H3 will be used to index them
//
        // invariant: node has two or more children
        while (!s.empty()) {
            n_type node = s.top();
            s.pop();
            
            auto v = get<0>(node);
            auto l_child = get<1>(node);
            auto r_child = get<2>(node);
            auto first = get<3>(node);
            auto depth = get<5>(node);
            if ( depth > max_depth )
                max_depth = depth;
            
            if (first) {  // first half
                R_size += size(get<4>(node));
                if ( is_real_node(v, l_child, r_child) ){
                    data.push_back(std::set<uint64_t>());
                }
                if ( size(get<4>(node)) > 0 ){
                    if ( depth+1 != data.size() ){
                        cout<<"ERROR!"<<endl;
                        return;
                    }
                    for (size_t i=get<4>(node).first; i<=get<4>(node).second; ++i){
                        data[depth].insert(dup2[i]);
                    }
                }
                // recurse down
                get<3>(node) = false;
                auto lb = cst.lb(cst.select_child(v, l_child));
                auto rb = cst.rb(cst.select_child(v, r_child));
                if ( count++ < 10 ){
                    bool b = is_real_node(v, l_child, r_child);
                    cout<<"..["<<lb<<","<<rb<<"] "<<b;
                    cout<<" _"<<l_child<<"."<<r_child<<"_";
                    if (size(get<4>(node)))
                        cout <<" ["<< get<4>(node).first<<","
                                   <<get<4>(node).second<<"]";
                    cout<<endl;
                }
                s.push(node);
                auto mid = l_child + (r_child - l_child) / 2;
                add_child(v, mid+1, r_child, depth);
                add_child(v, l_child, mid, depth);
            } else {  // second half
                if ( is_real_node(v, l_child, r_child) ){
                    auto in_order_id = map_node_to_dup2.id(v, l_child, r_child);
                    if ( reps[in_order_id].size() != 0 ){
                        cout<<"ERROR reps"<<endl;
                        return;
                    }
                    for(auto x : data.back()){
                        reps[in_order_id].push_back(x);
                    }
                    R_red_size += data.back().size();
                    data.pop_back();
                }
            }
            
        }
        cout<<"max_depth="<<max_depth<<endl;
        cout<<"R_size="<<R_size<<endl;
        cout<<"R_red_size="<<R_red_size<<endl;
        cout<<"real_nodes="<<real_nodes<<endl;
        cout<<"virtual_nodes="<<virtual_nodes<<endl;
        cout<<"cst.nodes()="<<cst.nodes()<<endl;
        cout<<"cst.size()="<<cst.size()<<endl;
        cout<<"cst inner nodes="<<cst.nodes()-cst.size()<<endl;
        {
            int_vector_buffer<1> h3(cache_file_name(surf::KEY_H3, cc), ios::out);
            int_vector_buffer<>  dup3(cache_file_name(surf::KEY_DUP3, cc), 
                                      ios::out, 1024*1024, dup2.width());
            for(size_t i=1; i<cst.size();++i){
                for(size_t j=0; j<reps[i].size(); ++j){
                    h3.push_back(0);
                    dup3.push_back(reps[i][j]);
                }
                h3.push_back(1);
            }
        }
        bit_vector h3;
        load_from_cache(h3, surf::KEY_H3, cc);
        rrr_vector<63> h3rrr(h3);
        store_to_cache(h3rrr, surf::KEY_H3, cc, true);
        rrr_vector<63>::select_1_type h3_select(&h3rrr);
        store_to_cache(h3_select, surf::KEY_H3_SELECT, cc, true);
    }
    cout<<"...WTD"<<endl;
    if (!cache_file_exists<t_wtd>(surf::KEY_WTD, cc) ){
        construct_darray<t_csa::alphabet_type::int_width>(cc, false);
        t_wtd wtd;
        construct(wtd, cache_file_name(surf::KEY_DARRAY, cc), cc);
        cout << "wtd.size() = " << wtd.size() << endl;
        cout << "wtd.sigma = " << wtd.sigma << endl;
        store_to_cache(wtd, surf::KEY_WTD, cc, true);
    }
// P corresponds to up-pointers
// W to the weight of the element
    cout<<"...P and W"<<endl;
    {
        construct_max_doc_len<t_csa::alphabet_type::int_width>(cc);    
        uint64_t max_len = 0;
        load_from_cache(max_len, surf::KEY_MAXDOCLEN, cc);

        int_vector<> dup3;
        load_from_cache(dup3, surf::KEY_DUP3, cc);
        cout<<"dup3.size()="<<dup3.size()<<endl;
        int_vector<> P(dup3.size(), 0, sdsl::bits::hi(max_depth)+1);

        int_vector<> W(dup3.size(), 0, bits::hi(max_len)+1);
        t_wtd wtd;
        load_from_cache(wtd, surf::KEY_WTD, cc, true);
        
        
        rrr_vector<63> h3rrr;
        load_from_cache(h3rrr, surf::KEY_H3, cc, true);
        rrr_vector<63>::select_1_type h3_select;
        load_from_cache(h3_select, surf::KEY_H3_SELECT, cc, true);
        h3_select.set_vector(&h3rrr);
        cst_type cst;
        load_from_file(cst, cache_file_name<cst_type>(surf::KEY_TMPCST, cc));
        map_node_to_dup3_type<cst_type, rrr_vector<63>::select_1_type> map_node_to_dup3(&h3_select, &cst);

        uint64_t doc_cnt = 1;
        load_from_cache(doc_cnt, KEY_DOCCNT, cc);
        typedef stack<uint32_t,vector<uint32_t>> t_stack;
//  HELPER to build the pointer structure
        vector<t_stack> depths(doc_cnt, t_stack(vector<uint32_t>(1,0)));// doc_cnt stack for last depth
        uint64_t depth=0;
        // DFS traversal of CST
        for (auto it=cst.begin(); it!=cst.end(); ++it) {
            auto v = *it; // get the node by dereferencing the iterator
            if ( !cst.is_leaf(v) ) {
                if (it.visit() == 1) {  // node visited the first time
                    ++depth;
                    range_type r = map_node_to_dup3(v);
                    if ( !empty(r) ){
                        for(size_t i=r.first; i<=r.second; ++i){
                            P[i] = depths[dup3[i]].top();
                            depths[dup3[i]].push(depth);
                            W[i] = wtd.rank(cst.rb(v)+1, dup3[i]) - 
                                   wtd.rank(cst.lb(v),dup3[i]);
                            if ( W[i] == 0 ){
                                cout << "ERROR: W["<<i<<"]=0"<< endl;
                                return;
                            }
                            W[i] = W[i] - 1; // store weight-1
                        }
                    }
                } else { // node visited the second time 
                    range_type r = map_node_to_dup3(v);
                    if ( !empty(r) ){
                         for(size_t i=r.first; i<=r.second; ++i){
                            depths[dup3[i]].pop();
                        }                       
                    }
                    --depth;
                }
            }
        }
        store_to_cache(P, surf::KEY_P, cc);
        store_to_cache(W, surf::KEY_WEIGHTS, cc);
        {
            wt_int<rrr_vector<63>> wtP;
            construct(wtP, cache_file_name(surf::KEY_P, cc));
            store_to_cache(wtP, surf::KEY_WTP, cc, true);
        }
        {
            wt_int<bit_vector,rank_support_v5<>,select_support_scan<>,select_support_scan<1>> wtP;
            construct(wtP, cache_file_name(surf::KEY_P, cc));
            store_to_cache(wtP, surf::KEY_WTP, cc, true);
       
        }
    }
    cout<<"...RMQ_C"<<endl;
    {
        int_vector<> C;
        load_from_cache(C, surf::KEY_C, cc);
        rmq_succinct_sct<> rmq_c(&C);
        store_to_cache(rmq_c, surf::KEY_RMQC, cc, true); 
    }
    cout<<"...W_AND_P"<<endl;
    {
        int_vector_buffer<> P_buf(cache_file_name(surf::KEY_P, cc));
        cout<<"P_buf.size()=" << P_buf.size() << endl;
        {
            int_vector<> id_v(P_buf.size(), 0,bits::hi(P_buf.size())+1);
            util::set_to_id(id_v);
            store_to_file(id_v, cache_file_name(surf::KEY_W_AND_P, cc)+".x");
        }
        {
            int_vector<> P;
            load_from_cache(P, surf::KEY_P, cc);
            store_to_file(P, cache_file_name(surf::KEY_W_AND_P,cc)+".y");
        }
        {
            int_vector<> W;
            load_from_cache(W, surf::KEY_WEIGHTS, cc);
            store_to_file(W, cache_file_name(surf::KEY_W_AND_P,cc)+".w");
        }
        cout<<"build k2treap"<<endl;
        k2treap_type k2treap;
        construct(k2treap, cache_file_name(surf::KEY_W_AND_P,cc));
        store_to_cache(k2treap, surf::KEY_W_AND_P, cc, true);
    }
/*    
    cout<<"...WEIGHT"<<endl;
    {
        int_vector<> W;
        load_from_cache(W, surf::KEY_WEIGHTS, cc);
        {
            //dac_vector<3,rrr_vector<63>> v(W);
            weight_type v(W);
            store_to_cache(v, surf::KEY_WEIGHTS, cc, true);
        }
        uint64_t n = W.size();
        uint64_t levels = sdsl::bits::hi(max_depth)+1;
        W.resize((levels-1)*W.size());
        for(size_t i=1; i<levels-1; ++i){
            for(size_t j=0; j<n; ++j){
                W[i*n+j] = W[j];
            }
        }
        rmq_succinct_sct<> rmq_w(&W);
        store_to_cache(rmq_w, surf::KEY_RMQW, cc, true); 
    }
*/
}

} // end namespace surf

#endif
