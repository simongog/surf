#ifndef DF_SADA_HPP
#define DF_SADA_HPP

#include "config.hpp"
#include "construct_doc_cnt.hpp"
#include "construct_doc_border.hpp"
#include "construct_darray.hpp"
#include <sdsl/bit_vectors.hpp>
#include <sdsl/suffix_trees.hpp>
#include <tuple>
#include <string>
#include <algorithm>
#include <unordered_set>

using std::string;

namespace surf{

template<typename t_alphabet>
struct df_sada_trait{
    typedef sdsl::cst_sct3<sdsl::csa_wt<wt_int<rrr_vector<64>>>,sdsl::lcp_dac<>,sdsl::bp_support_sada<>,sdsl::bit_vector,sdsl::rank_support_v<>,sdsl::select_support_mcl<>> cst_type;
};

template<>
struct df_sada_trait<sdsl::byte_alphabet_tag>{
    typedef sdsl::cst_sct3<sdsl::csa_wt<wt_huff<rrr_vector<64>>>,sdsl::lcp_dac<>,sdsl::bp_support_sada<>,sdsl::bit_vector,sdsl::rank_support_v<>,sdsl::select_support_mcl<>> cst_type;
};

//! Constant time and 2n+o(n) size structure for document frequency
/*! 
 * \tparam t_bv       Bitvector type.
 * \tparam t_sel      Select structure for ones
 * \tparam t_alphabet The alphabet category. 
 *
 * This data structure was described in [1].
 *
 * \par Reference
 *  [1] Kunihiko Sadakane: ,,Succinct data structures for flexible
 *      text retrieval systems'', JDA 2007.
 *  [2] Simon Gog and Matthias Petri: TODO 2014.
 */
template<typename t_bv=sdsl::bit_vector,
         typename t_sel=typename t_bv::select_1_type,
         typename t_alphabet=sdsl::int_alphabet_tag>
class df_sada{
    public:
        typedef typename sdsl::int_vector<>::size_type size_type;
        typedef t_bv  bit_vector_type;
        typedef t_sel select_type;
        typedef t_alphabet alphabet_category;

        typedef typename df_sada_trait<t_alphabet>::cst_type cst_type;
        typedef sdsl::wt_int<sdsl::bit_vector,
                             sdsl::rank_support_v<>,
                             sdsl::select_support_scan<1>,
                             sdsl::select_support_scan<0>>   wtc_type;
    private:
        bit_vector_type m_bv;
        select_type     m_sel;

    public:    

        df_sada()=default;

        //! Constructor
        /*! \param cc cache_config which should contain the
         *         following files:
         *           - 
         */
        df_sada(sdsl::cache_config& cc){
            using namespace sdsl;
            auto event = memory_monitor::event("construct df_sada");

            if (cache_file_exists(KEY_H, cc)){
                bit_vector h;
                load_from_cache(h, KEY_H, cc);
                store_to_cache(h, KEY_H, cc);
                // convert to proper bv type
                m_bv = bit_vector_type(h);
                m_sel = select_type(&m_bv);
                return;
            }

            cst_type temp_cst;
            wtc_type wtc;

            load_from_file(temp_cst, cache_file_name<cst_type>(surf::KEY_TMPCST, cc));
            load_from_file(wtc, cache_file_name<wtc_type>(surf::KEY_WTC, cc));

            // int_vector_buffer which will contain the positions of the duplicates in the
            // C array after this scope
            int_vector_buffer<> temp_dup(cache_file_name(KEY_TMPDUP, cc), std::ios::out,
                                         1024*1024, bits::hi(wtc.size())+1);

            string d_file = cache_file_name(surf::KEY_DARRAY, cc);
            int_vector_buffer<> D(d_file);

            uint64_t doc_cnt = 0;
            load_from_cache(doc_cnt, KEY_DOCCNT, cc);
            cout << "doc_cnt = " << doc_cnt << endl;

            cout<<"begin calc splits"<<endl;
            bit_vector D_split(D.size()+1, 0);
            {
                std::vector<int64_t> last_seen(doc_cnt+1, -2);
                int64_t last_border = -1;
                for (size_t i=0; i<D.size(); ++i){
                    uint64_t x = D[i];
                    if ( last_seen[x] >= last_border ){
                        D_split[i-1] = 1;
                        last_border = i;
                    }
                    last_seen[x] = i;
                }
            }
            cout<<"end calc splits"<<endl;
            rank_support_v<> D_split_rank(&D_split);
            cout << "D_split.size()="<<D_split.size()<<endl;
            cout << "D_split_rank(D_split.size())="<<D_split_rank(D_split.size())<<endl;
            cout << "avg dist="<<D_split.size()/(D_split_rank(D_split.size())+1.0)<<endl;
            

            // construct the bv
            bit_vector h(2 * D.size(), 0);
            util::set_to_value(h,0);
            size_t h_idx = 0, dup_idx = 0;
            uint64_t max_depth = 0;
            //                      cst node, l_child, r_child, first, real, 
            using n_type = std::tuple<cst_sct3<>::node_type, size_t, size_t, bool, bool>;
            std::stack<n_type> s;
            s.emplace(temp_cst.root(), 1, temp_cst.degree(temp_cst.root()), true, true);
            // invariant: node has two children
            while (!s.empty()) {
                n_type node = s.top();
                s.pop();
                auto v = std::get<0>(node);
                auto l_child = std::get<1>(node);
                auto r_child = std::get<2>(node);
                auto first = std::get<3>(node);
                auto real  = std::get<4>(node);
                if (first) {  // first half
                    // recurse down
                    std::get<3>(node) = false;
                    auto lb = temp_cst.lb(temp_cst.select_child(v, l_child));
                    auto rb = temp_cst.rb(temp_cst.select_child(v, r_child));
                    if ( D_split_rank(rb) > D_split_rank(lb) ) {
                        s.push(node);
                        if (r_child == l_child + 1) {
                            auto w = temp_cst.select_child(v, l_child);
                            if (!temp_cst.is_leaf(w)) s.emplace(w, 1, temp_cst.degree(w), true, true);
                        } else {
                            auto mid = l_child + (r_child - l_child) / 2;
                            s.emplace(v, l_child, mid, true, false);
                        }
                    } else {
                        for (size_t i=0; i<rb-lb; ++i) {
                            h[h_idx++] = 1;
                        }
                    }
                } else {  // second half
                    if ( real and v != temp_cst.root() ){
                        uint64_t depth = temp_cst.depth(v);
                        max_depth = std::max(depth, max_depth);
                        std::stack<std::array<uint64_t,2>> s_child;
                        s_child.push({l_child, r_child});
                        while ( !s_child.empty() ){
                            auto child = s_child.top();
                            s_child.pop(); 
                            auto lb = temp_cst.lb(temp_cst.select_child(v, child[0]));
                            auto rb = temp_cst.rb(temp_cst.select_child(v, child[1]));
                            auto mid = child[0] + (child[1] - child[0]) / 2;
                            size_t dup_elements = 0;
                            if (lb + 1 == rb) {
                                dup_elements = (wtc[rb] == lb);
                                if (dup_elements) {
                                    temp_dup[dup_idx++] = lb;
                                }
                            } else {
                                auto mid_rb = temp_cst.rb(temp_cst.select_child(v, mid));
                                auto mid_lb = mid_rb + 1;
                                auto dup_info = restricted_unique_range_values(wtc, mid_lb, rb, lb, mid_rb);
                                dup_elements = dup_info.size();
                                for (auto& dup : dup_info) {
                                    temp_dup[dup_idx++] = dup;
                                }
                            }
                            h_idx += dup_elements;
                            if ( mid+1 < child[1])
                                s_child.push({mid+1,child[1]});
                            if ( child[0] < mid )
                                s_child.push({child[0],mid});
                        }
                    }
                    h[h_idx++] = 1;
                    auto mid = l_child + (r_child - l_child) / 2;
                    if (mid + 1 == r_child) {
                        auto w = temp_cst.select_child(v, r_child);
                        if (!temp_cst.is_leaf(w)) s.emplace(w, 1, temp_cst.degree(w), true, true);
                    } else {
                        s.emplace(v, mid + 1, r_child, true, false);
                    }
                }
            }
            std::cerr<<"max_depth="<<max_depth<<std::endl;
            store_to_cache(max_depth, surf::KEY_MAXCSTDEPTH, cc);
            std::cerr<<"h_idx="<<h_idx<<std::endl;
            std::cerr<<"dup_idx="<<dup_idx<<std::endl;
            h.resize(h_idx);
            store_to_cache(h, KEY_TMPH, cc);
            util::clear(temp_cst);
            // convert to proper bv type
            m_bv = bit_vector_type(h);
            if (m_bv.size()<40){
                std::cerr<<"m_bv="<<m_bv<<std::endl;
            }
            util::clear(wtc);
            m_sel = select_type(&m_bv);
        }


        //! Get the document frequency and range in the duplication array
        /*! \param sp Left bound of term interval.
         *  \param ep Right bound of term interval.
         *  \return A triple (term frequency, dub_sp, dub_ep), where
         *          [dup_sp,dup_ep] is the corresponding interval in the
         *          duplication array.
         */
        std::tuple<uint64_t,uint64_t,uint64_t>
        operator()(uint64_t sp, uint64_t ep) const{
            uint64_t dup = 0;
            uint64_t y = m_sel(ep);
            uint64_t dup_end = (y + 1) - ep - 1; // # of 0 left to y - 1
            uint64_t dup_begin = 0;          // # of 0 left to x
            if (0 == sp) {
                dup = (y + 1) - ep;  // (#all elements)-#1=#0
            } else {
                uint64_t x = m_sel(sp);
                dup = (y + 1) - ep - ((x + 1) - sp);
                dup_begin = (x+1) - sp;
            }
            return std::make_tuple(ep - sp + 1 - dup, dup_begin, dup_end);
        }

        size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = NULL, string name = "") const {
            using namespace sdsl;
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += m_bv.serialize(out, child, "m_bv");
            written_bytes += m_sel.serialize(out, child, "m_sel");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        void load(std::istream &in){
            m_bv.load(in);
            m_sel.load(in);
            m_sel.set_vector(&m_bv);
        }
};

template<typename t_bv, typename t_sel, typename t_alphabet>
void construct(df_sada<t_bv,t_sel,t_alphabet> &idx, const string& file,
               sdsl::cache_config& cc, uint8_t){
    using namespace sdsl;
    typedef df_sada<t_bv, t_sel, t_alphabet> df_sada_type;

    cout << "construct(df_sada)"<< endl;
    if (!cache_file_exists(conf::KEY_SA, cc)) {
        construct_sa<t_alphabet::WIDTH>(cc);
    }
    register_cache_file(conf::KEY_SA, cc);
    cout << "sa constructed"<< endl;

    if (!cache_file_exists(conf::KEY_LCP, cc)) {
        if (t_alphabet::WIDTH == 8) {
            cout<< "byte lcp construct"<<endl;
            construct_lcp_semi_extern_PHI(cc);
        } else {
            cout<< "int lcp construct"<<endl;
            construct_lcp_PHI<t_alphabet::WIDTH>(cc);
        }
    }
    register_cache_file(conf::KEY_LCP, cc);
    using cst_type = typename df_sada<t_bv,t_sel,t_alphabet>::cst_type;
    if (!cache_file_exists<cst_type>(KEY_TMPCST, cc)) {
        auto event = memory_monitor::event("construct temp_cst");
//        cst_type temp_cst = cst_type(cc, true);
        cst_type temp_cst = cst_type(cc);
        store_to_file(temp_cst, cache_file_name<cst_type>(surf::KEY_TMPCST, cc));
    } 

    construct_doc_cnt<t_alphabet::WIDTH>(cc);
    uint64_t doc_cnt = 0;
    load_from_cache(doc_cnt, KEY_DOCCNT, cc);

    cout << "doc_cnt = " << doc_cnt << endl;

    construct_darray<t_alphabet::WIDTH>(cc);


    string d_file = cache_file_name(surf::KEY_DARRAY, cc);
    int_vector_buffer<> D(d_file);
    cout<<"n="<<D.size()<<endl;
    if (!cache_file_exists(surf::KEY_C, cc)){
        auto event = memory_monitor::event("construct c");
        int_vector<> C(D.size(), 0, bits::hi(D.size()) + 1);
        int_vector<> last_occ(doc_cnt, D.size(), bits::hi(D.size()) + 1);
        for (size_t i = 0; i < D.size(); ++i) {
            uint64_t d = D[i];
            C[i] = last_occ[d];
            last_occ[d] = i;
        }
        
        if (D.size() < 20){
            cout<<"D=";
            for(size_t i=0; i<D.size(); ++i){
                cout<<" "<<D[i];
            }
        }
        cout<<endl;
        
        util::bit_compress(C);
        store_to_file(C, cache_file_name(surf::KEY_C, cc));
    }

    using wtc_type  = typename df_sada_type::wtc_type;
    if (!cache_file_exists<wtc_type>(surf::KEY_WTC, cc)) {
        wtc_type wtc;
        auto event = memory_monitor::event("construct wt_c");
        construct(wtc, cache_file_name(surf::KEY_C, cc), cc, 0);
        store_to_cache(wtc, surf::KEY_WTC, cc, true);
        cout<<"wtc.size()="<<wtc.size() << endl;
        cout<<"wtc.sigma()="<<wtc.sigma << endl;
    }
    cout << "call df_sada_type construct" << endl;

    if ( !cache_file_exists<df_sada_type>(surf::KEY_SADADF, cc) ) {
        df_sada_type tmp_sadadf(cc);
        store_to_cache(tmp_sadadf, surf::KEY_SADADF,cc, true);
    }

    if (!cache_file_exists(surf::KEY_DUP, cc)){
        cout<<"construct dup"<<endl;
        auto event = memory_monitor::event("construct dup");
        int_vector_buffer<> tmpdup(cache_file_name(surf::KEY_TMPDUP,cc));
        string dup_file = cache_file_name(surf::KEY_DUP, cc);
        string H_file = cache_file_name(surf::KEY_H, cc);
        int_vector_buffer<1> h_buf(H_file, std::ios::out, 1024*1024);
        cout << "tmpdup.size()="<<tmpdup.size()<<endl;
        cout << "D.width()="<<(int)D.width()<<endl;
        {
            int_vector_buffer<> dup_buf(dup_file, std::ios::out,
                                             1024*1024, D.width());
            int_vector<> D_array;
            load_from_file(D_array, d_file);
            for (size_t i = 0; i < tmpdup.size(); ++i){
                dup_buf[i] = D_array[tmpdup[i]];
            }
        }
        int_vector_buffer<1> hv_tmph(cache_file_name(surf::KEY_TMPH,cc));
        int_vector_buffer<>  dup_buf(dup_file, std::ios::in); 

        uint64_t doc_cnt = 0;
        load_from_cache(doc_cnt, KEY_DOCCNT, cc);
        cout << "doc_cnt = " << doc_cnt << endl;
        int_vector<> dup_vec(dup_buf.size(), 0, bits::hi(doc_cnt)+1);
        size_t k=0;
        for (size_t i=0, j=0; i < hv_tmph.size(); ++i){
            if (hv_tmph[i] == 0 ){
                std::map<uint64_t,uint64_t> dup_in_node;
                while ( i < hv_tmph.size() and hv_tmph[i] == 0 ){
                    dup_in_node[dup_buf[j++]]++;
                    ++i;
                }
                for (auto x : dup_in_node){
                    dup_vec[k++] = x.first;
                    h_buf.push_back(0);
                }
                if ( i < hv_tmph.size() ) {
                    h_buf.push_back(1);
                }
            } else {
                h_buf.push_back(1);
            }
        }
        dup_buf.close();
        h_buf.close();
        dup_vec.resize(k);
        cout << "dup_vec.size() = " << dup_vec.size() << endl;
        cout << "h_buf.size() = " << h_buf.size() << endl;
        store_to_cache(dup_vec, surf::KEY_DUP, cc);
    }
    load_from_cache(idx, surf::KEY_SADADF, cc, true);
}


}

#endif
