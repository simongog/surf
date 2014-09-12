#ifndef DF_SADA_HPP
#define DF_SADA_HPP

#include "config.hpp"
#include "construct_doc_cnt.hpp"
#include "construct_doc_border.hpp"
#include "construct_darray.hpp"
#include "surf/construct_max_doc_len.hpp"
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
    typedef sdsl::cst_sct3<sdsl::csa_wt<sdsl::wt_int<sdsl::rrr_vector<64>>>,sdsl::lcp_dac<>,sdsl::bp_support_sada<>,sdsl::bit_vector,sdsl::rank_support_v<>,sdsl::select_support_mcl<>> cst_type;
};

template<>
struct df_sada_trait<sdsl::byte_alphabet_tag>{
    typedef sdsl::cst_sct3<sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<64>>>,sdsl::lcp_dac<>,sdsl::bp_support_sada<>,sdsl::bit_vector,sdsl::rank_support_v<>,sdsl::select_support_mcl<>> cst_type;
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


            construct_max_doc_len<alphabet_category::WIDTH>(cc);    
            uint64_t max_len = 0;
            load_from_cache(max_len, surf::KEY_MAXDOCLEN, cc);

            cst_type cst;

            load_from_file(cst, cache_file_name<cst_type>(surf::KEY_TMPCST, cc));

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

            std::string DUP_file = cache_file_name(surf::KEY_DUP, cc);
            std::string W_file = cache_file_name(surf::KEY_WEIGHTS, cc);

            int_vector_buffer<> dup_buf(DUP_file, std::ios::out, 1<<20, sdsl::bits::hi(doc_cnt)+1);
            int_vector_buffer<> weight_buf(W_file, std::ios::out, 1<<20, sdsl::bits::hi(max_len)+1);


            typedef WTD_TYPE t_wtd;
            t_wtd wtd;
            load_from_cache(wtd, surf::KEY_WTD, cc, true);

            // construct the bv
            bit_vector h(2 * D.size(), 0);
            util::set_to_value(h,0);
            size_t h_idx = 0, dup_idx = 0;
            size_t last_io_id = 0;
            uint64_t max_depth = 0;

            using node_type = typename cst_type::node_type;
            using n_type = std::tuple<node_type, bool>;
            std::stack<n_type> s;
            std::vector<node_type> child_vec;
            std::vector<range_type> range_vec;
            auto s_push = [&s,&cst,&D_split_rank](node_type v){
                if (D_split_rank(cst.rb(v)) > D_split_rank(cst.lb(v))){
                    s.emplace(v, true);
                }
            };
            s_push(cst.root());
            // invariant: node has two children
            while (!s.empty()) {
                n_type node = s.top();
                s.pop();
                auto v = std::get<0>(node);
                auto first = std::get<1>(node);
                if ( first ) {  // first half
                    // recurse down
                    std::get<1>(node) = false;
                    uint64_t depth = cst.depth(v);
                    max_depth = std::max(depth, max_depth);
                    s.push(node);
                    s_push(cst.select_child(v, 1));
                } else {  // second half
                    for (auto& child : cst.children(v)){
                        child_vec.push_back(child);
                    }
                    uint64_t node_io_id = cst.rb(cst.select_child(v, 1));
                    while (last_io_id+1 < node_io_id){
                        ++last_io_id;
                        h[h_idx++] = 1;
                    }
                    if ( v != cst.root() ){
                        for(auto& child : child_vec){
                            range_vec.emplace_back(cst.lb(child), cst.rb(child));
                        }
                        auto dups = intersect(wtd, range_vec, 2);
                        range_vec.clear();
                        for (auto &duplicate : dups) {
                            dup_buf[dup_idx] = duplicate.first;
                            weight_buf[dup_idx] = duplicate.second-1; 
                            ++dup_idx;
                            ++h_idx;
                        }
                    }
                    h[h_idx++] = 1;
                    last_io_id = node_io_id;
                    while ( child_vec.size() > 1 ){
                        s_push(child_vec.back());
                        child_vec.pop_back();
                    }
                    child_vec.pop_back();
                }
            }
            std::cerr<<"done last_io_id="<<last_io_id<<std::endl;
            while (last_io_id < wtd.size()){
                ++last_io_id;
                h[h_idx++] = 1;
            }
            std::cerr<<"max_depth="<<max_depth<<std::endl;
            store_to_cache(max_depth, surf::KEY_MAXCSTDEPTH, cc);
            std::cerr<<"h_idx="<<h_idx<<std::endl;
            std::cerr<<"dup_idx="<<dup_idx<<std::endl;
            h.resize(h_idx);
            store_to_cache(h, KEY_H, cc);
            util::clear(cst);
            // convert to proper bv type
            m_bv = bit_vector_type(h);
            if (m_bv.size()<40){
                std::cerr<<"m_bv="<<m_bv<<std::endl;
            }
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
        auto event = memory_monitor::event("construct cst");
        cst_type cst = cst_type(cc);
        store_to_file(cst, cache_file_name<cst_type>(surf::KEY_TMPCST, cc));
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
    typedef WTD_TYPE t_wtd;
    if (!cache_file_exists<t_wtd>(surf::KEY_WTD, cc) ){
        construct_darray<t_alphabet::WIDTH>(cc, false);
        t_wtd wtd;
        construct(wtd, cache_file_name(surf::KEY_DARRAY, cc), cc);
        cout << "wtd.size() = " << wtd.size() << endl;
        cout << "wtd.sigma = " << wtd.sigma << endl;
        store_to_cache(wtd, surf::KEY_WTD, cc, true);
    }

    cout << "call df_sada_type construct" << endl;
    if ( !cache_file_exists<df_sada_type>(surf::KEY_SADADF, cc) ) {
        df_sada_type tmp_sadadf(cc);
        store_to_cache(tmp_sadadf, surf::KEY_SADADF,cc, true);
    }
    load_from_cache(idx, surf::KEY_SADADF, cc, true);
}


}

#endif
