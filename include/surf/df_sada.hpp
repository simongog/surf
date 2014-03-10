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

using std::string;

namespace surf{

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

        typedef sdsl::cst_sct3<> cst_type;
        typedef sdsl::wt_int<>   wtc_type;
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

            cst_type temp_cst;
            wtc_type wtc;

            load_from_file(temp_cst, cache_file_name(surf::KEY_TMPCST, cc));
            load_from_file(wtc, cache_file_name(surf::KEY_WTC, cc));

            // int_vector_buffer which will contain the positions of the duplicates in the
            // C array after this scope
            int_vector_buffer<> temp_dup(cache_file_name(KEY_TMPDUP, cc), std::ios::out,
                                         1024*1024, bits::hi(wtc.size())+1);

            string d_file = cache_file_name(surf::KEY_DARRAY, cc);
            int_vector_buffer<> D(d_file);

            // construct the bv
            bit_vector h(2 * D.size(), 0);
            size_t h_idx = 0, dup_idx = 0;
            using n_type = std::tuple<cst_sct3<>::node_type, size_t, size_t, bool>;
            std::stack<n_type> s;
            s.emplace(temp_cst.root(), 1, temp_cst.degree(temp_cst.root()), true);
            // invariant: node has two children
            while (!s.empty()) {
                n_type node = s.top();
                s.pop();
                auto v = std::get<0>(node);
                auto l_child = std::get<1>(node);
                auto r_child = std::get<2>(node);
                auto first = std::get<3>(node);
                if (first) {  // first half
                    // recurse down
                    std::get<3>(node) = false;
                    s.push(node);
                    if (r_child == l_child + 1) {
                        auto w = temp_cst.select_child(v, l_child);
                        if (!temp_cst.is_leaf(w)) s.emplace(w, 1, temp_cst.degree(w), true);
                    } else {
                        auto mid = l_child + (r_child - l_child) / 2;
                        s.emplace(v, l_child, mid, true);
                    }
                } else {  // second half
                    auto lb = temp_cst.lb(temp_cst.select_child(v, l_child));
                    auto rb = temp_cst.rb(temp_cst.select_child(v, r_child));
                    auto mid = l_child + (r_child - l_child) / 2;
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
                    h[h_idx++] = 1;
                    if (mid + 1 == r_child) {
                        auto w = temp_cst.select_child(v, r_child);
                        if (!temp_cst.is_leaf(w)) s.emplace(w, 1, temp_cst.degree(w), true);
                    } else {
                        s.emplace(v, mid + 1, r_child, true);
                    }
                }
            }
            std::cerr<<"h_idx="<<h_idx<<std::endl;
            std::cerr<<"dup_idx="<<dup_idx<<std::endl;
            h.resize(h_idx);
            util::clear(temp_cst);
            // convert to proper bv type
            m_bv = bit_vector_type(h);
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
               sdsl::cache_config& cc, uint8_t num_bytes){
    using namespace sdsl;
    typedef df_sada<t_bv, t_sel, t_alphabet> df_sada_type;

    cout << "construct(df_sada)"<< endl;
    register_cache_file(sdsl::conf::KEY_TEXT_INT, cc);
    cout << "cache_file_name: "<< cache_file_name(sdsl::conf::KEY_TEXT_INT, cc) << endl;

    if (!cache_file_exists(sdsl::conf::KEY_SA, cc)) {
        construct_sa<t_alphabet::WIDTH>(cc);
    }
    register_cache_file(sdsl::conf::KEY_SA, cc);
    cout << "sa constructed"<< endl;

    if (!cache_file_exists(sdsl::conf::KEY_LCP, cc)) {
        if (t_alphabet::WIDTH == 8) {
            construct_lcp_semi_extern_PHI(cc);
        } else {
            construct_lcp_PHI<t_alphabet::WIDTH>(cc);
        }
    }
    register_cache_file(sdsl::conf::KEY_LCP, cc);
    using cst_type = typename df_sada<t_bv,t_sel,t_alphabet>::cst_type;
    cst_type temp_cst;
    if (!cache_file_exists(KEY_TMPCST, cc)) {
        auto event = memory_monitor::event("construct temp_cst");
        temp_cst = cst_type(cc, true);
        store_to_file(temp_cst, cache_file_name(surf::KEY_TMPCST, cc));
    } 

    if (!cache_file_exists(KEY_DOCCNT, cc)) {
        construct_doc_cnt<t_alphabet::WIDTH>(cc);
    }
    register_cache_file(KEY_DOCCNT, cc);
    uint64_t doc_cnt = 0;
    load_from_cache(doc_cnt, KEY_DOCCNT, cc);

    cout << "doc_cnt = " << doc_cnt << endl;

    if (!cache_file_exists(surf::KEY_DOCBORDER, cc)){
        construct_doc_border<t_alphabet::WIDTH>(cc);
    }
    if (!cache_file_exists(surf::KEY_DARRAY, cc)){
        construct_darray<t_alphabet::WIDTH>(cc);
    }

    typename df_sada_type::wtc_type wtc;
    string d_file = cache_file_name(surf::KEY_DARRAY, cc);
    int_vector_buffer<> D(d_file);
    cout<<"n="<<D.size()<<endl;
    if (!cache_file_exists(surf::KEY_WTC, cc)) {
        {
            auto event = memory_monitor::event("construct c");
            int_vector<> C(D.size(), 0, bits::hi(D.size()) + 1);
            int_vector<> last_occ(doc_cnt, D.size(), bits::hi(D.size()) + 1);
            for (size_t i = 0; i < D.size(); ++i) {
                uint64_t d = D[i];
                C[i] = last_occ[d];
                last_occ[d] = i;
            }
            util::bit_compress(C);
            store_to_file(C, cache_file_name(surf::KEY_C, cc));
        }
        auto event = memory_monitor::event("construct wt_c");
        construct(wtc, cache_file_name(surf::KEY_C, cc), cc, 0);
//        sdsl::remove(cache_file_name(surf::KEY_C, cc));
        store_to_file(wtc, cache_file_name(surf::KEY_WTC, cc));
    }
    cout << "call df_sada_type construct" << endl;

    if ( !cache_file_exists<df_sada_type>(surf::KEY_SADADF, cc) ) {
        df_sada_type tmp_sadadf(cc);
        store_to_cache(tmp_sadadf, surf::KEY_SADADF,cc, true);
    }

    if (!cache_file_exists(surf::KEY_WTDUP, cc)){
        auto event = memory_monitor::event("construct wtdup");
        int_vector<> D_array;
        load_from_file(D_array, d_file);
        int_vector_buffer<> wt_tmpdup(cache_file_name(surf::KEY_TMPDUP,cc));
        string wtdup_file = cache_file_name(surf::KEY_WTDUP, cc);
        int_vector_buffer<> wt_dup(wtdup_file, std::ios::out,
                                         1024*1024, D.width());
        for (size_t i = 0; i < wt_tmpdup.size(); ++i){
            wt_dup[i] = D_array[wt_tmpdup[i]];
        }
    }

    load_from_cache(idx, surf::KEY_SADADF, cc, true);
}


}

#endif
