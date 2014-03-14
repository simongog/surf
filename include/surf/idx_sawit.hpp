#ifndef SURF_IDX_SAWIT_HPP
#define SURF_IDX_SAWIT_HPP

#include "sdsl/suffix_trees.hpp"
#include "surf/df_sada.hpp"
#include "surf/rank_functions.hpp"
#include <algorithm>

namespace surf{

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
    typedef t_df    df_type;
    typedef t_wtdup wtdup_type;
    typedef rank_bm25<> ranker_type;
private:
    csa_type    m_csa;
    wtd_type    m_wtd;
    df_type     m_df;
    wtdup_type  m_wtdup; 
    ranker_type m_r;
public:
    result_t search(std::vector<uint64_t> qry,size_t k) {
        result_t res;
//        typedef std::tuple<uint64_t, size_type, size_type, >
        auto qry_frq = unique_and_freq(qry.begin(), qry.end());

        // get multiplicity
        // unique
        for (size_t i=0; i<qry_frq.size(); ++i){
            size_type sp=1, ep=0;
            if ( backward_search(m_csa, 0, m_csa.size()-1, qry_frq[i].first, sp, ep) > 0 ) {
                cout << "interval of " << qry[i] << " ["
                     << sp << "," << ep << "]" << endl;
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
