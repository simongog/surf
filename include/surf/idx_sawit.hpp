#ifndef SURF_IDX_SAWIT_HPP
#define SURF_IDX_SAWIT_HPP

#include "sdsl/suffix_trees.hpp"
#include "surf/df_sada.hpp"

namespace surf{

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
