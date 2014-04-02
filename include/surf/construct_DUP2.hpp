#ifndef SURF_CONSTRUCT_DUP2_HPP
#define SURF_CONSTRUCT_DUP2_HPP

#include <sdsl/int_vector.hpp>

namespace surf{

// generate the DUP2 array (= R in the paper) and
// the KEY_DUPMARK bitvector
template<typename t_df>
void construct_dup2(sdsl::cache_config& cc)
{
    using namespace sdsl;
    using namespace std;

    string dup2_file = cache_file_name(surf::KEY_DUP2,cc);
    if (!cache_file_exists(surf::KEY_DUP2,cc)){
        cout<<"......dup2 does not exist. Generate it..."<<endl;
        {
            t_df df;
            construct(df, "", cc, 0); // make sure that surf::KEY_DUP was generated 
        }
        int_vector_buffer<> dup(cache_file_name(surf::KEY_DUP, cc));
        cout<<".........dup.size()="<<dup.size()<<endl;
        cout<<".........dup.width()="<<(int)dup.width()<<endl;
        {
            bit_vector dup_mark(dup.size(), 1);
            store_to_cache(dup_mark, surf::KEY_DUPMARK, cc);
        }
        int_vector_buffer<1> dup_mark(cache_file_name(surf::KEY_DUPMARK, cc));
        int_vector_buffer<> dup2(dup2_file, std::ios::out, 1024*1024,
                                 dup.width());
        cout<<".........load df"<<endl;
        t_df df;
        load_from_cache(df, surf::KEY_SADADF, cc, true);
        cout<<".........load cst"<<endl;
        using cst_type =  typename t_df::cst_type;
        cst_type cst;
        load_from_file(cst, cache_file_name<cst_type>(surf::KEY_TMPCST, cc));
        cout<<".........cst.size()="<<cst.size()<<endl;
        auto root = cst.root();
        auto left_most_leaf = cst.select_leaf(1);
        uint64_t next_idx = 0;
        for (auto& v : cst.children(root)){
            if ( v == left_most_leaf )
                continue;
            auto lb = cst.lb(v);
            auto rb = cst.rb(v);
            auto df_info = df(lb, rb);
            std::vector<uint64_t> buf;
            for (uint64_t i = std::get<1>(df_info); i <= std::get<2>(df_info); ++i) {
                buf.push_back(dup[i]); 
            }
            for (uint64_t i = next_idx; i < std::get<1>(df_info); ++i){
                dup_mark[i] = 0;
            }
            next_idx = std::get<2>(df_info)+1;
            for (size_t i=0; i < buf.size(); ++i){
                dup2.push_back(buf[i]);
            }
        }
        for (uint64_t i = next_idx; i < dup_mark.size(); ++i){
            dup_mark[i]=0;
        }
    }
}

}// end namespace

#endif
