#ifndef SURF_CONSTRUCT_U_HPP
#define SURF_CONSTRUCT_U_HPP

#include <sdsl/int_vector.hpp>
#include <type_traits>

namespace surf{

// generate the U array (= D^1 in the paper) and
// the KEY_UMARK bitvector
template<typename t_df>
void construct_u(sdsl::cache_config& cc)
{
    using namespace sdsl;
    using namespace std;
    static_assert(std::is_same<typename t_df::cst_type::index_category, sdsl::cst_tag>::value, "CST class expected");

    string u_file = cache_file_name(surf::KEY_U,cc);
    if (!cache_file_exists(surf::KEY_U,cc)){
        cout<<"......U does not exist. Generate it..."<<endl;
        {
            t_df df;
            construct(df, "", cc, 0); // make sure that the cst was generated
        }
        int_vector_buffer<> D_array(cache_file_name(surf::KEY_DUP, cc));
        cout<<".........D.size()="<<D_array.size()<<endl;
        cout<<".........D.width()="<<(int)D_array.width()<<endl;
        cout<<".........load cst"<<endl;
        using t_cst = typename t_df::cst_type;
        t_cst cst;
        load_from_file(cst, cache_file_name<t_cst>(surf::KEY_TMPCST, cc));
        cout<<".........cst.size()="<<cst.size()<<endl;
        string u_file = cache_file_name(surf::KEY_U, cc);
        int_vector_buffer<> U(u_file, std::ios::out,
                                   1024*1024, D_array.width());
        string umark_file = cache_file_name(surf::KEY_UMARK, cc);
        int_vector_buffer<1> Umark(umark_file, std::ios::out);

        uint64_t doc_cnt = 0;
        load_from_cache(doc_cnt, KEY_DOCCNT, cc);
        cout << ".........doc_cnt = " << doc_cnt << endl;

        std::vector<int64_t> last_occ(doc_cnt+1, -1);

        auto root = cst.root();
        for (auto& v : cst.children(root)){
            auto lb = cst.lb(v);
            auto rb = cst.rb(v);
            std::vector<uint64_t> buf;
            for (auto i = lb; i<=rb; ++i){
                auto x = D_array[i];
                if ( last_occ[x] < (int64_t)lb ){
                    buf.push_back(x);
                }
                last_occ[x] = i;
            }
            std::sort(buf.begin(), buf.end());
            for (size_t i=0; i<buf.size();++i){
                U.push_back(buf[i]);
                Umark.push_back(1);
            }
            for (size_t i=0; i < rb-lb+1-buf.size(); ++i){
                Umark.push_back(0);
            }
        }
    }
    cout << "U and Umark generated" << endl;
}

}// end namespace

#endif
