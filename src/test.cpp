#include "surf/df_sada.hpp"
#include "sdsl/construct.hpp"

using namespace surf;
using namespace sdsl;

int main(){
//    csa_wt<wt_int<>> csa;
    df_sada<> df;
    cache_config cc(false, "../collections/wikishort/index", "SURF");
    construct(df, "", cc, 0);
}
