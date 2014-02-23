#include "surf/df_sada.hpp"
#include "sdsl/construct.hpp"

using namespace surf;
using namespace sdsl;

int main(){
    
    
    csa_wt<wt_int<>> csa;

    df_sada<> df;
    construct(df, "text.txt");
}
