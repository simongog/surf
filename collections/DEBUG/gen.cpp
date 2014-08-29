#include <iostream>
#include <string>
#include <sdsl/int_vector.hpp>

using namespace std;
using namespace sdsl;

int main(){
    string s = "LOL\1OLLL\1OOL\1";
    int_vector<8> v(s.size()+1,0);
    cout<<v.size()<<endl;
    for (auto i=0; i<s.size();++i)
        v[i]=s[i];
    store_to_file(v,"text_SURF.sdsl");
}
