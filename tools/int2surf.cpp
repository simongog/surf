#include <sdsl/int_vector.hpp>
#include <iostream>

using namespace sdsl;
using namespace std;

int main(int argc, char* argv[]){
    if ( argc < 2 ){
        cout<<"Usage: "<<argv[0]<<" file_name"<<endl;
        cout<<" Takes an int_vector<> and appends the zero symbol; the results is"<<endl;
        cout<<" stored in `text_SURF.sdsl`"<<endl;
    }
    int_vector<> text;
    load_from_file(text, argv[1]);
    text.resize(text.size()+1);
    text[text.size()-1]=0;
    store_to_file(text, "text_int_SURF.sdsl");
}
