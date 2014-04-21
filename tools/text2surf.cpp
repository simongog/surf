#include <sdsl/int_vector.hpp>
#include <iostream>

using namespace sdsl;
using namespace std;

int main(int argc, char* argv[]){
    if ( argc < 2 ){
        cout<<"Usage: "<<argv[0]<<" file_name"<<endl;
        cout<<" Takes byte character text and generates int_vector<8>"<<endl;
        cout<<" stored in `text_SURF.sdsl`"<<endl;
    }
    int_vector<8> text;
    load_vector_from_file(text, argv[1], 1);
    text.resize(text.size()+1);
    text[text.size()-1]=0;
    store_to_file(text, "text_SURF.sdsl");
}
