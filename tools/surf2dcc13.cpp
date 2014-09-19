#include <sdsl/int_vector.hpp>
#include <iostream>
#include <string>

using namespace sdsl;
using namespace std;

int main(int argc, char* argv[]){
    if ( argc < 4 ){
        cout << "Usage: ./" << argv[0] <<" surf_TEXT.sdsl" <<" existing_output_dir file_prefix" << endl;
        return 1;
    }
    string prefix = argv[3];
    string dir = argv[2];
    int_vector<8> text;
    load_from_file(text, argv[1]);
    size_t doc_id = 0;
    ofstream out(dir+"/"+prefix+to_string(doc_id)+".txt");
    for (size_t i=0; i<text.size(); ++i){
        if ( text[i]==1 ){
            cerr<<"Close file"<<endl;
            out.close();
            ++doc_id;
            continue;
        } else if (text[i-1]==1 and text[i]!=0){
            string name = dir+"/"+prefix+to_string(doc_id)+".txt";
            cerr<<"Open file "<<name<<endl;
            out.open(name);
        }
        if (text[i]!=0){
            if (text[i]==255){
                char x = '\1';
                out.write(&x, 1);
            } else {
                char x = text[i];
                out.write(&x,1); 
            }
        }
    }
}
