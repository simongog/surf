#include <iostream>
#include "sdsl/int_vector_buffer.hpp"

using namespace sdsl;
using namespace std;

int main(int argc, char* argv[]){
    if ( argc < 2 ){
        cout << "./" << argv[0] << " file" << endl;
        cout << "file has to contain a serialized sdsl::int_vector<>" << endl;
        cout << "Program outputs the size of elements and the width per element" << endl;
        return 1;
    }
    int_vector_buffer<> ivb(argv[1]);
    cout << ivb.size() << " " << (int) ivb.width() << endl;
}
