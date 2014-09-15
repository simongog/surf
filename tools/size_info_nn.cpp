#include <iostream>
#include <string>
#include <map>

using namespace std;

int main(){
    string s;
    while ( cin >> s ) {
        if ( s == "spaceJSON" ){
            cin >> s; // read = 
            cin >> s; // read {
            break;
        }
    }
    string json="{";
    size_t depth = 1;
    map<string,string> key_val;
    while ( cin >> s ) {
        if ( s[0] == '{' )
            ++depth;
        else if ( s[0] == '}' )
            --depth;
        else {
            size_t found = s.find(":");
            if ( found != std::string::npos ) {
                if ( depth == 1 ) {
                    string key = s.substr(1, found-2);
                    string val = s.substr(found+2);
                    val = val.substr(0, val.find('"'));
                    if ( key != "children" ) {
                        key_val[key+"1"] = val;
                    }
                }
                if ( depth == 2 ) {
                    string key = s.substr(1, found-2);
                    string val = s.substr(found+2);
                    val = val.substr(0, val.find('"'));
                    if ( key != "children" )
                        key_val[key] = val;
                    if ( key == "size" ){
                        for (auto x : key_val){
                            cout << x.second << ";";
                        }
                        cout << endl;
                    }
                }
            }
        }
//        json += s;
        if ( depth == 0 )
            break;
    }
}
