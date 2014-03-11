// extracts from an indri index a monoton sequence of integers in sdsl format
// which represent the parsed text collection.
#include <iostream>

#include "surf/config.hpp"
#include "surf/util.hpp"
#include "sdsl/int_vector_buffer.hpp"


int main( int argc, char** argv ) {
    if(argc != 3) {
        std::cout << "USAGE: " << argv[0] 
                  << " <string with # separated docs> <surf collection folder>" << std::endl;
        return EXIT_FAILURE;
    }
    std::string test_str = argv[1];
    std::string dir = argv[2];

    // setup collection directory
    if(surf::directory_exists(dir)) {
        std::cerr << "ERROR: collection directory already exists." << std::endl;
        return EXIT_FAILURE;
    }
    surf::create_directory(dir);

    if( test_str.back() != '#' ) {
        std::cerr << "ERROR: test string must end with doc seperator '#'" << std::endl;
        return EXIT_FAILURE;
    }
    
    // write collection string
    std::map<std::string::value_type,sdsl::int_vector<>::value_type> existing_syms;
    std::map<sdsl::int_vector<>::value_type,std::string::value_type> sym_mapping;
    sdsl::int_vector<> text_col(test_str.size()+1);
    size_t j=0;
    size_t num_docs = 0;
    for(const auto& sym : test_str) {
        if(sym == '#') {
            text_col[j++] = 1;
            num_docs++;
        } else {
            auto itr = existing_syms.find(sym);
            if(itr != existing_syms.end()) {
                text_col[j++] = itr->second;
            } else {
                sdsl::int_vector<>::value_type new_sym = existing_syms.size()+2;
                existing_syms[sym] = new_sym;
                sym_mapping[new_sym] = sym;
                text_col[j++] = new_sym;
            }
        }
    }
    text_col[j] = 0;
    std::ofstream ofs(dir+"/"+surf::TEXT_FILENAME);
    if(ofs.is_open()) {
        text_col.serialize(ofs);
    } else {
        std::cerr << "ERROR: could not write collection file." << std::endl;
        return EXIT_FAILURE;
    }

    // write the dict
    std::ofstream dict_ofs(dir+"/"+surf::DICT_FILENAME);
    if(dict_ofs.is_open()) {
        for(const auto& mapping : sym_mapping) {
            dict_ofs << mapping.second << " " << mapping.first << std::endl;
        }
    } else {
        std::cerr << "ERROR: could not write dictionary file." << std::endl;
        return EXIT_FAILURE;
    }

    // write docnames file
    std::ofstream docnames_ofs(dir+"/"+surf::DOCNAMES_FILENAME);
    if(docnames_ofs.is_open()) {
        for(size_t i=1;i<=num_docs;i++) {
            docnames_ofs << "DOCUMENT " << i << std::endl;
        }
    } else {
        std::cerr << "ERROR: could not write docnames file." << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "Created surf collection for string '" << test_str << "'" << std::endl;
    std::cout << "Found " << num_docs << " documents." << std::endl;
    std::cout << "Document delimiter = " << 1 << std::endl;
    std::cout << surf::TEXT_FILENAME << ": ";
    for(const auto& sym : text_col) {
        std::cout << sym << " ";
    }
    std::cout << std::endl;
    std::cout << "Mapping: ";
    for(const auto& mapping : sym_mapping) {
        std::cout << mapping.second << " -> " << mapping.first << "; ";
    }
    std::cout << std::endl;
    std::cout << "Document Names: ";
    for(size_t i=1;i<=num_docs;i++) {
        std::cout << "'DOCUMENT " << i << "'; ";
    }
    std::cout << std::endl;
}


