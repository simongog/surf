
#include <algorithm>
#include <vector>
#include <iostream>
#include <sstream>

#include "indri/KrovetzStemmer.hpp"

int main( int argc, char** argv ) {
    if(argc != 1) {
        std::cout << "USAGE: " << argv[0] << " < <input> > <output> " << std::endl;
        return EXIT_FAILURE;
    }

    using stemmer_t = indri::parse::KrovetzStemmer;
    stemmer_t ks;
    for(std::string line; std::getline(std::cin,line);) {
        auto id_sep_pos = line.find(';');
        auto qryid_str = line.substr(0,id_sep_pos);
        auto qry_id = std::stoull(qryid_str);
        std::istringstream qry_content_stream(line.substr(id_sep_pos+1));
        std::vector<std::string> stemmed_qry;
        for(std::string qry_token; std::getline(qry_content_stream,qry_token,' ');) {
            char stem_buf[stemmer_t::MAX_WORD_LENGTH+1] = {0};
            char original_word[stemmer_t::MAX_WORD_LENGTH+1] = {0};
            std::replace(qry_token.begin(),qry_token.end(),'-',' ');
            qry_token.erase(std::remove(qry_token.begin(),qry_token.end(),'\''),qry_token.end());
            qry_token.erase(std::remove(qry_token.begin(),qry_token.end(),'.'),qry_token.end());
            std::transform(qry_token.begin(), qry_token.end(), qry_token.begin(), ::tolower);
            std::copy(qry_token.begin(),qry_token.end(),std::begin(original_word));
            auto ret = ks.kstem_stem_tobuffer(original_word,stem_buf);
            if (ret > 0) {
                std::string tmp(stem_buf);
                stemmed_qry.push_back(tmp);
            } else {
                stemmed_qry.push_back(qry_token);
            }
        }
        std::cout << qry_id << ";";
        for(size_t i=0;i<stemmed_qry.size()-1;i++) {
            std::cout << stemmed_qry[i] << " ";
        }
        std::cout << stemmed_qry.back() << std::endl;
    }
}


