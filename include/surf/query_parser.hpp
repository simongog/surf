#ifndef SURF_QUERY_PARSER_HPP
#define SURF_QUERY_PARSER_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>

#include "surf/config.hpp"
#include "surf/query.hpp"

namespace surf{

struct query_parser {
    query_parser() = delete;
    using mapping_t = std::pair<std::unordered_map<std::string,uint64_t>,
                     std::unordered_map<uint64_t,std::string>
                     >;

    static mapping_t
         load_dictionary(const std::string& collection_dir)
    {
        std::unordered_map<std::string,uint64_t> id_mapping;
        std::unordered_map<uint64_t,std::string> reverse_id_mapping;
        {
            auto dict_file = collection_dir + "/" + surf::DICT_FILENAME;
            std::ifstream dfs(dict_file);
            if(!dfs.is_open()) {
                std::cerr << "cannot load dictionary file.";
                exit(EXIT_FAILURE);
            }
            std::string term_mapping;
            while( std::getline(dfs,term_mapping) ) {
                auto sep_pos = term_mapping.find(' ');
                auto term = term_mapping.substr(0,sep_pos);
                auto idstr = term_mapping.substr(sep_pos+1);
                uint64_t id = std::stoull(idstr);
                id_mapping[term] = id;
                reverse_id_mapping[id] = term;
            }
        }
        return {id_mapping,reverse_id_mapping};
    }

    static std::tuple<bool,std::vector<uint64_t>> 
        map_to_ids(const std::unordered_map<std::string,uint64_t>& id_mapping,
                   std::string query_str,bool only_complete)
    {
        std::vector<uint64_t> ids;
        std::istringstream qry_content_stream(query_str);
        for(std::string qry_token; std::getline(qry_content_stream,qry_token,' ');) {
            auto id_itr = id_mapping.find(qry_token);
            if(id_itr != id_mapping.end()) {
                ids.push_back(id_itr->second);
            } else {
                std::cerr << "ERROR: could not find '" << qry_token << "' in the dictionary." << std::endl;
                if(only_complete) {
                    return std::make_tuple(false,ids);
                }
            }
        }
        return std::make_tuple(true,ids);
    }

    static std::pair<bool,query_t> parse_query(const mapping_t& mapping,
                const std::string& query_str,bool only_complete = false)
    {
        auto id_sep_pos = query_str.find(';');
        auto qryid_str = query_str.substr(0,id_sep_pos);
        auto qry_id = std::stoull(qryid_str);
        auto qry_content = query_str.substr(id_sep_pos+1);
        const auto& id_mapping = mapping.first;
        const auto& reverse_mapping = mapping.second;
        auto mapped_qry = map_to_ids(id_mapping,qry_content,only_complete);

        bool parse_ok = std::get<0>(mapped_qry);
        if(parse_ok) {
            std::unordered_map<uint64_t,uint64_t> qry_set;
            const auto& qids = std::get<1>(mapped_qry);
            for(const auto& qid : qids) {
                qry_set[qid] += 1;
            }
            std::vector<query_token> query_tokens;
            for(const auto& qry_tok : qry_set) {
                std::vector<uint64_t> term;
                term.push_back(qry_tok.first);
                auto rmitr = reverse_mapping.find(qry_tok.first);
                std::string qry_str = rmitr->second;
                std::vector<std::string> term_str;
                term_str.push_back(qry_str);
                query_tokens.emplace_back(term,term_str,qry_tok.second);
            }
            query_t q(qry_id,query_tokens);
            return {true,q};
        }

        // error
        query_t q;
        return {false,q};
    }

    static std::vector<query_t> parse_queries(const std::string& collection_dir,
                                            const std::string& query_file,bool only_complete = false) 
    {
        std::vector<query_t> queries;

        /* load the mapping */
        auto mapping = load_dictionary(collection_dir);

        /* parse queries */
        std::ifstream qfs(query_file); 
        if(!qfs.is_open()) {
            std::cerr << "cannot load query file.";
            exit(EXIT_FAILURE);
        }

        std::string query_str;
        while( std::getline(qfs,query_str) ) {
            auto parsed_qry = parse_query(mapping,query_str,only_complete);
            if(parsed_qry.first) {
                queries.emplace_back(parsed_qry.second);
            }
        }

        return queries;
    }
};

}// end namespace

#endif
