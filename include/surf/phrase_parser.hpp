#ifndef SURF_PHRASE_PARSER_HPP
#define SURF_PHRASE_PARSER_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <ratio>

#include "surf/config.hpp"
#include "surf/query.hpp"

namespace surf{

template<class t_thres = std::ratio<1,2>>
struct phrase_parser {
    phrase_parser() = delete;

    template<class t_csa>
    static query_t phrase_segmentation(t_csa& csa,
    						const std::vector<uint64_t>& query_ids,
    						const std::unordered_map<uint64_t,std::string> reverse_mapping)
    {
    	double threshold = t_thres::num / t_thres::den;

    	//compute single term probabilities
    	std::vector<double> P_single;
    	for(size_t i=0;i<query_ids.size();i++) {
    		auto cnt = sdsl::count(csa,query_ids.begin()+i,query_ids.begin()+i+1);
    		double prob = (double)cnt / (double)csa.size();
    		P_single.push_back(prob);
    	}

    	//compute all probabilities
    	std::vector<std::vector<uint64_t>> phrases;
    	size_t start = 0;
    	size_t stop = query_ids.size();
    	while(start < stop) {
    		bool phrase_found = false;
    		bool phrase_added = false;
    		for(size_t i=start+1;i<stop;i++) {
    			auto cnt = sdsl::count(csa,query_ids.begin()+start,query_ids.begin()+i+1);
    			double prob = (double)cnt / (double)csa.size();

    			// single
    			double single = P_single[i];
    			for(size_t l=start;l<=i;l++) single *= P_single[l];

    			// calc ratio
    			double assoc_ratio = log(prob)-log(single);

    			if(assoc_ratio < threshold) {
    				// not a phrase. if the prev one was a phrase we use it
    				if(phrase_found) {
    					std::vector<uint64_t> phrase;
    					for(size_t j=start;j<i;j++) {
    						phrase.push_back(query_ids[j]);
    					}
    					phrases.push_back(phrase);
    				}
    				phrase_added = true;
    				start = i;
    				break;
    			} else {
    				// still a phrase. continue!
    				phrase_found = true;
    			}
    		}
    		if(!phrase_added) {
    			if(phrase_found) {
    				// we found a phrase that goes to the end of the id list
    				std::vector<uint64_t> phrase;
    				for(size_t i=start;i<stop;i++) {
    					phrase.push_back(query_ids[i]);
    					phrases.push_back(phrase);
    				}
    				start = stop;
    			} else {
    				// for this term we have not found any phrase 
    				// with it. add it as a single
    				std::vector<uint64_t> single;
    				single.push_back(query_ids[start]);
    				phrases.push_back(single);
    				start++;
    			}
    		}
    	}

    	// check if all phrases are uniq
    	query_t q;
    	auto itr = phrases.begin();
    	while(itr != phrases.end()) {
    		auto cur_list = *itr;
    		uint64_t num_equal = 0;
    		auto next = itr+1;
    		while(next != phrases.end()) {
    			auto next_list = *next;
    			if(std::equal(cur_list.begin(),cur_list.end(),next_list.begin())) {
    				num_equal++;
    				next = phrases.erase(next);
    			} else {
    				next++;
    			}
    		}

    		/* get the string representation */
    		std::vector<std::string> qry_str;
    		for(const auto& id : cur_list) {
    			auto rmitr = reverse_mapping.find(id);
    			qry_str.push_back(rmitr->second);
    		}
    		std::get<1>(q).emplace_back(*itr,qry_str,num_equal);
    		itr++;
    	}
    	return q;
    }

    static std::vector<query_t> parse_queries(sdsl::cache_config& cc,
    										  const std::string& collection_dir,
                                              const std::string& query_file) 
    {
        std::vector<query_t> queries;

        /* load the mapping */
        auto mapping = query_parser::load_dictionary(collection_dir);
        const auto& id_mapping = mapping.first;
        const auto& reverse_mapping = mapping.second;

        /* load csa */
        using csa_type = sdsl::csa_wt<sdsl::wt_int<sdsl::rrr_vector<63>>,1000000,1000000>;
        csa_type csa;
        load_from_cache(csa, surf::KEY_CSA, cc, true);

        /* parse queries */
        std::ifstream qfs(query_file); 
        if(!qfs.is_open()) {
            std::cerr << "cannot load query file.";
            exit(EXIT_FAILURE);
        }

        std::string query_str;
        while( std::getline(qfs,query_str) ) {
        	auto qry_mapping = query_parser::map_to_ids(id_mapping,query_str,true);
        	if( std::get<0>(qry_mapping)) {
        		auto qry_ids = std::get<1>(qry_mapping);
        		auto parsed_qry = phrase_segmentation(csa,qry_ids,reverse_mapping);
        		queries.emplace_back(parsed_qry);
        	}
        }

        return queries;
    }
};

}// end namespace

#endif
