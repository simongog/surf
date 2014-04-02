#ifndef SURF_PHRASE_PARSER_HPP
#define SURF_PHRASE_PARSER_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <ratio>
#include <chrono>

#include "surf/config.hpp"
#include "surf/query.hpp"

namespace surf{

struct phrase_parser {
    phrase_parser() = delete;

    template<class t_csa>
    static query_t phrase_segmentation(t_csa& csa,
    						const std::vector<uint64_t>& query_ids,
    						const std::unordered_map<uint64_t,std::string>& reverse_mapping,
                            double threshold)
    {
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
                // if we start at a very frequent word, a phrase can't start 
                // there.
                auto single_cnt = sdsl::count(csa,query_ids.begin()+start,query_ids.begin()+start+1);
                if( single_cnt * 100 > csa.size() ) {
                    break;
                }

    			auto cnt = sdsl::count(csa,query_ids.begin()+start,query_ids.begin()+i+1);
    			double prob = (double)cnt / (double)csa.size();

    			// single
    			double single = P_single[i];
    			for(size_t l=start;l<=i;l++) single *= P_single[l];

    			// calc ratio
    			double assoc_ratio = log(prob)-log(single);

                // debug
                /*
                std::cout << "SCORE(";            
                for(size_t l=start;l<=i;l++) {
                    auto id = query_ids[l];
                    auto stritr = reverse_mapping.find(id);
                    std::cout << stritr->second << " ";
                }
                std::cout << ") -> " << assoc_ratio << std::endl;
                */

    			if(assoc_ratio < threshold) {
    				// not a phrase. if the prev one was a phrase we use it
    				if(phrase_found) {
    					std::vector<uint64_t> phrase;
    					for(size_t j=start;j<i;j++) {
    						phrase.push_back(query_ids[j]);
    					}
    					phrases.push_back(phrase);
    				    phrase_added = true;
    				    start = i;
    				    break;
    				}
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
    				}
    				phrases.push_back(phrase);
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
};
}// end namespace

#endif
