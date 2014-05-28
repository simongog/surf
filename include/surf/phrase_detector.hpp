#ifndef SURF_PHRASE_DETECTOR_HPP
#define SURF_PHRASE_DETECTOR_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <ratio>
#include <chrono>
#include <tuple>

#include "surf/config.hpp"
#include "surf/query.hpp"

using parsed_qry = std::vector<std::vector<uint64_t>>;

namespace surf{

struct phrase_detector {
    phrase_detector() = delete;

    template<class t_csa>
    static parsed_qry parse_greedy_lr(t_csa& csa,
    								  const std::vector<uint64_t>& qry,
    								  double threshold)
    {
    	//compute single term probabilities
    	std::vector<double> P_single;
    	for(size_t i=0;i<qry.size();i++) {
    		auto cnt = sdsl::count(csa,qry.begin()+i,qry.begin()+i+1);
    		double prob = (double)cnt / (double)csa.size();
    		P_single.push_back(prob);
    	}

    	//compute all probabilities
    	parsed_qry phrases;
    	size_t start = 0;
    	size_t stop = qry.size();
    	while(start < stop) {
    		bool phrase_found = false;
    		bool phrase_added = false;
    		for(size_t i=start+1;i<stop;i++) {
    			auto cnt = sdsl::count(csa,qry.begin()+start,qry.begin()+i+1);
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
    						phrase.push_back(qry[j]);
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
    					phrase.push_back(qry[i]);
    				}
    				phrases.push_back(phrase);
    				start = stop;
    			} else {
    				// for this term we have not found any phrase 
    				// with it. add it as a single
    				std::vector<uint64_t> single;
    				single.push_back(qry[start]);
    				phrases.push_back(single);
    				start++;
    			}
    		}
    	}
    	return phrases;
    }


    template<class t_csa>
    static parsed_qry parse_greedy_paul(t_csa& csa,
    								    const std::vector<uint64_t>& qry,
    								    double threshold)
    {
    	parsed_qry phrases;

    	//compute single term probabilities
    	std::vector<double> P_single;
    	for(size_t i=0;i<qry.size();i++) {
    		auto cnt = sdsl::count(csa,qry.begin()+i,qry.begin()+i+1);
    		double prob = (double)cnt / (double)csa.size();
    		P_single.push_back(prob);
    	}

    	// compute adjacent pair probabilities
    	std::priority_queue<std::tuple<double,uint64_t,uint64_t>> assoc_pairs;
    	for(size_t i=0;i<qry.size()-1;i++) {
    		auto cnt = sdsl::count(csa,qry.begin()+i,qry.begin()+i+2);
    		double prob = (double)cnt / (double)csa.size();
			// single
			double single = P_single[i] * P_single[i+1];
			// calc ratio
			double assoc_ratio = log(prob)-log(single);
            //std::cout << "AR(" << i << "," << i+1 << ") = " << assoc_ratio << std::endl;
    		assoc_pairs.push(std::make_tuple(assoc_ratio,i,i+1));
    	}

    	std::unordered_set<uint64_t> used;
    	while(!assoc_pairs.empty()) {
    		auto cur = assoc_pairs.top(); assoc_pairs.pop();
            //std::cout << "dequeue " << std::get<0>(cur) << std::endl;
    		if( std::get<0>(cur) > threshold ) {
    			// check if we use one of the terms already
    			if( used.count( std::get<1>(cur)) == 0 && 
    				used.count(std::get<2>(cur)) == 0) 
  				{
  					used.insert(std::get<1>(cur));
  					used.insert(std::get<2>(cur));
  					std::vector<uint64_t> new_phrase;
  					new_phrase.push_back(qry[std::get<1>(cur)]);
  					new_phrase.push_back(qry[std::get<2>(cur)]);
  					phrases.push_back(new_phrase);
    			}
    		} else {
    			// no more phrases above threshold
    			break;
    		}
    	}

    	// add all singletons we have not used
    	for(size_t i=0;i<qry.size();i++) {
			if(used.count(i) == 0) {
				std::vector<uint64_t> singleton;
				singleton.push_back(qry[i]);
				phrases.push_back(singleton);
			}	
    	}

    	return phrases;
    }


    template<class t_csa>
    static parsed_qry parse_dp(t_csa& csa,
    								  const std::vector<uint64_t>& qry,
    								  double threshold)
    {
    	parsed_qry phrases;

    	return phrases;
    }

    template<class t_csa>
    static parsed_qry parse_none(t_csa& csa,
    								  const std::vector<uint64_t>& qry,
    								  double threshold)
    {
    	parsed_qry phrases;
    	for(const auto& id : qry) {
    		std::vector<uint64_t> single;
    		single.push_back(id);
    		phrases.push_back(single);
    	}
    	return phrases;
    }
};
}// end namespace

#endif
