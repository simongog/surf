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
using result_t = std::vector<std::pair<double,std::vector<uint64_t>>>;

namespace surf{

double prob(double x){
    return x==0 ? -99999999 : log(x);
}

template<size_t t_min_freq = 0,bool t_length_norm = false>
struct phrase_detector_sa_greedy {
    phrase_detector_sa_greedy() = delete;

    static std::string name() { 
        if(t_length_norm) return "SA-GREEDY-LN";
        return "SA-GREEDY"; 
    }

    template<class t_index>
    static result_t parse(t_index& index,const std::vector<uint64_t>& qry,double threshold)
    {
        result_t phrases;

    	//compute single term probabilities
    	std::vector<double> P_single;
    	for(size_t i=0;i<qry.size();i++) {
    		auto cnt = index.csa_count(qry.begin()+i,qry.begin()+i+1);
    		double single = (double)cnt / (double)index.csa_size();
    		P_single.push_back(single);
    	}

        // compute all phrases
        for (auto start = qry.begin(); start < qry.end(); ){
            double single = prob(P_single[start-qry.begin()]);
            auto end = start;
            double assoc_ratio = 0;
            do {
                ++end;
                if ( end == qry.end() )
                    break;
                single += prob(P_single[end-qry.begin()]);
    			auto cnt = index.csa_count(start, end+1);
                if(cnt < t_min_freq) {
                    assoc_ratio = 0;
                } else {
                    double joint = prob((double)cnt/index.csa_size());
                    assoc_ratio = joint-single;
                }
            } while ( assoc_ratio >= threshold );
            phrases.emplace_back(assoc_ratio,std::vector<uint64_t>(start, end+1));
            start = end;
        }

        return phrases;
    }
};

template<class t_index,class t_itr>
static double compute_x2(const std::vector<double>& freq_single,t_index& index,size_t i,t_itr begin,t_itr end)
{
    auto square = [](double a){ return a*a; };
    double N = index.csa_size();
    // expected freq of q_i,q_i+1
    double freq_qiqi1 = index.csa_count(begin,end);
    double freq_qinotqi1 = freq_single[i] - freq_qiqi1;
    double freq_notqiqi1 = freq_single[i+1] - freq_qiqi1;
    double freq_notqinotqi1 = N - freq_qinotqi1;

    double exp_freq_qiqi1 = (freq_single[i]*freq_single[i+1])/N;
    double exp_freq_qinotqi1 = (freq_single[i]*(freq_notqiqi1+freq_notqinotqi1))/N;
    double exp_freq_notqiqi1 = (freq_single[i+1]*(freq_notqiqi1+freq_notqinotqi1))/N;
    double exp_freq_notqinotqi1 = ((freq_qinotqi1+freq_notqinotqi1)*
                                   (freq_notqiqi1+freq_notqinotqi1))/N;

    double x2 = square(freq_qiqi1-exp_freq_qiqi1)/exp_freq_qiqi1 +
                square(freq_qinotqi1-exp_freq_qinotqi1)/exp_freq_qinotqi1 +
                square(freq_notqiqi1-exp_freq_notqiqi1)/exp_freq_notqiqi1 +
                square(freq_notqinotqi1-exp_freq_notqinotqi1)/exp_freq_notqinotqi1;
    return x2; 
}

template<size_t t_min_freq = 0>
struct phrase_detector_x2 {
    phrase_detector_x2() = delete;

    static std::string name() { return "X2"; }

    template<class t_index>
    static result_t parse(t_index& index,const std::vector<uint64_t>& qry,double threshold)
    {
        result_t phrases;

        // find stopwords
        std::vector<bool> b(qry.size(),0);
        for (auto begin = qry.begin(); begin != qry.end(); ++begin){
            auto scores = index.max_sim_scores(begin, begin+1);
            if(!scores.empty())
                b[begin-qry.begin()] = scores[0] < 1;
        }

        std::vector<double> freq_single;
        for(size_t i=0;i<qry.size();i++) {
            double freq = index.csa_count(qry.begin()+i,qry.begin()+i+1);
            freq_single.push_back(freq);
        }

        // compute adjacent pair probabilities
        for(size_t i=0;i<qry.size()-1;i++) {
            double x2 = compute_x2(freq_single,index,i,qry.begin()+i,qry.begin()+i+2);
            if ( !b[i] and !b[i+1] ){
                double freq = index.csa_count(qry.begin()+i,qry.begin()+i+2);
                if(freq >= t_min_freq) {
                    phrases.emplace_back(x2,std::vector<uint64_t>(qry.begin()+i,qry.begin()+i+2));
                }
            }
        }

        return phrases;
    }
};

template<size_t t_min_freq = 10,bool t_length_norm = false>
struct phrase_detector_x2_greedy {
    phrase_detector_x2_greedy() = delete;

    static std::string name() { 
        if(t_length_norm) return "X2-GREEDY-LN";
        return "X2-GREEDY"; 
    }

    template<class t_index>
    static result_t parse(t_index& index,const std::vector<uint64_t>& qry,double threshold)
    {
        result_t phrases;

        std::vector<bool> b(qry.size(),0);
        for (auto begin = qry.begin(); begin != qry.end(); ++begin){
            auto scores = index.max_sim_scores(begin, begin+1);
            if(!scores.empty())
                b[begin-qry.begin()] = scores[0] < 1;
        }

        std::vector<double> freq_single;
        for(size_t i=0;i<qry.size();i++) {
            double freq = index.csa_count(qry.begin()+i,qry.begin()+i+1);
            freq_single.push_back(freq);
        }
        // compute adjacent pair probabilities
        std::vector<double> x2_pairs;
        for(size_t i=0;i<qry.size()-1;i++) {
            double x2 = compute_x2(freq_single,index,i,qry.begin()+i,qry.begin()+i+2);
            if(t_length_norm) x2 *= log10(2);
            if ( !b[i] and !b[i+1] ){
                double freq = index.csa_count(qry.begin()+i,qry.begin()+i+2);
                if(freq >= t_min_freq) {
                    phrases.emplace_back(x2,std::vector<uint64_t>(qry.begin()+i,qry.begin()+i+2));
                }
            }
            x2_pairs.push_back(x2);
        }

        // add tuples by combining high adjacent pairs
        for (auto begin = qry.begin(); begin != qry.end(); ++begin){
            for (auto end = begin+1; end != qry.end(); ++end){

                double freq = index.csa_count(begin,end+1);
                if(freq < t_min_freq) continue;


                if ( !b[begin-qry.begin()] and !b[end-qry.begin()] ){
                    auto start = std::distance(qry.begin(),begin);
                    auto stop = std::distance(qry.begin(),end);
                    auto len = stop-start+1;
                    double tuple_score = x2_pairs[start];
                    for(size_t i=start;i<=stop;i++) {
                        tuple_score = std::min(tuple_score,x2_pairs[i]);
                    }
                    if(t_length_norm) tuple_score *= log10(len);
                    phrases.emplace_back(tuple_score,std::vector<uint64_t>(begin,end+1));
                }
            }
        }

        return phrases;
    }
};

template<size_t t_min_freq = 0>
struct phrase_detector_bm25 {
    phrase_detector_bm25() = delete;

    static std::string name() { 
        return "BM25"; 
    }

    template<class t_index>
    static result_t parse(t_index& index,const std::vector<uint64_t>& qry,double threshold)
    {
        result_t phrases;

        // find the stop words
        std::vector<bool> b(qry.size(),0);
        for (auto begin = qry.begin(); begin != qry.end(); ++begin){
            auto scores = index.max_sim_scores(begin, begin+1);
            if(!scores.empty())
                b[begin-qry.begin()] = scores[0] < 1;
        }

        // process
        for (auto begin = qry.begin(); begin != qry.end(); ++begin){
            for (auto end = begin+1; end != qry.end(); ++end){

                double freq = index.csa_count(begin,end+1);
                if(freq < t_min_freq) continue;

                if ( !b[begin-qry.begin()] and !b[end-qry.begin()] ){
                    auto scores = index.max_sim_scores(begin, end+1);
                    if(!scores.empty()) {
                        double score = scores[0];
                        phrases.emplace_back(score,std::vector<uint64_t>(begin,end+1));
                    }
                }
            }
        }
        return phrases;
    }
};

template<size_t t_min_freq = 0,bool t_length_norm = false>
struct phrase_detector_exist_prob {
    phrase_detector_exist_prob() = delete;

    static std::string name() { 
        if(t_length_norm) return "EXIST-PROB-LN";
        return "EXIST-PROB"; 
    }

    template<class t_index>
    static result_t parse(t_index& index,const std::vector<uint64_t>& qry,double threshold)
    {
        result_t phrases;

        // find stopwords
        std::vector<bool> b(qry.size(),0);
        for (auto begin = qry.begin(); begin != qry.end(); ++begin){
            auto scores = index.max_sim_scores(begin, begin+1);
            if(!scores.empty())
                b[begin-qry.begin()] = scores[0] < 1;
        }

        // find tuples
        for (auto begin = qry.begin(); begin != qry.end(); ++begin){
            for (auto end = begin+1; end != qry.end(); ++end){

                double freq = index.csa_count(begin,end+1);
                if(freq < t_min_freq) continue;

                if ( !b[begin-qry.begin()] and !b[end-qry.begin()] ){
                    auto prob = index.phrase_prob(begin,end+1);
                    size_t len = std::distance(begin,end+1);
                    if(t_length_norm) prob *= log10(len);
                    phrases.emplace_back(prob,std::vector<uint64_t>(begin,end+1));
                }
            }
        }

        return phrases;
    }
};


}// end namespace

#endif
