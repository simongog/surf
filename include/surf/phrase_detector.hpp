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
using heap_t = std::priority_queue<std::pair<double,std::vector<uint64_t>>>;

namespace surf{

double prob(double x){
    return x==0 ? -99999999 : log(x);
}

struct phrase_detector {
    phrase_detector() = delete;

    template<class t_index>
    static void parse_greedy_lr(t_index& index,
    								  const std::vector<uint64_t>& qry,
    								  double threshold,
                                      heap_t& heap)
    {
    	//compute single term probabilities
    	std::vector<double> P_single;
    	for(size_t i=0;i<qry.size();i++) {
    		auto cnt = index.csa_count(qry.begin()+i,qry.begin()+i+1);
    		double single = (double)cnt / (double)index.csa_size();
    		P_single.push_back(single);
    	}

        // compute all phrases
    	parsed_qry phrases;
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
                double joint = prob((double)cnt/index.csa_size());
                assoc_ratio = joint-single;
            } while ( assoc_ratio >= threshold );
            heap.emplace(assoc_ratio,std::vector<uint64_t>(start, end));
            start = end;
        }
    }

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

    template<class t_index>
    static void parse_x2(t_index& index,
                                        const std::vector<uint64_t>& qry,
                                        double threshold,
                                        heap_t& heap)
    {
        
        std::vector<double> freq_single;
        for(size_t i=0;i<qry.size();i++) {
            double freq = index.csa_count(qry.begin()+i,qry.begin()+i+1);
            freq_single.push_back(freq);
        }
        // compute adjacent pair probabilities
        for(size_t i=0;i<qry.size()-1;i++) {
           double x2 = compute_x2(freq_single,index,i,qry.begin()+i,qry.begin()+i+2);
           heap.emplace(x2,std::vector<uint64_t>(qry.begin()+i,qry.begin()+i+2));
        }
    }


    template<class t_index>
    static void parse_greedy_paul(t_index& index,
    								    const std::vector<uint64_t>& qry,
    								    double threshold,
                                        heap_t& heap)
    {
    	//compute single term probabilities
    	std::vector<double> P_single;
    	for(size_t i=0;i<qry.size();i++) {
    		auto cnt = index.csa_count(qry.begin()+i,qry.begin()+i+1);
    		double single = (double)cnt / (double)index.csa_size();
    		P_single.push_back(prob(single));
    	}

    	for(size_t i=0;i<qry.size()-1;i++) {
    		auto cnt = index.csa_count(qry.begin()+i,qry.begin()+i+2);
    		double joint = (double)cnt / (double)index.csa_size();
			// single
			double single = P_single[i]+P_single[i+1];
			// calc ratio
			double assoc_ratio = prob(joint)-single;
            heap.emplace(assoc_ratio,
                std::vector<uint64_t>(qry.begin()+i,qry.begin()+i+2));
    	}
    }

    typedef std::vector<bool> tVB;
    typedef std::vector<tVB> tVVB;

    typedef std::vector<size_t> tVI;
    typedef std::vector<tVI> tVVI;
    typedef std::pair<size_t,size_t> tPII;
    typedef std::vector<tPII> tVPII;

    static size_t max_phrase(size_t l, size_t r, const tVVB& is_phrase, tVVI& mem){
        if ( l > r )
            return 0;
        if ( mem[l][r] !=  (size_t)-1 ) {
            return mem[l][r];
        }
        mem[l][r] = 0;
        if (is_phrase[l][r]){
            mem[l][r] = 1;
        }
        for (size_t i=l+1; i<=r; ++i){
            size_t res = max_phrase(l, i, is_phrase, mem) + max_phrase(i+1, r, is_phrase, mem);
            if (res > mem[l][r]){
                mem[l][r] = res;
            }
        }
        return mem[l][r];
    }

    static size_t get_res(size_t l, size_t r, tVPII &res, const tVVB& is_phrase, tVVI& mem){
        size_t maxi = max_phrase(l, r, is_phrase, mem);
//        std::cout<<"get_res("<<l<<","<<r<<","<<maxi<<std::endl;
        if ( maxi > 0 ){
             for (size_t i=l+1; i<r; ++i){
                size_t m = max_phrase(l, i,is_phrase, mem) + max_phrase(i+1, r, is_phrase,mem);
                if (maxi == m) {
                    get_res(l,i,res, is_phrase,mem);
                    get_res(i+1, r, res, is_phrase, mem);
                    return maxi;
                }
            }       
            res.emplace_back(l,r);
        } 
        return maxi;
    }

    template<class t_index>
    static void parse_dp(t_index& index,
    								  const std::vector<uint64_t>& qry,
    								  double threshold,
                                      heap_t& heap)
    {
    	//compute single term probabilities
    	std::vector<double> P_single;
    	for(size_t i=0;i<qry.size();i++) {
    		auto cnt = index.csa_count(qry.begin()+i,qry.begin()+i+1);
    		double single = (double)cnt / (double)index.csa_size();
    		P_single.push_back(single);
    	}

    	//compute single term probabilities
        tVVB is_phrase(qry.size(), tVB(qry.size(), false));
        for (auto start = qry.begin(); start < qry.end(); ){
            double single = prob(P_single[start-qry.begin()]);
            auto end = start;
            double assoc_ratio = 0;
            while ( ++end != qry.end() ) {
                single += prob(P_single[end-qry.begin()]);
    			auto cnt = index.csa_count(start, end+1);
                double joint = prob((double)cnt/index.csa_size());
                assoc_ratio = joint-single;
                is_phrase[start-qry.begin()][end-qry.begin()] = assoc_ratio >= threshold;
            };
            start = end;
        }
        for (size_t i=0; i<qry.size(); ++i){
            for(size_t j=0; j<qry.size(); ++j)
                std::cout << is_phrase[i][j] << " ";
            std::cout<<std::endl;
        }
    	parsed_qry phrases;
        tVVI mem(qry.size(),tVI(qry.size(),(size_t)-1));

        max_phrase(0, qry.size()-1, is_phrase, mem);
//        std::cout << "max_p="<<max_p<<std::endl;
        tVPII res;
        get_res(0, qry.size()-1, res, is_phrase, mem);
//        std::cout<<"res.size()="<<res.size()<<std::endl;
        std::sort(res.begin(), res.end());
        for(size_t i=0; i<res.size(); ++i){
//            std::cout<<"...["<<i<<"]="<<res[i].first<<","<<res[i].second<<std::endl;
//            std::cout<<"___"<< is_phrase[res[i].first][res[i].second]<<std::endl;
        }

        for(size_t i=0,j=0; i<qry.size(); ){
            if (j >= res.size() or res[j].first > i ){
//                std::cout<<"qry["<<i<<"]="<<qry[i]<<std::endl;
                phrases.emplace_back(std::vector<uint64_t>(1, qry[i]));
                ++i;
            } else {
//                std::cout<<"res["<<j<<"]="<<res[j].first<<","<<res[j].second<<std::endl;
                phrases.emplace_back(std::vector<uint64_t>(qry.begin()+res[j].first, qry.begin()+res[j].second+1));
                i = res[j].second+1;
                ++j;
            }
        }
    }

    template<class t_index>
    static void parse_bm25(t_index& index,
                                const std::vector<uint64_t>& qry,
                                double threshold,
                                heap_t& heap)
    {
        std::vector<bool> b(qry.size(),0);
        for (auto begin = qry.begin(); begin != qry.end(); ++begin){
            auto scores = index.max_sim_scores(begin, begin+1);
            if(!scores.empty())
                b[begin-qry.begin()] = scores[0] < 1;
        }
        for (auto begin = qry.begin(); begin != qry.end(); ++begin){
            for (auto end = begin+1; end != qry.end(); ++end){
                if ( !b[begin-qry.begin()] and !b[end-qry.begin()] ){
                    auto scores = index.max_sim_scores(begin, end+1);
                    if(!scores.empty())
                        heap.emplace(scores[0],std::vector<uint64_t>(begin,end+1));
                }
            }
        }
    }

    template<class t_index>
    static void parse_exist_prob(t_index& index,
                                const std::vector<uint64_t>& qry,
                                double threshold,
                                heap_t& heap)
    {
        for (auto begin = qry.begin(); begin != qry.end(); ++begin) {
            for (auto end = begin+1; end != qry.end(); ++end){
                auto prob = index.phrase_prob(begin,end+1);
                heap.emplace(prob,std::vector<uint64_t>(begin,end+1));
            }
        }
    }
};
}// end namespace

#endif
