#include "surf/config.hpp"
#include "surf/indexes.hpp"
#include <iostream>
#include <chrono>
#include <algorithm>
#include <sstream>

using namespace std;
using namespace sdsl;
using namespace surf;

using idx_type = INDEX_TYPE;

const size_t buf_size=1024*128;
char   buffer[buf_size];

template<typename X>
struct myline {
    static string parse(char* str) {
        return string(str);
    }
};

template<>
struct myline<sdsl::int_alphabet_tag> {
    static vector<uint64_t> parse(char* str) {
        vector<uint64_t> res;
        stringstream ss(str);
        uint64_t x;
        while (ss >> x) {
            res.push_back(x);
        }
        return res;
    }
};

int main(int argc, char* argv[])
{
    if (argc < 4) {
        cout << "Usage: " << argv[0] << " collection_dir pattern_file k <debug>" << endl;
        cout << " Process all queries with the index." << endl;
        return 1;
    }
    string collection_dir = string(argv[1]);
    string pattern_file = string(argv[2]);
    size_t k = stoull(argv[3]);
    bool debug = argc > 4;
    idx_type idx;

    using timer = chrono::high_resolution_clock;

    if (!debug) cout<<"# index_file = "<< collection_dir <<endl;
    auto cc = surf::parse_collection<idx_type::alphabet_category>(collection_dir);
    idx.load(cc);
    if (!debug) {
        cout<<"# pattern_file = "<<pattern_file<<endl;
        cout<<"# doc_cnt = "<<idx.doc_cnt()<<endl;
        cout<<"# word_cnt = "<<idx.word_cnt()<<endl;
        cout<<"# k = "<<k <<endl;
    }
    ifstream in(pattern_file);
    if (!in) {
        cerr << "Could not load pattern file" << endl;
        return 1;
    }

    using timer = chrono::high_resolution_clock;
    size_t q_len = 0;
    size_t q_cnt = 0;
    size_t sum = 0;
    size_t sum_fdt = 0;
    bool tle = false; // flag: time limit exceeded
    auto start = timer::now();
    while (!tle and in.getline(buffer, buf_size)) {
        auto q_start = timer::now();
        auto query = myline<idx_type::alphabet_category>::parse(buffer);
//        cout<<query<<endl;
        q_len += query.size();
        ++q_cnt;
        size_t x = 0;
        auto res_it = idx.topk(query.begin(), query.end());
        while ( x < k and res_it ){
            ++x;
            sum_fdt += (*res_it).second;
            if ( debug ) {
                cout<<q_cnt<<";"<<x<<";"<<(*res_it).first<< ";"<<(*res_it).second << endl;
            }
            ++res_it;
        }
//        cout<<"x="<<x<<endl;
        sum += x;
        auto q_time = timer::now()-q_start;
        // single query should not take more then 5 seconds
        if (chrono::duration_cast<chrono::seconds>(q_time).count() > 5) {
            tle = true;
        }
    }
    auto stop = timer::now();
    auto elapsed = stop-start;
    if ( !debug ){
        cout<<"# TLE = " << tle << endl;
        cout<<"# query_len = "<<q_len/q_cnt<<endl;
        cout<<"# queries = " <<q_cnt <<endl;
        cout<<"# time_per_query = "<<chrono::duration_cast<chrono::microseconds>(elapsed).count()/q_cnt <<endl;
        cout<<"# check_sum = "<<sum<<endl;
        cout<<"# check_sum_fdt = "<<sum_fdt<<endl;
    }
}
