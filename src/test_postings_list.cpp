
#include <vector>
#include <iostream>

#include "surf/invfile_postings_list.hpp"

int main( int argc, char** argv ) {
    using plist_type = surf::postings_list<surf::compression_codec::optpfor,128>;

    // test small uncompressed lists 
    for(size_t i=0;i<500;i++) {
        size_t n = 1 + rand()%20;
        std::vector< std::pair<uint64_t,uint64_t> > A;
        uint64_t cur_id = rand()%5000;
        for(size_t j=0;j<n;j++) {
            cur_id += rand()%5000;
            uint64_t cur_freq = 1 + rand() % 50;
            A.emplace_back(cur_id,cur_freq);
        }
        plist_type pl(A);

        auto itr = pl.begin();
        auto end = pl.end();
        size_t j=0;
        while( itr != end) {
            auto id = itr.docid();
            auto freq = itr.freq();
            if(id != A[j].first && freq != A[j].second) {
                std::cerr << "ERROR: uncompressed list";
            }
            j++;
            ++itr;
        }
    }

    // test larger compressed lists 
    for(size_t i=0;i<500;i++) {
        size_t n = 1 + rand()%20000;
        std::vector< std::pair<uint64_t,uint64_t> > A;
        uint64_t cur_id = rand()%500;
        for(size_t j=0;j<n;j++) {
            cur_id += rand()%500;
            uint64_t cur_freq = 1 + rand() % 50;
            A.emplace_back(cur_id,cur_freq);
        }
        plist_type pl(A);

        auto itr = pl.begin();
        auto end = pl.end();
        size_t j=0;
        while( itr != end) {
            auto id = itr.docid();
            auto freq = itr.freq();
            if(id != A[j].first && freq != A[j].second) {
                std::cerr << "ERROR: uncompressed list";
            }
            j++;
            ++itr;
        }
    }

}


