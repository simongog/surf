
#ifndef SURF_IDX_DUMMY_HPP
#define SURF_IDX_DUMMY_HPP

#include "sdsl/int_vector.hpp"
#include "sdsl/io.hpp"
#include "surf/query.hpp"

namespace surf
{

class idx_dummy
{
    public:
        using size_type = sdsl::int_vector<>::size_type;
    public:

        idx_dummy() = default;

        idx_dummy(sdsl::cache_config&) { }

        result_t search(std::vector<query_token> qry,size_t k) {
        	result_t res;

        	return res;
        }

        void load(const sdsl::cache_config&){}
};

inline void construct(idx_dummy &, const std::string&,
               sdsl::cache_config&, uint8_t){
}



}

#endif
