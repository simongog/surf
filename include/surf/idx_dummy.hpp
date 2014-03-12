
#ifndef SURF_IDX_DUMMY_HPP
#define SURF_IDX_DUMMY_HPP

#include "sdsl/int_vector.hpp"

namespace surf
{

class idx_dummy
{
    public:
        using size_type = sdsl::int_vector<>::size_type;
    public:

        idx_dummy() = default;

        idx_dummy(sdsl::cache_config&) { }

        result_t search(std::vector<uint64_t> qry,size_t k) {
        	result_t res;

        	return res;
        }
};

inline void construct(idx_dummy &, const std::string&,
               sdsl::cache_config&, uint8_t){
}



}

#endif
