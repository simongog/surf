
#ifndef SURF_IDX_DUMMY_HPP
#define SURF_IDX_DUMMY_HPP

#include "sdsl/int_vector.hpp"

using namespace sdsl;

namespace surf
{

class idx_dummy
{
    public:
        using size_type = sdsl::int_vector<>::size_type;
    public:
        idx_dummy(std::istream&) {
        }
        idx_dummy(cache_config&) {

        }
        auto serialize(std::ostream& out, structure_tree_node* v=NULL, std::string name="") const -> size_type {
            size_type written_bytes = 0;
            return written_bytes;
        }
};

}

#endif
