
#ifndef SURF_IDX_DUMMY_HPP
#define SURF_IDX_DUMMY_HPP

#include "sdsl/int_vector.hpp"
#include "sdsl/construct.hpp"

namespace surf
{

class idx_dummy
{
    public:
        using size_type = sdsl::int_vector<>::size_type;
    public:

        idx_dummy() = default;

        idx_dummy(sdsl::cache_config&) { }

        auto serialize(std::ostream& out, sdsl::structure_tree_node* v=NULL, std::string name="") const -> size_type {
            size_type written_bytes = 0;
            return written_bytes;
        }

        void load(std::istream &){ }


};

inline void construct(idx_dummy &, const std::string&,
               sdsl::cache_config&, uint8_t){
}



}

#endif
