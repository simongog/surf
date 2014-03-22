#ifndef SURF_DOC_PERM_HPP
#define SURF_DOC_PERM_HPP

#include <sdsl/int_vector.hpp>
#include <string>

namespace surf{

struct doc_perm{
    typedef typename sdsl::int_vector<>::size_type size_type;
    sdsl::int_vector<> id2len; // doc id to length ordered id
    sdsl::int_vector<> len2id; // length ordered id to doc id

    inline size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = NULL, std::string name = "") const {
        using namespace sdsl;
        structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
        size_type written_bytes = 0;
        written_bytes += id2len.serialize(out, child, "id2len");
        written_bytes += len2id.serialize(out, child, "len2id");
        structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    inline void load(std::istream &in){
        id2len.load(in);
        len2id.load(in);
    }
};

}

#endif
