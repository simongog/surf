

#include <unistd.h>
#include <stdlib.h>
#include <algorithm>
#include <vector>
#include <iostream>
#include <sstream>

#include "surf/config.hpp"
#include "surf/util.hpp"
#include "sdsl/int_vector.hpp"
#include "sdsl/select_support_mcl.hpp"

typedef struct cmdargs {
    std::string collection_dir;
    uint64_t doc_id;
} cmdargs_t;

void
print_usage(char* program)
{
    fprintf(stdout,"%s -c <collection directory> -d <docid>",program);
    fprintf(stdout,"where\n");
    fprintf(stdout,"  -c <collection directory>  : the directory the collection is stored.\n");
    fprintf(stdout,"  -d <docid>  : the document to output\n");
};

cmdargs_t
parse_args(int argc,char* const argv[])
{
    cmdargs_t args;
    int op;
    args.collection_dir = "";
    int64_t doc_id = -1;
    while ((op=getopt(argc,argv,"c:d:")) != -1) {
        switch (op) {
            case 'c':
                args.collection_dir = optarg;
                break;
            case 'd':
                doc_id = std::strtoll(optarg,NULL,10);
                break;
            case '?':
            default:
                print_usage(argv[0]);
        }
    }
    if (args.collection_dir==""||doc_id<0) {
        std::cerr << "Missing command line parameters.\n";
        print_usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    args.doc_id = (uint64_t) doc_id;

    return args;
}

int main( int argc, char** argv ) {

    /* parse command line */
    cmdargs_t args = parse_args(argc,argv);

    /* parse repo */
    auto cc = surf::parse_collection(args.collection_dir);

    /* load doc border bv and build select structure */
    sdsl::bit_vector doc_border;
    sdsl::load_from_cache(doc_border, surf::KEY_DOCBORDER, cc);
    sdsl::bit_vector::select_1_type doc_border_select(&doc_border);

    /* load dictionary and create mapping */
    std::unordered_map<uint64_t,std::string> id_mapping;
    {
        auto dict_file = args.collection_dir + "/" + surf::DICT_FILENAME;
        std::ifstream dfs(dict_file);
        if(!dfs.is_open()) {
            std::cerr << "cannot load dictionary file.";
            exit(EXIT_FAILURE);
        }
        std::string term_mapping;
        while( std::getline(dfs,term_mapping) ) {
            auto sep_pos = term_mapping.find(' ');
            auto term = term_mapping.substr(0,sep_pos);
            auto idstr = term_mapping.substr(sep_pos+1);
            uint64_t id = std::stoull(idstr);
            id_mapping[id] = term;
        }
    }

    auto text_file = args.collection_dir + "/" + surf::TEXT_FILENAME;
    sdsl::int_vector_buffer<> T(text_file);
    uint64_t doc_id = args.doc_id;
    size_t doc_start = 0;
    if(doc_id != 0) {
      	doc_start = doc_border_select(doc_id) + 1;
    }
    auto doc_stop = doc_border_select(doc_id+1) - 1;

    std::cout << "document length = " << doc_stop - doc_start + 1 << std::endl;
    std::cout << "document content  = '";
    for(size_t i=doc_start;i<=doc_stop;i++) {
      	std::cout << id_mapping[T[i]] << " ";
    }
    std::cout << "'" << std::endl;
}


