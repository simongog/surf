

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
    std::string surf_file;
} cmdargs_t;

void
print_usage(char* program)
{
    fprintf(stdout,"%s -c <collection directory> -q <query file> -r <surf_res.csv>",program);
    fprintf(stdout,"where\n");
    fprintf(stdout,"  -c <collection directory>  : the directory the collection is stored.\n");
    fprintf(stdout,"  -r <surf results csv>  : the results file produced by surf.\n");
};

cmdargs_t
parse_args(int argc,char* const argv[])
{
    cmdargs_t args;
    int op;
    args.collection_dir = "";
    args.surf_file = "";
    while ((op=getopt(argc,argv,"c:r:")) != -1) {
        switch (op) {
            case 'c':
                args.collection_dir = optarg;
                break;
            case 'r':
                args.surf_file = optarg;
                break;
            case '?':
            default:
                print_usage(argv[0]);
        }
    }
    if (args.collection_dir==""||args.surf_file=="") {
        std::cerr << "Missing command line parameters.\n";
        print_usage(argv[0]);
        exit(EXIT_FAILURE);
    }
    return args;
}

std::vector<std::string>
tokenize(std::string line) {
    std::vector<std::string> tokens;
    size_t pos = 0;
    std::string token;
    while ((pos = line.find(";")) != std::string::npos) {
        token = line.substr(0, pos);
        tokens.push_back(token);
        line.erase(0, pos + 1);
    }
    tokens.push_back(line);
    return tokens;
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
    std::ifstream surfres_fs(args.surf_file);
    bool first = true;
    for(std::string line; std::getline(surfres_fs,line);) {
        if(first) {
            first = false;
            continue;
        }
        auto tokens = tokenize(line);
        auto qry_id = std::strtoul(tokens[0].c_str(),NULL,10);
        auto rank = std::strtoul(tokens[1].c_str(),NULL,10);
        auto doc_id = std::strtoul(tokens[2].c_str(),NULL,10);
        auto doc_score = std::strtod(tokens[3].c_str(),NULL);

        std::cout << "=====================================================================================\n";
        std::cout << "[Q]=" << qry_id << " rank=" << rank << " docid=" 
        		  << doc_id << " score=" << doc_score <<  std::endl;

        size_t doc_start = 0;
        if(doc_id != 0) {
        	doc_start = doc_border_select(doc_id) + 1;
        }
        auto doc_stop = doc_border_select(doc_id+1) - 1;

        std::cout << "document length = " << doc_stop - doc_start + 1 << std::endl;
        std::cout << "document content  = '";
        for(size_t i=doc_start;i<=doc_stop;i++) {
        	std::cout << "<" << id_mapping[T[i]] << ">";
        }
        std::cout << std::endl;
    }
}


