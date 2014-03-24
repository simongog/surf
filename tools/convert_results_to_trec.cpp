
#include <unistd.h>
#include <stdlib.h>
#include <algorithm>
#include <vector>
#include <iostream>
#include <sstream>

#include "surf/config.hpp"
#include "surf/util.hpp"

typedef struct cmdargs {
    std::string collection_dir;
    std::string surf_file;
    std::string trec_file;
} cmdargs_t;

void
print_usage(char* program)
{
    fprintf(stdout,"%s -c <collection directory> -q <query file> -r <surf_res.csv> -o <res.trec>\n",program);
    fprintf(stdout,"where\n");
    fprintf(stdout,"  -c <collection directory>  : the directory the collection is stored.\n");
    fprintf(stdout,"  -r <surf results csv>  : the results file produced by surf.\n");
    fprintf(stdout,"  -o <trec file>  : results converted to trec format.\n");
};

cmdargs_t
parse_args(int argc,char* const argv[])
{
    cmdargs_t args;
    int op;
    args.collection_dir = "";
    args.surf_file = "";
    args.trec_file = "";
    while ((op=getopt(argc,argv,"c:r:o:")) != -1) {
        switch (op) {
            case 'c':
                args.collection_dir = optarg;
                break;
            case 'r':
                args.surf_file = optarg;
                break;
            case 'o':
                args.trec_file = optarg;
                break;
            case '?':
            default:
                print_usage(argv[0]);
        }
    }
    if (args.collection_dir==""||args.surf_file==""||args.trec_file=="") {
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
    return tokens;
}

int main( int argc, char** argv ) {

    /* parse command line */
    cmdargs_t args = parse_args(argc,argv);

    /* parse repo */
    surf::parse_collection(args.collection_dir);

    /* load the docnames map */
    std::unordered_map<uint64_t,std::string> id_mapping;
    auto docnames_file = args.collection_dir + surf::DOCNAMES_FILENAME;
    std::ifstream dfs(docnames_file);
    size_t j=0;
    std::string name_mapping;
    while( std::getline(dfs,name_mapping) ) {
        id_mapping[j] = name_mapping;
        j++;
    }

    std::ofstream trec_out(args.trec_file);
    std::ifstream surfres_fs(args.surf_file);
    for(std::string line; std::getline(surfres_fs,line);) {
        auto tokens = tokenize(line);
        auto qry_id = std::strtoul(tokens[0].c_str(),NULL,10);
        auto rank = std::strtoul(tokens[0].c_str(),NULL,10);
        auto doc_id = std::strtoul(tokens[0].c_str(),NULL,10);
        auto doc_score = std::strtod(tokens[0].c_str(),NULL);

        trec_out 
            << qry_id << "\t"
            << "Q0" << "\t"
            << id_mapping[doc_id] << "\t"
            << rank               << "\t"
            << doc_score << "\t"
            << "SURF" << std::endl;
    }
}


