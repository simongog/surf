// extracts from an indri index a monoton sequence of integers in sdsl format
// which represent the parsed text collection.
#include <iostream>
#include <unistd.h>
#include <stdlib.h>

#include "surf/query.hpp"
#include "sdsl/config.hpp"
#include "surf/query_parser.hpp"
#include "surf/util.hpp"

typedef struct cmdargs {
    std::string collection_dir;
    std::string query_file;
    std::string output_file;
    uint64_t num_qrys;
} cmdargs_t;

void
print_usage(char* program)
{
    fprintf(stdout,"%s -c <collection directory> -q <query file> -n <top-k> -o <output.qry>\n",program);
    fprintf(stdout,"where\n");
    fprintf(stdout,"  -c <collection directory>  : the directory the collection is stored.\n");
    fprintf(stdout,"  -q <query file>  : the queries to be processed.\n");
    fprintf(stdout,"  -n <num-qrys>  : the number of queries to be selected.\n");
    fprintf(stdout,"  -o <output.qry>  : selected queries.\n");
};

cmdargs_t
parse_args(int argc,char* const argv[])
{
    cmdargs_t args;
    int op;
    args.collection_dir = "";
    args.query_file = "";
    args.output_file = "";
    args.num_qrys = 1000;
    while ((op=getopt(argc,argv,"c:q:n:o:")) != -1) {
        switch (op) {
            case 'c':
                args.collection_dir = optarg;
                break;
            case 'q':
                args.query_file = optarg;
                break;
            case 'o':
                args.output_file = optarg;
                break;
            case 'n':
                args.num_qrys = std::strtoul(optarg,NULL,10);
                break;
            case '?':
            default:
                print_usage(argv[0]);
        }
    }
    if (args.collection_dir==""||args.query_file=="") {
        std::cerr << "Missing command line parameters.\n";
        print_usage(argv[0]);
        exit(EXIT_FAILURE);
    }
    return args;
}
int main( int argc, char** argv ) {
    if(argc < 3) {
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }

    /* parse command line */
    cmdargs_t args = parse_args(argc,argv);

    /* parse repo */
    auto cc = surf::parse_collection(args.collection_dir);

    /* parse queries */
    std::cout << "Parsing query file '" << args.query_file << "'" << std::endl;
    auto queries = surf::query_parser::parse_queries(args.collection_dir,args.query_file);
    std::cout << "Found " << queries.size() << " queries." << std::endl;

    /* select num_queries random ones */
    std::mt19937 gen(4711);
    std::shuffle(queries.begin(), queries.end(), gen);
    auto id_sort = [](const surf::query_t& a,const surf::query_t& b) {
        return std::get<0>(a) < std::get<0>(b);
    };
    std::sort(queries.begin(),queries.begin()+args.num_qrys,id_sort);
    
    /* output */
    std::ofstream selected_fs(args.output_file);
    if(selected_fs.is_open()) {
        for(size_t i=0;i<args.num_qrys;i++) {
            selected_fs << std::get<0>(queries[i]) << ";";
            const auto& raw_tokens = std::get<2>(queries[i]);
            for(size_t j=0;j<raw_tokens.size()-1;j++) {
                selected_fs << raw_tokens[j] << " ";
            }
            selected_fs << raw_tokens.back() << std::endl;
        }
    } else {
        perror("could not open output file.");
    }
}


