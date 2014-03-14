#include <unistd.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "sdsl/config.hpp"
#include "surf/indexes.hpp"
#include "surf/query_parser.hpp"

typedef struct cmdargs {
    std::string collection_dir;
    std::string query_file;
    uint64_t k;
} cmdargs_t;

void
print_usage(char* program)
{
    fprintf(stdout,"%s -c <collection directory> -q <query file> -k <top-k>\n",program);
    fprintf(stdout,"where\n");
    fprintf(stdout,"  -c <collection directory>  : the directory the collection is stored.\n");
    fprintf(stdout,"  -q <query file>  : the queries to be performed.\n");
    fprintf(stdout,"  -k <top-k>  : the top-k documents to be retrieved for each query.\n");
};

cmdargs_t
parse_args(int argc,char* const argv[])
{
    cmdargs_t args;
    int op;
    args.collection_dir = "";
    args.query_file = "";
    args.k = 10;
    while ((op=getopt(argc,argv,"c:q:k:")) != -1) {
        switch (op) {
            case 'c':
                args.collection_dir = optarg;
                break;
            case 'q':
                args.query_file = optarg;
                break;
            case 'k':
                args.k = std::strtoul(optarg,NULL,10);
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

int main(int argc,char* const argv[])
{
    using clock = std::chrono::high_resolution_clock;
    /* parse command line */
    cmdargs_t args = parse_args(argc,argv);

    /* parse repo */
    auto cc = surf::parse_collection(args.collection_dir);

    /* parse queries */
    std::cout << "Parsing query file '" << args.query_file << "'" << std::endl;
    auto queries = surf::query_parser::parse_queries(args.collection_dir,args.query_file);
    std::cout << "Found " << queries.size() << " queries." << std::endl;

    /* define types */
    using surf_index_t = INDEX_TYPE;
    std::string index_name = IDXNAME;

    /* load the index */
    surf_index_t index;
    auto load_start = clock::now();
    construct(index, "", cc, 0);
    index.load(cc);
    auto load_stop = clock::now();
    auto load_time_sec = std::chrono::duration_cast<std::chrono::seconds>(load_stop-load_start);
    std::cout << "Index loaded in " << load_time_sec.count() << " seconds." << std::endl;

    /* process the queries */
    std::vector< std::tuple<uint64_t,result_t,std::chrono::microseconds> > timings;
    size_t num_runs = 3;
    for(size_t i=0;i<num_runs;i++) {
        for(const auto& query: queries) {
            auto id = query.first;
            auto qry_tokens = query.second;
            std::cout << "[" << id << "] |Q|=" << qry_tokens.size(); std::cout.flush();

            // run the query
            auto qry_start = clock::now();
            auto results = index.search(qry_tokens,args.k);
            auto qry_stop = clock::now();

            auto query_time = std::chrono::duration_cast<std::chrono::microseconds>(qry_stop-qry_start);
            std::cout << " TIME = " << std::setprecision(5)
                      << query_time.count() / 1000.0 
                      << " ms" << std::endl;
            timings.emplace_back(id,results,query_time);
        }
    }

    return EXIT_SUCCESS;
}
