#include <unistd.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <ctime>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "surf/query.hpp"
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
    fprintf(stdout,"%s -c <collection directory> -k <top-k>\n",program);
    fprintf(stdout,"where\n");
    fprintf(stdout,"  -c <collection directory>  : the directory the collection is stored.\n");
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
    while ((op=getopt(argc,argv,"c:k:")) != -1) {
        switch (op) {
            case 'c':
                args.collection_dir = optarg;
                break;
            case 'k':
                args.k = std::strtoul(optarg,NULL,10);
                break;
            case '?':
            default:
                print_usage(argv[0]);
        }
    }
//    if (args.collection_dir==""||args.query_file=="") {
    if (args.collection_dir=="") {
        std::cerr << "Missing command line parameters.\n";
        print_usage(argv[0]);
        exit(EXIT_FAILURE);
    }
    return args;
}

int main(int argc,char* const argv[])
{
    /* define types */
    using surf_index_t = INDEX_TYPE;
    using csa_t = CSA_TYPE;
    std::string index_name = IDXNAME;


    using clock = std::chrono::high_resolution_clock;
    /* parse command line */
    cmdargs_t args = parse_args(argc,argv);

    /* parse repo */
    auto cc = surf::parse_collection<csa_t::alphabet_category>(args.collection_dir);

    /* load the index */
    surf_index_t index;
    auto load_start = clock::now();
//    construct(index, "", cc, 0);
    index.load(cc);
    auto load_stop = clock::now();
    auto load_time_sec = std::chrono::duration_cast<std::chrono::seconds>(load_stop-load_start);
    std::cout << "Index loaded in " << load_time_sec.count() << " seconds." << std::endl;

    /* process the queries */
//    std::map<uint64_t,std::chrono::microseconds> query_times;
//    std::map<uint64_t,surf::result> query_results;
//    std::map<uint64_t,uint64_t> query_lengths;
    string qry;
    while ( std::cin >> qry ){
        auto res_iter = index.topk(qry.begin(), qry.end());
        uint64_t k=0;
        while ( res_iter and k < args.k ){
            auto docid_weight = *res_iter;
            uint64_t doc_id = docid_weight.first;
            std::cout<<"DOC_ID="<<doc_id<<" WEIGHT="<<docid_weight.second<<std::endl;
            string doc = index.doc(doc_id);
            auto found = doc.find(qry);
            uint64_t check_weight = 0;
            while (found!=std::string::npos) {
                std::cout << "found '"<<qry<<"' at position" << found << '\n';
                found = doc.find(qry, found+1);
                ++check_weight;
            }
            if ( check_weight != docid_weight.second ){
                std::cerr<<"ERROR: for query"<<qry<<" : ";
                std::cerr<<"check_weight="<<check_weight;
                std::cerr<<" != "<<docid_weight.second<<"=weight"<<std::endl;
                return 1;
            }
            ++res_iter;
            ++k;
        }
    }

    return EXIT_SUCCESS;
}
