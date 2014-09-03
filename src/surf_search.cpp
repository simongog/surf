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
    fprintf(stdout,"%s -c <collection directory> -q <query file> -k <top-k> -o <output.csv>\n",program);
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

    /* define types */
    using surf_index_t = INDEX_TYPE;
    std::string index_name = IDXNAME;

    /* parse repo */
    auto cc = surf::parse_collection<surf_index_t::alphabet_category>(args.collection_dir);

    /* parse queries */
    std::cout << "Parsing query file '" << args.query_file << "'" << std::endl;
    std::vector<surf::query_t> queries;
    queries = surf::query_parser::parse_queries<surf_index_t::alphabet_category>(args.collection_dir,args.query_file);
    std::cout << "Found " << queries.size() << " queries." << std::endl;

    /* load the index */
    surf_index_t index;
    auto load_start = clock::now();
//    construct(index, "", cc, 0);
    index.load(cc);
    auto load_stop = clock::now();
    auto load_time_sec = std::chrono::duration_cast<std::chrono::seconds>(load_stop-load_start);
    std::cout << "Index loaded in " << load_time_sec.count() << " seconds." << std::endl;

    /* process the queries */
    std::map<uint64_t,std::chrono::microseconds> query_times;
    std::map<uint64_t,surf::result> query_results;
    std::map<uint64_t,uint64_t> query_lengths;

    size_t num_runs = 1;
    for(size_t i=0;i<num_runs;i++) {
        for(const auto& query: queries) {
            auto id = std::get<0>(query);
            auto qry_tokens = std::get<1>(query);
            std::cout << "[" << id << "] |Q|=";
            if (qry_tokens.size()==0) {
                std::cout << qry_tokens.size();
            } else {
                std::cout << qry_tokens[0].token_ids.size();
            }
            std::cout.flush();
            // run the query
            auto qry_start = clock::now();
            auto results = index.search(qry_tokens,args.k);
            auto qry_stop = clock::now();

            auto query_time = std::chrono::duration_cast<std::chrono::microseconds>(qry_stop-qry_start);
            std::cout << " TIME = " << std::setprecision(5)
                      << query_time.count() / 1000.0 
                      << " ms" << std::endl;

            auto itr = query_times.find(id);
            if(itr != query_times.end()) {
                itr->second += query_time;
            } else {
                query_times[id] = query_time;
            }

            if(i==0) {
                query_results[id] = results;
                query_lengths[id] = qry_tokens.size();
            }
        }
    }

    /* output results to csv */
    char time_buffer [80] = {0};
    std::time_t t = std::time(NULL);
    auto timeinfo = localtime (&t);
    strftime (time_buffer,80,"%F-%H:%M:%S",timeinfo);
    std::string time_output_file = args.collection_dir + "/results/" 
                   + "surf-timings-" + index_name + "-k" + std::to_string(args.k) 
                   + "-" + std::string(time_buffer) + ".csv";
    std::string res_output_file = args.collection_dir + "/results/" 
                   + "surf-results-" + index_name + "-k" + std::to_string(args.k) 
                   + "-" + std::string(time_buffer) + ".csv";

    /* calc average */
    for(auto& timing : query_times) {
        timing.second = timing.second / num_runs;
    }

    /* output */
    {
        std::cout << "Writing timing results to '" << time_output_file << "'" << std::endl;     
        std::ofstream resfs(time_output_file);
        if(resfs.is_open()) {
            resfs << "id;index;k;num_terms;time_ms" << std::endl;
            for(const auto& timing: query_times) {
                auto qry_id = timing.first;
                auto qry_time = timing.second;
                resfs << qry_id << ";" << index_name << ";" << args.k << ";"
                          << query_lengths[qry_id] << ";"
                          << qry_time.count() / 1000.0  << "\n"; 
            }
        } else {
            perror("could not output results to file.");
        }
        std::cout << "Writing result listing to '" << res_output_file << "'" << std::endl;
        std::ofstream res_outfs(res_output_file);
        if(res_outfs.is_open()) {
            res_outfs << "id;rank;docid;score" << std::endl;
            for(const auto& result: query_results) {
                auto qry_id = result.first;
                auto qry_res = result.second.list;
                for(size_t i=1;i<=qry_res.size();i++) {
                    res_outfs << qry_id << ";" 
                              << i  << ";" 
                              << qry_res[i-1].doc_id << ";" 
                              << qry_res[i-1].score << "\n"; 
                }
            }
        } else {
            perror("could not output results to file.");
        }
    }


    return EXIT_SUCCESS;
}
