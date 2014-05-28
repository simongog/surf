#include <iostream>
#include <cstdlib>
#include <unistd.h>

#include "sdsl/suffix_arrays.hpp"
#include "surf/query_parser.hpp"
#include "surf/phrase_detector.hpp"
#include "surf/util.hpp"

typedef struct cmdargs {
    std::string collection_dir;
    std::string query_file;
    double threshold;
} cmdargs_t;

void
print_usage(char* program)
{
    fprintf(stdout,"%s -c <collection directory> -q <query file> -t <threshold>\n",program);
    fprintf(stdout,"where\n");
    fprintf(stdout,"  -c <collection directory>  : the directory the collection is stored.\n");
    fprintf(stdout,"  -q <query file>  : the queries to be performed.\n");
    fprintf(stdout,"  -t <threshold>  : threshold for phrase detector.\n");
};

cmdargs_t
parse_args(int argc,char* const argv[])
{
    cmdargs_t args;
    int op;
    args.collection_dir = "";
    args.query_file = "";
    args.threshold = 0.0;
    while ((op=getopt(argc,argv,"c:q:t:")) != -1) {
        switch (op) {
            case 'c':
                args.collection_dir = optarg;
                break;
            case 'q':
                args.query_file = optarg;
                break;
            case 't':
                args.threshold = std::strtod(optarg,NULL);
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

void
output_qry(uint64_t qid,
		   const std::string& method,
		   const std::vector<std::vector<uint64_t>> parsed_qry,
		   surf::query_parser::mapping_t& tm,
		   double threshold)
{
	std::cout << qid << ";";
	std::cout << method << ";";
	std::cout << threshold << ";";
	const auto& reverse_mapping = tm.second;
	for(const auto& token : parsed_qry) {
		std::cout << "[";
		bool first = true;
		for(const auto& id : token) {
			auto rmitr = reverse_mapping.find(id);
            if(rmitr != reverse_mapping.end()) {
            	if(!first) {
            		std::cout << " ";
            	}
            	std::cout << rmitr->second;
            }
            first = false;
		}
		std::cout << "]";
	}
	std::cout << std::endl;
}

int main(int argc,char* const argv[])
{
    /* parse command line */
    cmdargs_t args = parse_args(argc,argv);

    /* parse repo */
    auto cc = surf::parse_collection(args.collection_dir);

    /* parse queries */
    surf::query_parser::mapping_t term_map;
    std::cerr << "Loading dictionary and creating term map." << std::endl;
    term_map = surf::query_parser::load_dictionary(args.collection_dir);

    /* parse queries */
    std::cerr << "Parsing query file '" << args.query_file << "'" << std::endl;
    auto queries = surf::query_parser::parse_to_ids(args.collection_dir,args.query_file,term_map);
    std::cerr << "Found " << queries.size() << " queries." << std::endl;

    /* load or construct the csa */
    using csa_type = sdsl::csa_wt<sdsl::wt_int<sdsl::rrr_vector<63>>,1000000,1000000>;
    csa_type csa;
    if ( !sdsl::cache_file_exists<csa_type>(surf::KEY_CSA, cc) ) {
    	std::cerr << "Constructing CSA " << std::endl;	
        sdsl::construct(csa, "", cc, 0);
        std::cerr << "Storing CSA " << std::endl;	
        sdsl::store_to_cache(csa, surf::KEY_CSA, cc, true);
    } else {
	    std::cerr << "Loading CSA " << std::endl;	
    	sdsl::load_from_cache(csa, surf::KEY_CSA, cc, true);
    }
    
    /* phrase parsing */
    std::cerr << "Parsing Queries... " << std::endl;
    std::cout << "qryid;method;threshold;parsing" << std::endl;
    for(const auto& qry : queries) {
		auto parsed_query_none = surf::phrase_detector::parse_none(csa,qry.second,args.threshold);
    	output_qry(qry.first,"NO-PARSE",parsed_query_none,term_map,args.threshold);

    	auto parsed_query_greedy = surf::phrase_detector::parse_greedy_lr(csa,qry.second,args.threshold);
    	output_qry(qry.first,"GREEDY",parsed_query_greedy,term_map,args.threshold);

    	auto parsed_query_gpaul = surf::phrase_detector::parse_greedy_paul(csa,qry.second,args.threshold);
    	output_qry(qry.first,"GREEDY-PAUL",parsed_query_gpaul,term_map,args.threshold);

    	auto parsed_query_dp = surf::phrase_detector::parse_dp(csa,qry.second,args.threshold);
    	output_qry(qry.first,"DP",parsed_query_dp,term_map,args.threshold);
    }

    return EXIT_SUCCESS;
}
