#include <unistd.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <ctime>

#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "surf/query.hpp"
#include "sdsl/config.hpp"
#include "surf/indexes.hpp"
#include "surf/query_parser.hpp"

#include "zmq.hpp"

typedef struct cmdargs {
    std::string collection_dir;
} cmdargs_t;

void
print_usage(char* program)
{
    fprintf(stdout,"%s -c <collection directory> -q <query file> -k <top-k> -o <output.csv>\n",program);
    fprintf(stdout,"where\n");
    fprintf(stdout,"  -c <collection directory>  : the directory the collection is stored.\n");
};

cmdargs_t
parse_args(int argc,char* const argv[])
{
    cmdargs_t args;
    int op;
    args.collection_dir = "";
    while ((op=getopt(argc,argv,"c:q:k:")) != -1) {
        switch (op) {
            case 'c':
                args.collection_dir = optarg;
                break;
            case '?':
            default:
                print_usage(argv[0]);
        }
    }
    if (args.collection_dir=="") {
        std::cerr << "Missing command line parameters.\n";
        print_usage(argv[0]);
        exit(EXIT_FAILURE);
    }
    return args;
}

std::tuple<uint64_t,uint64_t,std::string>
parse_request(const void* req_data,size_t n)
{
	std::string req(static_cast<const char*>(req_data),n);
	std::cout << "processing request '" << req << "'" << std::endl;

    auto id_sep_pos = req.find(';');
    auto reqid_str = req.substr(0,id_sep_pos);
    uint64_t req_id = std::stoull(reqid_str);

    auto k_sep_pos = req.find(';',id_sep_pos+1);
    auto k_str = req.substr(id_sep_pos+1,k_sep_pos-(id_sep_pos+1));
    uint64_t k = std::stoull(k_str);

    std::string qry_str = req.substr(k_sep_pos+1);

	return make_tuple(req_id,k,qry_str);
}

int main(int argc,char* const argv[])
{
    using clock = std::chrono::high_resolution_clock;
    /* parse command line */
    cmdargs_t args = parse_args(argc,argv);

    /* parse repo */
    auto cc = surf::parse_collection(args.collection_dir);

    /* parse queries */
    std::cout << "Loading dictionary and creating term map.";
    auto term_map = surf::query_parser::load_dictionary(args.collection_dir);

    /* define types */
    using surf_index_t = INDEX_TYPE;
    std::string index_name = IDXNAME;

    /* load the index */
    std::cout << "Loading index." << std::endl;
    surf_index_t index;
    auto load_start = clock::now();
    construct(index, "", cc, 0);
    index.load(cc);
    auto load_stop = clock::now();
    auto load_time_sec = std::chrono::duration_cast<std::chrono::seconds>(load_stop-load_start);
    std::cout << "Index loaded in " << load_time_sec.count() << " seconds." << std::endl;


    /* daemon mode */
    {
    	std::cout << "Starting daemon mode on port 12345" << std::endl;
    	zmq::context_t context(1);
    	zmq::socket_t server(context, ZMQ_REP);
    	server.bind("tcp://*:12345");

    	while(true) {
    		zmq::message_t request;
    		/* wait for msg */
    		server.recv(&request);

    		auto parsed_request = parse_request(request.data(),request.size());
    		auto k = std::get<1>(parsed_request);
    		auto query_str = std::get<2>(parsed_request);
    		auto req_id = std::get<0>(parsed_request);

    		/* perform query */
    		auto qry_start = clock::now();

    		/* (1) parse qry terms */
    		auto parsed_query = surf::query_parser::parse_query(term_map,query_str);
    		if(!parsed_query.first) {
    			// error
    			zmq::message_t reply (5);
    			memcpy(reply.data(),"ERROR",5);
    			server.send (reply);
    		} else {
	    		/* (2) query the index */
	            auto qry_id = std::get<0>(parsed_query.second);
	            auto qry_tokens = std::get<1>(parsed_query.second);
	            auto search_start = clock::now();
	            auto results = index.search(qry_tokens,k);
	            auto search_stop = clock::now();
	            auto search_time = std::chrono::duration_cast<std::chrono::microseconds>(search_stop-search_start);

	    		auto qry_stop = clock::now();
	    		auto query_time = std::chrono::duration_cast<std::chrono::microseconds>(qry_stop-qry_start);

	    		/* (3) create answer */
	    		std::string response = std::to_string(req_id) + ";" +
	    							   std::to_string(k) + ";" +
	    							   std::to_string(qry_id) + ";" +
	    							   std::to_string(query_time.count()) + ";" +
	    							   std::to_string(search_time.count());

	    		zmq::message_t reply (response.size());
	    		memcpy(reply.data(),response.data(),response.size());
	    		server.send (reply);
	    	}
    	}
    }


    return EXIT_SUCCESS;
}
