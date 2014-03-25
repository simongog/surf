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
#include "surf/comm.hpp"

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
    	server.bind("tcp://127.0.0.1:12345");

    	while(true) {
    		zmq::message_t request;
    		/* wait for msg */
    		server.recv(&request);
            surf_qry_request* surf_req = (surf_qry_request*) request.data();

            if(surf_req->req_type == REQ_TYPE_QUIT) {
                std::cout << "Quitting..." << std::endl;
                break;
            }

    		/* perform query */
    		auto qry_start = clock::now();

    		/* (1) parse qry terms */
    		auto parsed_query = surf::query_parser::parse_query(term_map,std::string(surf_req->qry_str));
    		if(!parsed_query.first) {
                // error parsing the qry. send back error
                surf_time_resp surf_resp;
                surf_resp.status = REQ_PARSE_ERROR;
                surf_resp.req_id = surf_req->req_id;
    			zmq::message_t reply (sizeof(surf_time_resp));
    			memcpy(reply.data(),&surf_resp,sizeof(surf_time_resp));
    			server.send(reply);
    		} else {
              
	    		/* (2) query the index */
	            auto qry_id = std::get<0>(parsed_query.second);
	            auto qry_tokens = std::get<1>(parsed_query.second);
	            auto search_start = clock::now();
	            auto results = index.search(qry_tokens,surf_req->k);
	            auto search_stop = clock::now();
	            auto search_time = std::chrono::duration_cast<std::chrono::microseconds>(search_stop-search_start);

	    		auto qry_stop = clock::now();
	    		auto query_time = std::chrono::duration_cast<std::chrono::microseconds>(qry_stop-qry_start);

                /* (3a) output to qry to console */
                auto qry_strtokens = std::get<2>(parsed_query.second);
                std::cout << "REQ=" << std::left << std::setw(10) << surf_req->req_id << " " 
                          << " k="  << std::setw(5) << surf_req->k 
                          << " QID=" << std::setw(5) << qry_id 
                          << " TIME=" << std::setw(7) << query_time.count()/1000.0;
                std::cout << " [";
                for(const auto& qst : qry_strtokens) {
                    std::cout << qst << " ";
                }
                std::cout << "]" << std::endl;

	    		/* (3) create answer and send */
                surf_time_resp surf_resp;
                surf_resp.status = REQ_RESPONE_OK;
                surf_resp.req_id = surf_req->req_id;
                surf_resp.k = surf_req->k;
                surf_resp.qry_id = qry_id;
                surf_resp.qry_len = qry_tokens.size();
                surf_resp.result_size = results.size();
                surf_resp.qry_time = query_time.count();
                surf_resp.search_time = search_time.count();

	    		zmq::message_t reply (sizeof(surf_time_resp));
	    		memcpy(reply.data(),&surf_resp,sizeof(surf_time_resp));
	    		server.send (reply);
	    	}
    	}
    }


    return EXIT_SUCCESS;
}
