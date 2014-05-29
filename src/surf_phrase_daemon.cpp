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
#include "surf/phrase_parser.hpp"
#include "surf/rank_functions.hpp"

#include "zmq.hpp"

typedef struct cmdargs {
    std::string collection_dir;
    std::string port;
} cmdargs_t;

void
print_usage(char* program)
{
    fprintf(stdout,"%s -c <collection directory> -p <port> -r\n",program);
    fprintf(stdout,"where\n");
    fprintf(stdout,"  -c <collection directory>  : the directory the collection is stored.\n");
    fprintf(stdout,"  -p <port>  : the port the daemon is running on.\n");
};

cmdargs_t
parse_args(int argc,char* const argv[])
{
    cmdargs_t args;
    int op;
    args.collection_dir = "";
    args.port = std::to_string(12345);
    while ((op=getopt(argc,argv,"c:p")) != -1) {
        switch (op) {
            case 'c':
                args.collection_dir = optarg;
                break;
            case 'p':
                args.port = optarg;
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

int main(int argc,char* const argv[])
{
    /* parse command line */
    cmdargs_t args = parse_args(argc,argv);

    /* parse repo */
    auto cc = surf::parse_collection(args.collection_dir);
    char tmp_str[256] = {0};
    strncpy(tmp_str,args.collection_dir.c_str(),256);
    std::string base_name = basename(tmp_str);

    /* load dict */
    surf::query_parser::mapping_t term_map;
    std::cout << "Loading dictionary and creating term map." << std::endl;
    term_map = surf::query_parser::load_dictionary(args.collection_dir);
    const auto& id_mapping = term_map.first;
    const auto& reverse_mapping = term_map.second;

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

    /* daemon mode */
    {
    	std::cout << "Starting daemon mode on port " << args.port << std::endl;
    	zmq::context_t context(1);
    	zmq::socket_t server(context, ZMQ_REP);
    	server.bind(std::string("tcp://*:"+args.port).c_str());

    	while(true) {
    		zmq::message_t request;
    		/* wait for msg */
    		server.recv(&request);
            surf_phrase_request* surf_req = (surf_phrase_request*) request.data();

            surf_phrase_resp surf_resp;
            if(surf_req->type == REQ_TYPE_TERM2ID) {
                auto qry_mapping = surf::query_parser::map_to_ids(id_mapping,
                        std::string(surf_req->qry_str),true,false);
                if(std::get<0>(qry_mapping)) {
                    auto qid = std::get<1>(qry_mapping);
                    auto qry_ids = std::get<2>(qry_mapping);
	                surf_resp.qid = qid;
	                surf_resp.nids = qry_ids.size();
	                for(size_t i=0;i<surf_resp.nids;i++) 
	                	surf_resp.ids[i] = qry_ids[i];
                }

            }
            if(surf_req->type == REQ_TYPE_ID2TERM) {
            	auto id = surf_req->qids[0];
            	auto revitr = reverse_mapping.find(id);
	            if(revitr != reverse_mapping.end()) {
	            	std::copy(std::begin(revitr->second),std::end(revitr->second),
	            			  std::begin(surf_resp.term_str));
	            }
	        }
            if(surf_req->type == REQ_TYPE_COUNT) {
            	auto cnt = sdsl::count(csa,std::begin(surf_req->qids),
            							   std::begin(surf_req->qids)+surf_req->nids);
            	surf_resp.count = cnt;
            	surf_resp.size = csa.size();
            }

    		zmq::message_t reply (sizeof(surf_phrase_resp));
    		memcpy(reply.data(),&surf_resp,sizeof(surf_phrase_resp));
    		server.send (reply);
        }
    }

    return EXIT_SUCCESS;
}
