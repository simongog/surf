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
#include "surf/comm.hpp"

#include "zmq.hpp"

typedef struct cmdargs {
    std::string query_file;
    uint64_t k;
} cmdargs_t;

void
print_usage(char* program)
{
    fprintf(stdout,"%s -q <query file> -k <top-k> -o <output.csv>\n",program);
    fprintf(stdout,"where\n");
    fprintf(stdout,"  -q <query file>  : the queries to be performed.\n");
    fprintf(stdout,"  -k <top-k>  : the top-k documents to be retrieved for each query.\n");
};

cmdargs_t
parse_args(int argc,char* const argv[])
{
    cmdargs_t args;
    int op;
    args.query_file = "";
    args.k = 10;
    while ((op=getopt(argc,argv,"q:k:")) != -1) {
        switch (op) {
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
    if (args.query_file=="") {
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

    /* load queries from disk */
    std::cerr << "Loading queries from disk." << std::endl;
    std::ifstream qfs(args.query_file);
    std::string qry_str;
    std::vector<std::string> queries;
    while(std::getline(qfs,qry_str)) {
        if(qry_str.size() < MAX_QRY_LEN) {
            queries.push_back(qry_str);
        }
    }

    /* zmq magic! */
    std::cerr << "Connecting to surf daemon." << std::endl;
    zmq::context_t context (1);
    zmq::socket_t socket (context, ZMQ_REQ);
    socket.connect ("tcp://127.0.0.1:12345");
    if(!socket.connected()) {
        std::cerr << "Error connecting to daemon." << std::endl;
    }

    /* process the queries */
    std::cerr << "Processing queries..." << std::endl;
    std::map<uint64_t,std::string> query_results;
    size_t num_runs = 3;
    for(size_t i=0;i<num_runs;i++) {
        for(const auto& query: queries) {

            auto req_start = clock::now();

            surf_qry_request surf_req;
            surf_req.req_type = REQ_TYPE_QRY;
            surf_req.req_id = rand();
            surf_req.k = args.k;
            memcpy(surf_req.qry_str,query.data(),query.size());

            zmq::message_t request(sizeof(surf_qry_request));
            memcpy ((void *) request.data (), &surf_req, sizeof(surf_qry_request));
            socket.send (request);

            /* wait for reply */
            zmq::message_t reply;
            socket.recv (&reply);
            surf_time_resp* surf_resp = static_cast<surf_time_resp*>(reply.data());

            auto req_stop = clock::now();
            auto req_time = std::chrono::duration_cast<std::chrono::microseconds>(req_stop-req_start);

            if(surf_resp->req_id != surf_req.req_id) {
                std::cerr << "ERROR: got response for wrong request id!" << std::endl;
            }

            /* output */
            std::cout << surf_resp->qry_id << ";" 
                      << surf_resp->k << ";"
                      << surf_resp->qry_len << ";"
                      << surf_resp->result_size << ";"
                      << surf_resp->qry_time << ";"
                      << surf_resp->search_time << ";"
                      << req_time.count() << std::endl;
        }
    }


    return EXIT_SUCCESS;
}
