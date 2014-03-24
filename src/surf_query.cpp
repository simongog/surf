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
    /* parse command line */
    cmdargs_t args = parse_args(argc,argv);

    /* load queries from disk */
    std::cout << "Loading queries from disk." << std::endl;
    std::ifstream qfs(args.query_file);
    std::string qry_str;
    std::vector<std::string> queries;
    while(std::getline(qfs,qry_str)) {
        queries.push_back(qry_str);
    }

    /* zmq magic! */
    std::cout << "Connecting to surf daemon." << std::endl;
    zmq::context_t context (1);
    zmq::socket_t socket (context, ZMQ_REQ);
    socket.connect ("tcp://127.0.0.1:12345");
    if(!socket.connected()) {
        std::cerr << "Error connecting to daemon." << std::endl;
    }

    /* process the queries */
    std::cout << "Processing queries..." << std::endl;
    std::map<uint64_t,std::string> query_results;
    size_t num_runs = 3;
    for(size_t i=0;i<num_runs;i++) {
        for(const auto& query: queries) {
            uint64_t req_id = rand();
            std::string req = std::to_string(req_id) + ";" + std::to_string(args.k) + ";" + query;

            zmq::message_t request(req.size());
            memcpy ((void *) request.data (), req.data(), req.size());
            socket.send (request);

            /* wait for reply */
            zmq::message_t reply;
            socket.recv (&reply);
            std::string response(static_cast<const char*>(reply.data()),reply.size());
            std::cout << "GOT '" << response << "'" << std::endl;
        }
    }


    return EXIT_SUCCESS;
}
