#include <unistd.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <chrono>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "surf/comm.hpp"
#include "surf/util.hpp"
#include "surf/query_parser.hpp"

#include "zmq.hpp"

typedef struct cmdargs {
    std::string host;
    std::string query_file;
    uint64_t k;
    uint64_t runs;
    bool profile;
    bool quit;
    bool ranked_and;
    bool phrases;
    double phrase_threshold;
    bool output_results;
    bool integer_mode;
    std::string collection_dir;
} cmdargs_t;

void
print_usage(char* program)
{
    fprintf(stdout,"%s -h <host> -q <query file> -k <top-k> -r <runs> -p -P <thres> -s -a -R -i <collection>\n",program);
    fprintf(stdout,"where\n");
    fprintf(stdout,"  -h <host>  : host of the daemon.\n");
    fprintf(stdout,"  -q <query file>  : the queries to be performed.\n");
    fprintf(stdout,"  -k <top-k>  : the top-k documents to be retrieved for each query.\n");
    fprintf(stdout,"  -r <runs>  : the number of runs.\n");
    fprintf(stdout,"  -R : output results only\n");
    fprintf(stdout,"  -p : run queries in profile mode.\n");
    fprintf(stdout,"  -P <thres> : run queries with phrase parsing enabled and threshold <thres>.\n");
    fprintf(stdout,"  -s : stop the daemon after queries are processed.\n");
    fprintf(stdout,"  -a : perform ranked AND instead of ranked OR.\n");
    fprintf(stdout,"  -i : perform dict lookup at the client from <collection>.\n");
};

cmdargs_t
parse_args(int argc,char* const argv[])
{
    cmdargs_t args;
    int op;
    args.host = "127.0.0.1:12345";
    args.query_file = "";
    args.k = 10;    
    args.runs = 3;
    args.profile = false;
    args.quit = false;
    args.ranked_and = false;
    args.phrases = false;
    args.phrase_threshold = 0.0f;
    args.output_results = false;
    args.integer_mode = false;
    while ((op=getopt(argc,argv,"r:h:q:k:psaP:Ri:")) != -1) {
        switch (op) {
            case 'r':
                args.runs = std::strtoul(optarg,NULL,10);
                break;
            case 'h':
                args.host = optarg;
                break;
            case 'p':
                args.profile = true;
                break;
            case 'P':
                args.phrases = true;
                args.phrase_threshold = std::strtod(optarg,NULL);
                break;
            case 's':
                args.quit = true;
                break;
            case 'a':
                args.ranked_and = true;
                break;
            case 'R':
                args.output_results = true;
                break;
            case 'q':
                args.query_file = optarg;
                break;
            case 'k':
                args.k = std::strtoul(optarg,NULL,10);
                break;
            case 'i':
                args.integer_mode = true;
                args.collection_dir = optarg;
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

    if(args.integer_mode) {
        surf::parse_collection(args.collection_dir); // makes sure dir is valid
        std::cout << "Loading dictionary and creating term map." << std::endl;
        auto term_map = surf::query_parser::load_dictionary(args.collection_dir);
        const auto& id_mapping = term_map.first;
        std::vector<std::string> mapped_queries;
        for(auto& query: queries) {
            auto qry_mapping = surf::query_parser::map_to_ids(id_mapping,query,true,false);
            if(std::get<0>(qry_mapping)) {
                auto qid = std::get<1>(qry_mapping);
                auto qry_ids = std::get<2>(qry_mapping);
                std::string new_qry_str;
                new_qry_str += std::to_string(qid) + ";";
                for(size_t i=0;i<qry_ids.size()-1;i++) {
                    new_qry_str += std::to_string(qry_ids[i]) + " ";
                }
                // last one
                new_qry_str += std::to_string(qry_ids.back());
                mapped_queries.push_back(new_qry_str);
            }
        }
        queries = mapped_queries; // copy!!
    }

    /* zmq magic! */
    std::cerr << "Connecting to surf daemon." << std::endl;
    zmq::context_t context (1);
    zmq::socket_t socket (context, ZMQ_REQ);
    socket.connect (std::string("tcp://"+args.host).c_str());
    if(!socket.connected()) {
        std::cerr << "Error connecting to daemon." << std::endl;
    }

    /* process the queries */
    std::cerr << "Processing queries..." << std::endl;
    size_t num_runs = args.runs;
    for(size_t i=0;i<num_runs;i++) {
        for(const auto& query: queries) {

            auto req_start = clock::now();

            surf_qry_request surf_req;
            surf_req.type = REQ_TYPE_QRY_OR;
            uint8_t qry_mode = 0;
            if(args.ranked_and) {
                surf_req.type = REQ_TYPE_QRY_AND;
                qry_mode = 1;
            }

            surf_req.int_qry = 0;
            if(args.integer_mode) {
                surf_req.int_qry = 1;
            }

            if(args.phrases) {
                surf_req.phrases = 1;
                surf_req.phrase_threshold =  args.phrase_threshold;
                qry_mode += 2;
            } else {
                surf_req.phrases = 0;
                surf_req.phrase_threshold = 0.0;
            }

            if(args.profile) {
                surf_req.mode = REQ_MODE_PROFILE;
            } else {
                surf_req.mode = REQ_MODE_TIME;
            }

            if(args.output_results) {
                surf_req.output_results = 1;
            } else {
                surf_req.output_results = 0;
            }

            surf_req.id = rand();
            surf_req.k = args.k;
            memcpy(surf_req.qry_str,query.data(),query.size());

            zmq::message_t request(sizeof(surf_qry_request));
            memcpy ((void *) request.data (), &surf_req, sizeof(surf_qry_request));
            socket.send (request);

            /* wait for reply */
            if(!args.output_results) {
                zmq::message_t reply;
                socket.recv (&reply);
                surf_time_resp* surf_resp = static_cast<surf_time_resp*>(reply.data());

                auto req_stop = clock::now();
                auto req_time = std::chrono::duration_cast<std::chrono::microseconds>(req_stop-req_start);

                if(surf_resp->req_id != surf_req.id) {
                    std::cerr << "ERROR: got response for wrong request id!" << std::endl;
                }

                if(surf_resp->status != REQ_PARSE_ERROR) {
                    /* output */
                    std::cout << surf_resp->qry_id << ";" 
                              << surf_resp->collection << ";"
                              << surf_resp->ranker << ";"
                              << surf_resp->index << ";"
                              << (int)qry_mode << ";"
                              << surf_resp->k << ";"
                              << surf_resp->qry_len << ";"
                              << surf_resp->result_size << ";"
                              << surf_resp->qry_time << ";"
                              << surf_resp->search_time << ";"
                              << surf_resp->wt_search_space << ";"
                              << surf_resp->wt_nodes << ";"
                              << surf_resp->postings_evaluated << ";"
                              << surf_resp->postings_total << ";"
                              << req_time.count() << std::endl;
                } else {
                    std::cerr << "Error processing query '" << query << "'" << std::endl;
                }
            } else {
                zmq::message_t output;
                socket.recv (&output);
                surf_results* sr = (surf_results*)(output.data());
                for(size_t j=0;j<sr->size;j++) {
                    std::cout << "(" << j+1 << ") : " 
                              << (uint64_t)sr->data[j*2]
                              << " - "
                              << sr->data[j*2+1] << std::endl;
                }
            }
        }
    }

    // stop the daemon
    if(args.quit) {
        surf_qry_request surf_req;
        surf_req.type = REQ_TYPE_QUIT;
        zmq::message_t request(sizeof(surf_qry_request));
        memcpy ((void *) request.data (), &surf_req, sizeof(surf_qry_request));
        socket.send (request);
    }


    return EXIT_SUCCESS;
}
