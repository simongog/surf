#include <iostream>
#include <cstdlib>
#include <unistd.h>

#include "sdsl/suffix_arrays.hpp"
#include "surf/query_parser.hpp"
#include "surf/phrase_detector.hpp"
#include "surf/util.hpp"
#include "surf/comm.hpp"

#include "zmq.hpp"

typedef struct cmdargs {
    std::string host;
    std::string query_file;
    double threshold;
} cmdargs_t;

void
print_usage(char* program)
{
    fprintf(stdout,"%s -h <host> -q <query file> -t <threshold>\n",program);
    fprintf(stdout,"where\n");
    fprintf(stdout,"  -h <host>  : host of the daemon.\n");
    fprintf(stdout,"  -q <query file>  : the queries to be performed.\n");
    fprintf(stdout,"  -t <threshold>  : threshold for phrase detector.\n");
};

cmdargs_t
parse_args(int argc,char* const argv[])
{
    cmdargs_t args;
    int op;
    args.host = "127.0.0.1:12345";
    args.query_file = "";
    args.threshold = 0.0;
    while ((op=getopt(argc,argv,"h:q:t:")) != -1) {
        switch (op) {
            case 'h':
                args.host = optarg;
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
    if (args.host==""||args.query_file=="") {
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
		   zmq::socket_t& socket,
		   double threshold)
{
	std::cout << qid << ";";
	std::cout << method << ";";
	std::cout << threshold << ";";
	for(const auto& token : parsed_qry) {
		std::cout << "[";
		bool first = true;
		for(const auto& id : token) {
            // lookup str
            surf_phrase_request surf_req;
            surf_req.type = REQ_TYPE_ID2TERM;
            surf_req.qids[0] = id;
            zmq::message_t request(sizeof(surf_phrase_request));
            memcpy ((void *) request.data (), &surf_req, sizeof(surf_phrase_request));
            socket.send (request);
            // get answer
            zmq::message_t reply;
            socket.recv (&reply);
            surf_phrase_resp* surf_resp = static_cast<surf_phrase_resp*>(reply.data());
        	if(!first) {
        		std::cout << " ";
        	}
            std::cout << surf_resp->term_str;
            first = false;
		}
		std::cout << "]";
	}
	std::cout << std::endl;
}


    struct zmq_csa {
        size_t m_remote_size = 0;
        zmq::socket_t& socket;
        zmq_csa(zmq::socket_t& s) : socket(s) {}
        template<class t_itr>
        size_t count(t_itr begin,t_itr end) {
            // send count req
            surf_phrase_request surf_req;
            surf_req.type = REQ_TYPE_COUNT;
            surf_req.nids = end-begin;
            std::copy(begin,end,std::begin(surf_req.qids));
            zmq::message_t request(sizeof(surf_phrase_request));
            memcpy ((void *) request.data (), &surf_req, sizeof(surf_phrase_request));
            socket.send (request);
            // get answer
            zmq::message_t reply;
            socket.recv (&reply);
            surf_phrase_resp* surf_resp = static_cast<surf_phrase_resp*>(reply.data());
            m_remote_size = surf_resp->size;
            return surf_resp->count;
        }
        size_t size() {
            if(m_remote_size==0) {
                surf_phrase_request surf_req;
                surf_req.type = REQ_TYPE_COUNT;
                surf_req.nids = 1;
                surf_req.qids[0] = 1;
                zmq::message_t request(sizeof(surf_phrase_request));
                memcpy ((void *) request.data (), 
                        &surf_req, sizeof(surf_phrase_request));
                socket.send (request);
                // get answer
                zmq::message_t reply;
                socket.recv (&reply);
                surf_phrase_resp* surf_resp = static_cast<surf_phrase_resp*>(reply.data());
                m_remote_size = surf_resp->size;
            }
            return m_remote_size;
        }
    };

int main(int argc,char* const argv[])
{
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
    socket.connect (std::string("tcp://"+args.host).c_str());
    if(!socket.connected()) {
        std::cerr << "Error connecting to daemon." << std::endl;
    }


    /* process the queries */
    std::cerr << "Processing queries..." << std::endl;
    zmq_csa csa(socket);
    for(const auto& query: queries) {
        // get ids
        surf_phrase_request surf_req;
        surf_req.type = REQ_TYPE_TERM2ID;
        memcpy(surf_req.qry_str,query.data(),query.size());
        zmq::message_t request(sizeof(surf_phrase_request));
        memcpy ((void *) request.data (), &surf_req, sizeof(surf_phrase_request));
        socket.send (request);
        zmq::message_t reply;
        socket.recv (&reply);
        surf_phrase_resp* surf_resp = static_cast<surf_phrase_resp*>(reply.data());
        std::vector<uint64_t> qry_ids(surf_resp->nids);
        std::copy(std::begin(surf_resp->ids),std::begin(surf_resp->ids)+surf_resp->nids,
                  qry_ids.begin());

        // perform phrase stuff
        auto parsed_query_none = surf::phrase_detector::parse_none(csa,qry_ids,args.threshold);
        output_qry(surf_resp->qid,"NO-PARSE",parsed_query_none,socket,args.threshold);

        auto parsed_query_greedy = surf::phrase_detector::parse_greedy_lr(csa,qry_ids,args.threshold);
        output_qry(surf_resp->qid,"GREEDY",parsed_query_greedy,socket,args.threshold);

        auto parsed_query_gpaul = surf::phrase_detector::parse_greedy_paul(csa,qry_ids,args.threshold);
        output_qry(surf_resp->qid,"GREEDY-PAUL",parsed_query_gpaul,socket,args.threshold);

        auto parsed_query_dp = surf::phrase_detector::parse_dp(csa,qry_ids,args.threshold);
        output_qry(surf_resp->qid,"DP",parsed_query_dp,socket,args.threshold);
    }


    return EXIT_SUCCESS;
}
