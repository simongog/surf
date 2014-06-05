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
output_topk(const std::string& method,zmq::socket_t& socket,
            std::priority_queue<std::pair<double,std::vector<uint64_t>>> heap,size_t k) 
{
    if(heap.empty()) {
        std::cout << method << " -> NO PHRASES\n";
        return;
    }
    auto max_score = heap.top().first;
    for(size_t i=0;i<k;i++) {
        if(heap.empty()) break;
        auto top = heap.top(); heap.pop();
        const auto& tokens = top.second;
        auto score = top.first;
        std::cout << std::setw(12) << method<< " " << std::setw(6) << i+1 << "  " << std::setw(12) << score/max_score <<" [";
        bool first = true;
        for(const auto& id : tokens) {
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
        std::cout << "]" << std::endl;
    }
}

struct zmq_index {
    size_t m_remote_size = 0;
    zmq::socket_t& socket;
    zmq_index(zmq::socket_t& s) : socket(s) {}
    template<class t_itr>
    size_t csa_count(t_itr begin,t_itr end) {
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
    size_t csa_size() {
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
    template<class t_itr>
    double phrase_prob(t_itr begin,t_itr end) {
        // send count req
        surf_phrase_request surf_req;
        surf_req.type = REQ_TYPE_PHRASEPROB;
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
        return surf_resp->phrase_prob;
    }
    template<class t_itr>
    std::vector<double> max_sim_scores(t_itr begin,t_itr end) {
        // send count req
        surf_phrase_request surf_req;
        surf_req.type = REQ_TYPE_MAXSCORE;
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
        size_t num_scores = surf_resp->nscores;
        std::vector<double> scores(num_scores);
        std::copy(std::begin(surf_resp->max_score),std::begin(surf_resp->max_score)+num_scores,scores.begin());
        return scores;
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
    std::priority_queue<std::pair<double,std::vector<uint64_t>>>  heap_greedy_lr;
    std::priority_queue<std::pair<double,std::vector<uint64_t>>>  heap_greedy_paul;
    std::priority_queue<std::pair<double,std::vector<uint64_t>>>  heap_x2;
    std::priority_queue<std::pair<double,std::vector<uint64_t>>>  heap_greedy_x2;
    std::priority_queue<std::pair<double,std::vector<uint64_t>>>  heap_bm25;
    std::priority_queue<std::pair<double,std::vector<uint64_t>>>  heap_exist_prob;

    std::cerr << "Processing queries..." << std::endl;
    zmq_index index(socket);
    for(const auto& query: queries) {
        std::cout << "Processing '" << query << "'\n";
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
        surf::phrase_detector::parse_greedy_lr(index,qry_ids,args.threshold,heap_greedy_lr);
//        surf::phrase_detector::parse_greedy_paul(index,qry_ids,args.threshold,heap_greedy_paul);
        surf::phrase_detector::parse_x2(index,qry_ids,args.threshold,heap_x2);
        //surf::phrase_detector::parse_greedy_x2(index,qry_ids,args.threshold,heap_greedy_x2);
        surf::phrase_detector::parse_bm25(index,qry_ids,args.threshold,heap_bm25);
        surf::phrase_detector::parse_exist_prob(index,qry_ids,args.threshold,heap_exist_prob);
    }

    output_topk("GREEDY-LR",socket,heap_greedy_lr,100);
    output_topk("GREEDY-PAUL",socket,heap_greedy_paul,100);
    output_topk("X2",socket,heap_x2,100);
    output_topk("GREEDY-X2",socket,heap_greedy_x2,100);
    output_topk("BM25",socket,heap_bm25,100);
    output_topk("EXIST-PROB",socket,heap_exist_prob,100);

    return EXIT_SUCCESS;
}
