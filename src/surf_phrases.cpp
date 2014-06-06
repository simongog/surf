#include <iostream>
#include <cstdlib>
#include <unistd.h>
#include <future>


#include "sdsl/suffix_arrays.hpp"
#include "surf/query_parser.hpp"
#include "surf/phrase_detector.hpp"
#include "surf/util.hpp"
#include "surf/comm.hpp"

#include "zmq.hpp"

namespace std
{
    template<>
    struct hash<std::vector<uint64_t>>
    {
        typedef std::vector<uint64_t> argument_type;
        typedef std::size_t value_type;
 
        value_type operator()(argument_type const& s) const
        {
            std::hash<uint64_t> hash_fn;
            value_type hash = 4711;
            for(const auto& id : s) {
                hash ^= hash_fn(id);
            }
            return hash;
        }
    };
}


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

using res_heap = std::priority_queue<std::pair<double,std::vector<uint64_t>>>;

struct res_heaps {
    res_heap greedy_lr;
    std::unordered_set<std::vector<uint64_t>> greedy_lr_set;
    res_heap bm25;
    std::unordered_set<std::vector<uint64_t>> bm25_set;
    res_heap x2;
    std::unordered_set<std::vector<uint64_t>> x2_set;
    res_heap exist_prob;
    std::unordered_set<std::vector<uint64_t>> exist_prob_set;
};

template<class t_itr>
res_heaps
process_queries(t_itr begin,t_itr end,std::string host,double threshold) {
    res_heaps heaps;

    /* zmq magic! */
    std::cerr << "Connecting to surf daemon." << std::endl;
    zmq::context_t context (1);
    zmq::socket_t socket (context, ZMQ_REQ);
    socket.connect (std::string("tcp://"+host).c_str());
    if(!socket.connected()) {
        std::cerr << "Error connecting to daemon." << std::endl;
    }

    std::cerr << "Processing queries..." << std::endl;
    zmq_index index(socket);
    auto itr = begin;
    while(itr != end) {
        auto query = *itr;
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
        if(surf_resp->nids != 0) {
            std::vector<uint64_t> qry_ids(surf_resp->nids);
            std::copy(std::begin(surf_resp->ids),std::begin(surf_resp->ids)+surf_resp->nids,
                  qry_ids.begin());

            // perform phrase stuff
            surf::phrase_detector::parse_greedy_lr(index,qry_ids,threshold,heaps.greedy_lr);
            surf::phrase_detector::parse_x2(index,qry_ids,threshold,heaps.x2);
            //surf::phrase_detector::parse_greedy_x2(index,qry_ids,args.threshold,heap_greedy_x2);
            surf::phrase_detector::parse_bm25(index,qry_ids,threshold,heaps.bm25);
            surf::phrase_detector::parse_exist_prob(index,qry_ids,threshold,heaps.exist_prob);
        }
        itr++;
    }

    return heaps;
}

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

    size_t num_threads = 10;
    size_t n = queries.size();
    off_t block_size = n / num_threads + 1;
    std::vector<std::future<res_heaps>> v;
    for (size_t i=0; i<n; i+=block_size) {
        v.push_back(std::async(std::launch::async,[&,i] {
            auto begin = queries.begin()+i;
            size_t end_offset = std::min((size_t)n-1,(size_t)(i+block_size-1));
            auto end = queries.begin()+end_offset+1;
            return process_queries(begin,end,args.host,args.threshold);
        }));
    }

    // merge results
    res_heaps heaps;
    for (auto& f : v) {
        auto partial_heaps = f.get();
        while(!partial_heaps.greedy_lr.empty()) {
            auto top = partial_heaps.greedy_lr.top(); partial_heaps.greedy_lr.pop();
            if(heaps.greedy_lr_set.count(top.second) == 0) {
                heaps.greedy_lr.push(top);
                heaps.greedy_lr_set.emplace(top.second);
            }
        }
        while(!partial_heaps.bm25.empty()) {
            auto top = partial_heaps.bm25.top(); partial_heaps.bm25.pop();
            if(heaps.bm25_set.count(top.second) == 0) {
                heaps.bm25.push(top);
                heaps.bm25_set.emplace(top.second);
            }
        }
        while(!partial_heaps.x2.empty()) {
            auto top = partial_heaps.x2.top(); partial_heaps.x2.pop();
            if(heaps.x2_set.count(top.second) == 0) {
                heaps.x2.push(top);
                heaps.x2_set.emplace(top.second);
            }
        }
        while(!partial_heaps.exist_prob.empty()) {
            auto top = partial_heaps.exist_prob.top(); partial_heaps.exist_prob.pop();
            if(heaps.exist_prob_set.count(top.second) == 0) {
                heaps.exist_prob.push(top);
                heaps.exist_prob_set.emplace(top.second);
            }
        }
    }

    // lookup
    zmq::context_t context (1);
    zmq::socket_t socket (context, ZMQ_REQ);
    socket.connect (std::string("tcp://"+args.host).c_str());
    if(!socket.connected()) {
        std::cerr << "Error connecting to daemon." << std::endl;
    }

    output_topk("GREEDY-LR",socket,heaps.greedy_lr,100);
    //output_topk("GREEDY-PAUL",socket,heap_greedy_paul,100);
    output_topk("X2",socket,heaps.x2,100);
    //output_topk("GREEDY-X2",socket,heap_greedy_x2,100);
    output_topk("BM25",socket,heaps.bm25,100);
    output_topk("EXIST-PROB",socket,heaps.exist_prob,100);

    return EXIT_SUCCESS;
}
