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
    uint64_t k;
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
    args.k = 1000;
    while ((op=getopt(argc,argv,"h:q:t:k:")) != -1) {
        switch (op) {
            case 'h':
                args.host = optarg;
                break;
            case 'q':
                args.query_file = optarg;
                break;
            case 'k':
                args.k = std::stoul(optarg);
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

template<class t_method,class t_itr>
std::vector<std::pair<double,std::vector<uint64_t>>>
process_queries(t_itr begin,t_itr end,std::string host,double threshold) 
{
    std::vector<std::pair<double,std::vector<uint64_t>>> results;
    /* zmq magic! */
    zmq::context_t context (1);
    zmq::socket_t socket (context, ZMQ_REQ);
    socket.connect (std::string("tcp://"+host).c_str());
    if(!socket.connected()) {
        std::cerr << "Error connecting to daemon." << std::endl;
    }
    zmq_index index(socket);
    auto itr = begin;
    while(itr != end) {
        auto query = *itr;
        auto qry_res = t_method::parse(index,query,threshold);
        results.insert(results.end(), qry_res.begin(), qry_res.end());
        itr++;
    }
    return results;
}

template<class t_itr>
void output_results(const std::string& method,t_itr itr,t_itr end,std::string host,size_t k) 
{
    zmq::context_t context (1);
    zmq::socket_t socket (context, ZMQ_REQ);
    socket.connect (std::string("tcp://"+host).c_str());
    if(!socket.connected()) {
        std::cerr << "Error connecting to daemon." << std::endl;
        return;
    }

    double max_score = itr->first;
    size_t i=0;
    while(itr != end) {
        if(i ==k) break;
        double score = itr->first;
        if(score == 0) break;
        const auto& ids = itr->second;

        std::cout << std::setw(12) << method << " " 
                  << std::setw(6) << i+1 << "  " 
                  << std::setw(12) << score/max_score <<" [";

        bool first = true;
        for(const auto& id : ids) {
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

        i++;
        ++itr;
    }
}


template<class t_method>
void process_all_queries(const std::vector<std::vector<uint64_t>>& queries,std::string host,double threshold,size_t k)
{

    size_t num_threads = 10;
    size_t n = queries.size();
    off_t block_size = n / num_threads + 1;
    std::vector<std::future<std::vector<std::pair<double,std::vector<uint64_t>>>>> v;
    for (size_t i=0; i<n; i+=block_size) {
        v.push_back(std::async(std::launch::async,[&,i] {
            auto begin = queries.begin()+i;
            size_t end_offset = std::min((size_t)n-1,(size_t)(i+block_size-1));
            auto end = queries.begin()+end_offset+1;
            return process_queries<t_method>(begin,end,host,threshold);
        }));
    }

    // merge results
    std::vector<std::pair<double,std::vector<uint64_t>>> scores;
    for (auto& f : v) {
        auto partial_scores = f.get();
        scores.insert(scores.end(), partial_scores.begin(), partial_scores.end());
    }

    // sort
    std::sort(scores.begin(),scores.end(),std::greater<std::pair<double,std::vector<uint64_t>>>());
    auto last = std::unique(scores.begin(),scores.end());

    output_results(t_method::name(),scores.begin(),last,host,k);
}

int main(int argc,char* const argv[])
{
    /* parse command line */
    cmdargs_t args = parse_args(argc,argv);

    // lookup
    zmq::context_t context (1);
    zmq::socket_t socket (context, ZMQ_REQ);
    socket.connect (std::string("tcp://"+args.host).c_str());
    if(!socket.connected()) {
        std::cerr << "Error connecting to daemon." << std::endl;
        return EXIT_SUCCESS;
    }

    /* load queries from disk */
    std::cerr << "Loading queries from disk." << std::endl;
    std::ifstream qfs(args.query_file);
    std::string qry_str;
    std::vector<std::string> qrystrs;
    while(std::getline(qfs,qry_str)) {
        if(qry_str.size() < MAX_QRY_LEN) {
            qrystrs.push_back(qry_str);
        }
    }

    std::cerr << "Translate queries to ids" << std::endl;
    std::vector<std::vector<uint64_t>> queries;
    for(const auto& qstr : qrystrs) {
        surf_phrase_request surf_req;
        surf_req.type = REQ_TYPE_TERM2ID;
        memcpy(surf_req.qry_str,qstr.data(),qstr.size());
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
            queries.push_back(qry_ids);
        }
    }

    process_all_queries<surf::phrase_detector_sa_greedy<false>>(queries,args.host,args.threshold,args.k);
    process_all_queries<surf::phrase_detector_x2>(queries,args.host,args.threshold,args.k);
    process_all_queries<surf::phrase_detector_bm25<false>>(queries,args.host,args.threshold,args.k);
    process_all_queries<surf::phrase_detector_exist_prob<false>>(queries,args.host,args.threshold,args.k);
    process_all_queries<surf::phrase_detector_bm25<true>>(queries,args.host,args.threshold,args.k);
    process_all_queries<surf::phrase_detector_exist_prob<true>>(queries,args.host,args.threshold,args.k);
    process_all_queries<surf::phrase_detector_x2_greedy<false>>(queries,args.host,args.threshold,args.k);
    process_all_queries<surf::phrase_detector_x2_greedy<true>>(queries,args.host,args.threshold,args.k);

    return EXIT_SUCCESS;
}
