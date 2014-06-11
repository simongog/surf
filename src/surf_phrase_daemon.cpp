#include <unistd.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <thread>
#include <mutex>

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

std::unordered_map<surf::query_token,surf::result> score_cache;
std::mutex score_mutex;

void add_to_score_cache(const surf::query_token& q,const surf::result& r) {
    std::lock_guard<std::mutex> lock(score_mutex);
    if( score_cache.count(q) == 0 ) {
        score_cache[q] = r;
    }
}

bool is_score_cached(const surf::query_token& q) {
    std::lock_guard<std::mutex> lock(score_mutex);
    return score_cache.count(q) != 0;
}

surf::result get_cached_score(const surf::query_token& q) {
    std::lock_guard<std::mutex> lock(score_mutex);
    return score_cache.find(q)->second;
}

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

std::unordered_map<std::vector<uint64_t>,std::pair<double,double>> prob_cache;
std::mutex prob_mutex;

void add_to_prob_cache(const std::vector<uint64_t>& q,std::pair<double,double> r) {
    std::lock_guard<std::mutex> lock(prob_mutex);
    if( prob_cache.count(q) == 0 ) {
        prob_cache[q] = r;
    }
}

bool is_prob_cached(const std::vector<uint64_t>& q) {
    std::lock_guard<std::mutex> lock(prob_mutex);
    return prob_cache.count(q) != 0;
}

std::pair<double,double> get_cached_prob(const std::vector<uint64_t>& q) {
    std::lock_guard<std::mutex> lock(prob_mutex);
    return prob_cache.find(q)->second;
}


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

template<class t_index>
void worker(const t_index* index,const surf::query_parser::mapping_t* tm,zmq::context_t* context)
{
    const auto& id_mapping = tm->first;
    const auto& reverse_mapping = tm->second;
    zmq::socket_t socket(*context, ZMQ_REP);
    socket.connect("inproc://workers");
    while(true) {
        zmq::message_t request;
        socket.recv(&request);

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
            auto cnt = sdsl::count(index->m_csa,std::begin(surf_req->qids),
                                       std::begin(surf_req->qids)+surf_req->nids);
            surf_resp.count = cnt;
        }

        if(surf_req->type == REQ_TYPE_MAXSCORE) {
            std::vector<surf::query_token> qry_tokens;
            surf::query_token q;
            q.token_ids.resize(surf_req->nids);
            std::copy(std::begin(surf_req->qids),std::begin(surf_req->qids)+surf_req->nids,
                      q.token_ids.begin());

            bool found = false;
            if(surf_req->nids == 1) {
                if(is_score_cached(q)) {
                    auto res = get_cached_score(q);
                    size_t n = std::min((size_t)100,res.list.size());
                    surf_resp.nscores = n;
                    for(size_t i=0;i<n;i++) surf_resp.max_score[i] = res.list[i].score;
                    found = true;
                }
            }

            if(!found) {
                qry_tokens.push_back(q);
                auto results = index->search(qry_tokens,1,false,false);
                surf_resp.max_score[0] = 0.0f;
                surf_resp.nscores = 0;
                if(results.list.size() != 0) {
                    size_t n = std::min((size_t)1,results.list.size());
                    surf_resp.nscores = n;
                    for(size_t i=0;i<n;i++) surf_resp.max_score[i] = results.list[i].score;
                    if( surf_req->nids == 1 ) {
                        add_to_score_cache(q,results);
                    }
                }
            }
        }

        if(surf_req->type == REQ_TYPE_PHRASEPROB) {
            std::vector<uint64_t> ids;
            ids.resize(surf_req->nids);
            std::copy(std::begin(surf_req->qids),std::begin(surf_req->qids)+surf_req->nids,
                      ids.begin());
            if(is_prob_cached(ids)) {
                auto p = get_cached_prob(ids);
                surf_resp.phrase_prob = std::get<0>(p);
                surf_resp.single_prob = std::get<1>(p);
            } else {
                auto p = index->phrase_prob(ids);
                surf_resp.phrase_prob = std::get<0>(p);
                surf_resp.single_prob = std::get<1>(p);
                add_to_prob_cache(ids,p);
            }
        }

        surf_resp.size = index->m_csa.size();
        zmq::message_t reply (sizeof(surf_phrase_resp));
        memcpy(reply.data(),&surf_resp,sizeof(surf_phrase_resp));
        socket.send (reply);
    }
}

int main(int argc,char* const argv[])
{
#ifdef PHRASE_SUPPORT
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

    /* define types */
    using surf_index_t = INDEX_TYPE;
    std::string index_name = IDXNAME;

    /* load the index */
    std::cout << "Loading index." << std::endl;
    surf_index_t index;
    construct(index, "", cc, 0);
    index.load(cc);
    std::cout << "Index loaded." << std::endl;

    /* daemon mode */
    {
    	std::cout << "Starting daemon mode on port " << args.port << std::endl;
    	zmq::context_t context(1);
        zmq::socket_t clients(context, ZMQ_ROUTER);
        clients.bind(std::string("tcp://*:"+args.port).c_str());
        zmq::socket_t workers(context, ZMQ_DEALER);
        zmq_bind(workers, "inproc://workers");
        // create workers
        std::thread threads[10];
        for(size_t i=0;i<10;i++) {
            threads[i] = std::thread(worker<surf_index_t>,&index,&term_map,&context);
        }
        zmq_device(ZMQ_QUEUE, (void*)clients, (void*)workers);
    }
#endif

    return EXIT_SUCCESS;
}
