#include <unistd.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>

#include "surf/query.hpp"
#include "sdsl/config.hpp"
#include "surf/indexes.hpp"
#include "surf/query_parser.hpp"
#include "surf/comm.hpp"
#include "surf/phrase_parser.hpp"

#include "zmq.hpp"
#include <json/json.h>

typedef struct cmdargs {
    std::string collection_dir;
    std::string port;
    bool load_dictionary;
} cmdargs_t;

void
print_usage(char* program)
{
    fprintf(stdout,"%s -c <collection directory> -p <port> -r\n",program);
    fprintf(stdout,"where\n");
    fprintf(stdout,"  -c <collection directory>  : the directory the collection is stored.\n");
    fprintf(stdout,"  -p <port>  : the port the daemon is running on.\n");
    fprintf(stdout,"  -r : do not load the dictionary.\n");
};

cmdargs_t
parse_args(int argc,char* const argv[])
{
    cmdargs_t args;
    int op;
    args.collection_dir = "";
    args.port = std::to_string(12345);
    args.load_dictionary = true;
    while ((op=getopt(argc,argv,"c:p:r")) != -1) {
        switch (op) {
            case 'c':
                args.collection_dir = optarg;
                break;
            case 'p':
                args.port = optarg;
                break;
            case 'r':
                args.load_dictionary = false;
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

std::string
trim(const std::string& str, const std::string& whitespace = " \t")
{
    const auto strBegin = str.find_first_not_of(whitespace);
    if (strBegin == std::string::npos)
        return ""; // no content

    const auto strEnd = str.find_last_not_of(whitespace);
    const auto strRange = strEnd - strBegin + 1;

    return str.substr(strBegin, strRange);
}

int main(int argc,char* const argv[])
{
    using clock = std::chrono::high_resolution_clock;
    /* parse command line */
    cmdargs_t args = parse_args(argc,argv);

    /* parse repo */
    auto cc = surf::parse_collection(args.collection_dir);
    char tmp_str[256] = {0};
    strncpy(tmp_str,args.collection_dir.c_str(),256);
    std::string base_name = basename(tmp_str);

    /* load int vector */
    sdsl::int_vector<> text_col;
    std::ifstream ifs(args.collection_dir+"/"+surf::TEXT_FILENAME);
    if(ifs.is_open()) {
        text_col.load(ifs);
    } else {
        std::cerr << "ERROR: could not load collection file." << std::endl;
        return EXIT_FAILURE;
    }

    /* load doc names */
    ifs.close();
    ifs.clear();
    std::vector<std::string> docs;
    ifs.open(args.collection_dir+"/"+surf::DOCNAMES_FILENAME);
    std::string line;
    if(ifs.is_open()) {
        while(getline(ifs,line)){
            docs.push_back(line);
        }
    } else {
        std::cerr << "ERROR: could not load doc names file." << std::endl;
        return EXIT_FAILURE;
    }

    //std::cout << docs.at(0) << std::endl;

    /* parse queries */
    surf::query_parser::mapping_t term_map;
    if(args.load_dictionary) {
        std::cout << "Loading dictionary and creating term map." << std::endl;
        term_map = surf::query_parser::load_dictionary(args.collection_dir);
    }

    uint64_t eos = 0; // end of full document string
    uint64_t separator = 1; // end of document string
    uint64_t space; // Fix this static whitespace identifier

    const std::unordered_map<std::string,uint64_t>& id_mapping = term_map.first;
    auto id_itr = id_mapping.find("\0");
    if(id_itr != id_mapping.end()) {
        space = id_itr->second;
    }

    //std::cout << "whitespace: " << space << std::endl;

    uint x = 0;
    std::cout << "text (strings): ";
    for (uint i : text_col) {
        if (i == space) {
            std::cout << x << ":" << " ";
        } else if (i == separator) {
            std::cout << x << ":" << "# ";
        } else {
            auto id_itr = term_map.second.find(i);
            if(id_itr != term_map.second.end()) {
                std::cout << x << ":" << id_itr->second << " ";
            } else {
                std::cerr << "ERROR: could not find '" << i << "' in the dictionary." << std::endl;
            }
        }
        x++;
    }
    std::cout << std::endl;

    x = 0;
    std::cout << "text (ints): ";
    for (uint i : text_col) {
        std::cout << x << ":" << i << " ";
        x++;
    }
    std::cout << std::endl;


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
        std::cout << "Starting daemon mode on port " << args.port << std::endl;
        zmq::context_t context(1);
        zmq::socket_t server(context, ZMQ_REP);
        server.bind(std::string("tcp://*:"+args.port).c_str());

        while(true) {
            zmq::message_t request;
            /* wait for msg */
            server.recv(&request);

            std::string rpl = std::string(static_cast<char*>(request.data()), request.size());
            std::cout << rpl << std::endl;

            Json::Value root;
            Json::Reader reader;

            bool parsingSuccessful = reader.parse(rpl, root);
            if(!parsingSuccessful) {
                std::cout << "ERROR IN JSON PARSING PROCESS. SKIPPING QUERY" << reader.getFormattedErrorMessages();

                Json::Value res;
                res["status"] = REQ_JSON_PARSE_ERROR;
                Json::StyledWriter writer;
                std::string res_str = writer.write(res);
                zmq::message_t reply(res_str.length());
                memcpy(reply.data(), res_str.c_str(), res_str.length());
                server.send(reply);

                return 16;
            }

            surf_qry_request surf_req;
            surf_req.type = root.get("type", 0).asInt();
            surf_req.mode = root.get("mode", 0).asInt();
            surf_req.phrases = root.get("phrases", 0).asInt();
            surf_req.phrase_threshold = root.get("phrase_threshold", 0.0).asDouble();
            surf_req.id = root.get("id", 0).asUInt64();
            surf_req.k = root.get("k", 0).asUInt64();
            surf_req.output_results = root.get("output_results", 0).asInt();
            surf_req.int_qry = root.get("int_qry", 0).asInt();
            surf_req.proximity_num = root.get("proximity_num", 50).asInt();
            std::string qry_str = root.get("qry_str", "").asString();

            if(surf_req.type == REQ_TYPE_QUIT) {
                std::cout << "Quitting..." << std::endl;
                break;
            }

            bool parse_ok = true;

            std::istringstream buf(qry_str);
            std::istream_iterator<std::string> beg(buf), end;
            std::vector<std::string> tokens(beg, end); // done!
            std::vector<surf::query_token> q_ts;
            for(auto& s: tokens) {
                std::vector<uint64_t> token_ids;
                std::vector<std::string> token_strs;
                std::reverse(s.begin(), s.end());
                std::cout << s << std::endl;
                for(std::string::iterator it = s.begin(); it != s.end(); ++it) {
                    const char c = *it;
                    std::unordered_map<std::string,uint64_t> tm = term_map.first;
                    std::string s(1, c);
                    auto id_itr = tm.find(s);
                    if(id_itr != tm.end()) {
                        if (isspace(c))
                            token_ids.push_back(space);
                        else
                            token_ids.push_back(id_itr->second);
                    } else {
                        std::cerr << "ERROR: could not find '" << c << "' in the dictionary." << std::endl;
                        parse_ok = false;
                    }
                }
                surf::query_token q_t(token_ids, token_strs, token_ids.size());
                q_ts.push_back(q_t);
            }

            if(!parse_ok) {
                // error parsing the qry. send back error
                Json::Value res;
                res["status"] = REQ_PARSE_ERROR;
                res["req_id"] = surf_req.id;

                Json::StyledWriter writer;
                std::string res_str = writer.write(res);
                zmq::message_t reply(res_str.length());
                memcpy(reply.data(), res_str.c_str(), res_str.length());
                server.send(reply);

                std::cout << "ERROR IN QUERY PARSING PROCESS. SKIPPING QUERY" << std::endl;

                continue;
            }

            /* (1) parse qry terms */
            bool profile = false;
            bool ranked_and = false;
            if(surf_req.type == REQ_TYPE_QRY_AND) {
                ranked_and = true;
            }

            auto search_start = clock::now();

            // search
            auto results = index.search(q_ts,surf_req.k,ranked_and,profile);

            auto search_stop = clock::now();

            auto search_time = std::chrono::duration_cast<std::chrono::microseconds>(search_stop-search_start).count();

            // autocomplete
            // TODO: improve
            std::vector<std::vector<uint64_t>> autocompletes = index.autocomplete(q_ts.back(), space);

            Json::Value res;
            Json::Value resultsArray;

            for(size_t i=0;i<results.list.size();i++) {
                Json::Value result;
                uint64_t doc_id = results.list[i].doc_id;
                result["doc_id"] = doc_id;
                if (doc_id >= 0 && doc_id < docs.size())
                    result["doc_path"] = docs[doc_id];
                std::cout << "doc_id: " << results.list[i].doc_id << std::endl;
                result["score"] = results.list[i].score;

                ////////
                // proximity
                ////////
                Json::Value proximityArray;
                for (surf::prox q : results.list[i].query_proximities) {
                    /*surf::term_info ti = q.ti;
                    std::cout << ti.t.size() << std::endl;
                    for (uint64_t t : ti.t) {
                        std::cout << t << " ";
                    }
                    std::cout << std::endl;*/

                    std::string proximity;
                    Json::Value prox;
                    uint64_t idx = q.idx;
                    //std::cout << "proximity: " << idx << std::endl;
                    uint64_t c = text_col[idx];

                    // add idx itself
                    if (idx >= 0 && idx < text_col.size()) {
                        c = text_col[idx];
                        auto id_itr = term_map.second.find(c);
                        if(id_itr != term_map.second.end()) {
                            proximity.append(id_itr->second);
                        } else {
                            std::cerr << "ERROR: could not find '" << c << "' in the dictionary." << std::endl;
                        }
                    }

                    // set n chars before idx
                    uint64_t pred_offset = 1;
                    while (idx-pred_offset >= 0 && idx >= pred_offset &&
                           (c = text_col[idx-pred_offset]) != separator &&
                           (pred_offset <= surf_req.proximity_num || text_col[idx-pred_offset] != space)) {
                        if (c == space) {
                            proximity.insert(0, " ");
                        } else {
                            auto id_itr = term_map.second.find(c);
                            if (id_itr != term_map.second.end()) {
                                proximity.insert(0, id_itr->second);
                            } else {
                                std::cerr << "ERROR: could not find '" << c << "' in the dictionary." << std::endl;
                            }
                        }
                        pred_offset++;
                    }
                    pred_offset--;
                    while(proximity.front() == ' ') {
                        proximity.erase(0, 1);
                        pred_offset--;
                    }

                    // set n chars after idx
                    uint64_t succ_pred_offset = 1;
                    while (idx+succ_pred_offset < text_col.size() &&
                           (c = text_col[idx+succ_pred_offset]) != separator &&
                           (succ_pred_offset <= surf_req.proximity_num || text_col[idx+succ_pred_offset] != space)) {
                        if (c == space) {
                            proximity.append(" ");
                        } else {
                            auto id_itr = term_map.second.find(c);
                            if (id_itr != term_map.second.end()) {
                                proximity.append(id_itr->second);
                            } else {
                                std::cerr << "ERROR: could not find '" << c << "' in the dictionary." << std::endl;
                            }
                        }
                        succ_pred_offset++;
                    }
                    succ_pred_offset--;
                    while(proximity.back() == ' ') {
                        proximity.pop_back();
                        succ_pred_offset--;
                    }

                    prox["proximity_org"] = proximity;
                    std::reverse(proximity.begin(), proximity.end());
                    prox["pred_offset"] = pred_offset;
                    prox["succ_pred_offset"] = succ_pred_offset;
                    prox["proximity"] = proximity;
                    prox["start"] = (int)(succ_pred_offset - q.ti.t.size() + 1);
                    prox["offset"] = (int)q.ti.t.size()-1;
                    proximityArray.append(prox);
                }
                result["proximities"] = proximityArray;
                resultsArray.append(result);
            }

            // process autocompletes
            Json::Value autocompletesArray;
            for (std::vector<uint64_t> v : autocompletes) {
                std::string autocomplete;
                for (uint64_t c : v) {
                    if (c == space) {
                        autocomplete.append(" ");
                    } else {
                        auto id_itr = term_map.second.find(c);
                        if (id_itr != term_map.second.end()) {
                            autocomplete.append(id_itr->second);

                        } else {
                            std::cerr << "ERROR: could not find '" << c << "' in the dictionary." << std::endl;
                        }
                    }
                }
                std::reverse(autocomplete.begin(), autocomplete.end());
                autocompletesArray.append(autocomplete);
            }
            res["autocompletes"] = autocompletesArray;

            res["status"] = REQ_RESPONE_OK;
            res["req_id"] = surf_req.id;
            res["results"] = resultsArray;
            res["search_time"] = search_time;

            Json::StyledWriter writer;
            std::string res_str = writer.write(res);

            std::cout << res_str << std::endl;

            zmq::message_t reply(res_str.length());
            memcpy(reply.data(), res_str.c_str(), res_str.length());
            server.send(reply);
        }
    }

    return EXIT_SUCCESS;
}
