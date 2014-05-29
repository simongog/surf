#ifndef SURF_COMM_H
#define SURF_COMM_H

#define REQ_PARSE_ERROR		0
#define REQ_RESPONE_OK		1

#define REQ_TYPE_QRY_OR		0
#define REQ_TYPE_QRY_AND	1
#define REQ_TYPE_QUIT		2

#define REQ_MODE_PROFILE	0
#define REQ_MODE_TIME		1

#define MAX_QRY_LEN		 1024

#define REQ_TYPE_TERM2ID     0
#define REQ_TYPE_ID2TERM     1
#define REQ_TYPE_COUNT       2

#pragma pack(push, 1)
struct surf_time_resp {
	uint8_t status;
    uint64_t req_id;
    uint64_t qry_id;
    uint64_t qry_len;
    uint64_t k;
    uint64_t result_size;
    uint64_t qry_time;
    uint64_t search_time;
    uint64_t wt_search_space;
    uint64_t wt_nodes;
    uint64_t postings_evaluated;
    uint64_t postings_total;
    char index[256];
    char collection[256];
    char ranker[256];
};

struct surf_phrase_request {
    uint8_t type;
    uint64_t term_id;
    uint64_t nids;
    char qry_str[1024] = {0};
    uint64_t qids[256] = {0};
};

struct surf_phrase_resp {
    uint64_t count;
    uint8_t term_str[1024] = {0};
    uint64_t ids[256] = {0};
    uint64_t nids;
    uint64_t qid;
};

struct surf_qry_request {
	uint8_t type;
	uint8_t mode;
    uint8_t phrases;
    double phrase_threshold;
	uint64_t id;
	uint64_t k;
    uint8_t output_results;
    uint8_t int_qry;
	char qry_str[MAX_QRY_LEN] = {0};
};

struct surf_results {
    uint64_t size;
    double data[0];
};
#pragma pack(pop)

#endif