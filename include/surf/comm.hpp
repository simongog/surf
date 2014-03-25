#ifndef SURF_COMM_H
#define SURF_COMM_H

#define REQ_PARSE_ERROR		0
#define REQ_RESPONE_OK		1

#define REQ_TYPE_QRY		0
#define REQ_TYPE_PROFILE	1
#define REQ_TYPE_QUIT		2

#define MAX_QRY_LEN		 1024

struct surf_time_resp {
	uint8_t status;
    uint64_t req_id;
    uint64_t qry_id;
    uint64_t qry_len;
    uint64_t k;
    uint64_t result_size;
    uint64_t qry_time;
    uint64_t search_time;
};

struct surf_result_resp {
	uint64_t k;
	uint8_t* data;
};

struct surf_qry_request {
	uint8_t req_type;
	uint64_t req_id;
	uint64_t k;
	char qry_str[MAX_QRY_LEN] = {0};
};


#endif