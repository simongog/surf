
#include <unistd.h>
#include <stdlib.h>
#include <algorithm>
#include <vector>
#include <iostream>
#include <sstream>

#include "surf/config.hpp"
#include "surf/util.hpp"
#include "surf/construct_doc_cnt.hpp"

typedef struct cmdargs {
    std::string collection_dir;
    std::string surf_file;
    std::string trec_file;
} cmdargs_t;

void
print_usage(char* program)
{
    fprintf(stdout,"%s -c <collection directory>\n",program);
    fprintf(stdout,"where\n");
    fprintf(stdout,"  -c <collection directory>  : the directory the collection is stored.\n");
};

cmdargs_t
parse_args(int argc,char* const argv[])
{
    cmdargs_t args;
    int op;
    args.collection_dir = "";
    while ((op=getopt(argc,argv,"c:")) != -1) {
        switch (op) {
            case 'c':
                args.collection_dir = optarg;
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

int main( int argc, char** argv ) {
    /* parse command line */
    cmdargs_t args = parse_args(argc,argv);

    /* parse repo */
    auto cc = surf::parse_collection(args.collection_dir);
    sdsl::int_vector_buffer<> T(args.collection_dir+"/"+surf::TEXT_FILENAME);
    std::cout << "n = |T|= " << T.size() << std::endl;
    surf::construct_doc_cnt<sdsl::int_alphabet_tag::WIDTH>(cc);
    uint64_t doc_cnt = 0;
    load_from_cache(doc_cnt, surf::KEY_DOCCNT, cc);
    std::cout << "number of documents = N = " << doc_cnt << std::endl;
	std::ifstream dic_fs(args.collection_dir+"/"+surf::DICT_FILENAME);
	std::string line;
	size_t num_terms = 0;
	while( std::getline(dic_fs,line) ) {
		num_terms++;
	}
	std::cout << "number of terms = sigma = " << num_terms << std::endl;
	std::cout << "avg document length = " << T.size() / doc_cnt << std::endl;
}


