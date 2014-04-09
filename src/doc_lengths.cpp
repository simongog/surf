#include <unistd.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <ctime>

#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "surf/util.hpp"
#include "sdsl/config.hpp"
#include "surf/construct_doc_lengths.hpp"

typedef struct cmdargs {
    std::string collection_dir;
} cmdargs_t;

void
print_usage(char* program)
{
    fprintf(stdout,"%s -c <collection directory> -p <port> -r\n",program);
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

int main(int argc,char* const argv[])
{
    /* parse command line */
    cmdargs_t args = parse_args(argc,argv);

    /* parse repo */
    auto cc = surf::parse_collection(args.collection_dir);
    char tmp_str[256] = {0};
    strncpy(tmp_str,args.collection_dir.c_str(),256);
    std::string base_name = basename(tmp_str);

    sdsl::int_vector<> doc_lengths;
    if (!sdsl::cache_file_exists(surf::KEY_DOC_LENGTHS, cc)){
        surf::construct_doc_lengths<sdsl::int_alphabet_tag::WIDTH>(cc);
    }
    sdsl::load_from_cache(doc_lengths, surf::KEY_DOC_LENGTHS, cc);

    std::sort(doc_lengths.begin(),doc_lengths.end());

    std::cout << "count;len\n";
    auto cur = doc_lengths[0];
    size_t count = 1;
    for(size_t i=1;i<doc_lengths.size();i++) {
        if( doc_lengths[i] != cur) {
            std::cout << count << ";" << cur << "\n";
            cur = doc_lengths[i];
        }
        count++;
    }
    std::cout << count << ";" << cur << std::endl;

    return EXIT_SUCCESS;
}


