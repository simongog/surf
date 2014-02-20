#include <unistd.h>
#include <stdlib.h>
#include <iostream>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "sdsl/config.hpp"
#include "surf/indexes.hpp"

typedef struct cmdargs {
    std::string collection_dir;
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

int main(int argc,char* const argv[])
{
    using clock = std::chrono::high_resolution_clock;
    /* parse command line */
    cmdargs_t args = parse_args(argc,argv);

    /* define types */
    using surf_index_t = INDEX_TYPE;
    std::string index_name = IDXNAME;

    /* load the index */
    std::string index_file = args.collection_dir = "/index/" + index_name;
    std::ifstream ifs(index_file);
    if (ifs.is_open()) {
        auto load_start = clock::now();
        surf_index_t index(ifs);
        auto load_stop = clock::now();
        auto load_time_sec = std::chrono::duration_cast<std::chrono::seconds>(load_stop-load_start);
        std::cout << "Index loaded in " << load_time_sec.count() << " seconds." << std::endl;

        /* query stuff */

    } else {
        std::cerr << "Can not load index from file " << index_file << std::endl;
    }

    return EXIT_SUCCESS;
}
