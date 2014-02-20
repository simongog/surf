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

bool
directory_exists(std::string dir)
{
    struct stat sb;
    const char* pathname = dir.c_str();
    if (stat(pathname, &sb) == 0 && S_ISDIR(sb.st_mode)) {
        return true;
    }
    return false;
}

void
create_directory(std::string dir)
{
    if (!directory_exists(dir)) {
        if (mkdir(dir.c_str(),0777) == -1) {
            perror("could not create directory");
            exit(EXIT_FAILURE);
        }
    }
}

sdsl::cache_config
parse_collection(std::string collection_dir)
{
    sdsl::cache_config config;

    /* check if all the directories exist */
    if (! directory_exists(collection_dir)) {
        std::cerr << collection_dir << " is not a valid directory.\n";
        exit(EXIT_FAILURE);
    }
    create_directory(collection_dir+"/tmp/");
    create_directory(collection_dir+"/index/");

    /* populate cache config */


    return config;
}

int main(int argc,char* const argv[])
{
    using clock = std::chrono::high_resolution_clock;
    /* parse command line */
    cmdargs_t args = parse_args(argc,argv);

    /* parse repo */
    sdsl::cache_config config = parse_collection(args.collection_dir);

    /* define types */
    using surf_index_t = INDEX_TYPE;
    std::string index_name = IDXNAME;

    /* build the index */
    auto build_start = clock::now();
    surf_index_t index(config);
    auto build_stop = clock::now();
    auto build_time_sec = std::chrono::duration_cast<std::chrono::seconds>(build_stop-build_start);
    std::cout << "Index built in " << build_time_sec.count() << " seconds." << std::endl;

    // write index
    std::string output_file = args.collection_dir + "/index/" + index_name;
    std::cout << "Writing index to file " << output_file << std::endl;
    auto write_start = clock::now();
    std::ofstream ofile(output_file);
    index.serialize(ofile);
    ofile.close();
    auto write_stop = clock::now();
    auto write_time_sec = std::chrono::duration_cast<std::chrono::seconds>(write_stop-write_start);
    std::cout << "Index written to disk in " << write_time_sec.count() << " seconds." << std::endl;

    return EXIT_SUCCESS;
}
