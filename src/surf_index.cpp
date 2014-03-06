
#include "sdsl/config.hpp"
#include "surf/indexes.hpp"
#include "surf/util.hpp"

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


sdsl::cache_config
parse_collection(std::string collection_dir)
{
    /* check if all the directories exist */
    if( !surf::util::valid_collection(collection_dir) ) {
        exit(EXIT_FAILURE);
    }

    std::string index_directory = collection_dir+"/index/";
    surf::util::create_directory(index_directory);

    /* populate cache config */
    sdsl::cache_config config(false,collection_dir+"/index/","SURF");

    /* create symlink to text in index directory */
    std::string symlink_name = cache_file_name(sdsl::conf::KEY_TEXT_INT,config);
    if( ! surf::util::symlink_exists(cache_file_name(sdsl::conf::KEY_TEXT_INT,config)) ) {
        std::string collection_file = collection_dir+"/"+surf::TEXT_FILENAME;
        char* col_file_absolute = realpath(collection_file.c_str(), NULL);
        if( symlink(col_file_absolute,symlink_name.c_str()) != 0) {
            perror("cannot create symlink to collection file in index directory");
            exit(EXIT_FAILURE);
        }
        free(col_file_absolute);
    }

    /* register files that are present */
    for(const auto& key : surf::storage_keys) {
        cout<<"key="<<key<<endl;
        register_cache_file(key,config);
    }

    return config;
}

int main(int argc,char* const argv[])
{
    using clock = std::chrono::high_resolution_clock;
    /* parse command line */
    cmdargs_t args = parse_args(argc,argv);

    /* parse repo */
    sdsl::cache_config cc = parse_collection(args.collection_dir);
    std::cout<<"parse collections"<<std::endl;
    for(auto x : cc.file_map){
        std::cout<<x.first<<" "<<x.second<<std::endl;
    }

    /* define types */
    using surf_index_t = INDEX_TYPE;
    std::string index_name = IDXNAME;

    /* build the index */
    surf_index_t index;
    auto build_start = clock::now();
    construct(index, "", cc, 0);
    auto build_stop = clock::now();
    auto build_time_sec = std::chrono::duration_cast<std::chrono::seconds>(build_stop-build_start);
    std::cout << "Index built in " << build_time_sec.count() << " seconds." << std::endl;

    return EXIT_SUCCESS;
}
