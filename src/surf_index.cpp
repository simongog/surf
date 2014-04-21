
#include "sdsl/config.hpp"
#include "surf/indexes.hpp"
#include "surf/util.hpp"

typedef struct cmdargs {
    std::string collection_dir;
    bool print_memusage;
    bool byte_alphabet;
} cmdargs_t;

void
print_usage(char* program)
{
    fprintf(stdout,"%s -c <collection directory> -m\n",program);
    fprintf(stdout,"where\n");
    fprintf(stdout,"  -c <collection directory>  : the directory the collection is stored.\n");
    fprintf(stdout,"  -m : print memory usage.\n");
};

cmdargs_t
parse_args(int argc,char* const argv[])
{
    cmdargs_t args;
    int op;
    args.collection_dir = "";
    args.print_memusage = false;
    args.byte_alphabet  = false;
    while ((op=getopt(argc,argv,"c:m:b")) != -1) {
        switch (op) {
            case 'c':
                args.collection_dir = optarg;
                break;
            case 'm':
                args.print_memusage = true;
                break;
            case 'b':
                args.byte_alphabet = true;
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
    using surf_index_t = INDEX_TYPE;
    /* parse command line */
    cmdargs_t args = parse_args(argc,argv);

    /* parse repo */
    sdsl::cache_config cc = surf::parse_collection<typename surf_index_t::csa_type::alphabet_type>(args.collection_dir);
    std::cout<<"parse collections"<<std::endl;
    for(auto x : cc.file_map){
        std::cout<<x.first<<" "<<x.second<<std::endl;
    }

    /* define types */
    std::string index_name = IDXNAME;

    /* build the index */
    surf_index_t index;
    auto build_start = clock::now();
    construct(index, "", cc, args.byte_alphabet?1:0);
    auto build_stop = clock::now();
    auto build_time_sec = std::chrono::duration_cast<std::chrono::seconds>(build_stop-build_start);
    std::cout << "Index built in " << build_time_sec.count() << " seconds." << std::endl;

    /* visualize space usage */
    index.load(cc);
    std::cout<<"Write structure"<<std::endl;
    std::ofstream vofs(args.collection_dir+"/index/"+surf::SPACEUSAGE_FILENAME+"_"+IDXNAME+".html");
    write_structure<HTML_FORMAT>(index,vofs);

    /* print mem usage */
    if(args.print_memusage) {
        index.mem_info();
    }

    return EXIT_SUCCESS;
}
