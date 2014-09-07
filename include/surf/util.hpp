#ifndef SURF_UTIL_HPP
#define SURF_UTIL_HPP

#include "surf/config.hpp"
#include "sdsl/io.hpp"

#include <string>
#include <unistd.h>
#include <stdlib.h>
#include <iostream>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

namespace surf{

bool
directory_exists(std::string dir)
{
    struct stat sb;
    const char* pathname = dir.c_str();
    if (stat(pathname, &sb) == 0 && (S_IFDIR&sb.st_mode)) {
        return true;
    }
    return false;
}

bool
file_exists(std::string file_name)
{
    sdsl::isfstream in(file_name);
    if (in) {
        in.close();
        return true;
    }
    return false;
}

bool
symlink_exists(std::string file)
{
    struct stat sb;
    const char* filename = file.c_str();
    if (stat(filename, &sb) == 0 && (S_IFLNK&sb.st_mode) ) {
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

template<typename alphabet_tag=sdsl::int_alphabet_tag>
bool
valid_collection(std::string collection_dir)
{
    if (! surf::directory_exists(collection_dir)) {
        std::cerr << collection_dir << " is not a valid directory.\n";
        return false;
    } else {
        /* make sure the necessary files are present */
        if ( std::is_same<alphabet_tag, sdsl::int_alphabet_tag>::value ){
            if( ! surf::file_exists(collection_dir+"/"+surf::TEXT_FILENAME) ||
                ! surf::file_exists(collection_dir+"/"+surf::DICT_FILENAME) 
                )
            {
                std::cerr << collection_dir << " does not contain a valid surf collection.\n";
                std::cerr << "The files " << surf::TEXT_FILENAME << " , " << surf::DICT_FILENAME 
                          << " have to be present" << std::endl;
                return false;
            }
        } else {
            if( ! surf::file_exists(collection_dir+"/"+surf::TEXT_FILENAME_BYTE) )
            {
                std::cerr << collection_dir << " does not contain a valid surf collection.\n";
                std::cerr << "The files " << surf::TEXT_FILENAME_BYTE << " have to be present" << std::endl;
                return false;
            }
        }
    }
    return true;
}


template<typename alphabet_tag=sdsl::int_alphabet_tag>
sdsl::cache_config
parse_collection(std::string collection_dir)
{
    /* check if all the directories exist */
    if( !surf::valid_collection<alphabet_tag>(collection_dir) ) {
        exit(EXIT_FAILURE);
    }
    std::string index_directory = collection_dir+"/index/";
    surf::create_directory(index_directory);

    std::string results_directory = collection_dir+"/results/";
    surf::create_directory(results_directory);

    /* populate cache config */
    sdsl::cache_config config(false,collection_dir+"/index/","SURF");

    /* create symlink to text in index directory */
    if ( std::is_same<alphabet_tag, sdsl::int_alphabet_tag>::value ){
        std::string symlink_name = cache_file_name(sdsl::conf::KEY_TEXT_INT,config);
        if( ! surf::symlink_exists(cache_file_name(sdsl::conf::KEY_TEXT_INT,config)) ) {
            std::string collection_file = collection_dir+"/"+surf::TEXT_FILENAME;
            char* col_file_absolute = realpath(collection_file.c_str(), NULL);
            if( symlink(col_file_absolute,symlink_name.c_str()) != 0) {
                perror("cannot create symlink to collection file in index directory");
                exit(EXIT_FAILURE);
            }
            free(col_file_absolute);
        }
    } else {
         std::string symlink_name = cache_file_name(sdsl::conf::KEY_TEXT,config);
        if( ! surf::symlink_exists(cache_file_name(sdsl::conf::KEY_TEXT,config)) ) {
            std::string collection_file = collection_dir+"/"+surf::TEXT_FILENAME_BYTE;
            char* col_file_absolute = realpath(collection_file.c_str(), NULL);
            if( symlink(col_file_absolute,symlink_name.c_str()) != 0) {
                perror("cannot create symlink to collection file in index directory");
                exit(EXIT_FAILURE);
            }
            free(col_file_absolute);
        }   
    }

    /* register files that are present */
    for(const auto& key : surf::storage_keys) {
        register_cache_file(key,config);
    }

    return config;
}


} // end of surf namespace
#endif
