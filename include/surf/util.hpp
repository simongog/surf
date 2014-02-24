#ifndef SURF_UTIL_HPP
#define SURF_UTIL_HPP

#include <string>
#include <unistd.h>
#include <stdlib.h>
#include <iostream>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

namespace surf{

namespace util{


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
    isfstream in(file_name);
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

bool
valid_collection(std::string collection_dir)
{
    if (! surf::util::directory_exists(collection_dir)) {
        std::cerr << collection_dir << " is not a valid directory.\n";
        return false;
    } else {
        /* make sure the necessary files are present */
        if( ! surf::util::file_exists(collection_dir+"/"+surf::TEXT_FILENAME) ||
            ! surf::util::file_exists(collection_dir+"/"+surf::DICT_FILENAME) ||
            ! surf::util::file_exists(collection_dir+"/"+surf::DOCNAMES_FILENAME) )
        {
            std::cerr << collection_dir << " does not contain a valid surf collection.\n";
            std::cerr << "The files " << surf::TEXT_FILENAME << " , " << surf::DICT_FILENAME 
                      << " , " << surf::DOCNAMES_FILENAME << " have to be present" << std::endl;
            return false;
        }
    }
    return true;
}

} // end of util namespace

} // end of surf namespace
#endif
