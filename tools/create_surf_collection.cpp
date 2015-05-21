// extracts from an indri index a monoton sequence of integers in sdsl format
// which represent the parsed text collection.
#include <iostream>

#include "surf/config.hpp"
#include "surf/util.hpp"
#include "sdsl/int_vector_buffer.hpp"
#include <dirent.h>

typedef std::vector <std::string> DirListing_t;

void
GetDirListing( DirListing_t& result, const std::string& dirpath )
{
    DIR* dir = opendir( dirpath.c_str() );
    if (dir)
    {
        struct dirent* entry;
        while ((entry = readdir( dir )))
        {
            struct stat entryinfo;
            std::string entryname = entry->d_name;
            std::string entrypath = dirpath + "/" + entryname;
            if (!stat( entrypath.c_str(), &entryinfo ))
            {
                if (S_ISDIR( entryinfo.st_mode ))
                {
                    if      (entryname == "..");
                    else if (entryname == "." ) result.push_back( dirpath + "/" );
                    else                        GetDirListing( result, entrypath );
                }
                else
                {
                    result.push_back( entrypath );
                }
            }
        }
        closedir( dir );
    }
}

// TODO: redundant, fix
std::string
trim(const std::string& str, const std::string& whitespace = " \t")
{
    const auto strBegin = str.find_first_not_of(whitespace);
    if (strBegin == std::string::npos)
        return ""; // no content

    const auto strEnd = str.find_last_not_of(whitespace);
    const auto strRange = strEnd - strBegin + 1;

    return str.substr(strBegin, strRange);
}

int main( int argc, char** argv ) {
    if(argc != 3) {
        std::cout << "USAGE: " << argv[0] 
//                  << " <string with # separated docs> <surf collection folder>" << std::endl;
                    << " <input folder> <output folder>" << std::endl;
        return EXIT_FAILURE;
    }
    std::string input = argv[1];
    std::string output = argv[2];

    // setup collection directory
    if(surf::directory_exists(output)) {
        std::cerr << "ERROR: collection directory already exists." << std::endl;
        return EXIT_FAILURE;
    }
    surf::create_directory(output);

//    if( test_str.back() != '#' ) {
//        std::cerr << "ERROR: test string must end with doc seperator '#'" << std::endl;
//        return EXIT_FAILURE;
//    }

    DirListing_t dirtree;

    GetDirListing(dirtree, input);

    std::string global_str;

    for (unsigned n = 0; n < dirtree.size(); n++) {

        std::ifstream inFile;
        inFile.open(dirtree[n]);

        std::stringstream strStream;
        strStream << inFile.rdbuf();
        std::string current_str = strStream.str();

        // erase two or more whitespaces
        current_str.erase(std::unique(current_str.begin(), current_str.end(), [](char a, char b) { return a == ' ' && b == ' '; } ), current_str.end() );

        // reverse document
        std::reverse(current_str.begin(), current_str.end());

        // concat to global string
        global_str.append(" ");
        global_str.append(trim(current_str));
        global_str.append(" #");
    }

    // write collection string
    std::map<std::string::value_type,sdsl::int_vector<>::value_type> existing_syms;
    std::map<sdsl::int_vector<>::value_type,std::string::value_type> sym_mapping;
    sdsl::int_vector<> text_col(global_str.size()+1);
    size_t j=0;
    size_t num_docs = 0;
    for(const auto& sym : global_str) {

        std::cout << sym << std::endl;

        if(sym == '#') {
            text_col[j++] = 1;
            num_docs++;
        } else {
            auto itr = existing_syms.find(sym);
            if(itr != existing_syms.end()) {
                text_col[j++] = itr->second;
            } else {
                sdsl::int_vector<>::value_type new_sym = existing_syms.size()+2;
                existing_syms[sym] = new_sym;
                sym_mapping[new_sym] = sym;
                text_col[j++] = new_sym;
            }
        }
    }
    text_col[j] = 0;
    std::ofstream ofs(output+"/"+surf::TEXT_FILENAME);
    if(ofs.is_open()) {
        text_col.serialize(ofs);
    } else {
        std::cerr << "ERROR: could not write collection file." << std::endl;
        return EXIT_FAILURE;
    }

    // write the dict
    std::ofstream dict_ofs(output+"/"+surf::DICT_FILENAME);
    if(dict_ofs.is_open()) {
        for(const auto& mapping : sym_mapping) {
            dict_ofs << mapping.second << " " << mapping.first << std::endl;
        }
    } else {
        std::cerr << "ERROR: could not write dictionary file." << std::endl;
        return EXIT_FAILURE;
    }

    // write docnames file
    std::ofstream docnames_ofs(output+"/"+surf::DOCNAMES_FILENAME);
    if(docnames_ofs.is_open()) {
        for(size_t i=1;i<=num_docs;i++) {
            docnames_ofs << "DOCUMENT " << i << std::endl;
        }
    } else {
        std::cerr << "ERROR: could not write docnames file." << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "Created surf collection for string '" << global_str << "'" << std::endl;
    std::cout << "Found " << num_docs << " documents." << std::endl;
    std::cout << "Document delimiter = " << 1 << std::endl;
    std::cout << surf::TEXT_FILENAME << ": ";
    for(const auto& sym : text_col) {
        std::cout << sym << " ";
    }
    std::cout << std::endl;
    std::cout << "Mapping: ";
    for(const auto& mapping : sym_mapping) {
        std::cout << mapping.second << " -> " << mapping.first << "; ";
    }
    std::cout << std::endl;
    std::cout << "Document Names: ";
    for(size_t i=1;i<=num_docs;i++) {
        std::cout << "'DOCUMENT " << i << "'; ";
    }
    std::cout << std::endl;
}


