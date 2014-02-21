// extracts from an indri index a monoton sequence of integers in sdsl format
// which represent the parsed text collection.
#include <iostream>

#include "indri/Repository.hpp"
#include "indri/CompressedCollection.hpp"
#include "sdsl/int_vector_buffer.hpp"

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


int main( int argc, char** argv ) {
    if(argc != 3) {
        std::cout << "USAGE: " << argv[0] << " <indri repository> <surf collection folder>" << std::endl;
        return EXIT_FAILURE;
    }

    // parse cmd line
    std::string repository_name = argv[1];
    std::string surf_collection_folder = argv[2];
    create_directory(surf_collection_folder);
    std::string dict_file = surf_collection_folder + "/dict.txt";
    std::string doc_names_file = surf_collection_folder + "/doc_names.txt";
    std::string text_int_file = surf_collection_folder + "/text_int.sdsl";

    // load stuff
    indri::collection::Repository repo;
    repo.openRead( repository_name );

    // extract
    std::cout << "extracting sdsl integer file from indri index into file " << text_int_file << std::endl;
    std::vector<std::string> document_names;
    indri::collection::Repository::index_state state = repo.indexes();
    const auto& index = (*state)[0];
    uint64_t uniq_terms = index->uniqueTermCount();
    uniq_terms += 2; // we will shift all ids from idri by 2 so \0 and \1 is free
    uint8_t out_int_width = sdsl::bits::hi(uniq_terms)+1;
    sdsl::int_vector_buffer<> sdsl_col_file(text_int_file,std::ios::out,1024*1024,out_int_width,false);
    size_t written_term_ids = 0;
    indri::collection::CompressedCollection* collection = repo.collection();
    int64_t document_id = index->documentBase();
    indri::index::TermListFileIterator* iter = index->termListFileIterator();
    iter->startIteration();
    while( !iter->finished() ) {
        indri::index::TermList* list = iter->currentEntry();

        // find document name
        std::string doc_name = collection->retrieveMetadatum( document_id , "docno" );
        document_names.push_back(doc_name);

        if(document_id % 10000 == 0) {
            std::cout << ".";
            std::cout.flush();
        }

        // iterate over termlist
        for(const auto& term_id : list->terms()) {
            // we will shift all ids from idri by 2 so \0 and \1 is free
            if(term_id != 0) {
                sdsl_col_file[written_term_ids++] = term_id+2; 
            }
        }
        sdsl_col_file[written_term_ids++] = 1; // end of doc sep

        document_id++;
        iter->nextEntry();
    }
    std::cout << std::endl;
    sdsl_col_file[written_term_ids++] = 0; // end of collection sep

    // write document names
    {
        std::cout << "writing document names to " << doc_names_file << std::endl;
        std::ofstream of_doc_names(doc_names_file);
        for(const auto& doc_name : document_names) {
            of_doc_names << doc_name << std::endl;
        }
    }
    // write dictionary
    {
        std::cout << "writing dictionary to " << dict_file << std::endl;
        const auto& index = (*state)[0];
        std::ofstream of_dict(dict_file);
        for(size_t i=1;i<index->uniqueTermCount();i++) {
            auto term_str = index->term(i);
            of_dict << term_str << " " << i+2 << std::endl;
        }
    }
}


