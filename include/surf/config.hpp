#ifndef SURF_CONFIG_HPP
#define SURF_CONFIG_HPP

#include <string>
#include <vector>

namespace surf{

const std::string TEXT_FILENAME = "text_int.sdsl";
const std::string DICT_FILENAME = "dict.txt";
const std::string DOCNAMES_FILENAME = "doc_names.txt";


const std::string KEY_DOCWEIGHT = "docweights";
const std::string KEY_DARRAY = "darray";
const std::string KEY_DOCPERM = "docperm";
const std::string KEY_SADADF = "sadadf";
const std::string KEY_SADADFSEL = "sadadfsel";
const std::string KEY_WTD = "wtd";
const std::string KEY_C = "C";
const std::string KEY_WTC = "wtc";
const std::string KEY_TMPCST = "tempcst";
const std::string KEY_TMPDUP = "tmpdup";
const std::string KEY_WTDUP  = "wtdup";
const std::string KEY_DOCCNT  = "doccnt";

std::vector<std::string> storage_keys = {KEY_DOCCNT,
										 KEY_DARRAY,
										 KEY_DOCPERM,
										 KEY_SADADF,
										 KEY_SADADFSEL,
										 KEY_WTD,
										 KEY_C,
										 KEY_WTC,
										 KEY_TMPCST,
										 KEY_TMPDUP,
										 KEY_WTDUP,
										 KEY_DOCCNT};

} // end namespace
#endif
