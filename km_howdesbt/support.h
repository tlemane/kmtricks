#ifndef support_H
#define support_H

#include <string>
#include <vector>

//----------
//
// prototypes for functions in this module--
//
//----------

std::vector<std::string> parse_comma_list (const std::string& s);
std::vector<std::string> tokenize         (const std::string& s);
std::vector<std::string> quoted_tokenize  (const std::string& s);
void                     expand_filenames (const std::vector<std::string> filenames,
                                           int fileCount,
                                           std::vector<std::string>& newFilenames);

#endif // support_H
