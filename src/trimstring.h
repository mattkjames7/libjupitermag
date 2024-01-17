#ifndef __TRIMSTRING_H__
#define __TRIMSTRING_H__
#include <string>
#include <algorithm>
#include <cctype>


// Function to trim whitespace from the start of a string
std::string leftTrim(const std::string& str);

// Function to trim whitespace from the end of a string
std::string rightTrim(const std::string& str);

// Function to trim whitespace from both ends of a string
std::string trimString(const std::string& str);

#endif
