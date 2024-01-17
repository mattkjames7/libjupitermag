#include "trimstring.h"




// Function to trim whitespace from the start of a string
std::string leftTrim(const std::string& str) {
    std::string trimmed = str;
    trimmed.erase(trimmed.begin(), std::find_if(trimmed.begin(), trimmed.end(), [](unsigned char ch) {
        return !std::isspace(ch);
    }));
    return trimmed;
}

// Function to trim whitespace from the end of a string
std::string rightTrim(const std::string& str) {
    std::string trimmed = str;
    trimmed.erase(std::find_if(trimmed.rbegin(), trimmed.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), trimmed.end());
    return trimmed;
}

// Function to trim whitespace from both ends of a string
std::string trimString(const std::string& str) {
    return leftTrim(rightTrim(str));
}
