#include "generateheader.h"

std::vector<std::string> getVersion() {
    std::ifstream file("../VERSION");
    if (!file) {
        throw std::runtime_error("Unable to open VERSION file");
    }

    std::string line;
    std::getline(file, line);
    file.close();

    // Removing any potential trailing newline or whitespace
    line.erase(line.find_last_not_of(" \n\r\t") + 1);

    std::stringstream ss(line);
    std::string mj, mn, pa;
    std::getline(ss, mj, '.');
    std::getline(ss, mn, '.');
    std::getline(ss, pa);

    std::vector<std::string> out = {
        "#define LIBJUPITERMAG_VERSION_MAJOR " + mj + "\n",
        "#define LIBJUPITERMAG_VERSION_MINOR " + mn + "\n",
        "#define LIBJUPITERMAG_VERSION_PATCH " + pa + "\n",
    };

    return out;
}

std::string stringVectorToString(std::vector<std::string> strVec) {

    std::ostringstream out;
    for (auto &str : strVec) {
        if (str.find("typedef void (*modelFieldPtr)(double,double,double,double*,double*,double*);") == std::string::npos) {
            out << str;
        }
        
    }
    return out.str();
}



std::vector<std::string> removeDirectives(const std::vector<std::string>& lines) {
    std::vector<std::string> filteredLines;

    for (const auto& line : lines) {
        // Check if the line starts with '#' and skip it if it does
        if (line.find("#") != 0) {
            filteredLines.push_back(line);
        }
    }

    return filteredLines;
}

std::vector<std::string> readASCII(std::filesystem::path fileName) {
    std::ifstream file(fileName);
    std::string line;
    std::vector<std::string> out;
    while (std::getline(file,line)) {
        out.push_back(line+"\n");
    } 
    file.close();
    return out;
}


StrVecPair splitHeaderDefs(const std::vector<std::string>& lines) {
    std::vector<std::string> cCode;
    std::vector<std::string> cppCode;

    bool isC = false;
    for (const auto& line : lines) {
        std::string trimmedLine = trimString(line); // Add logic to trim whitespace if necessary
        if (isC && trimmedLine == "}") {
            isC = false;
        }
        if (isC) {
            cCode.push_back(line);
        } else {
            cppCode.push_back(line);
        }
        if (trimmedLine.find("extern \"C\"") != std::string::npos) {
            isC = true;
        }
    }

    return {cCode, cppCode};
}


std::string headerIncludes() {
    std::string a = R"(
#ifndef __LIBJUPITERMAG_H__
#define __LIBJUPITERMAG_H__
#define _USE_MATH_DEFINES

#ifdef __cplusplus
	#include <algorithm>
	#include <map>
	#include <math.h>
	#include <stdio.h>
	#include <stdlib.h>
	#include <string.h>
	#include <string>
	#include <vector>
#else
	#include <math.h>
	#include <stdio.h>
	#include <stdlib.h>
	#include <string.h>
	#include <stdbool.h>
#endif

typedef void (*modelFieldPtr)(double,double,double,double*,double*,double*);


)";

    std::vector<std::string> version = getVersion();
    std::stringstream b;
    for (const auto& line : version) {
        b << line;
    }

    std::string c = R"(

#ifdef __cplusplus
extern "C" {
#endif
)";

    return a + b.str() + c;
}



StrVecPair extractCandCC(std::filesystem::path fileName) {

    std::vector<std::string> lines = readASCII(fileName);
    StrVecPair code = splitHeaderDefs(lines);
    std::vector<std::string> cCode = code.first;
    std::vector<std::string> ccCode = code.second;
    std::vector<std::string> cOut = removeDirectives(cCode);
    std::vector<std::string> ccOut = removeDirectives(ccCode);

    return {cOut,ccOut};
}


std::string generateLibHeader() {

    std::ostringstream out;
    std::ostringstream cCode;
    std::ostringstream ccCode;
    StrVecPair codePair;
    std::vector<std::string> cVec, ccVec;

    /* add the includes and stuff at the top */
    out << headerIncludes();

    /* list of headers to read */
    std::vector<std::filesystem::path> headerFiles = {
        "../src/coordconv.h",
        "../src/footprint.h",
        "../lib/libspline/include/spline.h",
        "../src/interptraceclosestpos.h",
        "../lib/libcon2020/include/con2020.h",
        "../lib/libinternalfield/include/internalfield.h",
        "../src/model.h",
        "../src/trace.h",
        "../src/libjupitermag.h"
    };

    /* read each into the two code arrays */
    for (auto &headerFile : headerFiles) {
        codePair = extractCandCC(headerFile);
        cCode << stringVectorToString(codePair.first);
        ccCode << stringVectorToString(codePair.second);        
    }


    /* add the C-compatible bits to the extern C section */
    out << cCode.str();

    /* close the extern C bit and add the C++ only code */
    out << "#ifdef __cplusplus\n";
    out << "}\n\n";
    out << ccCode.str();
    out << "#endif\n";
    out << "#endif\n";

    return out.str();

}

void saveLibHeader(std::filesystem::path srcPath) {

    std::string headerCode = generateLibHeader();

    std::filesystem::path headerPath = srcPath.parent_path();
    headerPath /= "include/jupitermag.h";
    std::cout << "Saving library header: " << std::endl;
    std::cout << "\t" << headerPath << std::endl;
    std::ofstream file(headerPath);
    file << headerCode;
    file.close();
}

int main(int argc, char *argv[]) {


    std::filesystem::path srcPath = std::filesystem::current_path();

    /* generate the header for the library */
    saveLibHeader(srcPath);

    return 0;

}