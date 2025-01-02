#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <string>
#include <unordered_map>
#include <muParser.h>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <regex>

extern double eps_0;
extern double eps_r;
extern double EXACT_POINT_VALUE;
extern double EXACT_FLUX;
extern unsigned int NUM_PRELIMINARY_GLOBAL_REF;
extern unsigned int NUM_PRELIMINARY_REF;

extern bool MANUAL_LIFTING_ON;

class ConstantsParser {
public:
    explicit ConstantsParser(const std::string &filePath) {
        parseFile(filePath);
    }

    double getConstant(const std::string &name) const {
        auto it = constants.find(name);
        if (it != constants.end()) {
            return it->second;
        }
        throw std::runtime_error("Constant not found: " + name);
    }

    std::string getStringConstant(const std::string &name) const {
        auto it = stringConstants.find(name);
        if (it != stringConstants.end()) {
            return it->second;
        }
        throw std::runtime_error("String constant not found: " + name);
    }

private:
    std::unordered_map<std::string, double> constants;
    std::unordered_map<std::string, std::string> stringConstants;

    void parseFile(const std::string &filePath) {
        std::ifstream file(filePath);
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open file: " + filePath);
        }

        std::string line;
        std::regex constPattern(R"((\w+):\s*([^#]+))"); // YAML-style parsing
        std::smatch matches;

        mu::Parser parser;

        // First pass: parse constants and add them to the parser
        while (std::getline(file, line)) {
            if (std::regex_search(line, matches, constPattern)) {
                std::string name = matches[1];
                std::string valueExpr = matches[2];

                // Check if it's a numeric constant
                try {
                    parser.SetExpr(valueExpr);
                    double value = parser.Eval();
                    constants[name] = value;
                    parser.DefineVar(name, &constants[name]); // Add as variable for later
                } catch (mu::Parser::exception_type &e) {
                    // If parsing fails, treat it as a string constant
                    stringConstants[name] = valueExpr;
                }
            }
        }
    }

};

// Classe singleton per accedere globalmente alle costanti
class GlobalConstants {
public:
    static GlobalConstants& getInstance(const std::string &filePath = "") {
        static GlobalConstants instance(filePath);
        return instance;
    }

    // Get numeric constant
    double get(const std::string &name) {
        return parser.getConstant(name);
    }

    // Get string constant
    std::string getString(const std::string &name) {
        return parser.getStringConstant(name);
    }

private:
    GlobalConstants(const std::string &filePath) : parser(filePath.empty() ? "../constants.yaml" : filePath) {}

    ConstantsParser parser;

    // Evitare la copia e l'assegnazione
    GlobalConstants(const GlobalConstants&) = delete;
    GlobalConstants& operator=(const GlobalConstants&) = delete;
};

#endif // CONSTANTS_H
