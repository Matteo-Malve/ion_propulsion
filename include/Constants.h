#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <string>
#include <unordered_map>
#include <muParser.h>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <regex>

extern std::string OUTPUT_PATH;

extern double eps_0;
extern double eps_r;

extern double Ve;
extern double Vc;
extern std::string RHS_EXPRESSION;
//extern std::string UEX_EXPRESSION;


//extern unsigned int NUM_PRELIMINARY_GLOBAL_REF;
extern std::string PATH_TO_MESH;
extern unsigned int NUM_CONCENTRIC_REF;
extern unsigned int LOAD_FROM_SETUP;

extern bool MANUAL_LIFTING_ON;
extern unsigned int REFINEMENT_CRITERION;
extern unsigned int DUAL_FUNCTIONAL;

extern double EVALUATION_POINT_X;
extern double EVALUATION_POINT_Y;

extern double EXACT_POINT_VALUE;
extern double EXACT_FLUX;

//extern unsigned int NUM_PRELIMINARY_REF;

class ConstantsParser {
public:
    explicit ConstantsParser(const std::string &filePath) {
        parseFile(filePath);
    }

    double getConstant(const std::string &name, bool default_zero = false) const {
        auto it = constants.find(name);
        if (it != constants.end()) {
            return it->second;
        }
        if (default_zero)
            return 0.0;
        else
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

class GlobalConstants {
public:
    // Get the singleton instance, ensuring proper initialization
    static GlobalConstants& getInstance() {
        if (!instance) {
            throw std::runtime_error("GlobalConstants must be initialized before use.");
        }
        return *instance;
    }

    // Explicitly initialize the singleton with a configuration file
    static void initialize(const std::string &filePath) {
        delete instance; // Delete existing instance if any
        instance = new GlobalConstants(filePath);

        // Extract parent path and assign to OUTPUT_PATH
        std::filesystem::path configPath(filePath);
        OUTPUT_PATH = configPath.parent_path().string()+"/results";
    }

    // Get numeric constant
    double get(const std::string &name, bool default_zero = false)  {
        return parser.getConstant(name,default_zero);
    }

    // Get string constant
    std::string getString(const std::string &name) {
        return parser.getStringConstant(name);
    }

private:
    explicit GlobalConstants(const std::string &filePath) : parser(filePath.empty() ? "../config.yaml" : filePath) {}
    static GlobalConstants* instance;

    ConstantsParser parser;

    // Prevent copying and assignment
    GlobalConstants(const GlobalConstants&) = delete;
    GlobalConstants& operator=(const GlobalConstants&) = delete;
};

void useGlobalConstants();
void printParsedConstants();



#endif // CONSTANTS_H
