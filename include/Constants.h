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

extern unsigned int REFINEMENT_STRATEGY;
extern double TOP_FRACTION;
extern double BOTTOM_FRACTION;
extern unsigned int OPTIMIZE_ORDER;


//extern unsigned int NUM_PRELIMINARY_REF;

class ConstantsParser {
public:
    explicit ConstantsParser(const std::string &filePath) {
        parseFile(filePath);
    }

    double getConstant(const std::string &name, double default_value = std::numeric_limits<double>::lowest()) const {
        auto it = constants.find(name);
        if (it != constants.end()) {
            return it->second;
        }
        if (default_value > std::numeric_limits<double>::lowest() + 1e-10)  // Give some tolerance here
            return default_value;
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

        // Extract CONFIG_NAME
        std::string filename = configPath.filename().string(); // Get the file name: "config-4.yaml" or "config.yaml"
        size_t first_dash = filename.find('-');
        size_t second_dash = filename.find('-', first_dash + 1);

        std::string configName;
        if (first_dash != std::string::npos) {
            if (second_dash != std::string::npos) {
                configName = filename.substr(0, second_dash); // Extract up to the second dash
            } else {
                configName = filename.substr(0, filename.find('.')); // Extract up to the file extension
            }
        } else {
            configName = filename.substr(0, filename.find('.')); // Handle cases without dashes
        }

        // Create OUTPUT_PATH with configName
        OUTPUT_PATH = configPath.parent_path().string()+"/results/"+configName;

        // Check if the directory exists, and create it if not
        std::filesystem::path outputPath(OUTPUT_PATH);
        if (!std::filesystem::exists(outputPath)) {
            if (std::filesystem::create_directories(outputPath)) {
                std::cout << "Directory created: " << OUTPUT_PATH << std::endl;
            } else {
                std::cerr << "Failed to create directory: " << OUTPUT_PATH << std::endl;
            }
        } else {
            std::cout << "Directory already exists: " << OUTPUT_PATH << std::endl;
        }
    }

    // Get numeric constant
    double get(const std::string &name, double default_value = std::numeric_limits<double>::lowest())  {
        return parser.getConstant(name,default_value);
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
