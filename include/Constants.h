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

private:
    std::unordered_map<std::string, double> constants;

    void parseFile(const std::string &filePath) {
        std::ifstream file(filePath);
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open file: " + filePath);
        }

        std::string line;
        std::regex constPattern(R"((\w+):\s*([^#]+))"); // YAML-style parsing
        std::smatch matches;

        while (std::getline(file, line)) {
            if (std::regex_search(line, matches, constPattern)) {
                std::string name = matches[1];
                std::string valueExpr = matches[2];
                mu::Parser parser;
                try {
                    parser.SetExpr(valueExpr);
                    constants[name] = parser.Eval();
                } catch (mu::Parser::exception_type &e) {
                    std::cerr << "Error parsing expression for " << name << ": " << e.GetMsg() << std::endl;
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

    double get(const std::string &name) {
        return parser.getConstant(name);
    }

private:
    GlobalConstants(const std::string &filePath) : parser(filePath.empty() ? "../constants.yaml" : filePath) {}

    ConstantsParser parser;

    // Evitare la copia e l'assegnazione
    GlobalConstants(const GlobalConstants&) = delete;
    GlobalConstants& operator=(const GlobalConstants&) = delete;
};

#endif // CONSTANTS_H
