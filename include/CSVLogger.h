/* ------------------------------------------------------------------------
*
 * SPDX-License-Identifier: GPL-3-or-later
 * Copyright (C) 2023 - 2024 by Matteo Malvestiti
 *
 * This file is part of ion_propulsion.
 *
 * ion_propulsion is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * ion_propulsion is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with ion_propulsion; see the file COPYING.  If not, see
 * <https://www.gnu.org/licenses/>.
 *
 * ------------------------------------------------------------------------
 *
 * Author: Matteo Malvestiti, Politecnico di Milano, 2024
 *
 */

#ifndef CSVLOGGER_H
#define CSVLOGGER_H

#include <Constants.h>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <sstream>
#include <mutex>

#include <sstream>
#include <string>

inline std::string to_string_with_precision(const double value, const int precision = 6, const bool scientific = false) {
    std::ostringstream oss;
    if (scientific) {
        oss << std::scientific;
    } else {
        oss << std::fixed;
    }
    oss << std::setprecision(precision) << value;
    return oss.str();
}

class CSVLogger {
public:
    // Get the singleton instance
    static CSVLogger& getInstance(const std::string& filename = OUTPUT_PATH+"/"+"convergence_results.csv") {
        static CSVLogger instance(filename);
        return instance;
    }

    // Add or update a column value for the current row
    void addColumn(const std::string& column, const std::string& value) {
        std::lock_guard<std::mutex> lock(mutex_);
        if (columnOrder_.find(column) == columnOrder_.end()) {
            columnOrder_[column] = headers_.size();
            headers_.push_back(column);
        }
        currentRow_[column] = value;
    }

    // Write the current row to the file and reset
    void flushRow() {
        std::lock_guard<std::mutex> lock(mutex_);

        // Write header (if not already written)
        if (!headerWritten_) {
            writeHeader();
        }

        // Write the row data
        writeRow();
        currentRow_.clear(); // Reset the row buffer
    }

    // Disable copying and assignment
    CSVLogger(const CSVLogger&) = delete;
    CSVLogger& operator=(const CSVLogger&) = delete;

private:
    std::ofstream file_;
    std::unordered_map<std::string, std::string> currentRow_;
    std::vector<std::string> headers_;
    std::unordered_map<std::string, size_t> columnOrder_;
    bool headerWritten_ = false;
    std::mutex mutex_;

    CSVLogger(const std::string& filename) {
        file_.open(filename, std::ios::out | std::ios::trunc); // Clear file on creation
        if (!file_) {
            throw std::ios_base::failure("Failed to open file: " + filename);
        }
    }

    ~CSVLogger() {
        if (file_.is_open()) {
            file_.close();
        }
    }

    void writeHeader() {
        writeRowData(headers_);
        headerWritten_ = true;
    }

    void writeRow() {
        std::vector<std::string> rowData(headers_.size(), "");
        for (const auto& [column, value] : currentRow_) {
            size_t index = columnOrder_[column];
            rowData[index] = value;
        }
        writeRowData(rowData);
    }

    void writeRowData(const std::vector<std::string>& rowData) {
        std::ostringstream oss;
        for (size_t i = 0; i < rowData.size(); ++i) {
            if (i > 0) {
                oss << ",";
            }
            oss << rowData[i];
        }
        file_ << oss.str() << std::endl;
    }
};



#endif