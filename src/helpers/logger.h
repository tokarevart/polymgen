// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <string>
#include <iostream>
#include <vector>
#include <sstream>
#include <memory>


class Logger {
    enum class IomanipType {
        Width,
        Precision
    };

    template <IomanipType IomT>
    struct Iomanip;

public:
    static Logger::Iomanip<IomanipType::Width>     setw(std::streamsize new_width);
    static Logger::Iomanip<IomanipType::Precision> setprecision(std::streamsize new_precision);

    void open(std::ostream& stream);
    void flush();
    void clear();

    std::streamsize width() const;
    std::streamsize width(std::streamsize new_width);
    std::streamsize precision() const;
    std::streamsize precision(std::streamsize new_precision);
    std::ios::fmtflags setf(std::ios::fmtflags flags);
    std::ios::fmtflags setf(std::ios::fmtflags flags, std::ios::fmtflags mask);

    Logger& operator<<(short value);
    Logger& operator<<(unsigned short value);
    Logger& operator<<(int value);
    Logger& operator<<(unsigned int value);
    Logger& operator<<(long value);
    Logger& operator<<(unsigned long value);
    Logger& operator<<(long long value);
    Logger& operator<<(unsigned long long value);
    Logger& operator<<(float value);
    Logger& operator<<(double value);
    Logger& operator<<(long double value);
    Logger& operator<<(bool value);
    Logger& operator<<(const void* value);
    Logger& operator<<(const char* value);
    Logger& operator<<(const std::string& str);
    Logger& operator<<(std::ios_base& (*func)(std::ios_base&));
    Logger& operator<<(std::ostream& (*func)(std::ostream&));
    template <IomanipType IomT>
    Logger& operator<<(Logger::Iomanip<IomT> ioManip);

    Logger();
    Logger(std::ostream& stream);
    ~Logger();


private:
    template <IomanipType IomT>
    struct Iomanip {
        std::streamsize param;
    };

    std::ostream* m_stream = nullptr;

    std::vector<std::string> m_descriptions;
    std::vector<std::string> m_values;
    std::vector<std::string> m_endings;

    std::unique_ptr<std::stringstream> m_bufiss;
    int m_turnFlag = 0;

    void istreamTypeOperatorHelper();
};
