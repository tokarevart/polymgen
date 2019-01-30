#pragma once
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <memory>


namespace tva
{
    class Logger
    {
        struct Iomanip;

    public:
        enum IomanipType
        {
            WIDTH,
            PRECISION
        };


        static Logger::Iomanip setw(std::streamsize new_width);
        static Logger::Iomanip setprecision(std::streamsize new_precision);

        bool isOpen() const;
        void open(const std::string& filename, std::ios::openmode mode = std::ios::out);
        void close();
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
        Logger& operator<<(Logger::Iomanip ioManip);

        Logger();
        Logger(const std::string& filename, std::ios::openmode mode = std::ios::out);
        ~Logger();


    private:
        std::ofstream m_file;

        std::vector<std::string> m_descriptions;
        std::vector<std::string> m_values;
        std::vector<std::string> m_endings;

        std::unique_ptr<std::stringstream> m_bufiss;
        int m_turnFlag = 0;

        void flush();

        struct Iomanip
        {
            IomanipType manipType;
            std::streamsize param;
        };
    };
}