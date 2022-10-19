#pragma once
#include <plog/Record.h>
#include <plog/Util.h>
#include <iomanip>

//https://stackoverflow.com/questions/2616906/how-do-i-output-coloured-text-to-a-linux-terminal
namespace plog
{
    template<bool useUtcTime>
    class TxtFormatterImpl
    {
    public:
        static util::nstring header()
        {
            return util::nstring();
        }

        static util::nstring format(const Record& record)
        {
            tm t;
            useUtcTime ? util::gmtime_s(&t, &record.getTime().time) : util::localtime_s(&t, &record.getTime().time);

            util::nostringstream ss;
            ss << "\033[;35m[" << t.tm_year + 1900 << "-" << std::setfill(PLOG_NSTR('0')) << std::setw(2) << t.tm_mon + 1 << PLOG_NSTR("-") << std::setfill(PLOG_NSTR('0')) << std::setw(2) << t.tm_mday << PLOG_NSTR(" ");
            ss << std::setfill(PLOG_NSTR('0')) << std::setw(2) << t.tm_hour << PLOG_NSTR(":") << std::setfill(PLOG_NSTR('0')) << std::setw(2) << t.tm_min << PLOG_NSTR(":") << std::setfill(PLOG_NSTR('0')) << std::setw(2) << t.tm_sec << PLOG_NSTR(".") << std::setfill(PLOG_NSTR('0')) << std::setw(3) << static_cast<int> (record.getTime().millitm) << "]\033[0m" << PLOG_NSTR(" ");
            ss << std::setfill(PLOG_NSTR(' ')) << "\033[1;33m" << std::setw(5) << std::left << severityToString(record.getSeverity()) << "\033[0m" << PLOG_NSTR(" ");
//            ss << PLOG_NSTR("[") << record.getTid() << PLOG_NSTR("] ");
            ss << "\033[;36m" << PLOG_NSTR("[") << record.getFunc() << PLOG_NSTR("@") << record.getLine() << PLOG_NSTR("] ") << "\033[0m";
            ss << "\033[1;37m" << record.getMessage() << "\033[0m" << PLOG_NSTR("\n");

            return ss.str();
        }
    };

    class TxtFormatter : public TxtFormatterImpl<false> {};
    class TxtFormatterUtcTime : public TxtFormatterImpl<true> {};
}
