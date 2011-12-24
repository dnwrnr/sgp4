#include "Util.h"

#include <algorithm>
#include <locale>
#include <functional>

namespace Util
{
    namespace
    {
        struct IsDigit: std::unary_function<char, bool>
        {
            bool operator()(char c) const
            {
                return std::isdigit(c, std::locale::classic()) == 0;
            }
        };
    }

    void TrimLeft(std::string& s)
    {
        s.erase(s.begin(),
                std::find_if(s.begin(), s.end(), std::not1(IsDigit())));
    }

    void TrimRight(std::string& s)
    {
        s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(IsDigit())).base(),
                s.end());
    }

    void Trim(std::string& s)
    {
        TrimLeft(s);
        TrimRight(s);
    }

}
