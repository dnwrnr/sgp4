#include "Util.h"

#include <algorithm>
#include <locale>

namespace Util
{
    void TrimLeft(std::string& s)
    {
        s.erase(s.begin(),
                std::find_if(s.begin(), s.end(),
                    std::not1(std::ptr_fun<int, int>(std::isspace))));
    }

    void TrimRight(std::string& s)
    {
        s.erase(std::find_if(s.rbegin(), s.rend(),
                std::not1(std::ptr_fun<int, int>(std::isspace))).base(),
                s.end());
    }

    void Trim(std::string& s)
    {
        TrimLeft(s);
        TrimRight(s);
    }

}
