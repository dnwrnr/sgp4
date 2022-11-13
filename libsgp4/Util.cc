/*
 * Copyright 2013 Daniel Warner <contact@danrw.com>
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */


#include "Util.h"

#include <algorithm>
#include <locale>
#include <functional>

namespace libsgp4
{
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
} // namespace Util
} // namespace libsgp4
