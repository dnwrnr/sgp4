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

#include "CsvTleLoader.h"

#include <fstream>
#include <iostream>

namespace libsgp4
{

std::vector<Tle> LoadCsvTleFile(const std::string& filename)
{
    std::vector<Tle> tles;
    std::ifstream file(filename);

    if (!file.is_open())
    {
        return tles;
    }

    std::string line;
    bool headerSkipped = false;

    while (std::getline(file, line))
    {
        if (!headerSkipped)
        {
            headerSkipped = true;
            continue;
        }

        if (line.empty())
        {
            continue;
        }

        tles.push_back(Tle::FromCsv(line));
    }

    return tles;
}

} // namespace libsgp4
