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

#pragma once

#include "Tle.h"

#include <string>
#include <vector>

namespace libsgp4
{

/**
 * @brief Load Tle objects from a CelesTrak CSV file.
 *
 * The file must have a header row starting with "OBJECT_NAME" which is
 * skipped. Each subsequent line is parsed by Tle::FromCsv().
 * Empty lines are skipped.
 *
 * @param filename Path to the CSV file
 * @returns Vector of Tle objects
 * @throws TleException if a data line is malformed
 */
std::vector<Tle> LoadCsvTleFile(const std::string& filename);

} // namespace libsgp4
