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

#include <stdexcept>
#include <string>

namespace libsgp4
{

/**
 * @brief The exception that the Tle class throws on an error.
 *
 * The exception that the Tle decoder will throw on an error.
 */
class TleException : public std::runtime_error
{
public:
    /**
     * Constructor
     * @param message Exception message
     */
    explicit TleException(const char* message)
        : runtime_error(message)
    {
    }
};

} // namespace libsgp4
