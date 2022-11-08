/*
 * Copyright 2022 Andy Kirkham
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


#include <DateTime.h>

#if 0
#include <Tle.h>
#include <SGP4.h>
#include <Observer.h>
#include <CoordGeodetic.h>
#include <CoordTopocentric.h>
#endif

#include <list>
#include <string>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>


int RunTest(void)
{
    std::string dut("2022-11-08T06:14:56.037120");
    DateTime dt(dut);
    std::cout << dt.ToString();
    return 1;
}

int main()
{
    return RunTest();
}
