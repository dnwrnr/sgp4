# SGP4 library
## CMake support
### External
```cmake
# ./CMakeLists.txt

find_package(SGP4 1.1 CONFIG REQUIRED)
# ...
add_library(foo ...)
# ...
target_link_libraries(foo PRIVATE libsgp4::sgp4)
# target_link_libraries(foo PRIVATE libsgp4::sgp4-static) # or static
```

### Embedded
```cmake
# ./CMakeLists.txt

# this affects if the libsgp4::sgp4 target is shared or static
set(SGP4_BUILD_SHARED_LIBS ON)
# this affects if the libsgp4::sgp4-static target is available
set(SGP4_BUILD_STATIC_LIBS ON)

add_subdirectory(thirdparty/sgp4)
# ...
add_library(foo ...)
# ...
target_link_libraries(foo PRIVATE libsgp4::sgp4)
# target_link_libraries(foo PRIVATE libsgp4::sgp4-static) # or static
```

### Embedded (FetchContent)
```cmake
# ./CMakeLists.txt

include(FetchContent)

# this affects if the libsgp4::sgp4 target is shared or static
set(SGP4_BUILD_SHARED_LIBS ON)
# this affects if the libsgp4::sgp4-static target is available
set(SGP4_BUILD_STATIC_LIBS ON)

FetchContent_Declare(sgp4
  GIT_REPOSITORY https://github.com/dnwrnr/sgp4
  # TODO: update this with an actual commit
  GIT_TAG        bdff20d3eb6f010af95574990a662f8089db93f8 # https://github.com/dnwrnr/sgp4/tree/bdff20d3eb6f010af95574990a662f8089db93f8
)
FetchContent_MakeAvailable(sgp4)
# ...
add_library(foo ...)
# ...
target_link_libraries(foo PRIVATE libsgp4::sgp4)
# target_link_libraries(foo PRIVATE libsgp4::sgp4-static) # or static
```

### External with embedded fallback
```cmake
# ./CMakeLists.txt

# this affects if the libsgp4::sgp4 fallback target is shared or static
set(SGP4_BUILD_SHARED_LIBS ON)
# this affects if the libsgp4::sgp4-static fallback target is available
set(SGP4_BUILD_STATIC_LIBS ON)

# use include instead of add_subdirectory, otherwise targets won't get propagated
include(thirdparty/CMakeLists.txt)
add_subdirectory(src)
```

```cmake
# ./src/CMakeLists.txt

add_library(foo ...)
# ...
target_link_libraries(foo PRIVATE libsgp4::sgp4)
# target_link_libraries(foo PRIVATE libsgp4::sgp4-static) # or static
```

```cmake
# ./thirdparty/CMakeLists.txt

find_package(SGP4 1.1 CONFIG)
if(NOT SGP4_FOUND)
  message(STATUS "Using SGP4 submodule")
  add_subdirectory("${CMAKE_SOURCE_DIR}/thirdparty/sgp4")
endif()
```

## License

    Copyright 2017 Daniel Warner

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.
