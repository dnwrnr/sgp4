SGP4 Library
============

This repository contains a C++ SGP4 library plus a set of experimental acceleration components for batch propagation and conjunction-screening workflows.

The current branch builds around two main paths:

- A scalar SGP4 implementation in libsgp4::SGP4.
- A batch-oriented propagator in libsgp4::SGP4Batch with SoA storage, OpenMP parallelism, and SIMD-capable code paths with scalar and AVX2 fallbacks.

What Changed
------------

The repository was cleaned up to build on the current Linux toolchain without compile errors, and the documentation was aligned to the code that is actually present in this branch.

Build and code fixes applied:

- Restored the missing SGP4::SetAlongTrackBias API required by full_test.
- Applied the along-track bias correction in the active libsgp4 SGP4 implementation.
- Removed warning-producing local variable names in SGP4Batch that shadowed member fields.
- Removed an unused benchmark local variable.
- Made SLEEF linkage optional so the project still builds when that library is not installed.

Current Optimization Characteristics
------------------------------------

The codebase includes the following implemented optimization work:

- Structure-of-arrays storage in SGP4Batch for better cache locality.
- Batch propagation interfaces for large TLE sets.
- OpenMP parallel loops in batch and screening utilities.
- SIMD abstraction in SGP4Batch with AVX-512, AVX2, and scalar fallback code paths.
- Experimental runtime JIT kernel generation in JitPropagator.
- Experimental spatial filtering utilities such as LinearBVH, SpatialHash, SpatialPartition, and TemporalPruner.

Important caveats:

- Several performance-oriented utilities are experimental and not wired into the default library API.
- The JIT path currently generates AVX-512-oriented source and is hardware/toolchain dependent.
- Some standalone benchmark and research programs in the repository are environment-specific and are not part of the default build.
- Performance claims depend heavily on CPU ISA support, compiler behavior, and workload shape.

Build
-----

Requirements:

- CMake 3.13 or newer
- A C++17 compiler
- OpenMP support

Optional:

- SLEEF. If installed, CMake will link against it. If not installed, the library still builds.

Build commands:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

Validated on the current Linux environment with GCC 13 using AVX2/FMA flags from the project CMake configuration.

Targets
-------

Default configured targets include:

- libsgp4 static library
- libsgp4s shared library
- full_test
- benchmark
- sattrack
- runtest
- passpredict

Usage
-----

Scalar propagation:

```cpp
#include "libsgp4/SGP4.h"
#include "libsgp4/Tle.h"

using namespace libsgp4;

Tle tle("TEST",
   "1 25544U 98067A   24001.00000000  .00016717  00000-0  10270-3 0  9996",
   "2 25544  51.6417  25.1234 0005000 120.0000 240.0000 15.50000000000000");

SGP4 model(tle);
Eci state = model.FindPosition(60.0);
```

Batch propagation:

```cpp
#include "libsgp4/SGP4Batch.h"

using namespace libsgp4;

std::vector<Tle> tles = /* load TLEs */;
SGP4Batch batch(tles);
std::vector<Eci> results;

batch.Propagate(60.0, results, SGP4Batch::MathMode::Standard);
```

Experimental JIT propagation:

```cpp
#include "libsgp4/JitPropagator.h"

using namespace libsgp4;

SGP4Batch batch(tles);
JitPropagator jit(batch);
std::vector<Eci> results;

jit.Propagate(60.0, batch, results);
```

Repository Notes
----------------

- There are duplicate top-level source files and libsgp4 source files in this repository. The active CMake build uses the libsgp4 sources for the library.
- A committed binary artifact named run_5year_sim_vectorized is present in the repository root. It should not be committed as source control content if you want a clean Git history.

License
-------

Copyright 2017 Daniel Warner

Licensed under the Apache License, Version 2.0.
