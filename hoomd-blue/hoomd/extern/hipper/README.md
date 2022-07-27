# hipper

hipper is a compatability layer for the CUDA and HIP runtime APIs. It relaxes
dependencies for projects that would like to *optionally* use HIP to support
multiple GPU architectures, but do not want to require it. hipper only implements
the features that are common to both CUDA and HIP. These features must also be
supported across multiple versions of CUDA for maximum cross-compatibility.

## Quick start

hipper is meant to be a header-only drop in for `<cuda_runtime.h>`
or `<hip/hip_runtime.h>`:

```
#include <hipper/hipper_runtime.h>
```

One of the GPU runtime backends must be selected at compile time. These
can be either `HIPPER_CUDA` (to use the CUDA runtime) or `HIPPER_HIP`
(to use the HIP runtime). If `HIPPER_HIP` is selected, then HIP options
can still be used to select between the CUDA or HIP backends.

A few defines are additionally provided for interoperability in device code:

* `HIPPER_PLATFORM_NVCC` if targeting nvcc platform; otherwise undefined.
* `HIPPER_PLATFORM_HCC` if target hcc platform; otherwise undefined.
* `HIPPER_DEVICE_COMPILE` is 1 if compiling device code; otherwise undefined.

Refer to the HIP documentation for more detail about under which conditions
each of these flags is set, and how it can be used to replace CUDA-specific checks.

## Runtime functions

Most of the runtime API has been placed in the `hipper` namespace, and it is easy to
convert between CUDA or HIP commands to `hipper` commands by replacing either
`cuda` or `hip` with `hipper::` in the function signature.

```
cudaMemset -> hipper::memset
cudaMemcpyAsync -> hipper::memcpyAsync
```

See the [detailed API](doc/runtime.md) for details about what functions are supported.

Most of the CUDA or HIP data types or defines have also been spoofed in the `hipper`
namespace as enums or constants. You should work with the `hipper` types of these
parameters, which will usually be automatically converted as needed to pass to the
backend functions.

```
float* a;
const int N = 4;
hipper::error_t code = hipper::malloc(reinterpret_cast<void**>(&a), N*sizeof(float));
if (code != hipper::success)
    {
    std::cout << "GPU error: " << hipper::getErrorString(code) << std::endl;
    exit(code);
    }
```

The caller is responsible for only using properties of data types that are present in
the underlying API (e.g., members of `hipper::deviceProp_t` common to CUDA and HIP),
or a compile error may be generated. If certain features are required in CUDA builds,
use the `HIPPER_PLATFORM_NVCC` define.

## Kernel launch

Kernel launch is abstracted inside of a `KernelLauncher` class. The class has many different
constructors that allow different variations of specifying the launch dimensions and additional
resources (shared memory, streams). Once constructed, the launcher can be used to execute kernels:

```
hipper::KernelLauncher launcher(blocks, threadsPerBlock, sharedBytes, stream);
launcher(kernel, args...)
```

Within a kernel, a thread can be automatically mapped to a one-dimensional rank using the
`hipper::threadRank<GridDim,BlockDim>` function template.

```
int idx = threadRank<1,1>();
if (idx >= num_threads) return;
```

Currently, specializations are provided for `GridDim=1` and `BlockDim=1,2,3`, but additional
ones can be specified. See also the device functions for accessing the dimensions of the
launch configuration if more advanced control is required.

## Deprecated features

HIP has some deprecated features that are fully functional in CUDA. If you want to enable
these, you can opt in:

```
#define HIPPER_USE_DEPRECATED
#include <hipper/hipper_runtime.h>
```

## CMake setup & installation

hipper has a functional CMake build system if you would like to install it to a shared location
on your system. This installation can be detected using `find_package` in `CONFIG` mode.
Alternatively, hipper can be included directly in an existing project via `add_subdirectory`.
In either case, please link against the `hipper::hipper` target.

## Contributing

Support for additional features or defines is always welcome! Please submit a pull request in the
style of the existing code, including a new entry in the API documentation. If you identify any
issues with the code, please create an issue, and we will try to attend to it promptly.
