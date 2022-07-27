# Runtime functions supported by hipper

## **1. Device Management**

|   **CUDA**                            |   **HIP**                         |   **hipper**                          |
|---------------------------------------|-----------------------------------|---------------------------------------|
| `cudaDeviceProp`                      | `hipDeviceProp_t`                 | `hipper::deviceProp_t`                |
|                                       |                                   |                                       |
| `cudaChooseDevice`                    | `hipChooseDevice`                 | `hipper::chooseDevice`                |
| `cudaDeviceGetAttribute`              | `hipDeviceGetAttribute`           |                                       |
| `cudaDeviceGetByPCIBusId`             | `hipDeviceGetByPCIBusId`          | `hipper::deviceGetByPCIBusId`         |
| `cudaDeviceGetCacheConfig`            | `hipDeviceGetCacheConfig`         | `hipper::deviceGetCacheConfig`        |
| `cudaDeviceGetLimit`                  | `hipDeviceGetLimit`               | `hipper::deviceGetLimit`              |
| `cudaDeviceGetNvSciSyncAttributes`    |                                   |                                       |
| `cudaDeviceGetP2PAttribute`           |                                   |                                       |
| `cudaDeviceGetPCIBusId`               | `hipDeviceGetPCIBusId`            | `hipper:deviceGetPCIBusId`            |
| `cudaDeviceGetSharedMemConfig`        | `hipDeviceGetSharedMemConfig`     | `hipper:deviceGetSharedMemConfig`     |
| `cudaDeviceGetStreamPriorityRange`    | `hipDeviceGetStreamPriorityRange` | `hipper::deviceGetStreamPriorityRange`|
| `cudaDeviceReset`                     | `hipDeviceReset`                  | `hipper::deviceReset`                 |
| `cudaDeviceSetCacheConfig`            | `hipDeviceSetCacheConfig`         | `hipper::deviceSetCacheConfig`        |
| `cudaDeviceSetLimit`                  | `hipDeviceSetLimit`               | Not actually supported by HIP.        |
| `cudaDeviceSetSharedMemConfig`        | `hipDeviceSetSharedMemConfig`     | `hipper::deviceSetSharedMemConfig`    |
| `cudaDeviceSynchronize`               | `hipDeviceSynchronize`            | `hipper::deviceSynchronize`           |
| `cudaGetDevice`                       | `hipGetDevice`                    | `hipper::getDevice`                   |
| `cudaGetDeviceCount`                  | `hipGetDeviceCount`               | `hipper::getDeviceCount`              |
| `cudaGetDeviceFlags`                  | `hipCtxGetFlags`                  | `hipper::getDeviceFlags`\*            |
| `cudaGetDeviceProperties`             | `hipGetDeviceProperties`          | `hipper::getDeviceProperties`         |
| `cudaIpcCloseMemHandle`               | `hipIpcCloseMemHandle`            |                                       |
| `cudaIpcGetEventHandle`               | `hipIpcGetEventHandle`            |                                       |
| `cudaIpcGetMemHandle`                 | `hipIpcGetMemHandle`              |                                       |
| `cudaIpcOpenEventHandle`              | `hipIpcOpenEventHandle`           |                                       |
| `cudaIpcOpenMemHandle`                | `hipIpcOpenMemHandle`             |                                       |
| `cudaSetDevice`                       | `hipSetDevice`                    | `hipper::setDevice`                   |
| `cudaSetDeviceFlags`                  | `hipSetDeviceFlags`               | `hipper::setDeviceFlags`              |
| `cudaSetValidDevices`                 |                                   |                                       |

### Function Cache Configuration

|   **CUDA**                            |   **HIP**                         |   **hipper**                          |
|---------------------------------------|-----------------------------------|---------------------------------------|
| `cudaFuncCache`                       | `hipFuncCache_t`                  | `hipper::funcCache_t`                 |
|                                       |                                   |                                       |
| `cudaFuncCachePreferNone`             | `hipFuncCachePreferNone`          | `hipper::funcCachePreferNone`         |
| `cudaFuncCachePreferShared`           | `hipFuncCachePreferShared`        | `hipper::funcCachePreferShared`       |
| `cudaFuncCachePreferL1`               | `hipFuncCachePreferL1`            | `hipper::funcCachePreferL1`           |
| `cudaFuncCachePreferEqual`            | `hipFuncCachePreferEqual`         | `hipper::funcCachePreferEqual`        |

### Limits

|   **CUDA**                            |   **HIP**                         |   **hipper**                          |
|---------------------------------------|-----------------------------------|---------------------------------------|
| `cudaLimit`                           | `hipLimit_t`                      | `hipper::limit_t`                     |
|                                       |                                   |                                       |
| `cudaLimitStackSize`                  |                                   |                                       |
| `cudaLimitPrintfFifoSize`             |                                   |                                       |
| `cudaLimitMallocHeapSize`             | `hipLimitMallocHeapSize`          | `hipper::limitMallocHeapSize`         |
| `cudaLimitDevRuntimeSyncDepth`        |                                   |                                       |
|`cudaLimitDevRuntimePendingLaunchCount`|                                   |                                       |
| `cudaLimitMaxL2FetchGranularity`      |                                   |                                       |

### Shared Memory Configuration

|   **CUDA**                            |   **HIP**                         |   **hipper**                          |
|---------------------------------------|-----------------------------------|---------------------------------------|
| `cudaSharedMemConfig`                 | `hipSharedMemConfig_t`            | `hipper::sharedMemConfig_t`           |
|                                       |                                   |                                       |
| `cudaSharedMemBankSizeDefault`        | `hipSharedMemBankSizeDefault`     | `hipper::sharedMemBankSizeDefault`    |
| `cudaSharedMemBankSizeFourByte`       | `hipSharedMemBankSizeFourByte`    | `hipper::sharedMemBankSizeFourByte`   |
| `cudaSharedMemBankSizeEightByte`      | `hipSharedMemBankSizeEightByte`   | `hipper::sharedMemBankSizeEightByte`  |

### Constants

|   **CUDA**                            |   **HIP**                         |   **hipper**                          |
|---------------------------------------|-----------------------------------|---------------------------------------|
| `cudaDeviceScheduleAuto`              | `hipDeviceScheduleAuto`           | `hipper::deviceScheduleAuto`          |
| `cudaDeviceScheduleBlockingSync`      | `hipDeviceScheduleBlockingSync`   | `hipper::deviceScheduleBlockingSync`  |
| `cudaDeviceScheduleMask`              | `hipDeviceScheduleMask`           | `hipper::deviceScheduleMask`          |
| `cudaDeviceScheduleSpin`              | `hipDeviceScheduleSpin`           | `hipper::deviceScheduleSpin`          |
| `cudaDeviceScheduleYield`             | `hipDeviceScheduleYield`          | `hipper::deviceScheduleYield`         |
| `cudaDeviceLmemResizeToMax`           | `hipDeviceLmemResizeToMax`        | `hipper::deviceLmemResizeToMax`       |
| `cudaDeviceMapHost`                   | `hipDeviceMapHost`                | `hipper::deviceMapHost`               |

## **3. Error Handling**

|   **CUDA**                            |   **HIP**                         |   **hipper**                          |
|---------------------------------------|-----------------------------------|---------------------------------------|
| `cudaGetErrorName`                    | `hipGetErrorName`                 | `hipper::getErrorName`                |
| `cudaGetErrorString`                  | `hipGetErrorString`               | `hipper::getErrorString`              |
| `cudaGetLastError`                    | `hipGetLastError`                 | `hipper::getLastError`                |
| `cudaPeekAtLastError`                 | `hipPeekAtLastError`              | `hipper::peekAtLastError`             |

### Error Codes

|   **CUDA**                                |   **HIP**                             |   **hipper**                              |
|-------------------------------------------|---------------------------------------|-------------------------------------------|
| `cudaError_t`                             | `hipError_t`                          | `hipper::error_t`                         |
|                                           |                                       |                                           |
| `cudaSuccess`                             | `hipSuccess`                          | `hippe::success`                          |
| `cudaErrorInvalidValue`                   | `hipErrorInvalidValue`                | `hipper::errorInvalidValue`               |
| `cudaErrorMemoryAllocation`               | `hipErrorOutOfMemory`                 | `hipper::errorOutOfMemory`                |
| `cudaErrorInitializationError`            | `hipErrorNotInitialized`              | `hipper::errorNotInitialized`             |
| `cudaErrorCudartUnloading`                | `hipErrorDeinitialized`               | `hipper::errorDeinitialized`              |
| `cudaErrorProfilerDisabled`               | `hipErrorProfilerDisabled`            | `hipper::errorProfilerDisabled`           |
| `cudaErrorProfilerNotInitialized`         | `hipErrorProfilerNotInitialized`      | `hipper::errorProfilerNotInitialized`     |
| `cudaErrorProfilerAlreadyStarted`         | `hipErrorProfilerAlreadyStarted`      | `hipper::errorProfilerAlreadyStarted`     |
| `cudaErrorProfilerAlreadyStopped`         | `hipErrorProfilerAlreadyStopped`      | `hipper::errorProfilerAlreadyStopped`     |
| `cudaErrorInvalidConfiguration`           | `hipErrorInvalidConfiguration`        | `hipper::errorInvalidConfiguration`       |
| `cudaErrorInvalidPitchValue`              |                                       |                                           |
| `cudaErrorInvalidSymbol`                  | `hipErrorInvalidSymbol`               | `hipper::errorInvalidSymbol`              |
| `cudaErrorInvalidHostPointer`             |                                       |                                           |
| `cudaErrorInvalidDevicePointer`           | `hipErrorInvalidDevicePointer`        | `hipper::errorInvalidDevicePointer`       |
| `cudaErrorInvalidTexture`                 |                                       |                                           |
| `cudaErrorInvalidTextureBinding`          |                                       |                                           |
| `cudaErrorInvalidChannelDescriptor`       |                                       |                                           |
| `cudaErrorInvalidMemcpyDirection`         | `hipErrorInvalidMemcpyDirection`      | `hipper::errorInvalidMemcpyDirection`     |
| `cudaErrorAddressOfConstant`              |                                       |                                           |
| `cudaErrorTextureFetchFailed`             |                                       |                                           |
| `cudaErrorTextureNotBound`                |                                       |                                           |
| `cudaErrorSynchronizationError`           |                                       |                                           |
| `cudaErrorInvalidFilterSetting`           |                                       |                                           |
| `cudaErrorInvalidNormSetting`             |                                       |                                           |
| `cudaErrorMixedDeviceExecution`           |                                       |                                           |
| `cudaErrorNotYetImplemented`              |                                       |                                           |
| `cudaErrorMemoryValueTooLarge`            |                                       |                                           |
| `cudaErrorInsufficientDriver`             | `hipErrorInsufficientDriver`          | `hipper::errorInsufficientDriver`         |
| `cudaErrorInvalidSurface`                 |                                       |                                           |
| `cudaErrorDuplicateVariableName`          |                                       |                                           |
| `cudaErrorDuplicateTextureName`           |                                       |                                           |
| `cudaErrorDuplicateSurfaceName`           |                                       |                                           |
| `cudaErrorDevicesUnavailable`             |                                       |                                           |
| `cudaErrorIncompatibleDriverContext`      |                                       |                                           |
| `cudaErrorMissingConfiguration`           | `hipErrorMissingConfiguration`        | `hipper::errorMissingConfiguration`       |
| `cudaErrorPriorLaunchFailure`             | `hipErrorPriorLaunchFailure`          | `hipper::errorPriorLaunchFailure`         |
| `cudaErrorLaunchMaxDepthExceeded`         |                                       |                                           |
| `cudaErrorLaunchFileScopedTex`            |                                       |                                           |
| `cudaErrorLaunchFileScopedSurf`           |                                       |                                           |
| `cudaErrorSyncDepthExceeded`              |                                       |                                           |
| `cudaErrorLaunchPendingCountExceeded`     |                                       |                                           |
| `cudaErrorInvalidDeviceFunction`          | `hipErrorInvalidDeviceFunction`       | `hipper::errorInvalidDeviceFunction`      |
| `cudaErrorNoDevice`                       | `hipErrorNoDevice`                    | `hipper::errorNoDevice`                   |
| `cudaErrorInvalidDevice`                  | `hipErrorInvalidDevice`               | `hipper::errorInvalidDevice`              |
| `cudaErrorStartupFailure`                 |                                       |                                           |
| `cudaErrorInvalidKernelImage`             | `hipErrorInvalidImage`                | `hipper::errorInvalidImage`               |
| `cudaErrorDeviceUninitialized`            | `hipErrorInvalidContext`              | Not currently supported (CUDA 10.1).      |
| `cudaErrorMapBufferObjectFailed`          | `hipErrorMapFailed`                   | `hipper::errorMapFailed`                  |
| `cudaErrorUnmapBufferObjectFailed`        | `hipErrorUnmapFailed`                 | `hipper::errorUnmapFailed`                |
| `cudaErrorNoKernelImageForDevice`         | `hipErrorNoBinaryForGpu`              | `hipper::errorNoBinaryForGPU`             |
| `cudaErrorECCUncorrectable`               | `hipErrorECCNotCorrectable`           | `hipper::errorECCNotCorrectable`          |
| `cudaErrorUnsupportedLimit`               | `hipErrorUnsupportedLimit`            | `hipper::errorUnsupportedLimit`           |
| `cudaErrorDeviceAlreadyInUse`             |                                       |                                           |
| `cudaErrorPeerAccessUnsupported`          | `hipErrorPeerAccessUnsupported`       | `hippper::errorPeerAccessUnsupported`     |
| `cudaErrorInvalidPtx`                     | `hipErrorInvalidKernelFile`           | `hippper::errorInvalidKernelFile`         |
| `cudaErrorInvalidGraphicsContext`         | `hipErrorInvalidGraphicsContext`      | `hippper::errorInvalidGraphicsContext`    |
| `cudaErrorNvlinkUncorrectable`            |                                       |                                           |
| `cudaErrorJitCompilerNotFound`            |                                       |                                           |
| `cudaErrorInvalidSource`                  | `hipErrorInvalidSource`               | Not currently supported (CUDA 10.1).      |
| `cudaErrorFileNotFound`                   | `hipErrorFileNotFound`                | Not currently supported (CUDA 10.1).      |
| `cudaErrorSharedObjectSymbolNotFound`     | `hipErrorSharedObjectSymbolNotFound`  | `hipper::errorSharedObjectSymbolNotFound` |
| `cudaErrorSharedObjectInitFailed`         | `hipErrorSharedObjectInitFailed`      | `hipper::errorSharedObjectInitFailed`     |
| `cudaErrorOperatingSystem`                | `hipErrorOperatingSystem`             | `hipper::errorOperatingSystem`            |
| `cudaErrorInvalidResourceHandle`          | `hipErrorInvalidHandle`               | `hipper::errorInvalidHandle`              |
| `cudaErrorIllegalState`                   |                                       |                                           |
| `cudaErrorSymbolNotFound`                 | `hipErrorNotFound`                    | Not currently supported (CUDA 10.1).      |
| `cudaErrorNotReady`                       | `hipErrorNotReady`                    | `hipper::errorNotReady`                   |
| `cudaErrorIllegalAddress`                 | `hipErrorIllegalAddress`              | `hipper::errorIllegalAddress`             |
| `cudaErrorLaunchOutOfResources`           | `hipErrorLaunchOutOfResources`        | `hipper::errorLaunchOutOfResources`       |
| `cudaErrorLaunchTimeout`                  | `hipErrorLaunchTimeOut`               | `hipper::errorLaunchTimeOut`              |
| `cudaErrorLaunchIncompatibleTexturing`    |                                       |                                           |
| `cudaErrorPeerAccessAlreadyEnabled`       | `hipErrorPeerAccessAlreadyEnabled`    | `hipper::errorPeerAccessAlreadyEnabled`   |
| `cudaErrorPeerAccessNotEnabled`           | `hipErrorPeerAccessNotEnabled`        | `hipper::errorPeerAccessNotEnabled`       |
| `cudaErrorSetOnActiveProcess`             | `hipErrorSetOnActiveProcess`          | `hipper::errorSetOnActiveProcess`         |
| `cudaErrorContextIsDestroyed`             |                                       |                                           |
| `cudaErrorAssert`                         | `hipErrorAssert`                      | `hipper::errorAssert`                     |
| `cudaErrorTooManyPeers`                   |                                       |                                           |
| `cudaErrorHostMemoryAlreadyRegistered`    | `hipErrorHostMemoryAlreadyRegistered` | `hipper::errorHostMemoryAlreadyRegistered`|
| `cudaErrorHostMemoryNotRegistered`        | `hipErrorHostMemoryNotRegistered`     | `hipper::errorHostMemoryNotRegistered`    |
| `cudaErrorHardwareStackError`             |                                       |                                           |
| `cudaErrorIllegalInstruction`             |                                       |                                           |
| `cudaErrorMisalignedAddress`              |                                       |                                           |
| `cudaErrorInvalidAddressSpace`            |                                       |                                           |
| `cudaErrorInvalidPc`                      |                                       |                                           |
| `cudaErrorLaunchFailure`                  | `hipErrorLaunchFailure`               | `hipper::errorLaunchFailure`              |
| `cudaErrorCooperativeLaunchTooLarge`      | `hipErrorCooperativeLaunchTooLarge`   | Not currently supported (CUDA 9.0).       |
| `cudaErrorNotPermitted`                   |                                       |                                           |
| `cudaErrorNotSupported`                   | `hipErrorNotSupported`                | `hipper::errorNotSupported`               |
| `cudaErrorSystemNotReady`                 |                                       |                                           |
| `cudaErrorSystemDriverMismatch`           |                                       |                                           |
| `cudaErrorCompatNotSupportedOnDevice`     |                                       |                                           |
| `cudaErrorStreamCaptureUnsupported`       |                                       |                                           |
| `cudaErrorStreamCaptureInvalidated`       |                                       |                                           |
| `cudaErrorStreamCaptureMerge`             |                                       |                                           |
| `cudaErrorStreamCaptureUnmatched`         |                                       |                                           |
| `cudaErrorStreamCaptureUnjoined`          |                                       |                                           |
| `cudaErrorStreamCaptureIsolation`         |                                       |                                           |
| `cudaErrorStreamCaptureImplicit`          |                                       |                                           |
| `cudaErrorCapturedEvent`                  |                                       |                                           |
| `cudaErrorStreamCaptureWrongThread`       |                                       |                                           |
| `cudaErrorTimeout`                        |                                       |                                           |
| `cudaErrorGraphExecUpdateFailure`         |                                       |                                           |
| `cudaErrorUnknown`                        | `hipErrorUnknown`                     | `hipper::errorUnknown`                    |
| `cudaErrorApiFailureBase`                 |                                       |                                           |

## **4. Stream Management**

|   **CUDA**                            |   **HIP**                         |   **hipper**                          |
|---------------------------------------|-----------------------------------|---------------------------------------|
| `cudaStream_t`                        | `hipStream_t`                     | `hipper::stream_t`                    |
|                                       |                                   |                                       |
| `cudaStreamAddCallback`               | `hipStreamAddCallback`            | `hipper::streamAddCallback`           |
| `cudaStreamAttachMemAsync`            |                                   |                                       |
| `cudaStreamBeginCapture`              |                                   |                                       |
| `cudaStreamCreate`                    | `hipStreamCreate`                 | `hipper::streamCreate`                |
| `cudaStreamCreateWithFlags`           | `hipStreamCreateWithFlags`        | `hipper::streamCreateWithFlags`       |
| `cudaStreamCreateWithPriority`        | `hipStreamCreateWithPriority`     | `hipper::streamCreateWithPriority`    |
| `cudaStreamDestroy`                   | `hipStreamDestroy`                | `hipper::streamDestroy`               |
| `cudaStreamEndCapture`                |                                   |                                       |
| `cudaStreamGetCaptureInfo`            |                                   |                                       |
| `cudaStreamGetFlags`                  | `hipStreamGetFlags`               | `hipper::streamGetFlags`              |
| `cudaStreamGetPriority`               | `hipStreamGetPriority`            | `hipper::streamGetPriority`           |
| `cudaStreamIsCapturing`               |                                   |                                       |
| `cudaStreamQuery`                     | `hipStreamQuery`                  | `hipper::streamQuery`                 |
| `cudaStreamSynchronize`               | `hipStreamSynchronize`            | `hipper::streamSynchronize`           |
| `cudaStreamWaitEvent`                 | `hipStreamWaitEvent`              | `hipper::streamWaitEvent`             |
| `cudaThreadExchangeStreamCaptureMode` |                                   |                                       |

### Constants

|   **CUDA**                            |   **HIP**                         |   **hipper**                          |
|---------------------------------------|-----------------------------------|---------------------------------------|
| `cudaStreamDefault`                   | `hipStreamDefault`                | `hipper::streamDefault`               |
| `cudaStreamNonBlocking`               | `hipStreamNonBlocking`            | `hipper::streamNonBlocking`           |

## **5. Event Management**

|   **CUDA**                            |   **HIP**                         |   **hipper**                          |
|---------------------------------------|-----------------------------------|---------------------------------------|
| `cudaEvent_t`                         | `hipEvent_t`                      | `hipper::event_t`                     |
|                                       |                                   |                                       |
| `cudaEventCreate`                     | `hipEventCreate`                  | `hipper::eventCreate`                 |
| `cudaEventCreateWithFlags`            | `hipEventCreateWithFlags`         | `hipper::eventCreateWithFlags`        |
| `cudaEventDestroy`                    | `hipEventDestroy`                 | `hipper::eventDestroy`                |
| `cudaEventElapsedTime`                | `hipEventElapsedTime`             | `hipper::eventElapsedTime`            |
| `cudaEventQuery`                      | `hipEventQuery`                   | `hipper::eventQuery`                  |
| `cudaEventRecord`                     | `hipEventRecord`                  | `hipper::eventRecord`                 |
| `cudaEventSynchronize`                | `hipEventSynchronize`             | `hipper::eventSynchronize`            |

### Constants

|   **CUDA**                            |   **HIP**                         |   **hipper**                          |
|---------------------------------------|-----------------------------------|---------------------------------------|
| `cudaEventDefault`                    | `hipEventDefault`                 | `hipper::eventDefault`                |
| `cudaEventBlockingSync`               | `hipEventBlockingSync`            | `hipper::eventBlockingSync`           |
| `cudaEventDisableTiming`              | `hipEventDisableTiming`           | `hipper::eventDisableTiming`          |
| `cudaEventInterprocess`               | `hipEventInterprocess`            | `hipper::eventInterprocess`           |

## **7. Execution Control**

|   **CUDA**                                |   **HIP**                             |   **hipper**                          |
|-------------------------------------------|---------------------------------------|---------------------------------------|
| `cudaFuncAttributes`                      | `hipFuncAttributes`                   | `hipper::funcAttributes_t`            |
|                                           |                                       |                                       |
| `cudaFuncGetAttributes`                   | `hipFuncGetAttributes`                | `hipper::funcGetAttributes`           |
| `cudaFuncSetAttribute`                    |                                       |                                       |
| `cudaFuncSetCacheConfig`                  | `hipFuncSetCacheConfig`               | `hipper::funcSetCacheConfig`          |
| `cudaFuncSetSharedMemConfig`              |                                       |                                       |
| `cudaGetParameterBuffer`                  |                                       |                                       |
| `cudaGetParameterBufferV2`                |                                       |                                       |
| `cudaLaunchCooperativeKernel`             | `hipLaunchCooperativeKernel`          |                                       |
| `cudaLaunchCooperativeKernelMultiDevice`  |`hipLaunchCooperativeKernelMultiDevice`|                                       |
| `cudaLaunchHostFunc`                      |                                       |                                       |
| `cudaLaunchKernel`                        | `hipLaunchKernel`                     |                                       |
| `cudaSetDoubleForDevice`                  |                                       |                                       |
| `cudaSetDoubleForHost`                    |                                       |                                       |

### Grid dimensions (device)

|   **CUDA**                            |   **HIP**                         |   **hipper**                          |
|---------------------------------------|-----------------------------------|---------------------------------------|
| `threadIdx.x`                         | `hipThreadIdx_x`                  | `hipper::threadIndex().x`             |
| `blockDim.x`                          | `hipBlockDim_x`                   | `hipper::blockSize().x`               |
| `blockIdx.x`                          | `hipBlockIdx_x`                   | `hipper::blockIndex().x`              |
| `gridDim.x`                           | `hipGridDim_x`                    | `hipper::gridSize().x`                |

## **9. Memory Management**

|   **CUDA**                            |   **HIP**                         |   **hipper**                          |
|---------------------------------------|-----------------------------------|---------------------------------------|
| `cudaArrayGetInfo`                    |                                   |                                       |
| `cudaFree`                            | `hipFree`                         | `hipper::free`                        |
| `cudaFreeArray`                       | `hipFreeArray`                    |                                       |
| `cudaFreeHost`                        | `hipHostFree`                     | `hipper::hostFree`                    |
| `cudaFreeMipmappedArray`              |                                   |                                       |
| `cudaGetMipmappedArrayLevel`          |                                   |                                       |
| `cudaGetSymbolAddress`                | `hipGetSymbolAddress`             | `hipper::getSymbolAddress`            |
| `cudaGetSymbolSize`                   | `hipGetSymbolSize`                | `hipper::getSymbolSize`               |
| `cudaHostAlloc`                       | `hipHostMalloc`                   | `hipper::hostMalloc`                  |
| `cudaHostGetDevicePointer`            | `hipHostGetDevicePointer`         | `hipper::hostGetDevicePointer`        |
| `cudaHostGetFlags`                    | `hipHostGetFlags`                 | `hipper::hostGetFlags`                |
| `cudaHostRegister`                    | `hipHostRegister`                 | `hipper::hostRegister`                |
| `cudaHostUnregister`                  | `hipHostUnregister`               | `hipper::hostUnregister`              |
| `cudaMalloc`                          | `hipMalloc`                       | `hipper::malloc`                      |
| `cudaMalloc3D`                        | `hipMalloc3D`                     |                                       |
| `cudaMalloc3DArray`                   | `hipMalloc3DArray`                |                                       |
| `cudaMallocArray`                     | `hipMallocArray`                  |                                       |
| `cudaMallocHost`                      | `hipHostMalloc`                   | Use `hipper::hostMalloc`.             |
| `cudaMallocManaged`                   | `hipMallocManaged`                | `hipper::mallocManaged`               |
| `cudaMallocMipmappedArray`            |                                   |                                       |
| `cudaMallocPitch`                     |                                   |                                       |
| `cudaMemAdvise`                       |                                   |                                       |
| `cudaMemGetInfo`                      | `hipMemGetInfo`                   | `hipper::MemGetInfo`                  |
| `cudaMemPrefetchAsync`                |                                   |                                       |
| `cudaMemRangeGetAttribute`            |                                   |                                       |
| `cudaMemRangeGetAttributes`           |                                   |                                       |
| `cudaMemcpy`                          | `hipMemcpy`                       | `hipper::memcpy`                      |
| `cudaMemcpy2D`                        | `hipMemcpy2D`                     |                                       |
| `cudaMemcpy2DArrayToArray`            |                                   |                                       |
| `cudaMemcpy2DAsync`                   | `hipMemcpy2DAsync`                |                                       |
| `cudaMemcpy2DFromArray`               | `hipMemcpy2DFromArray`            |                                       |
| `cudaMemcpy2DFromArrayAsync`          | `hipMemcpy2DFromArrayAsync`       |                                       |
| `cudaMemcpy2DToArray`                 | `hipMemcpy2DToArray`              |                                       |
| `cudaMemcpy2DToArrayAsync`            |                                   |                                       |
| `cudaMemcpy3D`                        | `hipMemcpy3D`                     |                                       |
| `cudaMemcpy3DAsync`                   | `hipMemcpy3DAsync`                |                                       |
| `cudaMemcpy3DPeer`                    |                                   |                                       |
| `cudaMemcpy3DPeerAsync`               |                                   |                                       |
| `cudaMemcpyAsync`                     | `hipMemcpyAsync`                  | `hipper::memcpyAsync`                 |
| `cudaMemcpyFromSymbol`                | `hipMemcpyFromSymbol`             | `hipper::memcpyFromSymbol`            |
| `cudaMemcpyFromSymbolAsync`           | `hipMemcpyFromSymbolAsync`        | `hipper::memcpyFromSymbolAsync`       |
| `cudaMemcpyPeer`                      | `hipMemcpyPeer`                   | `hipper::memcpyPeer`                  |
| `cudaMemcpyPeerAsync`                 | `hipMemcpyPeerAsync`              | `hipper::memcpyPeerAsync`             |
| `cudaMemcpyToSymbol`                  | `hipMemcpyToSymbol`               | `hipper::memcpyToSymbol`              |
| `cudaMemcpyToSymbolAsync`             | `hipMemcpyToSymbolAsync`          | `hipper::memcpyToSymbolAsync`         |
| `cudaMemset`                          | `hipMemset`                       | `hipper::memset`                      |
| `cudaMemset2D`                        | `hipMemset2D`                     |                                       |
| `cudaMemset2DAsync`                   | `hipMemset2DAsync`                |                                       |
| `cudaMemset3D`                        | `hipMemset3D`                     |                                       |
| `cudaMemset3DAsync`                   | `hipMemset3DAsync`                |                                       |
| `cudaMemsetAsync`                     | `hipMemsetAsync`                  | `hipper::memsetAsync`                 |
| `make_cudaExtent`                     | `make_hipExtent`                  |                                       |
| `make_cudaPitchedPtr`                 | `make_hipPitchedPtr`              |                                       |
| `make_cudaPos`                        | `make_hipPos`                     |                                       |

### Memory Copy Kinds

|   **CUDA**                            |   **HIP**                         |   **hipper**                          |
|---------------------------------------|-----------------------------------|---------------------------------------|
| `cudaMemcpyKind`                      | `hipMemcpyKind`                   | `hipper::memcpyKind_t`                |
|                                       |                                   |                                       |
| `cudaMemcpyHostToHost`                | `hipMemcpyHostToHost`             | `hipper::memcpyHostToHost`            |
| `cudaMemcpyHostToDevice`              | `hipMemcpyHostToDevice`           | `hipper::memcpyHostToDevice`          |
| `cudaMemcpyDeviceToHost`              | `hipMemcpyDeviceToHost`           | `hipper::memcpyDeviceToHost`          |
| `cudaMemcpyDeviceToDevice`            | `hipMemcpyDeviceToDevice`         | `hipper::mempcyDeviceToDevice`        |
| `cudaMemcpyDefault`                   | `hipMemcpyDefault`                | `hipper::memcpyDefault`               |

### Constants

|   **CUDA**                            |   **HIP**                         |   **hipper**                          |
|---------------------------------------|-----------------------------------|---------------------------------------|
| `cudaHostAllocDefault`                | `hipHostMallocDefault`            | `hipper::hostMallocDefault`           |
| `cudaHostAllocMapped`                 | `hipHostMallocMapped`             | `hipper::hostMallocMapped`            |
| `cudaHostAllocPortable`               | `hipHostMallocPortable`           | `hipper::hostMallocPortable`          |
| `cudaHostAllocWriteCombined`          | `hipHostMallocWriteCombined`      | `hipper::hostMallocWriteCombined`     |
| `cudaHostRegisterDefault`             | `hipHostRegisterDefault`          | `hipper::hostRegisterDefault`         |
| `cudaHostRegisterIoMemory`            | `hipHostRegisterIoMemory`         | `hipper::hostRegisterIoMemory`        |
| `cudaHostRegisterMapped`              | `hipHostRegisterMapped`           | `hipper::hostRegisterMapped`          |
| `cudaHostRegisterPortable`            | `hipHostRegisterPortable`         | `hipper::hostRegisterPortable`        |
| `cudaMemAttachGlobal`                 | `hipMemAttachGlobal`              | `hipper::memAttachGlobal`             |
| `cudaMemAttachHost`                   | `hipMemAttachHost`                | `hipper::memAttachHost`               |

## **32. Profiler Control**

|   **CUDA**                            |   **HIP**                         |   **hipper**                          |
|---------------------------------------|-----------------------------------|---------------------------------------|
| `cudaProfilerInitialize`              |                                   |                                       |
| `cudaProfilerStart`                   | `hipProfilerStart`                | `hipper::profilerStart`\*             |
| `cudaProfilerStop`                    | `hipProfilerStop`                 | `hipper::profilerStop`\*              |

\*: feature only available with `HIPPER_USE_DEPRECATED`.

