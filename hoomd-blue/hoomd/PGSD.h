// Copyright (c) 2009-2023 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#pragma once

#include "hoomd/extern/pgsd.h"
#include <sstream>
#include <stdexcept>
#include <string>

namespace hoomd
    {
namespace detail
    {
/// Utility class to collect common parallel GSD file operations.
class PGSDUtils
    {
    public:
    /// Check and raise an exception if an error occurs
    static void checkError(int retval, const std::string& fname)
        {
        // checkError prints errors and then throws exceptions for common pgsd error codes
        if (retval == PGSD_ERROR_IO)
            {
            std::ostringstream s;
            s << "PGSD: " << strerror(errno) << " - " << fname;
            throw std::runtime_error(s.str());
            }
        else if (retval == PGSD_ERROR_INVALID_ARGUMENT)
            {
            std::ostringstream s;
            s << "PGSD: Invalid argument"
                 " - "
              << fname;
            throw std::invalid_argument(s.str());
            }
        else if (retval == PGSD_ERROR_NOT_A_PGSD_FILE)
            {
            std::ostringstream s;
            s << "PGSD: Not a PGSD file"
                 " - "
              << fname;
            throw std::runtime_error(s.str());
            }
        else if (retval == PGSD_ERROR_INVALID_PGSD_FILE_VERSION)
            {
            std::ostringstream s;
            s << "PGSD: Invalid PGSD file version"
                 " - "
              << fname;
            throw std::runtime_error(s.str());
            }
        else if (retval == PGSD_ERROR_FILE_CORRUPT)
            {
            std::ostringstream s;
            s << "PGSD: File corrupt"
                 " - "
              << fname;
            throw std::runtime_error(s.str());
            }
        else if (retval == PGSD_ERROR_MEMORY_ALLOCATION_FAILED)
            {
            std::ostringstream s;
            s << "PGSD: Memory allocation failed"
                 " - "
              << fname;
            throw std::runtime_error(s.str());
            }
        else if (retval == PGSD_ERROR_NAMELIST_FULL)
            {
            std::ostringstream s;
            s << "PGSD: Namelist full"
                 " - "
              << fname;
            throw std::runtime_error(s.str());
            }
        else if (retval == PGSD_ERROR_FILE_MUST_BE_WRITABLE)
            {
            std::ostringstream s;
            s << "PGSD: File must be writable"
                 " - "
              << fname;
            throw std::runtime_error(s.str());
            }
        else if (retval == PGSD_ERROR_FILE_MUST_BE_READABLE)
            {
            std::ostringstream s;
            s << "PGSD: File must be readable"
                 " - "
              << fname;
            throw std::runtime_error(s.str());
            }
        else if (retval != PGSD_SUCCESS)
            {
            std::ostringstream s;
            s << "PGSD: "
              << "Unknown error " << retval << " writing: " << fname;
            throw std::runtime_error(s.str());
            }
        }
    };
    } // namespace detail
    } // namespace hoomd
