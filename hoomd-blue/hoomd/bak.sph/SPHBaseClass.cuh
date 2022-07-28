// Copyright (c) 2009-2016 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

// Maintainer: schirwon

#include "hoomd/BoxDim.h"
#include <vector>

cudaError_t gpu_sphbase_apply_body_force(
                                 Scalar4 *d_force,
                                 Scalar4 *d_velocity,
                                 Scalar3 bforce,
                                 unsigned int *d_index_array,
                                 unsigned int group_size);