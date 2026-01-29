/* ---------------------------------------------------------
Copyright (c) 2025-2026 David Krach, Daniel Rostan.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
   may be used to endorse or promote products derived from this software without
   specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

maintainer: dkrach, david.krach@mib.uni-stuttgart.de
----------------------------------------------------------*/


/*! last bit 0 => type is SOLID
    last bit set => type is FLUID
    bit FLUID<<N set => type is fluidN (in this case the FLUID sould be set too)
*/

#ifndef __SolidFluidTypeBit_H__
#define __SolidFluidTypeBit_H__

namespace hoomd 
{
namespace sph 
{
struct SolidFluidTypeBit
    {
    enum Enum
        {
        NONE = 0,
        SOLID = 1<<1,
        FLUID = 1<<2,
        FLUID1 = FLUID<<1,
        FLUID2 = FLUID<<2,
        //FLUID3 = FLUID<<3,
        //FLUID4 = FLUID<<4,
        //FLUID5 = FLUID<<5,
        };
    };
    
//! Helper funciton to lookup type properties where the type id is storead as a `Scalar`
inline
bool checksolid(const unsigned int* type_props, Scalar mytype)
    {
    return type_props[__scalar_as_int(mytype)] & SolidFluidTypeBit::SOLID;
    }

inline
bool checkfluid(const unsigned int* type_props, Scalar mytype)
    {
    return type_props[__scalar_as_int(mytype)] & SolidFluidTypeBit::FLUID;
    }

inline
bool checkfluid1(const unsigned int* type_props, Scalar mytype)
    {
    return type_props[__scalar_as_int(mytype)] & SolidFluidTypeBit::FLUID1;
    }

inline
bool checkfluid2(const unsigned int* type_props, Scalar mytype)
    {
    return type_props[__scalar_as_int(mytype)] & SolidFluidTypeBit::FLUID2;
    }

} // end namespace sph 
} // end namespace hoomd 

#endif