/* ---------------------------------------------------------
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