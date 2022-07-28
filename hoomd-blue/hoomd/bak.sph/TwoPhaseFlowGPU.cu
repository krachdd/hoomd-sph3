// Copyright (c) 2009-2016 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

// Maintainer: schirwon

#include "TwoPhaseFlowGPU.cuh"
//#include "EquationOfState.h"
#include"GeneralGPUFunctions.cuh"
#include"GeneralGPUDeviceFunctions.cuh"

#include "SmoothingKernel.h"
#include "StateEquations.h"
#include "SolidFluidTypeBit.h"

// Constant memory for gridpoint weighting
#define CONSTANT_SIZE 2048

//! Helper function to compute particle pressures
template<StateEquationType SET_>
__global__ void gpu_tpf_compute_pressure_kernel(
                                            Scalar3 *d_dpe,
                                            unsigned int *d_index_array,
                                            unsigned int group_size,
                                            StateEquation<SET_> m_eos
                                               )
{
    // determine which particle this thread works on
    int group_idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (group_idx < group_size)
    {
        // Read particle index
        unsigned int i = d_index_array[group_idx];
        
        // Evaluate pressure
        d_dpe[i].y = m_eos.Pressure(d_dpe[i].x);
    }
}
template<StateEquationType SET_>
cudaError_t gpu_tpf_compute_pressure(
                                 Scalar3 *d_dpe,
                                 unsigned int *d_index_array,
                                 unsigned int group_size,
                                 StateEquation<SET_> *m_eos)
{
    // setup the grid to run the kernel
    unsigned int block_size=256;
    dim3 grid( (group_size / block_size) + 1, 1, 1);
    dim3 threads(block_size, 1, 1);
    gpu_tpf_compute_pressure_kernel<<< grid, threads >>>(
        d_dpe,
        d_index_array,
        group_size,
        *m_eos);
    
    return cudaSuccess;
}

/*! Helper function to compute particle number density
 * \post Particle number densities are stores in charge array
 */
template<SmoothingKernelType KT_>
__global__ void gpu_tpf_compute_ndensity_kernel(Scalar4 *d_pos,
                                            Scalar3 *d_dpe,
                                            Scalar4 *d_vel,
                                            Scalar *d_h,
                                            unsigned int *d_n_neigh,
                                            unsigned int *d_nlist,
                                            unsigned int *d_head_list,
                                            unsigned int N,
                                            const BoxDim box,
                                            Scalar m_ch,
                                            Scalar m_rcutsq,
                                            Scalar m_rho0,
                                            bool m_const_slength,
                                            SmoothingKernel<KT_> m_skernel,
                                            unsigned int *d_type_property_map,
                                            bool densinty_sum)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    
     if (i < N)
     {
         unsigned int myHead, size;
         Scalar w0 = m_skernel.w0(m_ch);
         // Access the particle's position
         Scalar3 pi;
         pi.x = d_pos[i].x;
         pi.y = d_pos[i].y;
         pi.z = d_pos[i].z;
         
        // Determine particle i type
        const unsigned int typei_props = d_type_property_map[__scalar_as_int(d_pos[i].w)];
        bool i_issolid = typei_props & SolidFluidTypeBit::SOLID;
        bool i_isfluid = typei_props & SolidFluidTypeBit::FLUID;

        // Do not compute number density based on summation
        if (!(!densinty_sum && !(i_issolid)) )
        {

            if ( i_issolid )
            {
                // If particle is a fictitious solid, compute fluid normalization constant
                // and store it in number density array
                d_dpe[i].x = 0;
            }
            else
            {
                // Initialize number density with self density of kernel
                d_dpe[i].x = m_const_slength ? w0 : m_skernel.w0(d_h[i]);
            }

            {
                // Loop over all of the neighbors of this particle
                myHead = d_head_list[i];
                size = (unsigned int)d_n_neigh[i];
                for (unsigned int j = 0; j < size; j++)
                {
                    // Index of neighbor
                    unsigned int k = d_nlist[myHead + j];
                    
                    // Read particle j type
                    bool j_issolid = d_type_property_map[__scalar_as_int(d_pos[k].w)] & SolidFluidTypeBit::SOLID;
                    
                    // If both particles are solid, continue with next neighbor in loop
                    if ( !(i_issolid && j_issolid) )
                    {

                        // Access neighbor position
                        Scalar3 pj;
                        pj.x = d_pos[k].x;
                        pj.y = d_pos[k].y;
                        pj.z = d_pos[k].z;
                        
                        // Compute distance vector
                        Scalar3 dx = pj - pi;
                        
                        // Apply periodic boundary conditions
                        dx = box.minImage(dx);
                        
                        // Calculate squared distance
                        Scalar rsq = dot(dx, dx);
                        
                        // If particle distance is too large, continue with next neighbor in loop
                        if (! (m_const_slength && rsq > m_rcutsq) )
                        {

                            // Calculate distance
                            Scalar r = sqrt(rsq);
                            
                            // If i_issolid and j_issolid - This part is not computed
                            // If i_issolid and j_isfluid - Add contribution to normalization constant
                            // If j_isfluid and j_issolid - Add contribution to particle number density
                            // If j_isfluid and j_isfluid - Add contribution to particle number density
                            d_dpe[i].x += m_const_slength ? m_skernel.wij(m_ch,r) : m_skernel.wij(Scalar(0.5)*(d_h[i]+d_h[k]),r);
                        }
                    }
                } // End of neighbor loop
            
                // Compute mass density from number density if particle i is a fluid particle
                // rho_i = m_i * \sum_j wij
                if ( densinty_sum && i_isfluid )
                    d_dpe[i].x= d_dpe[i].x * d_vel[i].w;

                // Do not compute number density based on summation
            }
        }

    } // End of particle loop
}

template<SmoothingKernelType KT_>
cudaError_t gpu_tpf_compute_ndensity(Scalar4 *d_pos,
                                 Scalar3 *d_dpe,
                                 Scalar4 *d_vel,
                                 Scalar *d_h,
                                 unsigned int *d_n_neigh,
                                 unsigned int *d_nlist,
                                 unsigned int *d_head_list,
                                 unsigned int N,
                                 const BoxDim& box,
                                 Scalar m_ch,
                                 Scalar m_rcutsq,
                                 Scalar m_rho0,
                                 bool m_const_slength,
                                 SmoothingKernel<KT_> *m_skernel,
                                 unsigned int *d_type_property_map,
                                 bool densinty_sum
)
{
    // setup the grid to run the kernel
    unsigned int block_size=256;
    dim3 grid( (N / block_size) + 1, 1, 1);
    dim3 threads(block_size, 1, 1);
    
    gpu_tpf_compute_ndensity_kernel<<< grid, threads >>>(
        d_pos,
        d_dpe,
        d_vel,
        d_h,
        d_n_neigh,
        d_nlist,
        d_head_list,
        N,
        box,
        m_ch,
        m_rcutsq,
        m_rho0,
        m_const_slength,
        *m_skernel,
        d_type_property_map,
        densinty_sum);
    
    return cudaSuccess;
}


/*! Helper function to compute renormalized number density
 * \post Fluid particle number density field is renormalized
 */
template<SmoothingKernelType KT_>
__global__ void gpu_tpf_compute_ndensityrenormalization_kernel(Scalar4 *d_pos,
                                                           Scalar3 *d_dpe,
                                                           Scalar4 *d_vel,
                                                           Scalar *d_h,
                                                           unsigned int *d_n_neigh,
                                                           unsigned int *d_nlist,
                                                           unsigned int *d_head_list,
                                                           unsigned int *d_index_array,
                                                           unsigned int group_size,
                                                           const BoxDim box,
                                                           Scalar m_rcutsq,
                                                           Scalar m_ch,
                                                           bool m_const_slength,
                                                           SmoothingKernel<KT_> m_skernel,
                                                           unsigned int *d_type_property_map
)
{
    
    // determine which particle this thread works on
    int group_idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (group_idx < group_size)
    {
        // Read particle index
        Scalar w0 = m_skernel.w0(m_ch);
        unsigned int i = d_index_array[group_idx];
        unsigned int myHead, size;
        // Access the particle's position
        Scalar3 pi;
        pi.x = d_pos[i].x;
        pi.y = d_pos[i].y;
        pi.z = d_pos[i].z;
        
        Scalar mi = d_vel[i].w;
        Scalar rhoi = d_dpe[i].x;

        // First compute renormalization factor
        // Initialize with self density of kernel
        Scalar normalization = m_const_slength ? w0 : m_skernel.w0(d_h[i]);
        normalization = normalization * ( mi / rhoi );

        // Check if fluid is of phase 1
        bool isfluid1 = d_type_property_map[__scalar_as_int(d_pos[i].w)] & SolidFluidTypeBit::FLUID1;

        // Loop over all of the neighbors of this particle
        // and compute normalization constant normwij = \sum_j wij*Vj
        myHead = d_head_list[i];
        size = (unsigned int)d_n_neigh[i];

        //bool skip(false);
        for (unsigned int j = 0; j < size; j++)
        {
            // Index of neighbor
            unsigned int k = d_nlist[myHead + j];
            
            // Only interpolate over fluid domain of same phase
            if ((bool)(d_type_property_map[__scalar_as_int(d_pos[k].w)] & SolidFluidTypeBit::FLUID1) != isfluid1)
                continue;
            
            // Access neighbor position
            Scalar3 pj;
            pj.x = d_pos[k].x;
            pj.y = d_pos[k].y;
            pj.z = d_pos[k].z;
                
            // Compute distance vector
            Scalar3 dx = pj - pi;
                
            // Apply periodic boundary conditions
            dx = box.minImage(dx);
                
            // Calculate squared distance
            Scalar rsq = dot(dx, dx);
                
            // If particle distance is too large, continue with next neighbor in loop
            if (! (m_const_slength && rsq > m_rcutsq))
            {
                // Calculate distance
                Scalar r = sqrt(rsq);
                    
                // Add contribution to renormalization
                Scalar Vj =  d_vel[k].w / d_dpe[k].x ;
                normalization += m_const_slength ? Vj*m_skernel.wij(m_ch,r) : Vj*m_skernel.wij(Scalar(0.5)*(d_h[i]+d_h[k]),r);

            }

        } // End of neighbor loop

        normalization = Scalar(1.0)/normalization;

        // Initialize density with normalized kernel self density
        d_dpe[i].x = m_const_slength ? w0*(mi*normalization): m_skernel.w0(d_h[i])*(mi*normalization);
        // Loop over all of the neighbors of this particle
        // and compute renormalied density rho_i = \sum_j wij*mj / normwij
        myHead = d_head_list[i];
        size = (unsigned int)d_n_neigh[i];
        for (unsigned int j = 0; j < size; j++)
        {
            // Index of neighbor
            unsigned int k = d_nlist[myHead + j];

            // Only interpolate over fluid domain of same phase
            if ((bool)(d_type_property_map[__scalar_as_int(d_pos[k].w)] & SolidFluidTypeBit::FLUID1) != isfluid1)
                continue;

            // Access neighbor position
            Scalar3 pj;
            pj.x = d_pos[k].x;
            pj.y = d_pos[k].y;
            pj.z = d_pos[k].z;

            // Compute distance vector
            Scalar3 dx = pj - pi;

            // Apply periodic boundary conditions
            dx = box.minImage(dx);

            // Calculate squared distance
            Scalar rsq = dot(dx, dx);

            // If particle distance is too large, continue with next neighbor in loop
            if (! (m_const_slength && rsq > m_rcutsq) )
            {

                // Calculate distance
                Scalar r = sqrt(rsq);

                // Add contribution to normalized density interpolation
                Scalar factor =  d_vel[k].w * normalization ;
                d_dpe[i].x += m_const_slength ? factor*m_skernel.wij(m_ch,r) : factor*m_skernel.wij(Scalar(0.5)*(d_h[i]+d_h[k]),r);
            }
       }

    } // End of particle loop
}

template<SmoothingKernelType KT_>
cudaError_t gpu_tpf_compute_ndensityrenormalization(Scalar4 *d_pos,
                                                Scalar3 *d_dpe,
                                                Scalar4 *d_vel,
                                                Scalar *d_h,
                                                unsigned int *d_n_neigh,
                                                unsigned int *d_nlist,
                                                unsigned int *d_head_list,
                                                unsigned int *d_index_array,
                                                unsigned int group_size,
                                                const BoxDim& box,
                                                Scalar m_rcutsq,
                                                Scalar m_ch,
                                                bool m_const_slength,
                                                SmoothingKernel<KT_> *m_skernel,
                                                unsigned int *d_type_property_map
)
{
    // setup the grid to run the kernel
    unsigned int block_size=256;
    dim3 grid( (group_size / block_size) + 1, 1, 1);
    dim3 threads(block_size, 1, 1);
    
    gpu_tpf_compute_ndensityrenormalization_kernel<<< grid, threads >>>(
        d_pos,
        d_dpe,
        d_vel,
        d_h,
        d_n_neigh,
        d_nlist,
        d_head_list,
        d_index_array,
        group_size,
        box,
        m_rcutsq,
        m_ch,
        m_const_slength,
        *m_skernel,
        d_type_property_map);
    
    return cudaSuccess;
}

/*! Helper function to compute fictitious solid particle properties (pressures and velocities)
 * \pre Ghost particle number densities (i.e. charge array) must be up-to-date
 * \pre Normalization constant \sum_fluids Wij must be computed and stored in orientation array
 * \post Fictitious particle properties are computed and stored in orentiation array
 */
template<SmoothingKernelType KT_, StateEquationType SET_>
__global__ void gpu_tpf_compute_noslip_kernel(Scalar3 *d_vf,
                                          Scalar4 *d_pos,
                                          Scalar3 *d_dpe,
                                          Scalar4 *d_velocity,
                                          Scalar3 *d_accel,
                                          Scalar *d_h,
                                          unsigned int *d_n_neigh,
                                          unsigned int *d_nlist,
                                          unsigned int *d_head_list,
                                          unsigned int *d_index_array,
                                          unsigned int group_size,
                                          const BoxDim box,
                                          Scalar m_rcutsq,
                                          Scalar m_ch,
                                          Scalar m_rho0,
                                          bool m_const_slength,
                                          Scalar3 g_accel,
                                          StateEquation<SET_> m_eos,
                                          SmoothingKernel<KT_> m_skernel,
                                          unsigned int *d_type_property_map,
                                          bool m_body_acceleration
                                             )
{
    // determine which particle this thread works on
    int group_idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (group_idx < group_size)
    {
        //Read particle index
        unsigned int i = d_index_array[group_idx];
        unsigned int fluidneighbors;

        // Access the particle's position, velocity, mass and type
        Scalar3 pi;
        pi.x = d_pos[i].x;
        pi.y = d_pos[i].y;
        pi.z = d_pos[i].z;

        // Read normalization constant of solid particle i
        Scalar norm_constant = 1./d_dpe[i].x;

        // Read acceleration of solid particle i if content is not NaN
        Scalar3 accel_i = make_scalar3(0,0,0);
        if ( d_accel[i].x != d_accel[i].x ||
             d_accel[i].y != d_accel[i].y ||
             d_accel[i].z != d_accel[i].z )
            {
            }
        else
            {
            accel_i.x = d_accel[i].x;
            accel_i.y = d_accel[i].y;
            accel_i.z = d_accel[i].z;
            }

        // Initialize fictitious solid velocity vector
        Scalar3 uf_c0 = make_scalar3(0, 0, 0);

        // Initialize fictitious solid pressure scalar
        Scalar pf_c0= Scalar(0);

        // Initialize hydrostatic pressure contribution
        Scalar3 pd_c0 = make_scalar3(0, 0, 0);

        // Loop over all of the neighbors of this particle
        fluidneighbors = 0;
        const unsigned int myHead = d_head_list[i];
        const unsigned int size = (unsigned int)d_n_neigh[i];
        for (unsigned int j = 0; j < size; j++)
            {
            // Index of neighbor (MEM TRANSFER: 1 scalar)
            unsigned int k = d_nlist[myHead + j];

            // If neighbor particle is solid, continue with next element in loop
            // i.e. interpolations only apply to fluid particles
            if (d_type_property_map[__scalar_as_int(d_pos[k].w)] & SolidFluidTypeBit::SOLID)
                continue;
            else
                fluidneighbors += 1;

            // Access neighbor position
            Scalar3 pj;
            pj.x = d_pos[k].x;
            pj.y = d_pos[k].y;
            pj.z = d_pos[k].z;

            // Compute distance vector (FLOPS: 3)
            Scalar3 dx = pi - pj;

            // Apply periodic boundary conditions (FLOPS: 9)
            dx = box.minImage(dx);

            // Calculate squared distance (FLOPS: 5)
            Scalar rsq = dot(dx, dx);

            // If particle distance is too large, skip this loop
            if ( m_const_slength && rsq > m_rcutsq )
                continue;

            // Access neighbor velocity and mass
            Scalar3 vj;
            vj.x = d_velocity[k].x;
            vj.y = d_velocity[k].y;
            vj.z = d_velocity[k].z;

            // Read particle k pressure
            Scalar Pj = d_dpe[k].y;

            // Calculate absolute and normalized distance
            Scalar r = sqrt(rsq);

            // Evaluate kernel function
            Scalar wij = m_const_slength ? m_skernel.wij(m_ch,r) : m_skernel.wij(Scalar(0.5)*(d_h[i]+d_h[k]),r);

            // Add contribution to solid fictitious velocity
            uf_c0.x += vj.x*wij;
            uf_c0.y += vj.y*wij;
            uf_c0.z += vj.z*wij;

            // Add contribution to solid fictitious pressure
            pf_c0 += Pj*wij;

            // Add contribution to hydrostatic pressure term
            if ( m_body_acceleration )
                {
                pd_c0.x += d_dpe[k].x * dx.x * wij;
                pd_c0.y += d_dpe[k].x * dx.y * wij;
                pd_c0.z += d_dpe[k].x * dx.z * wij;
                }
            } // End neighbor loop

        // Store fictitious solid particle velocity
        if (fluidneighbors > 0)
            {
            d_vf[i].x = 2 * d_velocity[i].x - norm_constant * uf_c0.x;
            d_vf[i].y = 2 * d_velocity[i].y - norm_constant * uf_c0.y;
            d_vf[i].z = 2 * d_velocity[i].z - norm_constant * uf_c0.z;
            }
        else
            {
            d_vf[i].x = 0;
            d_vf[i].y = 0;
            d_vf[i].z = 0;
            }

        // Store fictitious solid particle pressure
        if (fluidneighbors > 0)
            d_dpe[i].y = norm_constant * pf_c0 + dot( g_accel - accel_i , norm_constant * pd_c0);
        else
            d_dpe[i].y = m_eos.getBackgroundPressure();

        // Compute solid densities by inverting equation of state
        d_dpe[i].x = m_eos.Density(d_dpe[i].y);

    }
}


template<SmoothingKernelType KT_, StateEquationType SET_>
cudaError_t gpu_tpf_compute_noslip(Scalar3 *d_vf,
                               Scalar4 *d_pos,
                               Scalar3 *d_dpe,
                               Scalar4 *d_velocity,
                               Scalar3 *d_accel,
                               Scalar *d_h,
                               unsigned int *d_n_neigh,
                               unsigned int *d_nlist,
                               unsigned int *d_head_list,
                               unsigned int *d_index_array,
                               unsigned int group_size,
                               const BoxDim& box,
                               Scalar m_rcutsq,
                               Scalar m_ch,
                               Scalar m_rho0,
                               bool m_const_slength,
                               Scalar3 g_accel,
                               StateEquation<SET_> *m_eos,
                               SmoothingKernel<KT_> *m_skernel,
                               unsigned int *d_type_property_map,
                               bool m_body_acceleration
                                  )
{
    // setup the grid to run the kernel
    unsigned int block_size=256;
    dim3 grid( (group_size / block_size) + 1, 1, 1);
    dim3 threads(block_size, 1, 1);
    
    gpu_tpf_compute_noslip_kernel<<< grid, threads >>>(
        d_vf,
        d_pos,
        d_dpe,
        d_velocity,
        d_accel,
        d_h,
        d_n_neigh,
        d_nlist,
        d_head_list,
        d_index_array,
        group_size,
        box,
        m_rcutsq,
        m_ch,
        m_rho0,
        m_const_slength,
        g_accel,
        *m_eos,
        *m_skernel,
        d_type_property_map,
        m_body_acceleration);
    
    return cudaSuccess;
}



//! Kernel to propagate the positions and velocities, second half of NPT update
template<SmoothingKernelType KT_>
__global__ void gpu_tpf_forcecomputation_kernel(Scalar4 *d_pos,
                                                   Scalar4 *d_velocity,
                                                   Scalar3 *d_dpe,
                                                   Scalar3 *d_vf,
                                                   Scalar *d_h,
                                                   Scalar4 *d_force,
                                                   Scalar4 *d_ratedpe,
                                                   Scalar3 *d_sf,
                                                   unsigned int *d_n_neigh,
                                                   unsigned int *d_nlist,
                                                   unsigned int *d_head_list,
                                                   unsigned int *d_index_array,
                                                   unsigned int group_size,
                                                   const BoxDim box,
                                                   Scalar m_rcutsq,
                                                   Scalar m_ch,
                                                   Scalar m_mu1,
                                                   Scalar m_mu2,
                                                   Scalar m_rho01,
                                                   Scalar m_rho02,
                                                   Scalar m_c1,
                                                   Scalar m_c2,
                                                   bool m_const_slength,
                                                   Scalar m_artificial_viscosity,
                                                   Scalar m_avalpha,
                                                   Scalar m_cmax,
                                                   Scalar m_avbeta,
                                                   Scalar m_ddiff,
                                                   SmoothingKernel<KT_> m_skernel,
                                                   unsigned int *d_type_property_map,
                                                   bool density_sum,
                                                   bool density_cont,
                                                   bool m_compute_solid_forces,
                                                   bool m_density_diffusion
                                               )
{
    // determine which particle this thread works on
    int group_idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (group_idx < group_size)
    {

         // Read particle index
         unsigned int i = d_index_array[group_idx];
         

        // Access the particle's position, velocity, mass and type
        Scalar3 pi;
        pi.x = d_pos[i].x;
        pi.y = d_pos[i].y;
        pi.z = d_pos[i].z;

        Scalar3 vi;
        vi.x = d_velocity[i].x;
        vi.y = d_velocity[i].y;
        vi.z = d_velocity[i].z;
        Scalar mi = d_velocity[i].w;

        // Read particle i pressure
        Scalar Pi = d_dpe[i].y;

        // Read particle i density and volume
        Scalar rhoi = d_dpe[i].x;
        Scalar Vi   = mi / rhoi;

        // Read particle i type, viscosity, speed of sound and rest density
        bool i_isfluid1 = d_type_property_map[__scalar_as_int(d_pos[i].w)] & SolidFluidTypeBit::FLUID1;
        Scalar mui   = i_isfluid1 ? m_mu1 : m_mu2;
        Scalar rho0i = i_isfluid1 ? m_rho01 : m_rho02;
        Scalar ci    = i_isfluid1 ? m_c1 : m_c2;

        // Loop over all of the neighbors of this particle
        const unsigned int myHead = d_head_list[i];
        const unsigned int size = (unsigned int)d_n_neigh[i];
        for (unsigned int j = 0; j < size; j++)
            {
            // Index of neighbor (MEM TRANSFER: 1 scalar)
            unsigned int k = d_nlist[myHead + j];

            // Access neighbor position
            Scalar3 pj;
            pj.x = d_pos[k].x;
            pj.y = d_pos[k].y;
            pj.z = d_pos[k].z;

            // Determine neighbor type
            const unsigned int typej_props = d_type_property_map[__scalar_as_int(d_pos[k].w)];
            bool j_issolid  = typej_props & SolidFluidTypeBit::SOLID;
            bool j_isfluid1 = typej_props & SolidFluidTypeBit::FLUID1;

            // Read particle j viscosity, speed of sound and rest density
            Scalar muj   = j_isfluid1 ? m_mu1 : m_mu2;
            Scalar rho0j = j_isfluid1 ? m_rho01 : m_rho02;
            Scalar cj    = j_isfluid1 ? m_c1 : m_c2;
            // If particle j is solid, set parameters equal to those of particle i
            muj   = j_issolid ? mui : muj;
            rho0j = j_issolid ? rho0i : rho0j;
            cj    = j_issolid ? ci : cj;

            // Compute distance vector (FLOPS: 3)
            Scalar3 dx = pi - pj;

            // Apply periodic boundary conditions (FLOPS: 9)
            dx = box.minImage(dx);

            // Calculate squared distance (FLOPS: 5)
            Scalar rsq = dot(dx, dx);

            // If particle distance is too large, skip this loop
            if ( m_const_slength && rsq > m_rcutsq )
                continue;

            // Access neighbor velocity; depends on fluid or fictitious solid particle
            Scalar3 vj  = make_scalar3(0.0, 0.0, 0.0);
            Scalar mj   = d_velocity[k].w;
            if ( j_issolid )
                {
                vj.x = d_vf[k].x;
                vj.y = d_vf[k].y;
                vj.z = d_vf[k].z;
                }
            else
                {
                vj.x = d_velocity[k].x;
                vj.y = d_velocity[k].y;
                vj.z = d_velocity[k].z;
                }
            Scalar rhoj = d_dpe[k].x;
            Scalar Vj   = mj / rhoj;

            // Read particle k pressure
            Scalar Pj = d_dpe[k].y;

            // Compute velocity difference
            Scalar3 dv = vi - vj;

            // Calculate absolute and normalized distance
            Scalar r = sqrt(rsq);

            // Mean smoothing length and denominator modifier
            Scalar meanh  = m_const_slength ? m_ch : Scalar(0.5)*(d_h[i]+d_h[k]);
            Scalar eps    = Scalar(0.1)*meanh;
            Scalar epssqr = eps*eps;

            // Kernel function derivative evaluation
            Scalar dwdr   = m_skernel.dwijdr(meanh,r);
            Scalar dwdr_r = dwdr/(r+eps);

            // Evaluate inter-particle pressure forces
            //temp0 = -((mi*mj)/(rhoj*rhoi))*(Pi+Pj);
            //temp0 = -Vi*Vj*( Pi + Pj );
            //temp0 = -mi*mj*(Pi+Pj)/(rhoi*rhoj);
            //temp0 = -mi*mj*( Pi/(rhoi*rhoj) + Pj/(rhoj*rhoj) );
            Scalar temp0;
            if ( density_sum )
                temp0 = -(Vi*Vi+Vj*Vj)*((rhoj*Pi+rhoi*Pj)/(rhoi+rhoj));
            else
                temp0 = -mi*mj*(Pi+Pj)/(rhoi*rhoj);


            // Optionally add artificial viscosity
            // Monaghan (1983) J. Comput. Phys. 52 (2) 374–389
            if ( m_artificial_viscosity && !j_issolid )
                {
                Scalar dotdvdx = dot(dv,dx);
                if ( dotdvdx < Scalar(0) )
                    {
                    Scalar muij    = meanh*dotdvdx/(rsq+epssqr);
                    Scalar meanrho = Scalar(0.5)*(rhoi+rhoj);
                    temp0 += mi*mj*(m_avalpha*m_cmax*muij+m_avbeta*muij*muij)/meanrho;
                    }
                }

            // Add contribution to fluid particle
            atomicAdd(&d_force[i].x,temp0*dwdr_r*dx.x);
            atomicAdd(&d_force[i].y,temp0*dwdr_r*dx.y);
            atomicAdd(&d_force[i].z,temp0*dwdr_r*dx.z);

            // Add contribution to solid particle
            if ( j_issolid && m_compute_solid_forces )
                {
                atomicAdd(&d_force[k].x,-(mj/mi)*temp0*dwdr_r*dx.x);
                atomicAdd(&d_force[k].y,-(mj/mi)*temp0*dwdr_r*dx.y);
                atomicAdd(&d_force[k].z,-(mj/mi)*temp0*dwdr_r*dx.z);
                }

            // Evaluate viscous interaction forces
            temp0 = ((Scalar(2)*(mui*muj))/(mui+muj)) * (Vi*Vi+Vj*Vj) * dwdr_r;
            atomicAdd(&d_force[i].x,temp0*dv.x);
            atomicAdd(&d_force[i].y,temp0*dv.y);
            atomicAdd(&d_force[i].z,temp0*dv.z);

            // Add contribution to solid particle
            if ( j_issolid && m_compute_solid_forces )
                {
                atomicAdd(&d_force[k].x,-(mj/mi)*temp0*dv.x);
                atomicAdd(&d_force[k].y,-(mj/mi)*temp0*dv.y);
                atomicAdd(&d_force[k].z,-(mj/mi)*temp0*dv.z);
                }

            // Evaluate rate of change of density if CONTINUITY approach is used
            if ( density_cont )
                {
                if ( j_issolid )
                    {
                    // Use physical advection velocity rather than fictitious velocity here
                    vj.x = d_velocity[k].x;
                    vj.y = d_velocity[k].y;
                    vj.z = d_velocity[k].z;

                    // Recompute velocity difference
                    dv = vi - vj;

                    //Vj = mj / m_rho0;
                    }

                // Compute density rate of change
                d_ratedpe[i].x += rhoi*Vj*dot(dv,dwdr_r*dx);
                //d_ratedpe[i].x += mj*dot(dv,dwdr_r*dx);

                // Add density diffusion if requested
                // Molteni and Colagrossi, Computer Physics Communications 180 (2009) 861–872
                if ( !j_issolid && m_density_diffusion )
                    d_ratedpe[i].x -= (Scalar(2)*m_ddiff*meanh*m_cmax*mj*(rhoi/rhoj-Scalar(1))*dot(dx,dwdr_r*dx))/(rsq+epssqr);
                }

            } // Closing Neighbor Loop
        // Add surface force
        atomicAdd(&d_force[i].x,d_sf[i].x);
        atomicAdd(&d_force[i].y,d_sf[i].y);
        atomicAdd(&d_force[i].z,d_sf[i].z);
    }
}

/*! \param d_vel array of particle velocities
 * \ *param d_accel array of particle accelerations
 * \param d_group_members Device array listing the indicies of the mebers of the group to integrate
 * \param group_size Number of members in the group
 * \param mat_exp_v Matrix exponential for velocity update
 * \param d_net_force Net force on each particle
 * \param deltaT Time to move forward in one whole step
 * 
 * This is just a kernel driver for gpu_SuspendedObject_step_kernel(). See it for more details.
 */
template<SmoothingKernelType KT_>
cudaError_t gpu_tpf_forcecomputation(Scalar4 *d_pos,
                                     Scalar4 *d_velocity,
                                     Scalar3 *d_dpe,
                                     Scalar3 *d_vf,
                                     Scalar *d_h,
                                     Scalar4 *d_force,
                                     Scalar4 *d_ratedpe,
                                     Scalar3 *d_sf,
                                     unsigned int *d_n_neigh,
                                     unsigned int *d_nlist,
                                     unsigned int *d_head_list,
                                     unsigned int *d_group_members,
                                     unsigned int group_size,
                                     const BoxDim& box,
                                     Scalar m_rcutsq,
                                     Scalar m_h,
                                     Scalar m_mu1,
                                     Scalar m_mu2,
                                     Scalar m_rho01,
                                     Scalar m_rho02,
                                     Scalar m_c1,
                                     Scalar m_c2,
                                     bool m_const_slength,
                                     Scalar m_artificial_viscosity,
                                     Scalar m_avalpha,
                                     Scalar m_cmax,
                                     Scalar m_avbeta,
                                     Scalar m_ddiff,
                                     SmoothingKernel<KT_> *m_skernel,
                                     unsigned int *d_type_property_map,
                                     bool density_sum,
                                     bool density_cont,
                                     bool m_compute_solid_forces,
                                     bool m_density_diffusion
)
{
    // setup the grid to run the kernel
    unsigned int block_size=256;
    dim3 grid( (group_size / block_size) + 1, 1, 1);
    dim3 threads(block_size, 1, 1);
    
    // run the kernel
    gpu_tpf_forcecomputation_kernel<<< grid, threads >>>(d_pos,
                                                         d_velocity,
                                                         d_dpe,
                                                         d_vf,
                                                         d_h,
                                                         d_force,
                                                         d_ratedpe,
                                                         d_sf,
                                                         d_n_neigh,
                                                         d_nlist,
                                                         d_head_list,
                                                         d_group_members,
                                                         group_size,
                                                         box,
                                                         m_rcutsq,
                                                         m_h,
                                                         m_mu1,
                                                         m_mu2,
                                                         m_rho01,
                                                         m_rho02,
                                                         m_c1,
                                                         m_c2,
                                                         m_const_slength,
                                                         m_artificial_viscosity,
                                                         m_avalpha,
                                                         m_cmax,
                                                         m_avbeta,
                                                         m_ddiff,
                                                         *m_skernel,
                                                         d_type_property_map,
                                                         density_sum,
                                                         density_cont,
                                                         m_compute_solid_forces,
                                                         m_density_diffusion
    );
    
    return cudaSuccess;
}



//! Kernel to propagate the positions and velocities, second half of NPT update
template<SmoothingKernelType KT_>
__global__ void gpu_tpf_colorgradients_kernel(Scalar4 *d_pos,
                                    Scalar4 *d_velocity,
                                    Scalar3 *d_dpe,
                                    Scalar *d_h,
                                    Scalar3 *d_sn,
                                    Scalar3 *d_fn,
                                    unsigned int *d_n_neigh,
                                    unsigned int *d_nlist,
                                    unsigned int *d_head_list,
                                    unsigned int N,
                                    const BoxDim box,
                                    Scalar m_rcutsq,
                                    bool m_const_slength,
                                    Scalar m_artificial_viscosity,
                                    Scalar m_avalpha,
                                    Scalar m_ch,
                                    Scalar m_avbeta,
                                    Scalar m_ddiff,
                                    SmoothingKernel<KT_> m_skernel,
                                    unsigned int *d_type_property_map,
                                    bool m_compute_solid_forces,
                                    bool m_density_diffusion
                                               )
{
    // determine which particle this thread works on
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (i < N)
    {
        // Access the particle's position, mass and type
        Scalar3 pi;
        pi.x = d_pos[i].x;
        pi.y = d_pos[i].y;
        pi.z = d_pos[i].z;
        Scalar mi = d_velocity[i].w;

        // Read particle i density and volume
        Scalar rhoi = d_dpe[i].x;
        Scalar Vi   = mi / rhoi;

        // Detect particle i type
        const unsigned int typei_props = d_type_property_map[__scalar_as_int(d_pos[i].w)];
        bool i_issolid  = typei_props & SolidFluidTypeBit::SOLID;
        bool i_isfluid1 = typei_props & SolidFluidTypeBit::FLUID1;

        // Loop over all of the neighbors of this particle
        const unsigned int myHead = d_head_list[i];
        const unsigned int size = (unsigned int)d_n_neigh[i];
        for (unsigned int j = 0; j < size; j++)
            {
            // Index of neighbor (MEM TRANSFER: 1 scalar)
            unsigned int k = d_nlist[myHead + j];

            // Access neighbor position
            Scalar3 pj;
            pj.x = d_pos[k].x;
            pj.y = d_pos[k].y;
            pj.z = d_pos[k].z;

            // Determine neighbor type
            const unsigned int typej_props = d_type_property_map[__scalar_as_int(d_pos[k].w)];
            bool j_issolid  = typej_props & SolidFluidTypeBit::SOLID;
            bool j_isfluid1 = typej_props & SolidFluidTypeBit::FLUID1;

            // Skip color gradient computation if both particles belong to same phase
            if ( (i_issolid==true && j_issolid==true) || i_isfluid1 == j_isfluid1 )
                continue;

            // Compute distance vector (FLOPS: 3)
            Scalar3 dx = pi - pj;

            // Apply periodic boundary conditions (FLOPS: 9)
            dx = box.minImage(dx);

            // Calculate squared distance (FLOPS: 5)
            Scalar rsq = dot(dx, dx);

            // If particle distance is too large, skip this loop
            if ( m_const_slength && rsq > m_rcutsq )
                continue;

            // Access neighbor mass and density
            Scalar mj   = d_velocity[k].w;
            Scalar rhoj = d_dpe[k].x;
            Scalar Vj   = mj / rhoj;

            // Calculate absolute and normalized distance
            Scalar r = sqrt(rsq);

            // Mean smoothing length and denominator modifier
            Scalar meanh  = m_const_slength ? m_ch : Scalar(0.5)*(d_h[i]+d_h[k]);
            Scalar eps    = Scalar(0.1)*meanh;

            // Kernel function derivative evaluation
            Scalar dwdr   = m_skernel.dwijdr(meanh,r);
            Scalar dwdr_r = dwdr/(r+eps);

            // Evaluate interpolation contribution
            //Scalar temp0 = (Vi*Vi+Vj*Vj)*(rhoi/(rhoi+rhoj))*(mi/rhoj);
            Scalar temp0 = (Vj*Vj/Vi);

            // If either on of the particle is a solid, interface must be solid-fluid
            if ( i_issolid || j_issolid )
                {
                    d_sn[i].x += temp0*dwdr_r*dx.x;
                    d_sn[i].y += temp0*dwdr_r*dx.y;
                    d_sn[i].z += temp0*dwdr_r*dx.z;
                }
            // Otherwise, interface must be fluid-fluid
            else
                {
                    d_fn[i].x += temp0*dwdr_r*dx.x;
                    d_fn[i].y += temp0*dwdr_r*dx.y;
                    d_fn[i].z += temp0*dwdr_r*dx.z;
                }

            } // Closing Neighbor Loop

        // Make sure that color gradients point from solid to fluid
        // and from fluid 1 to fluid 2 (affects sign of normals)
        if ( i_issolid )
            {
                d_sn[i].x = -d_sn[i].x;
                d_sn[i].y = -d_sn[i].y;
                d_sn[i].z = -d_sn[i].z;
            }
        if ( i_isfluid1 )
            {
                d_fn[i].x = -d_fn[i].x;
                d_fn[i].y = -d_fn[i].y;
                d_fn[i].z = -d_fn[i].z;
            }

    }
}

/*! \param d_vel array of particle velocities
 * \ *param d_accel array of particle accelerations
 * \param d_group_members Device array listing the indicies of the mebers of the group to integrate
 * \param group_size Number of members in the group
 * \param mat_exp_v Matrix exponential for velocity update
 * \param d_net_force Net force on each particle
 * \param deltaT Time to move forward in one whole step
 * 
 * This is just a kernel driver for gpu_SuspendedObject_step_kernel(). See it for more details.
 */
template<SmoothingKernelType KT_>
cudaError_t gpu_tpf_colorgradients(Scalar4 *d_pos,
                                    Scalar4 *d_velocity,
                                    Scalar3 *d_dpe,
                                    Scalar *d_h,
                                    Scalar3 *d_sn,
                                    Scalar3 *d_fn,
                                    unsigned int *d_n_neigh,
                                    unsigned int *d_nlist,
                                    unsigned int *d_head_list,
                                    unsigned int N,
                                    const BoxDim& box,
                                    Scalar m_rcutsq,
                                    bool m_const_slength,
                                    Scalar m_artificial_viscosity,
                                    Scalar m_avalpha,
                                    Scalar m_ch,
                                    Scalar m_avbeta,
                                    Scalar m_ddiff,
                                    SmoothingKernel<KT_> *m_skernel,
                                    unsigned int* d_type_property_map,
                                    bool m_compute_solid_forces,
                                    bool m_density_diffusion
)
{
    // setup the grid to run the kernel
    unsigned int block_size=256;
    dim3 grid( (N / block_size) + 1, 1, 1);
    dim3 threads(block_size, 1, 1);
    
    // run the kernel
    gpu_tpf_colorgradients_kernel<<< grid, threads >>>(d_pos,
                                    d_velocity,
                                    d_dpe,
                                    d_h,
                                    d_sn,
                                    d_fn,
                                    d_n_neigh,
                                    d_nlist,
                                    d_head_list,
                                    N,
                                    box,
                                    m_rcutsq,
                                    m_const_slength,
                                    m_artificial_viscosity,
                                    m_avalpha,
                                    m_ch,
                                    m_avbeta,
                                    m_ddiff,
                                    *m_skernel,
                                    d_type_property_map,
                                    m_compute_solid_forces,
                                    m_density_diffusion
    );
    
    return cudaSuccess;
}



//! Kernel to propagate the positions and velocities, second half of NPT update
template<SmoothingKernelType KT_>
__global__ void gpu_tpf_compute_surfaceforce_kernel(Scalar4 *d_pos,
                                         Scalar4 *d_velocity,
                                         Scalar3 *d_dpe,
                                         Scalar3 *d_sf,
                                         Scalar *d_h,
                                         Scalar3 *d_sn,
                                         Scalar3 *d_fn,
                                         unsigned int *d_n_neigh,
                                         unsigned int *d_nlist,
                                         unsigned int *d_head_list,
                                         unsigned int *d_index_array,
                                         unsigned int group_size,
                                         const BoxDim box,
                                         Scalar m_rcutsq,
                                         bool m_const_slength,
                                         Scalar m_sigma01,
                                         Scalar m_sigma02,
                                         Scalar m_sigma12,
                                         Scalar m_ch,
                                         SmoothingKernel<KT_> m_skernel,
                                         unsigned int *d_type_property_map,
                                         bool m_compute_solid_forces,
                                         bool m_density_diffusion
                                               )
{
    // determine which particle this thread works on
    int group_idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (group_idx < group_size)
    {
        // Read particle index
        unsigned int i = d_index_array[group_idx];

        // Access the particle's position and type
        Scalar3 pi;
        pi.x = d_pos[i].x;
        pi.y = d_pos[i].y;
        pi.z = d_pos[i].z;
        bool i_isfluid1 = d_type_property_map[__scalar_as_int(d_pos[i].w)] & SolidFluidTypeBit::FLUID1;

        // Check if there is any fluid particle near the current particle, if not continue
        // This makes sure that only particle near a fluid interface experience an interfacial force.
        // In other words, fluid particles only near solid interfaces are omitted.
        bool nearfluidinterface = false;
        // Loop over all of the neighbors of this particle
        unsigned int myHead = d_head_list[i];
        unsigned int size = (unsigned int)d_n_neigh[i];
        for (unsigned int j = 0; j < size; j++)
            {
            // Index of neighbor (MEM TRANSFER: 1 scalar)
            unsigned int k = d_nlist[myHead + j];
            const unsigned int typej_props = d_type_property_map[__scalar_as_int(d_pos[k].w)];
            bool j_issolid  = typej_props & SolidFluidTypeBit::SOLID;
            bool j_isfluid1 = typej_props & SolidFluidTypeBit::FLUID1;
            if ( !(j_issolid) && i_isfluid1 != j_isfluid1 )
                {
                    nearfluidinterface = true;
                    break;
                }

            }
        if ( nearfluidinterface )
            {

            // Access the particle's mass
            Scalar mi = d_velocity[i].w;

            // Read particle i density and volume
            Scalar rhoi = d_dpe[i].x;
            Scalar Vi   = mi / rhoi;

            // Read particle i color gradients
            Scalar3 sni;
            sni.x = d_sn[i].x;
            sni.y = d_sn[i].y;
            sni.z = d_sn[i].z;
            Scalar normsni = sqrt(dot(sni,sni));
            Scalar3 fni;
            fni.x = d_fn[i].x;
            fni.y = d_fn[i].y;
            fni.z = d_fn[i].z;
            Scalar normfni = sqrt(dot(fni,fni));

            // Evaluate particle i interfacial stress tensor
            Scalar istress[6] = {0};
            // Fluid phase 1 - Fluid phase 2 interface
            if ( m_sigma12 > 0 && normfni > 0 )
                {
                Scalar temp0 = (m_sigma12/normfni);
                istress[0] +=  temp0*(normfni*normfni-fni.x*fni.x); // xx
                istress[1] +=  temp0*(normfni*normfni-fni.y*fni.y); // yy
                istress[2] +=  temp0*(normfni*normfni-fni.z*fni.z); // zz
                istress[3] += -temp0*(fni.x*fni.y);                 // xy yx
                istress[4] += -temp0*(fni.x*fni.z);                 // xz zx
                istress[5] += -temp0*(fni.y*fni.z);                 // yz zy
                }
            // Fluid phase 1 - Solid interface
            if ( i_isfluid1 && m_sigma01 > 0 && normsni > 0 )
                {
                Scalar temp0 = (m_sigma01/normsni);
                istress[0] +=  temp0*(normsni*normsni-sni.x*sni.x); // xx
                istress[1] +=  temp0*(normsni*normsni-sni.y*sni.y); // yy
                istress[2] +=  temp0*(normsni*normsni-sni.z*sni.z); // zz
                istress[3] += -temp0*(sni.x*sni.y);                 // xy yx
                istress[4] += -temp0*(sni.x*sni.z);                 // xz zx
                istress[5] += -temp0*(sni.y*sni.z);                 // yz zy
                }
            // Fluid phase 2 - Solid interface
            if ( !i_isfluid1 && m_sigma02 > 0 && normsni > 0 )
                {
                Scalar temp0 = (m_sigma02/normsni);
                istress[0] +=  temp0*(normsni*normsni-sni.x*sni.x); // xx
                istress[1] +=  temp0*(normsni*normsni-sni.y*sni.y); // yy
                istress[2] +=  temp0*(normsni*normsni-sni.z*sni.z); // zz
                istress[3] += -temp0*(sni.x*sni.y);                 // xy yx
                istress[4] += -temp0*(sni.x*sni.z);                 // xz zx
                istress[5] += -temp0*(sni.y*sni.z);                 // yz zy
                }

            // Loop over all of the neighbors of this particle
            myHead = d_head_list[i];
            size = (unsigned int)d_n_neigh[i];
            for (unsigned int j = 0; j < size; j++)
                {

                // Index of neighbor (MEM TRANSFER: 1 scalar)
                unsigned int k = d_nlist[myHead + j];


                // Access neighbor position
                Scalar3 pj;
                pj.x = d_pos[k].x;
                pj.y = d_pos[k].y;
                pj.z = d_pos[k].z;

                // Determine neighbor type
                const unsigned int typej_props = d_type_property_map[__scalar_as_int(d_pos[k].w)];
                bool j_issolid  = typej_props & SolidFluidTypeBit::SOLID;
                bool j_isfluid1 = typej_props & SolidFluidTypeBit::FLUID1;

                // Compute normalized color gradients
                Scalar3 snj;
                snj.x = d_sn[k].x;
                snj.y = d_sn[k].y;
                snj.z = d_sn[k].z;
                Scalar normsnj = sqrt(dot(snj,snj));
                Scalar3 fnj;
                fnj.x = d_fn[k].x;
                fnj.y = d_fn[k].y;
                fnj.z = d_fn[k].z;
                Scalar normfnj = sqrt(dot(fnj,fnj));

                // Compute distance vector (FLOPS: 3)
                Scalar3 dx = pi - pj;

                // Apply periodic boundary conditions (FLOPS: 9)
                dx = box.minImage(dx);

                // Calculate squared distance (FLOPS: 5)
                Scalar rsq = dot(dx, dx);

                // If particle distance is too large, skip this loop
                if ( m_const_slength && rsq > m_rcutsq )
                    continue;

                // Calculate absolute and normalized distance
                Scalar r = sqrt(rsq);

                // Access neighbor mass and density
                Scalar mj   = d_velocity[k].w;
                Scalar rhoj = d_dpe[k].x;
                Scalar Vj   = mj / rhoj;

                // Mean smoothing length and denominator modifier
                Scalar meanh  = m_const_slength ? m_ch : Scalar(0.5)*(d_h[i]+d_h[k]);
                Scalar eps    = Scalar(0.1)*meanh;

                // Kernel function derivative evaluation
                Scalar dwdr   = m_skernel.dwijdr(meanh,r);
                Scalar dwdr_r = dwdr/(r+eps);

                // Evaluate particle i interfacial stress tensor
                Scalar jstress[6] = {0};
                // Fluid phase 1 - Fluid phase 2 interface
                if ( !(j_issolid) && m_sigma12 > 0 && normfnj > 0 )
                    {
                    Scalar temp0 = (m_sigma12/normfnj);
                    jstress[0] +=  temp0*(normfnj*normfnj-fnj.x*fnj.x); // xx
                    jstress[1] +=  temp0*(normfnj*normfnj-fnj.y*fnj.y); // yy
                    jstress[2] +=  temp0*(normfnj*normfnj-fnj.z*fnj.z); // zz
                    jstress[3] += -temp0*(fnj.x*fnj.y);                 // xy yx
                    jstress[4] += -temp0*(fnj.x*fnj.z);                 // xz zx
                    jstress[5] += -temp0*(fnj.y*fnj.z);                 // yz zy
                    }
                // Fluid phase 1 - Solid interface
                if ( j_isfluid1 && m_sigma01 > 0 && normsnj > 0 )
                    {
                    Scalar temp0 = (m_sigma01/normsnj);
                    jstress[0] +=  temp0*(normsnj*normsnj-snj.x*snj.x); // xx
                    jstress[1] +=  temp0*(normsnj*normsnj-snj.y*snj.y); // yy
                    jstress[2] +=  temp0*(normsnj*normsnj-snj.z*snj.z); // zz
                    jstress[3] += -temp0*(snj.x*snj.y);                 // xy yx
                    jstress[4] += -temp0*(snj.x*snj.z);                 // xz zx
                    jstress[5] += -temp0*(snj.y*snj.z);                 // yz zy
                    }
                // Fluid phase 2 - Solid interface
                if ( !j_isfluid1 && m_sigma02 > 0 && normsnj > 0 )
                    {
                    Scalar temp0 = (m_sigma02/normsnj);
                    jstress[0] +=  temp0*(normsnj*normsnj-snj.x*snj.x); // xx
                    jstress[1] +=  temp0*(normsnj*normsnj-snj.y*snj.y); // yy
                    jstress[2] +=  temp0*(normsnj*normsnj-snj.z*snj.z); // zz
                    jstress[3] += -temp0*(snj.x*snj.y);                 // xy yx
                    jstress[4] += -temp0*(snj.x*snj.z);                 // xz zx
                    jstress[5] += -temp0*(snj.y*snj.z);                 // yz zy
                    }

                // Add contribution to surface force
                d_sf[i].x += dwdr_r*dx.x*(Vi*Vi*istress[0]+Vj*Vj*jstress[0])+
                                dwdr_r*dx.y*(Vi*Vi*istress[3]+Vj*Vj*jstress[3])+
                                dwdr_r*dx.z*(Vi*Vi*istress[4]+Vj*Vj*jstress[4]);
                d_sf[i].y += dwdr_r*dx.x*(Vi*Vi*istress[3]+Vj*Vj*jstress[3])+
                                dwdr_r*dx.y*(Vi*Vi*istress[1]+Vj*Vj*jstress[1])+
                                dwdr_r*dx.z*(Vi*Vi*istress[5]+Vj*Vj*jstress[5]);
                d_sf[i].z += dwdr_r*dx.x*(Vi*Vi*istress[4]+Vj*Vj*jstress[4])+
                                dwdr_r*dx.y*(Vi*Vi*istress[5]+Vj*Vj*jstress[5])+
                                dwdr_r*dx.z*(Vi*Vi*istress[2]+Vj*Vj*jstress[2]);

                } // End of neighbor loop

            // Set component normal to solid surface at solid interface to zero
            }

    }
}

/*! \param d_vel array of particle velocities
 * \ *param d_accel array of particle accelerations
 * \param d_group_members Device array listing the indicies of the mebers of the group to integrate
 * \param group_size Number of members in the group
 * \param mat_exp_v Matrix exponential for velocity update
 * \param d_net_force Net force on each particle
 * \param deltaT Time to move forward in one whole step
 * 
 * This is just a kernel driver for gpu_SuspendedObject_step_kernel(). See it for more details.
 */
template<SmoothingKernelType KT_>
cudaError_t gpu_tpf_compute_surfaceforce(Scalar4 *d_pos,
                                         Scalar4 *d_velocity,
                                         Scalar3 *d_dpe,
                                         Scalar3 *d_sf,
                                         Scalar *d_h,
                                         Scalar3 *d_sn,
                                         Scalar3 *d_fn,
                                         unsigned int *d_n_neigh,
                                         unsigned int *d_nlist,
                                         unsigned int *d_head_list,
                                         unsigned int *d_group_members,
                                         unsigned int group_size,
                                         const BoxDim& box,
                                         Scalar m_rcutsq,
                                         bool m_const_slength,
                                         Scalar m_sigma01,
                                         Scalar m_sigma02,
                                         Scalar m_sigma12,
                                         Scalar m_ch,
                                         SmoothingKernel<KT_> *m_skernel,
                                         unsigned int *d_type_property_map,
                                         bool m_compute_solid_forces,
                                         bool m_density_diffusion
)
{
    // setup the grid to run the kernel
    unsigned int block_size=256;
    dim3 grid( (group_size / block_size) + 1, 1, 1);
    dim3 threads(block_size, 1, 1);
    
    // run the kernel
    gpu_tpf_compute_surfaceforce_kernel<<< grid, threads >>>(d_pos,
                                         d_velocity,
                                         d_dpe,
                                         d_sf,
                                         d_h,
                                         d_sn,
                                         d_fn,
                                         d_n_neigh,
                                         d_nlist,
                                         d_head_list,
                                         d_group_members,
                                         group_size,
                                         box,
                                         m_rcutsq,
                                         m_const_slength,
                                         m_sigma01,
                                         m_sigma02,
                                         m_sigma12,
                                         m_ch,
                                         *m_skernel,
                                         d_type_property_map,
                                         m_compute_solid_forces,
                                         m_density_diffusion
);
    
    return cudaSuccess;
}



/*
 * explicit template instantiation
 */
#define INST_gpu_tpf_noslip(r,seq) template cudaError_t gpu_tpf_compute_noslip(Scalar3 *d_vf,Scalar4 *d_pos, Scalar3 *d_dpe,Scalar4 *d_velocity, \
    Scalar3 *d_accel,Scalar *d_h,unsigned int *d_n_neigh, unsigned int *d_nlist,unsigned int *d_head_list, unsigned int *d_index_array, \
    unsigned int group_size, const BoxDim& box,Scalar m_rcutsq,Scalar m_ch, Scalar m_rho0, bool m_const_slength, Scalar3 g_accel, \
    StateEquation<BOOST_PP_SEQ_ELEM(1,seq)> *m_eos, SmoothingKernel<BOOST_PP_SEQ_ELEM(0,seq)> *m_skernel, unsigned int *d_type_property_map, \
    bool m_body_acceleration);

#define INST_gpu_tpf_ndensity(r,data,k) template cudaError_t gpu_tpf_compute_ndensity(Scalar4 *d_pos, Scalar3 *d_dpe, Scalar4 *d_vel, Scalar *d_h, \
   unsigned int *d_n_neigh, unsigned int *d_nlist, unsigned int *d_head_list, unsigned int N, const BoxDim& box, Scalar m_ch, Scalar m_rcutsq, \
   Scalar m_rho0, bool m_const_slength, SmoothingKernel<k> *m_skernel,unsigned int *d_type_property_map, \
   bool densinty_sum);

#define INST_gpu_tpf_ndensitynorm(r,data,k) template cudaError_t gpu_tpf_compute_ndensityrenormalization(Scalar4 *d_pos, Scalar3 *d_dpe, Scalar4 *d_vel, \
   Scalar *d_h, unsigned int *d_n_neigh, unsigned int *d_nlist, unsigned int *d_head_list, unsigned int *d_index_array, unsigned int group_size, \
   const BoxDim& box, Scalar m_rcutsq, Scalar m_ch, bool m_const_slength, SmoothingKernel<k> *m_skernel, unsigned int *d_type_property_map);

#define INST_gpu_tpf_force(r,data,k) template cudaError_t gpu_tpf_forcecomputation(Scalar4 *d_pos, Scalar4 *d_velocity, Scalar3 *d_dpe, Scalar3 *d_vf, \
   Scalar *d_h, Scalar4 *d_force, Scalar4 *d_ratedpe, Scalar3 *d_sf, unsigned int *d_n_neigh, unsigned int *d_nlist, unsigned int *d_head_list, \
   unsigned int *d_group_members, unsigned int group_size, const BoxDim& box, Scalar m_rcutsq, Scalar m_h, Scalar m_mu1, Scalar m_mu2, Scalar m_rho01, \
   Scalar m_rho02, Scalar m_c1, Scalar m_c2, bool m_const_slength, Scalar m_artificial_viscosity, Scalar m_avalpha, Scalar m_cmax, Scalar m_avbeta, \
   Scalar m_ddiff, SmoothingKernel<k> *m_skernel, unsigned int *d_type_property_map, bool density_sum, \
   bool density_cont, bool m_compute_solid_forces, bool m_density_diffusion);

#define INST_gpu_tpf_grad(r,data,k) template cudaError_t gpu_tpf_colorgradients(Scalar4 *d_pos, Scalar4 *d_velocity, Scalar3 *d_dpe, Scalar *d_h, \
   Scalar3 *d_sn, Scalar3 *d_fn, unsigned int *d_n_neigh, unsigned int *d_nlist, unsigned int *d_head_list, unsigned int N, const BoxDim& box, \
   Scalar m_rcutsq, bool m_const_slength, Scalar m_artificial_viscosity, Scalar m_avalpha, Scalar m_ch, Scalar m_avbeta, Scalar m_ddiff, \
   SmoothingKernel<k> *m_skernel, unsigned int *d_type_property_map, bool m_compute_solid_forces, \
   bool m_density_diffusion);

#define INST_gpu_tpf_sfaceforce(r,data,k) template cudaError_t gpu_tpf_compute_surfaceforce(Scalar4 *d_pos, Scalar4 *d_velocity, Scalar3 *d_dpe, \
   Scalar3 *d_sf, Scalar *d_h, Scalar3 *d_sn, Scalar3 *d_fn, unsigned int *d_n_neigh, unsigned int *d_nlist, unsigned int *d_head_list, \
   unsigned int *d_group_members, unsigned int group_size, const BoxDim& box, Scalar m_rcutsq, bool m_const_slength, Scalar m_sigma01, Scalar m_sigma02, \
   Scalar m_sigma12, Scalar m_ch, SmoothingKernel<k> *m_skernel, unsigned int *d_type_property_map, \
   bool m_compute_solid_forces, bool m_density_diffusion);

BOOST_PP_SEQ_FOR_EACH_PRODUCT(INST_gpu_tpf_noslip, (KERNELTYPES)(SEQTYPES));
BOOST_PP_SEQ_FOR_EACH(INST_gpu_tpf_ndensity, ~, KERNELTYPES);
BOOST_PP_SEQ_FOR_EACH(INST_gpu_tpf_ndensitynorm, ~, KERNELTYPES);
BOOST_PP_SEQ_FOR_EACH(INST_gpu_tpf_force, ~, KERNELTYPES);
BOOST_PP_SEQ_FOR_EACH(INST_gpu_tpf_grad, ~, KERNELTYPES);
BOOST_PP_SEQ_FOR_EACH(INST_gpu_tpf_sfaceforce, ~, KERNELTYPES);

template cudaError_t gpu_tpf_compute_pressure(Scalar3 *d_dpe, unsigned int *d_index_array,
  unsigned int group_size, StateEquation<tait> *m_eos);

template cudaError_t gpu_tpf_compute_pressure(Scalar3 *d_dpe, unsigned int *d_index_array,
  unsigned int group_size, StateEquation<linear> *m_eos);

#undef INST_gpu_tpf_noslip
#undef INST_gpu_tpf_ndensity
#undef INST_gpu_tpf_ndensitynorm
#undef INST_gpu_tpf_force
#undef INST_gpu_tpf_grad
#undef INST_gpu_tpf_sfaceforce
