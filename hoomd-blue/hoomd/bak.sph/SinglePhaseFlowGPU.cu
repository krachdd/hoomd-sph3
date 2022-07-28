// Copyright (c) 2009-2016 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

// Maintainer: schirwon

#include "SinglePhaseFlowGPU.cuh"
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
__global__ void gpu_spf_compute_pressure_kernel(
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
cudaError_t gpu_spf_compute_pressure(
                                 Scalar3 *d_dpe,
                                 unsigned int *d_index_array,
                                 unsigned int group_size,
                                 StateEquation<SET_> *m_eos)
{
    // setup the grid to run the kernel
    unsigned int block_size=256;
    dim3 grid( (group_size / block_size) + 1, 1, 1);
    dim3 threads(block_size, 1, 1);
    gpu_spf_compute_pressure_kernel<<< grid, threads >>>(
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
__global__ void gpu_spf_compute_ndensity_kernel(Scalar4 *d_pos,
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
        unsigned int typei = __scalar_as_int(d_pos[i].w);
        bool i_issolid = d_type_property_map[typei] & SolidFluidTypeBit::SOLID;
        bool i_isfluid = d_type_property_map[typei] & SolidFluidTypeBit::FLUID;

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
            // Skip neighbor loop if this solid solid particle does not have fluid neighbors
            bool solid_w_fluid_neigh = false;
            if ( i_issolid )
                {
                myHead = d_head_list[i];
                size = (unsigned int)d_n_neigh[i];
                for (unsigned int j = 0; j < size; j++)
                    {
                        unsigned int k = d_nlist[myHead + j];
                        solid_w_fluid_neigh = d_type_property_map[__scalar_as_int(d_pos[k].w)] & SolidFluidTypeBit::FLUID;

                        if ( solid_w_fluid_neigh == true )
                            {
                              break;
                            }
                        }
                }
                bool skip(false);
                if ( i_issolid && !(solid_w_fluid_neigh) )
                {
                  d_dpe[i].x = m_rho0;
                  skip=true;
                }
            if(!skip)
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
cudaError_t gpu_spf_compute_ndensity(Scalar4 *d_pos,
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
    
    gpu_spf_compute_ndensity_kernel<<< grid, threads >>>(
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
__global__ void gpu_spf_compute_ndensityrenormalization_kernel(Scalar4 *d_pos,
                                                           Scalar3 *d_dpe,
                                                           Scalar3 *d_dpe_old,
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
                                                           SmoothingKernel<KT_> m_skernel
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

        // Loop over all of the neighbors of this particle
        // and compute normalization constant normwij = \sum_j wij*Vj
        myHead = d_head_list[i];
        size = (unsigned int)d_n_neigh[i];

        //bool skip(false);
        for (unsigned int j = 0; j < size; j++)
        {
            // Index of neighbor
            unsigned int k = d_nlist[myHead + j];
            
            
            
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
                Scalar Vj =  d_vel[k].w / d_dpe_old[k].x ;
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
cudaError_t gpu_spf_compute_ndensityrenormalization(Scalar4 *d_pos,
                                                Scalar3 *d_dpe,
                                                Scalar3 *d_dpe_old,
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
                                                SmoothingKernel<KT_> *m_skernel
)
{
    // setup the grid to run the kernel
    unsigned int block_size=256;
    dim3 grid( (group_size / block_size) + 1, 1, 1);
    dim3 threads(block_size, 1, 1);

    gpu_spf_compute_ndensityrenormalization_kernel<<< grid, threads >>>(
        d_pos,
        d_dpe,
        d_dpe_old,
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
        *m_skernel);
    
    return cudaSuccess;
}

/*! Helper function to compute fictitious solid particle properties (pressures and velocities)
 * \pre Ghost particle number densities (i.e. charge array) must be up-to-date
 * \pre Normalization constant \sum_fluids Wij must be computed and stored in orientation array
 * \post Fictitious particle properties are computed and stored in orentiation array
 */
template<SmoothingKernelType KT_, StateEquationType SET_>
__global__ void gpu_spf_compute_noslip_kernel(Scalar3 *d_vf,
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
        unsigned int myHead, size;
        //Access the particle's position, velocity, mass and type
        Scalar3 pi;
        pi.x = d_pos[i].x;
        pi.y = d_pos[i].y;
        pi.z = d_pos[i].z;

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
        // Count fluid neighbors before setting solid particle properties
        unsigned int fluidneighbors = 0;

        // Skip neighbor loop if this solid particle does not have fluid neighbors
        bool solid_w_fluid_neigh = false;
        myHead = d_head_list[i];
        size = (unsigned int)d_n_neigh[i];

        for (unsigned int j = 0; j < size; j++)
        {
            unsigned int k = d_nlist[myHead + j];
            solid_w_fluid_neigh = d_type_property_map[__scalar_as_int(d_pos[k].w)] & SolidFluidTypeBit::FLUID;
            

            if ( solid_w_fluid_neigh == true )
            {
                break;
            }
        }
        bool skip(false);
         if ( !(solid_w_fluid_neigh) )
            {
            // Solid particles which do not have fluid neighbors are marked
            // using energy=1 so that they can be deleted during simulation
            d_dpe[i].z = Scalar(1);
            // Set fictitious solid velocity to zero
            d_vf[i].x = 0;
            d_vf[i].y = 0;
            d_vf[i].z = 0;
            // If no fluid neighbors are present,
            // Set pressure to background pressure
            d_dpe[i].y = m_eos.getBackgroundPressure();
            // Density to rest density
            d_dpe[i].x = m_rho0;
            skip=true;
            }
        if(!skip)
        {

            myHead = d_head_list[i];
            size = (unsigned int)d_n_neigh[i];
            for (unsigned int j = 0; j < size; j++)
                {
                // Index of neighbor (MEM TRANSFER: 1 scalar)
                unsigned int k = d_nlist[myHead + j];



                // If neighbor particle is solid, continue with next element in loop
                // i.e. interpolations only apply to fluid particles
                bool solid_skip = d_type_property_map[__scalar_as_int(d_pos[k].w)] & SolidFluidTypeBit::SOLID;

                if ( solid_skip )
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
                if (! (m_const_slength && rsq > m_rcutsq) )
                {

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
                }
            } // End neighbor loop

            // Store fictitious solid particle velocity
            if (fluidneighbors > 0 && d_dpe[i].x > 0 )
                {
                //  Compute zeroth order normalization constant
                Scalar norm_constant = 1./d_dpe[i].x;
                // Set fictitious velocity
                d_vf[i].x = 2 * d_velocity[i].x - norm_constant * uf_c0.x;
                d_vf[i].y = 2 * d_velocity[i].y - norm_constant * uf_c0.y;
                d_vf[i].z = 2 * d_velocity[i].z - norm_constant * uf_c0.z;

                d_dpe[i].y = norm_constant * pf_c0 + dot( g_accel - accel_i , norm_constant * pd_c0);
                // Compute solid densities by inverting equation of state
                d_dpe[i].x = m_eos.Density(d_dpe[i].y);
                }
            else
                {
                // Set fictitious solid velocity to zero
                d_vf[i].x = 0;
                d_vf[i].y = 0;
                d_vf[i].z = 0;

                // If no fluid neighbors are present,
                // Set pressure to background pressure
                d_dpe[i].y = m_eos.getBackgroundPressure();
                // Density to rest density
                d_dpe[i].x = m_rho0;
                }

        }
    }
}


template<SmoothingKernelType KT_, StateEquationType SET_>
cudaError_t gpu_spf_compute_noslip(Scalar3 *d_vf,
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

    gpu_spf_compute_noslip_kernel<<< grid, threads >>>(
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



/*! Helper function to set communication flags and update ghosts orientations
 * \param timestep The time step
 * \post Ghost particle orientations are up-to-date
 */

cudaError_t gpu_spf_compute_properties(
    Scalar4 *d_vel,
    unsigned int *d_group_members,
    unsigned int group_size, Scalar3 *sum)
{
    return gpu_sum_xyz(d_vel, group_size, sum, d_group_members);
}


//! Kernel to propagate the positions and velocities, second half of NPT update
template<SmoothingKernelType KT_>
__global__ void gpu_spf_forcecomputation_kernel(Scalar4 *d_pos,
                                                Scalar4 *d_velocity,
                                                Scalar3 *d_dpe,
                                                Scalar3 *d_vf,
                                                Scalar *d_h,
                                                Scalar4 *d_force,
                                                Scalar4 *d_ratedpe,   
                                                unsigned int *d_n_neigh,
                                                unsigned int *d_nlist,
                                                unsigned int *d_head_list,
                                                unsigned int *d_index_array,
                                                unsigned int group_size,
                                                const BoxDim box,
                                                Scalar m_rcutsq,
                                                Scalar m_ch,
                                                Scalar m_mu,
                                                bool m_const_slength,
                                                Scalar m_artificial_viscosity,
                                                Scalar m_avalpha,
                                                Scalar m_c,
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
            bool issolid = d_type_property_map[__scalar_as_int(d_pos[k].w)] & SolidFluidTypeBit::SOLID;

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
            Scalar3 vj; //  = make_scalar3(0.0, 0.0, 0.0);
            Scalar mj   = d_velocity[k].w;
            if ( issolid )
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
            if ( m_artificial_viscosity && !issolid )
                {
                Scalar dotdvdx = dot(dv,dx);
                if ( dotdvdx < Scalar(0) )
                    {
                    Scalar muij    = meanh*dotdvdx/(rsq+epssqr);
                    Scalar meanrho = Scalar(0.5)*(rhoi+rhoj);
                    temp0 += mi*mj*(m_avalpha*m_c*muij+m_avbeta*muij*muij)/meanrho;
                    }
                }

            d_force[i].x += temp0*dwdr_r*dx.x;
            d_force[i].y += temp0*dwdr_r*dx.y;
            d_force[i].z += temp0*dwdr_r*dx.z;

            // Add contribution to solid particle
            if ( issolid && m_compute_solid_forces )
                {
                d_force[k].x -= (mj/mi)*temp0*dwdr_r*dx.x;
                d_force[k].y -= (mj/mi)*temp0*dwdr_r*dx.y;
                d_force[k].z -= (mj/mi)*temp0*dwdr_r*dx.z;
                }

            // Evaluate viscous interaction forces
            temp0 = m_mu * (Vi*Vi+Vj*Vj) * dwdr_r;
            d_force[i].x  += temp0*dv.x;
            d_force[i].y  += temp0*dv.y;
            d_force[i].z  += temp0*dv.z;

            // Add contribution to solid particle
            if ( issolid && m_compute_solid_forces )
                {
                d_force[k].x -= (mj/mi)*temp0*dv.x;
                d_force[k].y -= (mj/mi)*temp0*dv.y;
                d_force[k].z -= (mj/mi)*temp0*dv.z;
                }

            // Evaluate rate of change of density if CONTINUITY approach is used
            if ( density_cont )
                {
                if ( issolid )
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
                if ( !issolid && m_density_diffusion )
                    d_ratedpe[i].x -= (Scalar(2)*m_ddiff*meanh*m_c*mj*(rhoi/rhoj-Scalar(1))*dot(dx,dwdr_r*dx))/(rsq+epssqr);
                }

            } // Closing Neighbor Loop

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
cudaError_t gpu_spf_forcecomputation(Scalar4 *d_pos,
                                     Scalar4 *d_velocity,
                                     Scalar3 *d_dpe,
                                     Scalar3 *d_vf,
                                     Scalar *d_h,
                                     Scalar4 *d_force,
                                     Scalar4 *d_ratedpe,   
                                     unsigned int *d_n_neigh,
                                     unsigned int *d_nlist,
                                     unsigned int *d_head_list,
                                     unsigned int *d_group_members,
                                     unsigned int group_size,
                                     const BoxDim& box,
                                     Scalar m_rcutsq,
                                     Scalar m_h,
                                     Scalar m_mu,
                                     bool m_const_slength,
                                     Scalar m_artificial_viscosity,
                                     Scalar m_avalpha,
                                     Scalar m_c,
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
    gpu_spf_forcecomputation_kernel<<< grid, threads >>>(d_pos,
                                                         d_velocity,
                                                         d_dpe,
                                                         d_vf,
                                                         d_h,
                                                         d_force,
                                                         d_ratedpe,
                                                         d_n_neigh,
                                                         d_nlist,
                                                         d_head_list,
                                                         d_group_members,
                                                         group_size,
                                                         box,
                                                         m_rcutsq,
                                                         m_h,
                                                         m_mu,
                                                         m_const_slength,
                                                         m_artificial_viscosity,
                                                         m_avalpha,
                                                         m_c,
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


/*
 * explicit template instantiation
 */


#define INST_gpu_spf_noslip(r,seq) template cudaError_t gpu_spf_compute_noslip(Scalar3 *d_vf,Scalar4 *d_pos, Scalar3 *d_dpe,Scalar4 *d_velocity, \
    Scalar3 *d_accel,Scalar *d_h,unsigned int *d_n_neigh, unsigned int *d_nlist,unsigned int *d_head_list, unsigned int *d_index_array, \
    unsigned int group_size, const BoxDim& box,Scalar m_rcutsq,Scalar m_ch,Scalar m_rho0, bool m_const_slength, Scalar3 g_accel, \
    StateEquation<BOOST_PP_SEQ_ELEM(1,seq)> *m_eos, SmoothingKernel<BOOST_PP_SEQ_ELEM(0,seq)> *m_skernel, unsigned int *d_type_property_map, \
    bool m_body_acceleration);

#define INST_gpu_spf_ndensity(r,data,k) template cudaError_t gpu_spf_compute_ndensity(Scalar4 *d_pos, Scalar3 *d_dpe, Scalar4 *d_vel, Scalar *d_h, \
   unsigned int *d_n_neigh, unsigned int *d_nlist, unsigned int *d_head_list, unsigned int N, const BoxDim& box, Scalar m_ch, Scalar m_rcutsq, Scalar m_rho0, \
   bool m_const_slength, SmoothingKernel<k> *m_skernel, unsigned int *d_type_property_map, bool densinty_sum);

#define INST_gpu_spf_ndensitynorm(r,data,k) template cudaError_t gpu_spf_compute_ndensityrenormalization(Scalar4 *d_pos, Scalar3 *d_dpe, Scalar3 *d_dpe_alt,\
Scalar4 *d_vel, Scalar *d_h, unsigned int *d_n_neigh, unsigned int *d_nlist, unsigned int *d_head_list, unsigned int *d_index_array, unsigned int group_size, \
   const BoxDim& box, Scalar m_rcutsq, Scalar m_ch, bool m_const_slength, SmoothingKernel<k> *m_skernel);

#define INST_gpu_spf_force(r,data,k) template cudaError_t gpu_spf_forcecomputation(Scalar4 *d_pos, Scalar4 *d_velocity, Scalar3 *d_dpe, Scalar3 *d_vf, \
   Scalar *d_h, Scalar4 *d_force, Scalar4 *d_ratedpe, unsigned int *d_n_neigh, unsigned int *d_nlist, unsigned int *d_head_list, unsigned int *d_group_members, \
   unsigned int group_size, const BoxDim& box, Scalar m_rcutsq, Scalar m_h, Scalar m_mu, bool m_const_slength, Scalar m_artificial_viscosity, Scalar m_avalpha, \
   Scalar m_c, Scalar m_avbeta, Scalar m_ddiff, SmoothingKernel<k> *m_skernel, unsigned int* d_type_property_map, bool density_sum, bool densinty_cont, \
   bool m_compute_solid_forces, bool m_density_diffusion);

BOOST_PP_SEQ_FOR_EACH_PRODUCT(INST_gpu_spf_noslip, (KERNELTYPES)(SEQTYPES));
BOOST_PP_SEQ_FOR_EACH(INST_gpu_spf_ndensity, ~, KERNELTYPES);
BOOST_PP_SEQ_FOR_EACH(INST_gpu_spf_ndensitynorm, ~, KERNELTYPES);
BOOST_PP_SEQ_FOR_EACH(INST_gpu_spf_force, ~, KERNELTYPES);

template cudaError_t gpu_spf_compute_pressure(Scalar3 *d_dpe, unsigned int *d_index_array, unsigned int group_size, StateEquation<tait> *m_eos);
template cudaError_t gpu_spf_compute_pressure(Scalar3 *d_dpe, unsigned int *d_index_array, unsigned int group_size, StateEquation<linear> *m_eos);

#undef INST_gpu_spf_noslip
#undef INST_gpu_spf_ndensity
#undef INST_gpu_spf_ndensitynorm
#undef INST_gpu_spf_force