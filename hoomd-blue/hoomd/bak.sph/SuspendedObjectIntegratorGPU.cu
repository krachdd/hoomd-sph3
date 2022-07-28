// Copyright (c) 2009-2016 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

// Maintainer: schirwon

#include "SuspendedObjectIntegratorGPU.cuh"
#include "GeneralGPUFunctions.cuh"
#include "GeneralGPUDeviceFunctions.cuh"

#ifndef M_TWOPI
#define M_TWOPI 6.28318530717958647692
#endif

/*! \file SuspendedObjectIntegratorGPU.cu
 \ *brief Defines GPU kernel code for NPT integration on the GPU using the Martyna-Tobias-Klein update equations. Used by SuspendedObjectIntegratorGPU.
 */

//! Kernel to propagate the positions and velocities, first half of NPT update
__global__ void gpu_SuspendedOIntegrator_step_one_kernel(Scalar4 *d_pos,
                                                         Scalar4 *d_vel,
                                                         unsigned int *d_group_members,
                                                         unsigned int group_size,
                                                         const BoxDim box,
                                                         Scalar3 centerofmass,
                                                         vec3<Scalar> angularvel,
                                                         Scalar3 translationvel, Scalar deltaT)
{
    // determine which particle this thread works on
    int group_idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    // initialize eigenvectors
    if (group_idx < group_size)
    {
        unsigned int idx = d_group_members[group_idx];
        
        // Read particle position and compute relative position
        Scalar3 pos = make_scalar3(d_pos[idx].x, d_pos[idx].y, d_pos[idx].z);
        
        // Relative position to center of mass
        Scalar3 rdiff = pos - centerofmass;
        
        // Apply periodic boundary conditions
        rdiff = box.minImage(rdiff);
        
        // Relative position to center of mass
        vec3<Scalar> rdiffv(rdiff.x, rdiff.y, rdiff.z);
        
        
        
        // Rotational contribution
        vec3<Scalar> rot = cross(angularvel,rdiffv);
        
        // Update slave particle velocities
        d_vel[idx].x = translationvel.x + rot.x;
        d_vel[idx].y = translationvel.y + rot.y;
        d_vel[idx].z = translationvel.z + rot.z;
        
        // Update slave particle positions
        d_pos[idx].x += d_vel[idx].x*deltaT;
        d_pos[idx].y += d_vel[idx].y*deltaT;
        d_pos[idx].z += d_vel[idx].z*deltaT;
    }
}

/*! \param d_pos array of particle positions
 \ *param d_vel array of particle velocities
 \param d_accel array of particle accelerations
 \param d_group_members Device array listing the indicies of the mebers of the group to integrate
 \param group_size Number of members in the group
 \param exp_thermo_fac Update factor for thermostat
 \param mat_exp_v Matrix exponential for velocity update
 \param mat_exp_r Matrix exponential for position update
 \param mat_exp_r_int Integrated matrix exp for position update
 \param deltaT Time to advance (for one full step)
 \param deltaT Time to move forward in one whole step
 \param rescale_all True if all particles in the system should be rescaled at once
 
 This is just a kernel driver for gpu_SuspendedObject_step_one_kernel(). See it for more details.
 */
cudaError_t gpu_SuspendedOIntegrator_step_one(Scalar4 *d_pos,
                                              Scalar4 *d_vel,
                                              unsigned int *d_group_members,
                                              unsigned int group_size,
                                              const BoxDim& box,
                                              Scalar3 centerofmass,
                                              vec3<Scalar> angularvel,
                                              Scalar3 translationvel, Scalar deltaT)
{
    // setup the grid to run the kernel
    unsigned int block_size = 256;
    dim3 grid( (group_size / block_size) + 1, 1, 1);
    dim3 threads(block_size, 1, 1);
    
    // run the kernel
    gpu_SuspendedOIntegrator_step_one_kernel<<< grid, threads >>>(d_pos,
                                                                  d_vel,
                                                                  d_group_members,
                                                                  group_size,
                                                                  box,
                                                                  centerofmass,
                                                                  angularvel,
                                                                  translationvel, deltaT);
    
    return cudaSuccess;
}

//! Kernel to propagate the positions and velocities, second half of NPT update
__global__ void gpu_SuspendedOIntegrator_step_two_a_kernel(const Scalar4 *d_pos,
                                                           const Scalar4 *d_net_force,
                                                           unsigned int *d_group_members,
                                                           unsigned int group_size,
                                                           const BoxDim box,
                                                           Scalar3 centerofmass,
                                                           Scalar* totalforcetorque)
{
    // determine which particle this thread works on
    int group_idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    vec3<Scalar> net_force_j;
    Scalar3 torque_j;
    
    if (group_idx < group_size) {
        unsigned int idx = d_group_members[group_idx];
        
        // Net force acting on slave particle
        net_force_j = vec3<Scalar>(d_net_force[idx].x, d_net_force[idx].y, d_net_force[idx].z);
        
        // Read particle position and compute relative position
        Scalar3 pos = make_scalar3(d_pos[idx].x, d_pos[idx].y, d_pos[idx].z);
        
        // Relative position to center of mass
        Scalar3 rdiff = pos - centerofmass;
        
        // Apply periodic boundary conditions
        rdiff = box.minImage(rdiff);
        
        // Relative position to center of mass
        vec3<Scalar> rdiffv(rdiff.x, rdiff.y, rdiff.z);
        
        // Torque contribution
        Scalar3 torque_j = vec_to_scalar3(cross(rdiffv,net_force_j));
    } else {
        net_force_j = vec3<Scalar>(0, 0, 0);
        torque_j = make_scalar3(0, 0, 0);
    }

    __shared__ Scalar shared[REDUCTION_BLOCK_SIZE];

    // TODO too many __syncthreads? -> this can be replaced with an incomplete reduction op 24->12->6->atomicAdd
    //      check other _gpu_sums too
    // Compute total force
    _gpu_sum(shared, group_idx, net_force_j.x, &totalforcetorque[0]);
    _gpu_sum(shared, group_idx, net_force_j.y, &totalforcetorque[1]);
    _gpu_sum(shared, group_idx, net_force_j.z, &totalforcetorque[2]);

    // Compute total torque
    _gpu_sum(shared, group_idx, torque_j.x, &totalforcetorque[3]);
    _gpu_sum(shared, group_idx, torque_j.y, &totalforcetorque[4]);
    _gpu_sum(shared, group_idx, torque_j.z, &totalforcetorque[5]);
}

/*! \param d_vel array of particle velocities
 \ *param d_accel array of particle accelerations
 \param d_group_members Device array listing the indicies of the mebers of the group to integrate
 \param group_size Number of members in the group
 \param mat_exp_v Matrix exponential for velocity update
 \param d_net_force Net force on each particle
 
 \param deltaT Time to move forward in one whole step
 
 This is just a kernel driver for gpu_SuspendedObject_step_kernel(). See it for more details.
 */
cudaError_t gpu_SuspendedOIntegrator_step_two_a(Scalar4 *d_pos,
                                                Scalar4 *d_net_force,
                                                unsigned int *d_group_members,
                                                unsigned int group_size,
                                                const BoxDim& box,
                                                Scalar3 centerofmass,
                                                Scalar* totalforcetorque)
{
    // setup the grid to run the kernel
    unsigned int block_size = REDUCTION_BLOCK_SIZE;
    dim3 grid( (group_size / block_size) + 1, 1, 1);
    dim3 threads(block_size, 1, 1);
    
    // run the kernel
    gpu_SuspendedOIntegrator_step_two_a_kernel<<< grid, threads >>>(d_pos,
                                                                    d_net_force,
                                                                    d_group_members,
                                                                    group_size,
                                                                    box,
                                                                    centerofmass,
                                                                    totalforcetorque);
    
    return cudaSuccess;
}

//! Kernel to propagate the positions and velocities, second half of NPT update
__global__ void gpu_SuspendedOIntegrator_step_two_b_kernel(Scalar4 *d_vel,
                                                           Scalar4 *d_pos,
                                                           unsigned int *d_group_members,
                                                           unsigned int group_size,
                                                           const BoxDim box,
                                                           Scalar3 centerofmass,
                                                           Scalar3 translationvel,
                                                           vec3<Scalar> angularvel)
{
    // determine which particle this thread works on
    int group_idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (group_idx < group_size)
    {
        unsigned int idx = d_group_members[group_idx];
        
        Scalar3 pos = make_scalar3(d_pos[idx].x, d_pos[idx].y, d_pos[idx].z);
        
        // Relative position to center of mass
        Scalar3 rdiff = pos - centerofmass;
        
        // Apply periodic boundary conditions
        rdiff = box.minImage(rdiff);
        
        // Relative position to center of mass
        vec3<Scalar> rdiffv(rdiff.x, rdiff.y, rdiff.z);
        
        // Rotational contribution
        vec3<Scalar> rot = cross(angularvel,rdiffv);
        
        // Update slave particle velocities
        d_vel[idx].x = translationvel.x + rot.x;
        d_vel[idx].y = translationvel.y + rot.y;
        d_vel[idx].z = translationvel.z + rot.z;
    }
}

/*! \param d_vel array of particle velocities
 \ *param d_accel array of particle accelerations
 \param d_group_members Device array listing the indicies of the mebers of the group to integrate
 \param group_size Number of members in the group
 \param mat_exp_v Matrix exponential for velocity update
 \param d_net_force Net force on each particle
 \param deltaT Time to move forward in one whole step
 
 This is just a kernel driver for gpu_SuspendedObject_step_kernel(). See it for more details.
 */
cudaError_t gpu_SuspendedOIntegrator_step_two_b(Scalar4 *d_vel,
                                                Scalar4 *d_pos,
                                                unsigned int *d_group_members,
                                                unsigned int group_size,
                                                const BoxDim& box,
                                                Scalar3 centerofmass,
                                                Scalar3 translationvel,
                                                vec3<Scalar> angularvel)
{
    // setup the grid to run the kernel
    unsigned int block_size=256;
    dim3 grid( (group_size / block_size) + 1, 1, 1);
    dim3 threads(block_size, 1, 1);
    
    // run the kernel
    gpu_SuspendedOIntegrator_step_two_b_kernel<<< grid, threads >>>(d_vel,
                                                                    d_pos,
                                                                    d_group_members,
                                                                    group_size,
                                                                    box,
                                                                    centerofmass,
                                                                    translationvel,
                                                                    angularvel);
    
    return cudaSuccess;
}


__global__ void gpu_SuspendedOIntegrator_compute_center_of_mass_kernel(
    const Scalar4* d_pos, 
    const unsigned int* indices, unsigned int count, 
    Scalar* dst_angles,
    Scalar3 lo, Scalar3 Ld
) {
    unsigned int ix = blockIdx.x * blockDim.x + threadIdx.x;
    Scalar3 pos;
    if (ix < count) {
        unsigned int j = indices[ix];
        pos = make_scalar3(d_pos[j].x, d_pos[j].y, d_pos[j].z);
    } else {
        pos = make_scalar3(0, 0, 0);
    }

    __shared__ Scalar shared[REDUCTION_BLOCK_SIZE];

    // _gpu_sum will only consider values with tid < count
    // therefore it is ok to calcualte some trash values
    Scalar theta_x = ((pos.x - lo.x) / Ld.x) * Scalar(M_TWOPI);
    Scalar theta_y = ((pos.y - lo.y) / Ld.y) * Scalar(M_TWOPI);
    Scalar theta_z = ((pos.z - lo.z) / Ld.z) * Scalar(M_TWOPI);

    // TODO this has some more __syncthreads than neccessary. These can be reduced at the cost of more shared mem 
    _gpu_sum(shared, count, (Ld.x / Scalar(M_TWOPI))*cos(theta_x), &dst_angles[0]);
    _gpu_sum(shared, count, (Ld.x / Scalar(M_TWOPI))*sin(theta_x), &dst_angles[1]);
    _gpu_sum(shared, count, (Ld.y / Scalar(M_TWOPI))*cos(theta_y), &dst_angles[2]);
    _gpu_sum(shared, count, (Ld.y / Scalar(M_TWOPI))*sin(theta_y), &dst_angles[3]);
    _gpu_sum(shared, count, (Ld.z / Scalar(M_TWOPI))*cos(theta_z), &dst_angles[4]);
    _gpu_sum(shared, count, (Ld.z / Scalar(M_TWOPI))*sin(theta_z), &dst_angles[5]);
}


cudaError_t gpu_SuspendedOIntegrator_compute_center_of_mass(
    const Scalar4* pos, 
    const unsigned int* indices, unsigned int count, 
    Scalar* dst_angles,
    Scalar3 lo, Scalar3 Ld
) {
    // setup the grid to run the kernel
    unsigned int block_size = REDUCTION_BLOCK_SIZE;
    dim3 grid( (count / block_size) + 1, 1, 1);
    dim3 threads(block_size, 1, 1);

    // run the kernel
    gpu_SuspendedOIntegrator_compute_center_of_mass_kernel<<<grid, threads>>>(pos, indices, count, dst_angles, lo, Ld);

    return cudaSuccess;
}

__global__ void gpu_SuspendedOIntegrator_compute_translation_velocity_kernel(
    const Scalar4* d_vel,
    Scalar totalmass,
    const unsigned int* indices, unsigned int count, 
    Scalar* translationvel
) {
    unsigned int ix = blockIdx.x * blockDim.x + threadIdx.x;

    Scalar3 vel;
    Scalar m_M;

    if (ix < count) {
        unsigned int j = indices[ix];
        vel = make_scalar3(d_vel[j].x, d_vel[j].y, d_vel[j].z);
        m_M = d_vel[j].w/totalmass;
    } else {
        vel = make_scalar3(0, 0, 0);
        m_M = 0;
    }

    __shared__ Scalar shared[REDUCTION_BLOCK_SIZE];

    _gpu_sum(shared, count, m_M*vel.x, &translationvel[0]);
    _gpu_sum(shared, count, m_M*vel.y, &translationvel[1]);
    _gpu_sum(shared, count, m_M*vel.z, &translationvel[2]);
}

cudaError_t gpu_SuspendedOIntegrator_compute_translation_velocity(
    const Scalar4* vel,
    Scalar totalmass,
    const unsigned int* indices, unsigned int count, 
    Scalar* translationvel
) {
    // setup the grid to run the kernel
    unsigned int block_size = REDUCTION_BLOCK_SIZE;
    dim3 grid( (count / block_size) + 1, 1, 1);
    dim3 threads(block_size, 1, 1);

    // run the kernel
    gpu_SuspendedOIntegrator_compute_translation_velocity_kernel<<<grid, threads>>>(vel, totalmass, indices, count, translationvel);

    return cudaSuccess;
}

__global__ void gpu_SuspendedOIntegrator_compute_moment_of_intertia_kernel(
    const Scalar4* d_pos, const Scalar4* vel, 
    const unsigned int* indices, unsigned int count, 
    Scalar3 centerofmass, BoxDim box,
    Scalar* J
) {
    unsigned int ix = blockIdx.x * blockDim.x + threadIdx.x;

    Scalar mass;
    Scalar dotrdiff;
    Scalar3 rdiff;
    if (ix < count) {
        unsigned int j = indices[ix];

        // Read particle position and mass
        Scalar3 pos = make_scalar3(d_pos[j].x, d_pos[j].y, d_pos[j].z);
        mass = vel[j].w;

        // Relative position to center of mass
        Scalar3 rdiff = pos - centerofmass;

        // Apply periodic boundary conditions
        rdiff = box.minImage(rdiff);

        Scalar dotrdiff = dot(rdiff,rdiff);
    } else {
        mass = 0;
        rdiff = make_scalar3(0, 0, 0);
        dotrdiff = 0;
    }

    __shared__ Scalar shared[REDUCTION_BLOCK_SIZE];

    // Add contribution to moment of inertia tensor
    // TODO many __syncthreads
    _gpu_sum(shared, count, mass*(dotrdiff-rdiff.x*rdiff.x), &J[0]); // 00
    _gpu_sum(shared, count, mass*(        -rdiff.x*rdiff.y), &J[1]); // 01
    _gpu_sum(shared, count, mass*(        -rdiff.x*rdiff.z), &J[2]); // 02
    _gpu_sum(shared, count, mass*(        -rdiff.y*rdiff.x), &J[3]); // 10
    _gpu_sum(shared, count, mass*(dotrdiff-rdiff.y*rdiff.y), &J[4]); // 11
    _gpu_sum(shared, count, mass*(        -rdiff.y*rdiff.z), &J[5]); // 12
    _gpu_sum(shared, count, mass*(        -rdiff.z*rdiff.x), &J[6]); // 20
    _gpu_sum(shared, count, mass*(        -rdiff.z*rdiff.y), &J[7]); // 21
    _gpu_sum(shared, count, mass*(dotrdiff-rdiff.z*rdiff.z), &J[8]); // 22
}

cudaError_t gpu_SuspendedOIntegrator_compute_moment_of_intertia(
    const Scalar4* pos, const Scalar4* vel, 
    const unsigned int* indices, unsigned int count, 
    Scalar3 centerofmass, BoxDim box,
    Scalar* J
) {
    // setup the grid to run the kernel
    unsigned int block_size = REDUCTION_BLOCK_SIZE;
    dim3 grid( (count / block_size) + 1, 1, 1);
    dim3 threads(block_size, 1, 1);

    // run the kernel
    gpu_SuspendedOIntegrator_compute_moment_of_intertia_kernel<<<grid, threads>>>(pos, vel, indices, count, centerofmass, box, J);

    return cudaSuccess;
}

__global__ void gpu_suspendedOIntegrator_compute_angular_velocity_kernel(
    Scalar4* d_pos, Scalar4* d_vel, 
    Scalar3 centerofmass, Scalar3 translationvel, BoxDim box,
    const unsigned int* indices, unsigned int count, 
    Scalar* angmomentum
) {
    unsigned int ix = blockIdx.x * blockDim.x + threadIdx.x;

    Scalar3 angmomentum_j;
    if (ix < count) {
        unsigned int j = indices[ix];

        // Read particle position and mass
        Scalar3 pos = make_scalar3(d_pos[j].x, d_pos[j].y, d_pos[j].z);
        Scalar3 vel = make_scalar3(d_vel[j].x, d_vel[j].y, d_vel[j].z);
        Scalar mass = d_vel[j].w;

        // Relative position to center of mass
        Scalar3 rdiff = pos - centerofmass;

        // Apply periodic boundary conditions
        rdiff = box.minImage(rdiff);

        // Relative position to center of mass
        vec3<Scalar> rdiffm(rdiff.x*mass,rdiff.y*mass,rdiff.z*mass);
        vec3<Scalar> vdiff(vel.x-translationvel.x,vel.y-translationvel.y,vel.z-translationvel.z);

        // Angular momentum due to particle j
        angmomentum_j = vec_to_scalar3(cross(rdiffm,vdiff));
    } else {
        angmomentum_j = make_scalar3(0, 0, 0);
    }

    __shared__ Scalar shared[REDUCTION_BLOCK_SIZE];


    // Compute total angular momentum
    _gpu_sum(shared, count, angmomentum_j.x, &angmomentum[0]);
    _gpu_sum(shared, count, angmomentum_j.y, &angmomentum[1]);
    _gpu_sum(shared, count, angmomentum_j.z, &angmomentum[2]);
}

cudaError_t gpu_suspendedOIntegrator_compute_angular_velocity(
    Scalar4* pos, Scalar4* vel, 
    Scalar3 centerofmass, Scalar3 translationvel, BoxDim box,
    const unsigned int* indices, unsigned int count, 
    Scalar* angmomentum
) {
    // setup the grid to run the kernel
    unsigned int block_size = REDUCTION_BLOCK_SIZE;
    dim3 grid( (count / block_size) + 1, 1, 1);
    dim3 threads(block_size, 1, 1);

    // run the kernel
    gpu_suspendedOIntegrator_compute_angular_velocity_kernel<<<grid, threads>>>(
        pos, vel, 
        centerofmass, translationvel, box,
        indices, count,
        angmomentum
    );

    return cudaSuccess;
}