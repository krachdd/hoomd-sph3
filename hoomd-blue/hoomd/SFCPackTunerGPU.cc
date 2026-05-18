// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

/*! \file SFCPackTunerGPU.cc
    \brief Defines the SFCPackTunerGPU class
*/

#ifdef ENABLE_HIP

#include "SFCPackTunerGPU.h"
#include "SFCPackTunerGPU.cuh"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdexcept>

using namespace std;

namespace hoomd
    {
//! Constructor
/*! \param sysdef System to perform sorts on
 */
SFCPackTunerGPU::SFCPackTunerGPU(std::shared_ptr<SystemDefinition> sysdef,
                                 std::shared_ptr<Trigger> trigger)
    : SFCPackTuner(sysdef, trigger)
    {
    m_exec_conf->msg->notice(5) << "Constructing SFCPackTunerGPU" << endl;

    // perform lots of sanity checks
    assert(m_pdata);

    GPUArray<unsigned int> gpu_sort_order(m_pdata->getMaxN(), m_exec_conf);
    m_gpu_sort_order.swap(gpu_sort_order);

    GPUArray<unsigned int> gpu_particle_bins(m_pdata->getMaxN(), m_exec_conf);
    m_gpu_particle_bins.swap(gpu_particle_bins);
    }

/*! reallocate the internal arrays
 */
void SFCPackTunerGPU::reallocate()
    {
    m_gpu_sort_order.resize(m_pdata->getMaxN());
    m_gpu_particle_bins.resize(m_pdata->getMaxN());
    }

/*! Destructor
 */
SFCPackTunerGPU::~SFCPackTunerGPU()
    {
    m_exec_conf->msg->notice(5) << "Destroying SFCPackTunerGPU" << endl;
    }

void SFCPackTunerGPU::getSortedOrder2D()
    {
    // on the GPU, getSortedOrder3D handles both cases
    getSortedOrder3D();
    }

void SFCPackTunerGPU::getSortedOrder3D()
    {
    // start by checking the saneness of some member variables
    assert(m_pdata);
    assert(m_gpu_sort_order.getNumElements() >= m_pdata->getN());

    // make even bin dimensions
    const BoxDim& box = m_pdata->getBox();

    // reallocate memory arrays if m_grid changed
    // also regenerate the traversal order
    if ((m_last_grid != m_grid || m_last_dim != 3) && m_sysdef->getNDimensions() == 3)
        {
        if (m_grid > 256)
            {
            unsigned int mb = m_grid * m_grid * m_grid * 4 / 1024 / 1024;
            m_exec_conf->msg->warning()
                << "sorter is about to allocate a very large amount of memory (" << mb << "MB)"
                << " and may crash." << endl;
            m_exec_conf->msg->warning() << "            Reduce the amount of memory allocated to "
                                           "prevent this by decreasing the "
                                        << endl;
            m_exec_conf->msg->warning() << "            grid dimension (i.e. "
                                           "sorter.set_params(grid=128) ) or by disabling it "
                                        << endl;
            m_exec_conf->msg->warning()
                << "            ( sorter.disable() ) before beginning the run()." << endl;
            }

        // generate the traversal order
        GPUArray<unsigned int> traversal_order(m_grid * m_grid * m_grid, m_exec_conf);
        m_traversal_order.swap(traversal_order);

        vector<unsigned int> reverse_order(m_grid * m_grid * m_grid);
        reverse_order.clear();

        // we need to start the hilbert curve with a seed order 0,1,2,3,4,5,6,7
        unsigned int cell_order[8];
        for (unsigned int i = 0; i < 8; i++)
            cell_order[i] = i;
        generateTraversalOrder(0, 0, 0, m_grid, m_grid, cell_order, reverse_order);

        // access traversal order
        ArrayHandle<unsigned int> h_traversal_order(m_traversal_order,
                                                    access_location::host,
                                                    access_mode::overwrite);

        for (unsigned int i = 0; i < m_grid * m_grid * m_grid; i++)
            h_traversal_order.data[reverse_order[i]] = i;

        m_last_grid = m_grid;
        // store the last system dimension computed so we can be mindful if that ever changes
        m_last_dim = m_sysdef->getNDimensions();
        }

    // sanity checks
    assert(m_gpu_particle_bins.getNumElements() >= m_pdata->getN());

    // access arrays
    ArrayHandle<Scalar4> d_pos(m_pdata->getPositions(), access_location::device, access_mode::read);
    ArrayHandle<unsigned int> d_gpu_particle_bins(m_gpu_particle_bins,
                                                  access_location::device,
                                                  access_mode::overwrite);
    ArrayHandle<unsigned int> d_gpu_sort_order(m_gpu_sort_order,
                                               access_location::device,
                                               access_mode::overwrite);
    ArrayHandle<unsigned int> d_traversal_order(m_traversal_order,
                                                access_location::device,
                                                access_mode::read);

    // put the particles in the bins and sort
    kernel::gpu_generate_sorted_order(m_pdata->getN(),
                                      d_pos.data,
                                      d_gpu_particle_bins.data,
                                      d_traversal_order.data,
                                      m_grid,
                                      d_gpu_sort_order.data,
                                      box,
                                      m_sysdef->getNDimensions() == 2,
                                      m_exec_conf->getCachedAllocator());

    if (m_exec_conf->isCUDAErrorCheckingEnabled())
        CHECK_CUDA_ERROR();
    }

void SFCPackTunerGPU::applySortOrder()
    {
    assert(m_pdata);
    assert(m_gpu_sort_order.getNumElements() >= m_pdata->getN());

        {
        // access alternate arrays to write to
        ArrayHandle<Scalar4> d_pos_alt(m_pdata->getAltPositions(),
                                       access_location::device,
                                       access_mode::overwrite);
        ArrayHandle<Scalar4> d_vel_alt(m_pdata->getAltVelocities(),
                                       access_location::device,
                                       access_mode::overwrite);
        ArrayHandle<Scalar> d_density_alt(m_pdata->getAltDensities(),
                                          access_location::device,
                                          access_mode::overwrite);
        ArrayHandle<Scalar> d_pressure_alt(m_pdata->getAltPressures(),
                                           access_location::device,
                                           access_mode::overwrite);
        ArrayHandle<Scalar> d_energy_alt(m_pdata->getAltEnergies(),
                                         access_location::device,
                                         access_mode::overwrite);
        ArrayHandle<Scalar3> d_aux1_alt(m_pdata->getAltAuxiliaries1(),
                                        access_location::device,
                                        access_mode::overwrite);
        ArrayHandle<Scalar3> d_aux2_alt(m_pdata->getAltAuxiliaries2(),
                                        access_location::device,
                                        access_mode::overwrite);
        ArrayHandle<Scalar3> d_aux3_alt(m_pdata->getAltAuxiliaries3(),
                                        access_location::device,
                                        access_mode::overwrite);
        ArrayHandle<Scalar3> d_aux4_alt(m_pdata->getAltAuxiliaries4(),
                                        access_location::device,
                                        access_mode::overwrite);
        ArrayHandle<Scalar> d_slength_alt(m_pdata->getAltSlenghts(),
                                          access_location::device,
                                          access_mode::overwrite);
        ArrayHandle<Scalar3> d_accel_alt(m_pdata->getAltAccelerations(),
                                         access_location::device,
                                         access_mode::overwrite);
        ArrayHandle<Scalar3> d_dpedt_alt(m_pdata->getAltDPEdts(),
                                         access_location::device,
                                         access_mode::overwrite);
        ArrayHandle<int3> d_image_alt(m_pdata->getAltImages(),
                                      access_location::device,
                                      access_mode::overwrite);
        ArrayHandle<unsigned int> d_body_alt(m_pdata->getAltBodies(),
                                             access_location::device,
                                             access_mode::overwrite);
        ArrayHandle<unsigned int> d_tag_alt(m_pdata->getAltTags(),
                                            access_location::device,
                                            access_mode::overwrite);
        ArrayHandle<Scalar4> d_net_force_alt(m_pdata->getAltNetForce(),
                                             access_location::device,
                                             access_mode::overwrite);
        ArrayHandle<Scalar4> d_net_ratedpe_alt(m_pdata->getAltNetRateDPE(),
                                               access_location::device,
                                               access_mode::overwrite);

        // access live particle data to read from
        ArrayHandle<Scalar4> d_pos(m_pdata->getPositions(),
                                   access_location::device,
                                   access_mode::read);
        ArrayHandle<Scalar4> d_vel(m_pdata->getVelocities(),
                                   access_location::device,
                                   access_mode::read);
        ArrayHandle<Scalar> d_density(m_pdata->getDensities(),
                                      access_location::device,
                                      access_mode::read);
        ArrayHandle<Scalar> d_pressure(m_pdata->getPressures(),
                                       access_location::device,
                                       access_mode::read);
        ArrayHandle<Scalar> d_energy(m_pdata->getEnergies(),
                                     access_location::device,
                                     access_mode::read);
        ArrayHandle<Scalar3> d_aux1(m_pdata->getAuxiliaries1(),
                                    access_location::device,
                                    access_mode::read);
        ArrayHandle<Scalar3> d_aux2(m_pdata->getAuxiliaries2(),
                                    access_location::device,
                                    access_mode::read);
        ArrayHandle<Scalar3> d_aux3(m_pdata->getAuxiliaries3(),
                                    access_location::device,
                                    access_mode::read);
        ArrayHandle<Scalar3> d_aux4(m_pdata->getAuxiliaries4(),
                                    access_location::device,
                                    access_mode::read);
        ArrayHandle<Scalar> d_slength(m_pdata->getSlengths(),
                                      access_location::device,
                                      access_mode::read);
        ArrayHandle<Scalar3> d_accel(m_pdata->getAccelerations(),
                                     access_location::device,
                                     access_mode::read);
        ArrayHandle<Scalar3> d_dpedt(m_pdata->getDPEdts(),
                                     access_location::device,
                                     access_mode::read);
        ArrayHandle<int3> d_image(m_pdata->getImages(),
                                  access_location::device,
                                  access_mode::read);
        ArrayHandle<unsigned int> d_body(m_pdata->getBodies(),
                                         access_location::device,
                                         access_mode::read);
        ArrayHandle<unsigned int> d_tag(m_pdata->getTags(),
                                        access_location::device,
                                        access_mode::read);
        ArrayHandle<Scalar4> d_net_force(m_pdata->getNetForce(),
                                         access_location::device,
                                         access_mode::read);
        ArrayHandle<Scalar4> d_net_ratedpe(m_pdata->getNetRateDPE(),
                                           access_location::device,
                                           access_mode::read);

        // access rtags
        ArrayHandle<unsigned int> d_rtag(m_pdata->getRTags(),
                                         access_location::device,
                                         access_mode::readwrite);

        // access sorted order
        ArrayHandle<unsigned int> d_gpu_sort_order(m_gpu_sort_order,
                                                   access_location::device,
                                                   access_mode::read);

        // apply sorted order and re-build rtags
        kernel::gpu_apply_sorted_order(m_pdata->getN(),
                                       m_pdata->getNGhosts(),
                                       d_gpu_sort_order.data,
                                       d_pos.data,
                                       d_pos_alt.data,
                                       d_vel.data,
                                       d_vel_alt.data,
                                       d_density.data,
                                       d_density_alt.data,
                                       d_pressure.data,
                                       d_pressure_alt.data,
                                       d_energy.data,
                                       d_energy_alt.data,
                                       d_aux1.data,
                                       d_aux1_alt.data,
                                       d_aux2.data,
                                       d_aux2_alt.data,
                                       d_aux3.data,
                                       d_aux3_alt.data,
                                       d_aux4.data,
                                       d_aux4_alt.data,
                                       d_slength.data,
                                       d_slength_alt.data,
                                       d_accel.data,
                                       d_accel_alt.data,
                                       d_dpedt.data,
                                       d_dpedt_alt.data,
                                       d_image.data,
                                       d_image_alt.data,
                                       d_body.data,
                                       d_body_alt.data,
                                       d_tag.data,
                                       d_tag_alt.data,
                                       d_net_force.data,
                                       d_net_force_alt.data,
                                       d_net_ratedpe.data,
                                       d_net_ratedpe_alt.data,
                                       d_rtag.data);

        if (m_exec_conf->isCUDAErrorCheckingEnabled())
            CHECK_CUDA_ERROR();
        }

    // make alternate arrays current
    m_pdata->swapPositions();
    m_pdata->swapVelocities();
    m_pdata->swapDensities();
    m_pdata->swapPressures();
    m_pdata->swapEnergies();
    m_pdata->swapAuxiliaries1();
    m_pdata->swapAuxiliaries2();
    m_pdata->swapAuxiliaries3();
    m_pdata->swapAuxiliaries4();
    m_pdata->swapSlengths();
    m_pdata->swapAccelerations();
    m_pdata->swapDPEdts();
    m_pdata->swapImages();
    m_pdata->swapBodies();
    m_pdata->swapTags();
    m_pdata->swapNetForce();
    m_pdata->swapNetRateDPE();
    }

namespace detail
    {
void export_SFCPackTunerGPU(pybind11::module& m)
    {
    pybind11::class_<SFCPackTunerGPU, SFCPackTuner, std::shared_ptr<SFCPackTunerGPU>>(
        m,
        "SFCPackTunerGPU")
        .def(pybind11::init<std::shared_ptr<SystemDefinition>, std::shared_ptr<Trigger>>());
    }

    } // end namespace detail

    } // end namespace hoomd

#endif
