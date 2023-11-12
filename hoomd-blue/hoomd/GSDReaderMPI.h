// Copyright (c) 2009-2023 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#ifdef __HIPCC__
#error This header cannot be compiled by nvcc
#endif

#include "ParticleData.h"
#include "hoomd/extern/pgsd.h"
#include <string>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#ifndef __GSD_INITIALIZER_MPI_H__
#define __GSD_INITIALIZER_MPI_H__

namespace hoomd
    {
//! Forward declarations
template<class Real> struct SnapshotSystemData;

//! Reads a GSD input file
/*! Read an input GSD file and generate a system snapshot. GSDReaderMPI can read any frame from a GSD
    file into the snapshot. For information on the GSD specification, see http://gsd.readthedocs.io/

    \ingroup data_structs
*/
class PYBIND11_EXPORT GSDReaderMPI
    {
    public:
    //! Loads in the file and parses the data
    GSDReaderMPI(std::shared_ptr<const ExecutionConfiguration> exec_conf,
              const std::string& name,
              const uint64_t frame,
              bool from_end);

    //! Destructor
    ~GSDReaderMPI();

    //! Returns the timestep of the simulation
    uint64_t getTimeStep() const
        {
        uint64_t timestep = m_timestep;

// timestep is only read on the root, broadcast to the other nodes
#ifdef ENABLE_MPI
        const MPI_Comm mpi_comm = m_exec_conf->getMPICommunicator();
        bcast(timestep, 0, mpi_comm);
#endif

        return timestep;
        }

    //! initializes a snapshot with the particle data
    std::shared_ptr<SnapshotSystemData<float>> getSnapshot() const
        {
        return m_snapshot;
        }

    //! initializes a snapshot with the particle data
    uint64_t getFrame() const
        {
        return m_frame;
        }

    //! Helper function to read a quantity from the file
    bool readChunk(void* data,
                   uint64_t frame,
                   const char* name,
                   size_t expected_size,
                   unsigned int cur_n = 0, 
                   uint32_t *offset=nullptr);

    //! clears the snapshot object
    void clearSnapshot()
        {
        m_snapshot.reset();
        }

    //! get handle
    pgsd_handle getHandle(void) const
        {
        return m_handle;
        }

    pybind11::list readTypeShapesPy(uint64_t frame);

    private:
    std::shared_ptr<const ExecutionConfiguration> m_exec_conf; //!< The execution configuration
    uint64_t m_timestep;                                       //!< Timestep at the selected frame
    std::string m_name;                                        //!< Cached file name
    uint64_t m_frame;                                          //!< Cached frame
    std::shared_ptr<SnapshotSystemData<float>> m_snapshot;     //!< The snapshot to read
    pgsd_handle m_handle;                                       //!< Handle to the file
    std::vector<unsigned int> m_part_per_rank;

    //! Helper function to read a type list from the file
    std::vector<std::string> readTypes(uint64_t frame, const char* name);

    // helper functions to read sections of the file
    void readHeader();
    void readParticles();
    void readTopology();
    };

namespace detail
    {
/// Exports GSDReaderMPI to python
void export_GSDReaderMPI(pybind11::module& m);

    } // end namespace detail

    } // end namespace hoomd

#endif
