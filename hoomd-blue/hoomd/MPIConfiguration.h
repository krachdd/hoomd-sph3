// Copyright (c) 2009-2023 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#pragma once

// ensure that HOOMDMath.h is the first header included to work around broken mpi headers
#include "HOOMDMath.h"

#ifdef ENABLE_MPI
#include <mpi.h>
#endif

#include "ClockSource.h"

/*! \file MPIConfiguration.h
    \brief Declares MPIConfiguration, which initializes the MPI environment
*/

#ifdef __HIPCC__
#error This header cannot be compiled by nvcc
#endif

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace hoomd
    {
//! Defines the MPI configuration for the simulation
/*! \ingroup data_structs
    MPIConfiguration is class that stores the MPI communicator and splits it into partitions if
   needed.

    It is constructed *before* ExecutionConfiguration and other classes (Messenger) that need the
   MPI world communicator.
*/
class PYBIND11_EXPORT MPIConfiguration
    {
    public:
    //! Constructor with externally provided MPI communicator (only in MPI enabled builds)
    /*! \param mpi_comm world MPI communicator
     */
    MPIConfiguration(
#ifdef ENABLE_MPI
        MPI_Comm hoomd_world = MPI_COMM_WORLD
#endif
    );

    //! Destructor
    virtual ~MPIConfiguration() {};

#ifdef ENABLE_MPI
    MPI_Comm operator()() const
        {
        return getCommunicator();
        }

    //! Returns the MPI communicator
    MPI_Comm getCommunicator() const
        {
        return m_mpi_comm;
        }
    //! Returns the World MPI communicator
    MPI_Comm getHOOMDWorldCommunicator() const
        {
        return m_hoomd_world;
        }
#endif

    //!< Partition the communicator
    /*! \param nrank Number of ranks per partition
     */
    void splitPartitions(unsigned int nrank);

    //! Return the rank of this processor in the partition
    unsigned int getRank() const
        {
        return m_rank;
        }

    //! Return the global rank of this processor
    unsigned int getRankGlobal() const
        {
#ifdef ENABLE_MPI
        // get rank on world communicator
        int rank;
        MPI_Comm_rank(m_hoomd_world, &rank);
        return rank;
#else
        return 0;
#endif
        }

    //! Return the global communicator size
    unsigned int getNRanksGlobal() const
        {
#ifdef ENABLE_MPI
        int size;
        MPI_Comm_size(m_hoomd_world, &size);
        return size;
#else
        return 1;
#endif
        }

    //! Returns the partition number of this processor
    unsigned int getPartition() const
        {
        return getRankGlobal() / m_n_rank;
        }

    //! Returns the number of partitions
    unsigned int getNPartitions() const
        {
        return getNRanksGlobal() / m_n_rank;
        }

    //! Return the number of ranks in this partition
    unsigned int getNRanks() const;

    //! Returns true if this is the root processor
    bool isRoot() const
        {
        return getRank() == 0;
        }

    //! Perform a job-wide MPI barrier
    void barrier()
        {
#ifdef ENABLE_MPI
        MPI_Barrier(m_mpi_comm);
#endif
        }

    double getWalltime()
        {
        double walltime = static_cast<double>(m_clock.getTime()) / 1e9;
#ifdef ENABLE_MPI
        MPI_Bcast(&walltime, 1, MPI_DOUBLE, 0, m_mpi_comm);
#endif
        return walltime;
        }

    //! added dkrach 
    double bcast_double(double v)
    {
#ifdef ENABLE_MPI
        MPI_Bcast(&v, 1, MPI_DOUBLE, 0, m_mpi_comm);
#endif
    return v;
    }


//     //! Wrapper around MPI_Allgatherv
// std::vector<double> all_gather_dist(double in_value, std::vector<double> out_values, MPI_Comm mpi_comm)
//     {

//     MPI_Allgather(&in_value, 1, MPI_DOUBLE, &out_values, 1, MPI_DOUBLE, mpi_comm);
//     return out_values;
//     // int rank;
//     // int size;
//     // MPI_Comm_rank(mpi_comm, &rank);
//     // MPI_Comm_size(mpi_comm, &size);

//     // // serialize in_value
//     // std::stringstream s(std::ios_base::out | std::ios_base::binary);
//     // cereal::BinaryOutputArchive ar(s);

//     // ar << in_value;
//     // s.flush();

//     // // copy into send buffer
//     // std::string str = s.str();
//     // unsigned int send_count = (unsigned int)str.length();

//     // // allocate memory for buffer lengths
//     // out_values.resize(size);
//     // int* recv_counts = new int[size];
//     // int* displs = new int[size];

//     // // gather lengths of buffers
//     // MPI_Allgather(&send_count, 1, MPI_INT, recv_counts, 1, MPI_INT, mpi_comm);

//     // // allocate receiver buffer
//     // unsigned int len = 0;
//     // for (unsigned int i = 0; i < (unsigned int)size; i++)
//     //     {
//     //     displs[i] = (i > 0) ? displs[i - 1] + recv_counts[i - 1] : 0;
//     //     len += recv_counts[i];
//     //     }
//     // char* rbuf = new char[len];

//     // // now gather actual objects
//     // MPI_Allgatherv((void*)str.data(),
//     //                send_count,
//     //                MPI_BYTE,
//     //                rbuf,
//     //                recv_counts,
//     //                displs,
//     //                MPI_BYTE,
//     //                mpi_comm);

//     // // de-serialize data
//     // for (unsigned int i = 0; i < out_values.size(); i++)
//     //     {
//     //     std::stringstream s(std::string(rbuf + displs[i], recv_counts[i]),
//     //                         std::ios_base::in | std::ios_base::binary);
//     //     cereal::BinaryInputArchive ar(s);

//     //     ar >> out_values[i];
//     //     }

//     // delete[] displs;
//     // delete[] recv_counts;
//     // delete[] rbuf;
//     // return out_values;
//     }



    protected:
#ifdef ENABLE_MPI
    MPI_Comm m_mpi_comm;    //!< The MPI communicator
    MPI_Comm m_hoomd_world; //!< The HOOMD world communicator
#endif
    unsigned int m_rank;   //!< Rank of this processor (0 if running in single-processor mode)
    unsigned int m_n_rank; //!< Ranks per partition

    /// Clock to provide rank synchronized walltime.
    ClockSource m_clock;
    };

namespace detail
    {
//! Exports MPIConfiguration to python
#ifndef __HIPCC__
void export_MPIConfiguration(pybind11::module& m);
#endif

    } // end namespace detail

    } // end namespace hoomd
