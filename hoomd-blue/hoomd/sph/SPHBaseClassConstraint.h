/* ---------------------------------------------------------
maintainer: dkrach, david.krach@mib.uni-stuttgart.de
----------------------------------------------------------*/

#include "hoomd/Compute.h"
#include "hoomd/Index1D.h"
#include "hoomd/ParticleGroup.h"
//#include "hoomd/ForceCompute.h"
#include "hoomd/ForceConstraint.h"
#include "hoomd/nsearch/NeighborList.h"

#ifdef ENABLE_MPI
#include "hoomd/Communicator.h"
#include "hoomd/HOOMDMPI.h"
#include "hoomd/Autotuner.h"
#endif

#include <memory>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>

#include "SmoothingKernel.h"
#include "StateEquations.h"

#include "EvaluationMethodDefinition.h"

/*! \file SPHBaseClassConstraint.cc
    \brief Contains base class for any SPH Force compute. Takes care of
           storing SmoothingKernel and NeighborList class instances.
*/

#ifdef __HIPCC__
#error This header cannot be compiled by nvcc
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef __SPHBaseClassConstraint_H__
#define __SPHBaseClassConstraint_H__

const unsigned int NO_AGGREGATE = (unsigned int)0xffffffff;

namespace hoomd
{
namespace sph
{
// //! Enum for various density evaluation approaches
// enum DensityMethod
// {
//     DENSITYSUMMATION,    //!< Summation approach
//     DENSITYCONTINUITY,    //!< Continuity approach
// };

// //! Enum for various viscosity evaluation approaches
// enum ViscosityMethod
// {
//     HARMONICAVERAGE, //!< Viscosity operator based on inter-particle averaged shear stress
// };

template<SmoothingKernelType KT_, StateEquationType SET_>
class PYBIND11_EXPORT SPHBaseClassConstraint : public ForceConstraint
    {
    public:
        //! Constructor
        SPHBaseClassConstraint(std::shared_ptr<SystemDefinition> sysdef,
                     std::shared_ptr<SmoothingKernel<KT_> > skernel,
                     std::shared_ptr<StateEquation<SET_> > eos,
                     std::shared_ptr<nsearch::NeighborList> nlist);

        //! Destructor
        virtual ~SPHBaseClassConstraint();

        //! Return the number of DOF removed by this constraint
        virtual Scalar getNDOFRemoved(std::shared_ptr<ParticleGroup> query)
        {
        return 0;
        }

#ifdef ENABLE_MPI
    //! Get ghost particle fields requested by this pair potential
    virtual CommFlags getRequestedCommFlags(uint64_t timestep)
        {
        CommFlags flags = CommFlags(0);

        // request communication of tags
        flags[comm_flag::tag] = 1;

        flags |= ForceConstraint::getRequestedCommFlags(timestep);

        return flags;
        }
#endif


    //! Return aggregate index
    const Index2D& getAggregateIndexer()
        {
        checkParticlesSorted();

        return m_aggregate_indexer;
        }

    //! Return aggregate list
    const GlobalVector<unsigned int>& getAggregateList()
        {
        checkParticlesSorted();

        return m_aggregate_list;
        }

    //! Return aggregate lengths
    const GlobalVector<unsigned int>& getAggregateLengths()
        {
        checkParticlesSorted();

        return m_aggregate_length;
        }

    //! Return aggregate order
    const GlobalVector<unsigned int>& getAggregateOrder()
        {
        checkParticlesSorted();

        return m_aggregate_order;
        }

    //! Return reverse lookup array
    const GlobalVector<unsigned int>& getAggregateIndex()
        {
        checkParticlesSorted();

        return m_aggregate_idx;
        }

    /// Get the number of aggregates (global)
    unsigned int getNAggregatesGlobal() const
        {
        return m_n_aggregates_global;
        }

        /*! Helper function to compute available type ids for a given group of particles
         * \param pgroup Group of particles to construct type id vectors for
         */
        void constructTypeVectors(std::shared_ptr<ParticleGroup> const pgroup,
                                  std::vector<unsigned int> *global_typeids);

        /*! Helper function to apply external body force to a given group of particles
         * \param pgroup Group of particles to apply body force to
         */
        void applyBodyForce(uint64_t timestep, std::shared_ptr<ParticleGroup> pgroup);
        
// #ifdef ENABLE_HIP
//         void applyBodyForceGPU(uint64_t timestep, std::shared_ptr<ParticleGroup> pgroup);
// #endif

        /*! Set the volumetric acceleration
         * \param gx Volumetric acceleration in x-Direction
         * \param gy Volumetric acceleration in y-Direction
         * \param gz Volumetric acceleration in z-Direction
         * \param damp damping time in units of time steps during which body acceleration is smoothly applied
         */
        void setAcceleration(Scalar gx, Scalar gy, Scalar gz, unsigned int damptime);

        // Get the volumetric acceleration
        Scalar3 getAcceleration(uint64_t timestep);

        
// #ifdef ENABLE_MPI
//         /// The system's communicator.
//         std::shared_ptr<Communicator> m_comm;
//         //! Set the communicator to use
//         void setCommunicator(std::shared_ptr<Communicator> comm)
//             {
//             if (!m_comm)
//                 {
//                 assert(comm);
//                 }
//             Compute::setCommunicator(comm);
//             }
// #endif

    protected:
    std::shared_ptr<SmoothingKernel<KT_> > m_skernel; //!< The kernel function class this method is associated with
    std::shared_ptr<StateEquation<SET_> > m_eos; //!< The equation of state class this method is associated with
    std::shared_ptr<nsearch::NeighborList> m_nlist; //!< The neighbor list to use for the computation

    Index2D m_typpair_idx;        //!< Helper class for indexing per type pair arrays

    GlobalVector<unsigned int> m_aggregate_tag; //!< Aggregate tag per particle tag
    unsigned int m_n_aggregates_global;         //!< Global number of aggregates

    bool m_rebuild_aggregates; //!< True if we need to rebuild indices

    //! Helper function to check if particles have been sorted and rebuild indices if necessary
    virtual void checkParticlesSorted()
        {
        if (m_rebuild_aggregates)
            {
            // rebuild aggregate list
            initAggregates();
            m_rebuild_aggregates = false;
            }
        }

        DensityMethod m_densitymethod;
        ViscosityMethod m_viscositymethod;

        Scalar3 m_bodyforce; //!< Volumetric force
        unsigned int m_damptime; //!< Damping time
        bool m_body_acceleration; //!< True if body acceleration has been set and not null

    private:
    /// 2D Array of aggregate members. Use m_aggregate_indexer to index into this array. The data
    /// stored is
    /// m_aggregate_list[
    ///     m_aggregate_indexer(particle_aggregate_index, aggregate_index)
    /// ] == local_particle_index
    GlobalVector<unsigned int> m_aggregate_list;

    /// List of aggregate lengths
    GlobalVector<unsigned int> m_aggregate_length;

    /// Index of particle in a aggregate. Accessed through local particle index.
    GlobalVector<unsigned int> m_aggregate_order;

    /// Reverse-lookup into aggregate list, specifically
    /// m_aggregate_idx[particle_index] == / aggregate_index (note that this is the temporary
    /// particle index not the permanent particle tag).
    GlobalVector<unsigned int> m_aggregate_idx;

#ifdef ENABLE_HIP
    std::shared_ptr<Autotuner<1>>
        m_tuner_fill; //!< Autotuner for block size for filling the aggregate table
#endif

    /// Functor for indexing into a 1D array as if it were a 2-D array. Index is
    /// [constituent_number, aggregate_number].
    Index2D m_aggregate_indexer;

    void setRebuildAggregates()
        {
        m_rebuild_aggregates = true;
        }

    //! construct a list of local aggregates
    virtual void initAggregates();

#ifdef ENABLE_HIP
    //! construct a list of local aggregates on the GPU
    virtual void initAggregatesGPU();
#endif

#ifdef ENABLE_HIP
    GPUPartition m_gpu_partition; //!< Partition of the aggregates on GPUs
#endif

    // private:
    //     //! Connection to the signal notifying when number of particle types changes
    //     boost::signals2::connection m_num_type_change_connection;

    //     //! Connection to the signal notifying when number of particle types changes
    //     boost::signals2::connection m_particle_num_change_connection;
    };


namespace detail 
{

template<SmoothingKernelType KT_, StateEquationType SET_>
void export_SPHBaseClassConstraint(pybind11::module& m, std::string name);

void export_DensityMethod(pybind11::module& m);

void export_ViscosityMethod(pybind11::module& m);

} // end namespace detail

} // end namespace sph
} // end namespace hoomd

#endif // __SPHBaseClassConstraint_H__
