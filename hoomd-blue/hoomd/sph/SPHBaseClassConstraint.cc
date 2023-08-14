/* ---------------------------------------------------------
maintainer: dkrach, david.krach@mib.uni-stuttgart.de
----------------------------------------------------------*/

#include "SPHBaseClassConstraint.h"
#include "hoomd/nsearch/NeighborList.h"

#include "hoomd/Autotuner.h"
#include "hoomd/CachedAllocator.h"

#ifdef ENABLE_HIP
#include "SPHBaseClassConstraint.cuh"
#endif


#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>

#include <stdexcept>
#include <utility>
#include <set>
#include <algorithm>

#include <map>
#include <string.h>

using namespace std;

namespace hoomd
{
namespace sph
{
/*! Constructor
*/
template<SmoothingKernelType KT_, StateEquationType SET_>
SPHBaseClassConstraint<KT_, SET_>::SPHBaseClassConstraint(std::shared_ptr<SystemDefinition> sysdef,
                           std::shared_ptr<SmoothingKernel<KT_> > skernel,
                           std::shared_ptr<StateEquation<SET_> > eos,
                           std::shared_ptr<nsearch::NeighborList> nlist)
    : ForceConstraint(sysdef), m_skernel(skernel), m_eos(eos), m_nlist(nlist),
      m_typpair_idx(m_pdata->getNTypes()), m_aggregate_tag(m_exec_conf),
      m_n_aggregates_global(0), m_rebuild_aggregates(true), m_aggregate_list(m_exec_conf),
      m_aggregate_length(m_exec_conf), m_aggregate_order(m_exec_conf), m_aggregate_idx(m_exec_conf)
        {

        // connect to the ParticleData to receive notifications when particles change order in memory
        m_pdata->getParticleSortSignal()
            .connect<SPHBaseClassConstraint, &SPHBaseClassConstraint::setRebuildAggregates>(this);

        TAG_ALLOCATION(m_aggregate_tag);
        TAG_ALLOCATION(m_aggregate_list);
        TAG_ALLOCATION(m_aggregate_length);
        TAG_ALLOCATION(m_aggregate_order);
        TAG_ALLOCATION(m_aggregate_idx);

        m_exec_conf->msg->notice(5) << "Constructing SPHBaseClassConstraint" << std::endl;

        // Sanity checks
        assert(m_nlist);

        assert(m_skernel);
        assert(m_eos);

        // Set default variables
        m_bodyforce = make_scalar3(Scalar(0), Scalar(0), Scalar(0));
        m_damptime = 0;
        m_body_acceleration = false;

#ifdef ENABLE_HIP
    if (m_exec_conf->isCUDAEnabled())
        {
        // initialize autotuner
        std::vector<unsigned int> valid_params;
        unsigned int warp_size = m_exec_conf->dev_prop.warpSize;
        for (unsigned int block_size = warp_size; block_size <= 1024; block_size += warp_size)
            valid_params.push_back(block_size);

        m_tuner_fill.reset(new Autotuner<1>({AutotunerBase::makeBlockSizeRange(this->m_exec_conf)},
                                            this->m_exec_conf,
                                            "fill_aggregate_table"));
        this->m_autotuners.push_back(m_tuner_fill);
        }
#endif
      }

/*! Destructor
*/
template<SmoothingKernelType KT_, StateEquationType SET_>
SPHBaseClassConstraint<KT_, SET_>::~SPHBaseClassConstraint()
    {
    m_pdata->getParticleSortSignal()
        .disconnect<SPHBaseClassConstraint, &SPHBaseClassConstraint::setRebuildAggregates>(this);
    m_exec_conf->msg->notice(5) << "Destroying SPHBaseClassConstraint" << std::endl;
    }

#ifdef ENABLE_HIP
template<SmoothingKernelType KT_, StateEquationType SET_>
void SPHBaseClassConstraint<KT_, SET_>::initAggregatesGPU()
    {
    unsigned int nptl_local = m_pdata->getN() + m_pdata->getNGhosts();

    unsigned int n_local_aggregates = 0;

    // maximum aggregate length
    unsigned int nmax = 0;

    // number of local particles that are part of aggregates
    unsigned int n_local_ptls_in_aggregates = 0;

    // resize to maximum possible number of local aggregates
    m_aggregate_length.resize(nptl_local);
    m_aggregate_idx.resize(nptl_local);

    ScopedAllocation<unsigned int> d_idx_sorted_by_tag(m_exec_conf->getCachedAllocator(),
                                                       nptl_local);
    ScopedAllocation<unsigned int> d_local_aggregates_lowest_idx(m_exec_conf->getCachedAllocator(),
                                                                nptl_local);

        {
        ArrayHandle<unsigned int> d_aggregate_tag(m_aggregate_tag,
                                                 access_location::device,
                                                 access_mode::read);
        ArrayHandle<unsigned int> d_tag(m_pdata->getTags(),
                                        access_location::device,
                                        access_mode::read);
        ArrayHandle<unsigned int> d_aggregate_length(m_aggregate_length,
                                                    access_location::device,
                                                    access_mode::overwrite);
        ArrayHandle<unsigned int> d_aggregate_idx(m_aggregate_idx,
                                                 access_location::device,
                                                 access_mode::overwrite);

        // temporary buffers
        ScopedAllocation<unsigned int> d_local_aggregate_tags(m_exec_conf->getCachedAllocator(),
                                                             nptl_local);
        ScopedAllocation<unsigned int> d_local_unique_aggregate_tags(
            m_exec_conf->getCachedAllocator(),
            m_n_aggregates_global);
        ScopedAllocation<unsigned int> d_sorted_by_tag(m_exec_conf->getCachedAllocator(),
                                                       nptl_local);
        ScopedAllocation<unsigned int> d_idx_sorted_by_aggregate_and_tag(
            m_exec_conf->getCachedAllocator(),
            nptl_local);

        ScopedAllocation<unsigned int> d_lowest_idx(m_exec_conf->getCachedAllocator(),
                                                    m_n_aggregates_global);
        ScopedAllocation<unsigned int> d_lowest_idx_sort(m_exec_conf->getCachedAllocator(),
                                                         m_n_aggregates_global);
        ScopedAllocation<unsigned int> d_lowest_idx_in_aggregates(m_exec_conf->getCachedAllocator(),
                                                                 m_n_aggregates_global);
        ScopedAllocation<unsigned int> d_lowest_idx_by_aggregate_tag(
            m_exec_conf->getCachedAllocator(),
            m_aggregate_tag.getNumElements());

        kernel::gpu_sort_by_aggregate(nptl_local,
                                     d_tag.data,
                                     d_aggregate_tag.data,
                                     d_local_aggregate_tags.data,
                                     d_local_aggregates_lowest_idx.data,
                                     d_local_unique_aggregate_tags.data,
                                     d_aggregate_idx.data,
                                     d_sorted_by_tag.data,
                                     d_idx_sorted_by_tag.data,
                                     d_idx_sorted_by_aggregate_and_tag.data,
                                     d_lowest_idx.data,
                                     d_lowest_idx_sort.data,
                                     d_lowest_idx_in_aggregates.data,
                                     d_lowest_idx_by_aggregate_tag.data,
                                     d_aggregate_length.data,
                                     n_local_aggregates,
                                     nmax,
                                     n_local_ptls_in_aggregates,
                                     m_exec_conf->getCachedAllocator(),
                                     m_exec_conf->isCUDAErrorCheckingEnabled());

        if (m_exec_conf->isCUDAErrorCheckingEnabled())
            CHECK_CUDA_ERROR();
        }

    // set up indexer
    m_aggregate_indexer = Index2D(nmax, n_local_aggregates);

    m_exec_conf->msg->notice(7) << "SPHBaseClassConstraint: " << n_local_aggregates << " aggregates, "
                                << n_local_ptls_in_aggregates << " particles in aggregates "
                                << std::endl;

    // resize aggregate list
    m_aggregate_list.resize(m_aggregate_indexer.getNumElements());

    // resize aggregate lookup to size of local particle data
    m_aggregate_order.resize(m_pdata->getMaxN());

        {
        // write out aggregate list and order
        ArrayHandle<unsigned int> d_aggregate_list(m_aggregate_list,
                                                  access_location::device,
                                                  access_mode::overwrite);
        ArrayHandle<unsigned int> d_aggregate_order(m_aggregate_order,
                                                   access_location::device,
                                                   access_mode::overwrite);
        ArrayHandle<unsigned int> d_aggregate_idx(m_aggregate_idx,
                                                 access_location::device,
                                                 access_mode::read);

        m_tuner_fill->begin();
        unsigned int block_size = m_tuner_fill->getParam()[0];

        kernel::gpu_fill_aggregate_table(nptl_local,
                                        n_local_ptls_in_aggregates,
                                        m_aggregate_indexer,
                                        d_aggregate_idx.data,
                                        d_local_aggregates_lowest_idx.data,
                                        d_idx_sorted_by_tag.data,
                                        d_aggregate_list.data,
                                        d_aggregate_order.data,
                                        block_size,
                                        m_exec_conf->getCachedAllocator());

        if (m_exec_conf->isCUDAErrorCheckingEnabled())
            CHECK_CUDA_ERROR();

        m_tuner_fill->end();
        }

    // distribute aggregates evenly over GPUs
    // NOTE: going forward we could slave the GPU partition of the aggregates
    // to that of the local particles in the ParticleData
    m_gpu_partition = GPUPartition(m_exec_conf->getGPUIds());
    m_gpu_partition.setN(n_local_aggregates);

#if defined(ENABLE_HIP) && defined(__HIP_PLATFORM_NVCC__)
    if (m_exec_conf->allConcurrentManagedAccess())
        {
        auto gpu_map = m_exec_conf->getGPUIds();

        for (unsigned int idev = 0; idev < m_exec_conf->getNumActiveGPUs(); ++idev)
            {
            std::pair<unsigned int, unsigned int> range = m_gpu_partition.getRange(idev);
            unsigned int nelem = range.second - range.first;

            if (nelem == 0)
                continue;

            cudaMemAdvise(m_aggregate_length.get() + range.first,
                          sizeof(unsigned int) * nelem,
                          cudaMemAdviseSetPreferredLocation,
                          gpu_map[idev]);
            cudaMemPrefetchAsync(m_aggregate_length.get() + range.first,
                                 sizeof(unsigned int) * nelem,
                                 gpu_map[idev]);

            if (m_aggregate_indexer.getW() == 0)
                continue;

            cudaMemAdvise(m_aggregate_list.get() + m_aggregate_indexer(0, range.first),
                          sizeof(unsigned int) * nelem * m_aggregate_indexer.getW(),
                          cudaMemAdviseSetPreferredLocation,
                          gpu_map[idev]);
            cudaMemPrefetchAsync(m_aggregate_list.get() + m_aggregate_indexer(0, range.first),
                                 sizeof(unsigned int) * nelem * m_aggregate_indexer.getW(),
                                 gpu_map[idev]);
            }

        for (unsigned int idev = 0; idev < m_exec_conf->getNumActiveGPUs(); ++idev)
            {
            auto range = m_pdata->getGPUPartition().getRange(idev);
            unsigned int nelem = range.second - range.first;

            // skip if no hint set
            if (!nelem)
                continue;

            cudaMemAdvise(m_aggregate_idx.get() + range.first,
                          sizeof(unsigned int) * nelem,
                          cudaMemAdviseSetPreferredLocation,
                          gpu_map[idev]);
            cudaMemPrefetchAsync(m_aggregate_idx.get() + range.first,
                                 sizeof(unsigned int) * nelem,
                                 gpu_map[idev]);
            }

        CHECK_CUDA_ERROR();
        }
#endif
    }
#endif

template<SmoothingKernelType KT_, StateEquationType SET_>
void SPHBaseClassConstraint<KT_, SET_>::initAggregates()
    {
    // return early if no aggregates are defined
    if (!m_n_aggregates_global)
        {
        return;
        }

    m_exec_conf->msg->notice(7) << "SPHBaseClassConstraint initializing aggregate table" << std::endl;

#ifdef ENABLE_HIP
    if (m_exec_conf->isCUDAEnabled())
        {
        initAggregatesGPU();
        return;
        }
#endif

    // construct local aggregate table
    unsigned int nptl_local = m_pdata->getN() + m_pdata->getNGhosts();

    ArrayHandle<unsigned int> h_aggregate_tag(m_aggregate_tag,
                                             access_location::host,
                                             access_mode::read);
    ArrayHandle<unsigned int> h_tag(m_pdata->getTags(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_rtag(m_pdata->getRTags(), access_location::host, access_mode::read);

    std::set<unsigned int> local_aggregate_tags;

    unsigned int n_local_aggregates = 0;

    std::vector<unsigned int> local_aggregate_idx(nptl_local, NO_AGGREGATE);

    // keep track of particle with lowest tag within a aggregate. This is assumed/required to be the
    // central particle for the aggregate.
    std::map<unsigned int, unsigned int> lowest_tag_by_aggregate;

    for (unsigned int particle_index = 0; particle_index < nptl_local; ++particle_index)
        {
        unsigned int tag = h_tag.data[particle_index];
        assert(tag < m_aggregate_tag.getNumElements());

        unsigned int mol_tag = h_aggregate_tag.data[tag];
        if (mol_tag == NO_AGGREGATE)
            {
            continue;
            }

        auto it = lowest_tag_by_aggregate.find(mol_tag);
        unsigned int min_tag = tag;
        if (it != lowest_tag_by_aggregate.end())
            {
            min_tag = std::min(it->second, tag);
            }

        lowest_tag_by_aggregate[mol_tag] = min_tag;
        }

    // sort local aggregates by the index of the smallest particle tag in a aggregate, and sort within
    // the aggregate by particle tag.
    std::map<unsigned int, std::set<unsigned int>> local_aggregates_sorted;

    for (unsigned int particle_index = 0; particle_index < nptl_local; ++particle_index)
        {
        unsigned int tag = h_tag.data[particle_index];
        assert(tag < m_aggregate_tag.getNumElements());

        unsigned int mol_tag = h_aggregate_tag.data[tag];
        if (mol_tag == NO_AGGREGATE)
            {
            continue;
            }

        unsigned int lowest_tag = lowest_tag_by_aggregate[mol_tag];
        unsigned int lowest_idx = h_rtag.data[lowest_tag];
        assert(lowest_idx < m_pdata->getN() + m_pdata->getNGhosts());

        local_aggregates_sorted[lowest_idx].insert(tag);
        }

    n_local_aggregates = static_cast<unsigned int>(local_aggregates_sorted.size());

    m_exec_conf->msg->notice(7) << "SPHBaseClassConstraint: " << n_local_aggregates << " aggregates"
                                << std::endl;

    m_aggregate_length.resize(n_local_aggregates);

    ArrayHandle<unsigned int> h_aggregate_length(m_aggregate_length,
                                                access_location::host,
                                                access_mode::overwrite);

    // reset lengths
    for (unsigned int imol = 0; imol < n_local_aggregates; ++imol)
        {
        h_aggregate_length.data[imol] = 0;
        }

    // count aggregate lengths
    unsigned int i = 0;
    for (auto it = local_aggregates_sorted.begin(); it != local_aggregates_sorted.end(); ++it)
        {
        h_aggregate_length.data[i++] = (unsigned int)it->second.size();
        }

    // find maximum length
    unsigned nmax = 0;
    for (unsigned int imol = 0; imol < n_local_aggregates; ++imol)
        {
        if (h_aggregate_length.data[imol] > nmax)
            {
            nmax = h_aggregate_length.data[imol];
            }
        }

    // set up indexer
    m_aggregate_indexer = Index2D(nmax, n_local_aggregates);

    // resize aggregate list
    m_aggregate_list.resize(m_aggregate_indexer.getNumElements());

    // reset lengths again
    for (unsigned int imol = 0; imol < n_local_aggregates; ++imol)
        {
        h_aggregate_length.data[imol] = 0;
        }

    // resize and reset aggregate lookup to size of local particle data
    m_aggregate_order.resize(m_pdata->getMaxN());
    ArrayHandle<unsigned int> h_aggregate_order(m_aggregate_order,
                                               access_location::host,
                                               access_mode::overwrite);
    memset(h_aggregate_order.data,
           0,
           sizeof(unsigned int) * (m_pdata->getN() + m_pdata->getNGhosts()));

    // resize reverse-lookup
    m_aggregate_idx.resize(nptl_local);

    // fill aggregate list
    ArrayHandle<unsigned int> h_aggregate_list(m_aggregate_list,
                                              access_location::host,
                                              access_mode::overwrite);
    ArrayHandle<unsigned int> h_aggregate_idx(m_aggregate_idx,
                                             access_location::host,
                                             access_mode::overwrite);

    // reset reverse lookup
    memset(h_aggregate_idx.data, 0, sizeof(unsigned int) * nptl_local);

    unsigned int i_mol = 0;
    for (auto it_mol = local_aggregates_sorted.begin(); it_mol != local_aggregates_sorted.end();
         ++it_mol)
        {
        // Since the set is ordered by value, and this orders the particles within the aggregate by
        // tag, and types should have been validated by validateRigidBodies, then this ordering in
        // h_aggregate_order should preserve types even though it is indexed by particle index.
        for (std::set<unsigned int>::iterator it_tag = it_mol->second.begin();
             it_tag != it_mol->second.end();
             ++it_tag)
            {
            unsigned int particle_index = h_rtag.data[*it_tag];
            assert(particle_index < m_pdata->getN() + m_pdata->getNGhosts());
            // Gets the current aggregate index for the particle while incrementing the length of the
            // aggregate.
            unsigned int n = h_aggregate_length.data[i_mol]++;
            h_aggregate_list.data[m_aggregate_indexer(n, i_mol)] = particle_index;
            h_aggregate_idx.data[particle_index] = i_mol;
            h_aggregate_order.data[particle_index] = n;
            }
        i_mol++;
        }
    }


/*! \post Returns a vector of type IDs associated with a particle group
*/
template<SmoothingKernelType KT_, StateEquationType SET_>
void SPHBaseClassConstraint<KT_, SET_>::constructTypeVectors(std::shared_ptr<ParticleGroup> const pgroup,
                                        std::vector<unsigned int> *global_typeids)
    {
    // Clear the input vector
    global_typeids->clear();

    // Position and type read handle
    ArrayHandle<Scalar4> h_pos(m_pdata->getPositions(), access_location::host, access_mode::read);

    // Types vectors
    vector<unsigned int> local_typeids;

    // Loop through all local group particles
    unsigned int group_size = pgroup->getNumMembers();
    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
        {
        // Read particle index
        unsigned int i = pgroup->getMemberIndex(group_idx);

        // Read particle type
        unsigned int type = __scalar_as_int(h_pos.data[i].w);

        // Check if type is already present in vector. If not, push to vector
        if ( find(local_typeids.begin(), local_typeids.end(), type) != local_typeids.end() )
            continue;
        else
            local_typeids.push_back(type);
        }

    // Construct global list in parallel
#ifdef ENABLE_MPI
        // Number of ranks
        unsigned int n_ranks = m_exec_conf->getNRanks();

        // Gather results from all processors for fluid types
        vector< vector<unsigned int> > typeids_proc(n_ranks);
        all_gather_v(local_typeids, typeids_proc, m_exec_conf->getMPICommunicator());

        // Combine all types into an ordered set
        set<unsigned int> typeids_set;
        for (unsigned int irank = 0; irank < n_ranks; ++irank)
            {
            typeids_set.insert(typeids_proc[irank].begin(), typeids_proc[irank].end());
            }

        // Store in vector
        for (set<unsigned int>::iterator it=typeids_set.begin(); it!=typeids_set.end(); ++it)
            global_typeids->push_back(*it);

#else
        // Copy data
        for (unsigned int i=0;i<local_typeids.size();i++)
            global_typeids->push_back(local_typeids[i]);
#endif

    m_exec_conf->msg->notice(7) << "Available types of requested group: " << std::endl;
    for (unsigned int i=0;i<global_typeids->size();i++)
        m_exec_conf->msg->notice(7) << "Typenumber " << i << ": " << (*global_typeids)[i] << std::endl;
    }

/*! \post Return current body force
*/
template<SmoothingKernelType KT_, StateEquationType SET_>
Scalar3 SPHBaseClassConstraint<KT_, SET_>::getAcceleration(uint64_t timestep)
    {
    if ( m_damptime > 0 && timestep < m_damptime )
        return m_bodyforce * Scalar(0.5)*(sin(M_PI*(-Scalar(0.5)+Scalar(timestep)/Scalar(m_damptime)))+Scalar(1));
    else
        return m_bodyforce;
    }


/*! \post Body forces mass*body acceleration are applied
 * damp factor for 5000 steps
1 +-----------------------------------------------------------------------------------------+
      |                                                                               *****     |
      |                                                                            ****         |
      |                                                                         ****            |
      |                                                                       ***               |
  0.9 |                                                                     ***                 |
      |                                                                   ***                   |
      |                                                                 ***                     |
      |                                                                **                       |
  0.8 |                                                              ***                        |
      |                                                            ***                          |
      |                                                           **                            |
      |                                                          **                             |
      |                                                        ***                              |
  0.7 |                                                       **                                |
      |                                                      **                                 |
      |                                                    ***                                  |
      |                                                   **                                    |
  0.6 |                                                  **                                     |
      |                                                ***                                      |
      |                                               **                                        |
      |                                              **                                         |
      |                                             **                                          |
  0.5 |                                           ***                                           |
      |                                          **                                             |
      |                                         **                                              |
      |                                        **                                               |
      |                                      ***                                                |
  0.4 |                                     **                                                  |
      |                                    **                                                   |
      |                                  ***                                                    |
      |                                 **                                                      |
  0.3 |                                **                                                       |
      |                              ***                                                        |
      |                             **                                                          |
      |                            **                                                           |
      |                          ***                                                            |
  0.2 |                        ***                                                              |
      |                       **                                                                |
      |                     ***                                                                 |
      |                   ***                                                                   |
  0.1 |                 ***                                                                     |
      |               ***                                                                       |
      |            ****                                                                         |
      |         ****                                                                            |
      |     *****                                                                               |
    0 +-----------------------------------------------------------------------------------------+
      0       500      1000     1500     2000     2500     3000     3500     4000     4500     5000
*/



template<SmoothingKernelType KT_, StateEquationType SET_>
void SPHBaseClassConstraint<KT_, SET_>::applyBodyForce(uint64_t timestep, std::shared_ptr<ParticleGroup> pgroup)
    {
    if ( m_body_acceleration )
        {

        Scalar damp = Scalar(1);
        if ( m_damptime > 0 && timestep < m_damptime )
            damp = Scalar(0.5)*(sin(M_PI*(-Scalar(0.5)+Scalar(timestep)/Scalar(m_damptime)))+Scalar(1));
        Scalar3 bforce = m_bodyforce * damp;

        m_exec_conf->msg->notice(7) << "Computing SPHBaseClassConstraint::applyBodyForce with damp factor " << damp << std::endl;

        // Grab handles for particle data
        ArrayHandle<Scalar4> h_force(m_force,access_location::host, access_mode::readwrite);
        ArrayHandle<Scalar4> h_velocity(m_pdata->getVelocities(), access_location::host, access_mode::read);

        // for each particle in given group
        unsigned int group_size = pgroup->getNumMembers();
        for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
            {
            // Read particle index
            unsigned int i = pgroup->getMemberIndex(group_idx);

            // Read particle mass
            Scalar mi = h_velocity.data[i].w;

            // Add contribution to force
            h_force.data[i].x += bforce.x*mi;
            h_force.data[i].y += bforce.y*mi;
            h_force.data[i].z += bforce.z*mi;
            }
        }
    }

// #ifdef ENABLE_HIP
// template<SmoothingKernelType KT_,StateEquationType SET_>
// void SPHBaseClassConstraint<KT_, SET_>::applyBodyForceGPU(uint64_t timestep, std::shared_ptr<ParticleGroup> pgroup)
//     {
//     if ( this->m_body_acceleration )
//         {

//         Scalar damp = Scalar(1);
//         if ( this->m_damptime > 0 && timestep < this->m_damptime )
//             damp = Scalar(0.5)*(sin(M_PI*(-Scalar(0.5)+Scalar(timestep)/Scalar(this->m_damptime)))+Scalar(1));
//         Scalar3 bforce = this->m_bodyforce * damp;

//         this->m_exec_conf->msg->notice(7) << "Computing SPHBaseClassConstraintGPU::applyBodyForce with damp factor " << damp << std::endl;

//         // Grab handles for particle data
//         ArrayHandle<Scalar4> d_force(this->m_force,access_location::device, access_mode::readwrite);
//         ArrayHandle<Scalar4> d_velocity(this->m_pdata->getVelocities(), access_location::device, access_mode::read);

//         // for each particle in given group
//         ArrayHandle< unsigned int > d_index_array(pgroup->getIndexArray(), access_location::device, access_mode::read);
//         unsigned int group_size = pgroup->getNumMembers();

//         gpu_sphbase_apply_body_force(d_force.data,
//                                     d_velocity.data,
//                                     bforce,
//                                     d_index_array.data,
//                                     group_size
//         );
//         }
//     }
// #endif

/*! \post Set acceleration components
 */
template<SmoothingKernelType KT_, StateEquationType SET_>
void SPHBaseClassConstraint<KT_, SET_>::setAcceleration(Scalar gx, Scalar gy, Scalar gz, unsigned int damptime)
    {
    m_bodyforce = make_scalar3(gx,gy,gz);
    if ( gx != 0 || gy != 0 || gz != 0 )
        {
        m_body_acceleration = true;
        m_damptime = damptime;
        m_exec_conf->msg->notice(7) << "Non-Zero volumetric acceleration: gx:" << gx << " gy:" << gy << " gz:" << gz << std::endl;
        }
    else
        {
        m_exec_conf->msg->notice(7) << "No volumetric acceleration" << std::endl;
        }
    }

namespace detail {

template<SmoothingKernelType KT_, StateEquationType SET_>
void export_SPHBaseClassConstraint(pybind11::module& m, std::string name)
{   
    // using Class = SPHBaseClassConstraint<KT_, SET_>; 
    pybind11::class_<SPHBaseClassConstraint<KT_, SET_>, ForceCompute, std::shared_ptr<SPHBaseClassConstraint<KT_, SET_>>>
        sphbaseclass(m, name.c_str());
    sphbaseclass
        .def(pybind11::init<std::shared_ptr<SystemDefinition>,
                            std::shared_ptr<SmoothingKernel<KT_> >,
                           std::shared_ptr<StateEquation<SET_> >,
                           std::shared_ptr<nsearch::NeighborList>>())
        .def("constructTypeVectors", &SPHBaseClassConstraint<KT_, SET_>::constructTypeVectors)
        .def("getAcceleration", &SPHBaseClassConstraint<KT_, SET_>::getAcceleration)
        .def("applyBodyForce", &SPHBaseClassConstraint<KT_, SET_>::applyBodyForce)
        // .def("applyBodyForceGPU", &SPHBaseClassConstraint<KT_, SET_>::applyBodyForceGPU)
        .def("setAcceleration", &SPHBaseClassConstraint<KT_, SET_>::setAcceleration);
}


// void export_DensityMethodConstraint(pybind11::module& m)
// {
//     pybind11::enum_<DensityMethodConstraint>(m, "PhaseFlowDensityMethod")
//         .value("DENSITYSUMMATION", DensityMethod::DENSITYSUMMATION)
//         .value("DENSITYCONTINUITY", DensityMethod::DENSITYCONTINUITY)
//         ;
// }

// void export_ViscosityMethodConstraint(pybind11::module& m)
// {
//     pybind11::enum_<ViscosityMethodConstraint>(m, "PhaseFlowViscosityMethod")
//         .value("HARMONICAVERAGE", ViscosityMethod::HARMONICAVERAGE)
//         ;
// }



} // end namespace detail 

//! Explicit template instantiations
template class PYBIND11_EXPORT SPHBaseClassConstraint<wendlandc2, linear>;
template class PYBIND11_EXPORT SPHBaseClassConstraint<wendlandc2, tait>;
template class PYBIND11_EXPORT SPHBaseClassConstraint<wendlandc4, linear>;
template class PYBIND11_EXPORT SPHBaseClassConstraint<wendlandc4, tait>;
template class PYBIND11_EXPORT SPHBaseClassConstraint<wendlandc6, linear>;
template class PYBIND11_EXPORT SPHBaseClassConstraint<wendlandc6, tait>;
template class PYBIND11_EXPORT SPHBaseClassConstraint<quintic, linear>;
template class PYBIND11_EXPORT SPHBaseClassConstraint<quintic, tait>;
template class PYBIND11_EXPORT SPHBaseClassConstraint<cubicspline, linear>;
template class PYBIND11_EXPORT SPHBaseClassConstraint<cubicspline, tait>;


namespace detail
{

    template void export_SPHBaseClassConstraint<wendlandc2, linear>(pybind11::module& m, std::string name = "SPHBaseClassConstraint_WC2_L");
    template void export_SPHBaseClassConstraint<wendlandc2, tait>(pybind11::module& m, std::string name = "SPHBaseClassConstraint_WC2_T");
    template void export_SPHBaseClassConstraint<wendlandc4, linear>(pybind11::module& m, std::string name = "SPHBaseClassConstraint_WC4_L");
    template void export_SPHBaseClassConstraint<wendlandc4, tait>(pybind11::module& m, std::string name = "SPHBaseClassConstraint_WC4_T");
    template void export_SPHBaseClassConstraint<wendlandc6, linear>(pybind11::module& m, std::string name = "SPHBaseClassConstraint_WC6_L");
    template void export_SPHBaseClassConstraint<wendlandc6, tait>(pybind11::module& m, std::string name = "SPHBaseClassConstraint_WC6_T");
    template void export_SPHBaseClassConstraint<quintic, linear>(pybind11::module& m, std::string name = "SPHBaseClassConstraint_Q_L");
    template void export_SPHBaseClassConstraint<quintic, tait>(pybind11::module& m, std::string name = "SPHBaseClassConstraint_Q_T");
    template void export_SPHBaseClassConstraint<cubicspline, linear>(pybind11::module& m, std::string name = "SPHBaseClassConstraint_CS_L");
    template void export_SPHBaseClassConstraint<cubicspline, tait>(pybind11::module& m, std::string name = "SPHBaseClassConstraint_CS_T");

} // end namespace detail
} // end namespace sph
} // end namespace hoomd