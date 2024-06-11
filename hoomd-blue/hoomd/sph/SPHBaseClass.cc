/* ---------------------------------------------------------
maintainer: dkrach, david.krach@mib.uni-stuttgart.de
----------------------------------------------------------*/

#include "SPHBaseClass.h"
#include "hoomd/nsearch/NeighborList.h"

#ifdef ENABLE_HIP
#include "SPHBaseClass.cuh"
#endif


#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>

#include <stdexcept>
#include <utility>
#include <set>
#include <algorithm>

using namespace std;

namespace hoomd
{
namespace sph
{
/*! Constructor
*/
template<SmoothingKernelType KT_, StateEquationType SET_>
SPHBaseClass<KT_, SET_>::SPHBaseClass(std::shared_ptr<SystemDefinition> sysdef,
                           std::shared_ptr<SmoothingKernel<KT_> > skernel,
                           std::shared_ptr<StateEquation<SET_> > eos,
                           std::shared_ptr<nsearch::NeighborList> nlist)
    : ForceCompute(sysdef), m_skernel(skernel), m_eos(eos), m_nlist(nlist), m_typpair_idx(m_pdata->getNTypes())
      {
        m_exec_conf->msg->notice(5) << "Constructing SPHBaseClass" << std::endl;

        // Sanity checks
        assert(m_nlist);

        assert(m_skernel);
        assert(m_eos);

        // Set default variables
        m_bodyforce = make_scalar3(Scalar(0), Scalar(0), Scalar(0));
        m_damptime = 0;
        m_body_acceleration = false;
      }

/*! Destructor
*/
template<SmoothingKernelType KT_, StateEquationType SET_>
SPHBaseClass<KT_, SET_>::~SPHBaseClass()
    {
    m_exec_conf->msg->notice(5) << "Destroying SPHBaseClass" << std::endl;
    }

/*! \post Returns a vector of type IDs associated with a particle group
*/
template<SmoothingKernelType KT_, StateEquationType SET_>
void SPHBaseClass<KT_, SET_>::constructTypeVectors(std::shared_ptr<ParticleGroup> const pgroup,
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
Scalar3 SPHBaseClass<KT_, SET_>::getAcceleration(uint64_t timestep)
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
void SPHBaseClass<KT_, SET_>::applyBodyForce(uint64_t timestep, std::shared_ptr<ParticleGroup> pgroup)
    {
    if ( m_body_acceleration )
        {

        Scalar damp = Scalar(1);
        if ( m_damptime > 0 && timestep < m_damptime )
            damp = Scalar(0.5)*(sin(M_PI*(-Scalar(0.5)+Scalar(timestep)/Scalar(m_damptime)))+Scalar(1));
        Scalar3 bforce = m_bodyforce * damp;

        m_exec_conf->msg->notice(7) << "Computing SPHBaseClass::applyBodyForce with damp factor " << damp << std::endl;

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
// void SPHBaseClass<KT_, SET_>::applyBodyForceGPU(uint64_t timestep, std::shared_ptr<ParticleGroup> pgroup)
//     {
//     if ( this->m_body_acceleration )
//         {

//         Scalar damp = Scalar(1);
//         if ( this->m_damptime > 0 && timestep < this->m_damptime )
//             damp = Scalar(0.5)*(sin(M_PI*(-Scalar(0.5)+Scalar(timestep)/Scalar(this->m_damptime)))+Scalar(1));
//         Scalar3 bforce = this->m_bodyforce * damp;

//         this->m_exec_conf->msg->notice(7) << "Computing SPHBaseClassGPU::applyBodyForce with damp factor " << damp << std::endl;

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
void SPHBaseClass<KT_, SET_>::setAcceleration(Scalar gx, Scalar gy, Scalar gz, unsigned int damptime)
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
void export_SPHBaseClass(pybind11::module& m, std::string name)
{   
    // using Class = SPHBaseClass<KT_, SET_>; 
    pybind11::class_<SPHBaseClass<KT_, SET_>, ForceCompute, std::shared_ptr<SPHBaseClass<KT_, SET_>>>
        sphbaseclass(m, name.c_str());
    sphbaseclass
        .def(pybind11::init<std::shared_ptr<SystemDefinition>,
                            std::shared_ptr<SmoothingKernel<KT_> >,
                           std::shared_ptr<StateEquation<SET_> >,
                           std::shared_ptr<nsearch::NeighborList>>())
        .def("constructTypeVectors", &SPHBaseClass<KT_, SET_>::constructTypeVectors)
        .def("getAcceleration", &SPHBaseClass<KT_, SET_>::getAcceleration)
        .def("applyBodyForce", &SPHBaseClass<KT_, SET_>::applyBodyForce)
        // .def("applyBodyForceGPU", &SPHBaseClass<KT_, SET_>::applyBodyForceGPU)
        .def("setAcceleration", &SPHBaseClass<KT_, SET_>::setAcceleration);
}


void export_DensityMethod(pybind11::module& m)
{
    pybind11::enum_<DensityMethod>(m, "PhaseFlowDensityMethod")
        .value("DENSITYSUMMATION", DensityMethod::DENSITYSUMMATION)
        .value("DENSITYCONTINUITY", DensityMethod::DENSITYCONTINUITY)
        ;
}

void export_ViscosityMethod(pybind11::module& m)
{
    pybind11::enum_<ViscosityMethod>(m, "PhaseFlowViscosityMethod")
        .value("HARMONICAVERAGE", ViscosityMethod::HARMONICAVERAGE)
        ;
}

void export_MaterialModel(pybind11::module& m)
{
    pybind11::enum_<MaterialModel>(m, "PhaseFlowMaterialModel")
        .value("REGULARIZEDBINGHAM", MaterialModel::REGULARIZEDBINGHAM)
        .value("BIVISCOUS", MaterialModel::BIVISCOUS)
        ;
}


} // end namespace detail 

//! Explicit template instantiations
template class PYBIND11_EXPORT SPHBaseClass<wendlandc2, linear>;
template class PYBIND11_EXPORT SPHBaseClass<wendlandc2, tait>;
template class PYBIND11_EXPORT SPHBaseClass<wendlandc4, linear>;
template class PYBIND11_EXPORT SPHBaseClass<wendlandc4, tait>;
template class PYBIND11_EXPORT SPHBaseClass<wendlandc6, linear>;
template class PYBIND11_EXPORT SPHBaseClass<wendlandc6, tait>;
template class PYBIND11_EXPORT SPHBaseClass<quintic, linear>;
template class PYBIND11_EXPORT SPHBaseClass<quintic, tait>;
template class PYBIND11_EXPORT SPHBaseClass<cubicspline, linear>;
template class PYBIND11_EXPORT SPHBaseClass<cubicspline, tait>;


namespace detail
{

    template void export_SPHBaseClass<wendlandc2, linear>(pybind11::module& m, std::string name = "SPHBaseClass_WC2_L");
    template void export_SPHBaseClass<wendlandc2, tait>(pybind11::module& m, std::string name = "SPHBaseClass_WC2_T");
    template void export_SPHBaseClass<wendlandc4, linear>(pybind11::module& m, std::string name = "SPHBaseClass_WC4_L");
    template void export_SPHBaseClass<wendlandc4, tait>(pybind11::module& m, std::string name = "SPHBaseClass_WC4_T");
    template void export_SPHBaseClass<wendlandc6, linear>(pybind11::module& m, std::string name = "SPHBaseClass_WC6_L");
    template void export_SPHBaseClass<wendlandc6, tait>(pybind11::module& m, std::string name = "SPHBaseClass_WC6_T");
    template void export_SPHBaseClass<quintic, linear>(pybind11::module& m, std::string name = "SPHBaseClass_Q_L");
    template void export_SPHBaseClass<quintic, tait>(pybind11::module& m, std::string name = "SPHBaseClass_Q_T");
    template void export_SPHBaseClass<cubicspline, linear>(pybind11::module& m, std::string name = "SPHBaseClass_CS_L");
    template void export_SPHBaseClass<cubicspline, tait>(pybind11::module& m, std::string name = "SPHBaseClass_CS_T");

} // end namespace detail
} // end namespace sph
} // end namespace hoomd