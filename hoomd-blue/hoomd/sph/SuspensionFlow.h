/* ---------------------------------------------------------
maintainer: dkrach, david.krach@mib.uni-stuttgart.de
----------------------------------------------------------*/

#include "hoomd/Compute.h"
#include "hoomd/Index1D.h"
#include "hoomd/ParticleGroup.h"
#include "hoomd/ForceCompute.h"
#include "hoomd/Integrator.h"
#include "hoomd/nsearch/NeighborList.h"

#ifdef ENABLE_MPI
#include "hoomd/Communicator.h"
#include "hoomd/HOOMDMPI.h"
#endif

#include <memory>
#include <vector>
#include <stdexcept>
#include <utility>
#include <set>
#include <algorithm>

#ifdef ENABLE_HIP
#include "hoomd/GPUPartition.cuh"
#endif

#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>

#include "SmoothingKernel.h"
#include "StateEquations.h"
#include "SPHBaseClassConstraint.h"
#include "SolidFluidTypeBit.h"

#include "EvaluationMethodDefinition.h"


/*! \file SuspensionFlow.h
    \brief Contains code for the Quasi-incompressible Navier-Stokes solver
          for Single-phase flow with solid aggregates
*/

#ifdef __HIPCC__
#error This header cannot be compiled by nvcc
#endif

#ifndef __SuspensionFlow_H__
#define __SuspensionFlow_H__


namespace hoomd 
{
namespace sph
{

//! Computes SuspensionFlow forces on each particle
/*!
*/
template<SmoothingKernelType KT_,StateEquationType SET_>
class PYBIND11_EXPORT SuspensionFlow : public SPHBaseClassConstraint<KT_, SET_>
    {
    public:

        //! Constructor
        SuspensionFlow(std::shared_ptr<SystemDefinition> sysdef,
                        std::shared_ptr<SmoothingKernel<KT_> > skernel,
                        std::shared_ptr<StateEquation<SET_> > equationofstate,
                        std::shared_ptr<nsearch::NeighborList> nlist,
                        std::shared_ptr<ParticleGroup> fluidgroup,
                        std::shared_ptr<ParticleGroup> solidgroup,
                        std::shared_ptr<ParticleGroup> aggregategroup,
                        DensityMethod   mdensitymethod=DENSITYSUMMATION,
                        ViscosityMethod mviscositymethod=HARMONICAVERAGE);

        //! Destructor
        virtual ~SuspensionFlow();

        //! Set the rcut for a single type pair
        virtual void setRcut(unsigned int typ1, unsigned int typ2, Scalar rcut);

        /// Set the rcut for a single type pair using a tuple of strings
        virtual void setRCutPython(pybind11::tuple types, Scalar r_cut);

        /// Validate that types are within Ntypes
        void validateTypes(unsigned int typ1, unsigned int typ2, std::string action);

        //! Returns a list of log quantities this compute calculates
        virtual std::vector<double> getProvidedTimestepQuantities(uint64_t timestep);

        /*! Set the parameters
         * \param mu Dynamic viscosity
         */
        virtual void setParams(Scalar mu, Scalar rhoS, Scalar f0);

        //! Getter and Setter methods for density method
        DensityMethod getDensityMethod()
            {
            return m_density_method;
            }
        void setDensityMethod(DensityMethod densitymethod)
            {
            m_density_method = densitymethod;
            }

        //! Getter and Setter methods for viscosity method
        ViscosityMethod getViscosityMethod()
            {
            return m_viscosity_method;
            }
        void setViscosityMethod(ViscosityMethod viscositymethod)
            {
            m_viscosity_method = viscositymethod;
            }

        // Set constant smoothing length option to true for faster computation
        void setConstSmoothingLength(Scalar h)
            {
            m_const_slength = true;
            // constant slength in most cases
            m_ch = h;
            m_rcut = m_kappa * m_ch;
            // squared cutoff radius to compare with distance dot(dx, dx)
            m_rcutsq = m_rcut * m_rcut;  
            }

        /*! Set compute solid forces option to true. This is necessary if suspended object
         *  are present or if solid drag forces are to be evaluated.
         */
        void computeSolidForces()
            {
            m_compute_solid_forces = true;
            }

        // /*! Set compute wall forces option to true. This is necessary if 
        //  *  if there are walls (BC) defined, then solid body IDs start with 1 instead of 0
        //  */
        // void computeWallForces()
        //     {
        //     m_walls = true;
        //     }

        /*! Turn Monaghan type artificial viscosity option on.
         * \param alpha Volumetric diffusion coefficient for artificial viscosity operator
         * \param beta Shock diffusion coefficient for artificial viscosity operator
         */
        void activateArtificialViscosity(Scalar alpha, Scalar beta)
            {
            m_artificial_viscosity = true;
            m_avalpha = alpha;
            m_avbeta = beta;
            }

        /*! Turn Monaghan type artificial viscosity option off.
         */
        void deactivateArtificialViscosity()
            {
            m_artificial_viscosity = false;
            }

        /*! Turn Molteni type density diffusion option on.
         * \param ddiff Diffusion coefficient for artificial density diffusion operator
         */
        void activateDensityDiffusion(Scalar ddiff)
            {
            m_density_diffusion = true;
            m_ddiff = ddiff;
            }

        /*! Turn Molteni type density diffusion off.
         */
        void deactivateDensityDiffusion()
            {
            m_density_diffusion = false;
            }

        /*! Turn Shepard type density reinitialization on
         * \param shepardfreq Number of timesteps the renormalization is to be applied
         */
        void activateShepardRenormalization(unsigned int shepardfreq);

        /*! Turn Shepard type density reinitialization off.
         */
        void deactivateShepardRenormalization()
            {
            m_shepard_renormalization = false;
            }

        //! Computes forces
        void computeForces(uint64_t timestep);

        //! Set the coordinates for the template for a rigid body of type typeid
        /*! \param body_type The type of rigid body
         * \param type Types of the constituent particles
         * \param pos Relative positions of the constituent particles
         * \param orientation Orientations of the constituent particles
         */
        virtual void setRigidParams(unsigned int body_typeid,
                              std::vector<unsigned int>& type,
                              std::vector<Scalar3>& pos);
                              // std::vector<Scalar4>& orientation,
                              // std::vector<Scalar>& charge,
                              // std::vector<Scalar>& diameter);

        // TODO: Keine Ahnung ob noetig
        // //! Returns true because we compute the torque on the central particle
        // virtual bool isAnisotropic()
        //     {
        //     return true;
        //     }

        /// Update the constituent particles of a composite particle using the position, velocity
        /// and orientation of the central particle.
        virtual void updateCompositeParticles(uint64_t timestep);

        /// Validate rigid body constituent particles. The method purposely does not check
        /// positions or orientation.
        virtual void validateRigidBodies();

        //! Create rigid body constituent particles
        virtual void createRigidBodies();

        /// Construct from a Python dictionary
        void setBody(std::string typ, pybind11::object v)
            {
            if (v.is_none())
                {
                return;
                }
            pybind11::list types = v["constituent_types"];
            pybind11::list positions = v["positions"];
            // pybind11::list orientations = v["orientations"];
            // pybind11::list charges = v["charges"];
            // pybind11::list diameters = v["diameters"];
            auto N = pybind11::len(positions);
            // Ensure proper list lengths
            // for (const auto& list : {types, orientations, charges, diameters})
            //     {
            //     if (pybind11::len(list) != N)
            //         {
            //         throw std::runtime_error("All attributes of a rigid body must be the same length.");
            //         }
            //     }
            // for (const auto& list : {types})
            //     {
            //     if (pybind11::len(list) != N)
            //         {
            //         throw std::runtime_error("All attributes of a rigid body must be the same length.");
            //         }
            //     }

            // extract the data from the python lists
            std::vector<Scalar3> pos_vector;
            // std::vector<Scalar4> orientation_vector;
            // std::vector<Scalar> charge_vector;
            // std::vector<Scalar> diameter_vector;
            std::vector<unsigned int> type_vector;

            for (size_t i(0); i < N; ++i)
                {
                pybind11::tuple position_i(positions[i]);
                pos_vector.emplace_back(make_scalar3(position_i[0].cast<Scalar>(),
                                                     position_i[1].cast<Scalar>(),
                                                     position_i[2].cast<Scalar>()));

                // pybind11::tuple orientation_i(orientations[i]);
                // orientation_vector.emplace_back(make_scalar4(orientation_i[0].cast<Scalar>(),
                //                                              orientation_i[1].cast<Scalar>(),
                //                                              orientation_i[2].cast<Scalar>(),
                //                                              orientation_i[3].cast<Scalar>()));

                // charge_vector.emplace_back(charges[i].cast<Scalar>());
                // diameter_vector.emplace_back(diameters[i].cast<Scalar>());
                type_vector.emplace_back(m_pdata->getTypeByName(types[i].cast<std::string>()));
                }

            setRigidParams(m_pdata->getTypeByName(typ),
                     type_vector,
                     pos_vector);
                     // orientation_vector,
                     // charge_vector,
                     // diameter_vector);
            }

        /// Convert parameters to a python dictionary
        pybind11::object getBody(std::string body_type)
            {
            auto body_type_id = m_pdata->getTypeByName(body_type);
            ArrayHandle<unsigned int> h_body_len(m_body_len,
                                                 access_location::host,
                                                 access_mode::readwrite);
            unsigned int N = h_body_len.data[body_type_id];
            if (N == 0)
                {
                return pybind11::none();
                }
            ArrayHandle<Scalar3> h_body_pos(m_body_pos, access_location::host, access_mode::read);
            // ArrayHandle<Scalar4> h_body_orientation(m_body_orientation,
            //                                         access_location::host,
            //                                         access_mode::read);
            ArrayHandle<unsigned int> h_body_types(m_body_types,
                                                   access_location::host,
                                                   access_mode::read);

            pybind11::list positions;
            // pybind11::list orientations;
            pybind11::list types;
            // pybind11::list charges;
            // pybind11::list diameters;

            for (unsigned int i = 0; i < N; i++)
                {
                auto index = m_body_idx(body_type_id, i);
                positions.append(pybind11::make_tuple(static_cast<Scalar>(h_body_pos.data[index].x),
                                                      static_cast<Scalar>(h_body_pos.data[index].y),
                                                      static_cast<Scalar>(h_body_pos.data[index].z)));
                // orientations.append(
                //     pybind11::make_tuple(static_cast<Scalar>(h_body_orientation.data[index].x),
                //                          static_cast<Scalar>(h_body_orientation.data[index].y),
                //                          static_cast<Scalar>(h_body_orientation.data[index].z),
                //                          static_cast<Scalar>(h_body_orientation.data[index].w)));
                types.append(m_pdata->getNameByType(h_body_types.data[index]));
                // charges.append(m_body_charge[body_type_id][i]);
                // diameters.append(m_body_diameter[body_type_id][i]);
                }
            pybind11::dict v;
            v["constituent_types"] = types;
            v["positions"] = positions;
            // v["orientations"] = orientations;
            // v["charges"] = charges;
            // v["diameters"] = diameters;
            return std::move(v);
            }

        /// Get the number of free particles (global)
        unsigned int getNFreeParticlesGlobal()
            {
            return m_n_free_particles_global;
            }

    #ifdef ENABLE_MPI
        /// The system's communicator.
        std::shared_ptr<Communicator> m_comm;
    #endif

    #ifdef ENABLE_MPI
        //! Get requested ghost communication flags
        virtual CommFlags getRequestedCommFlags(uint64_t timestep)
            {
            // Request communication of all field required during ForceCompute
            CommFlags flags(0);
            flags[comm_flag::net_force] = 0;
            flags[comm_flag::position] = 1; // Stores position and type
            flags[comm_flag::velocity] = 1; // Stores velocity and mass
            flags[comm_flag::density] = 1; // Stores density 
            flags[comm_flag::pressure] = 1; // Stores pressure
            flags[comm_flag::energy] = 0; // Stores density and pressure
            flags[comm_flag::auxiliary1] = 1; // Stores fictitious velocity
            flags[comm_flag::auxiliary2] = 1; // Stores orientation of center particle
//            flags[comm_flag::auxiliary3] = 1; // Stores contact force
            flags[comm_flag::slength] = 1; // Stores smoothing length TODO is this needed
            // Add flags requested by base class
            flags |= ForceCompute::getRequestedCommFlags(timestep);
            return flags;
            }
    #endif

        //! Returns true because we compute dpe array content
        virtual bool ComputesDPE()
            {
            return true;
            }

    protected:

        // Shared pointers
        std::shared_ptr<ParticleGroup> m_fluidgroup; //!< Group of fluid particles
        std::shared_ptr<ParticleGroup> m_solidgroup; //!< Group of solid particles (BC)
        std::shared_ptr<ParticleGroup> m_aggregategroup; //!< Group of solid aggregate particles


        /// r_cut (not squared) given to the neighbor list
        std::shared_ptr<GlobalArray<Scalar>> m_r_cut_nlist;


        // Index for rcut pair info -> nlist
        Index2D m_typpair_idx;        //!< Helper class for indexing per type pair arrays

        // std::shared_ptr<nsearch::NeighborList> m_nlist; //!< the neighborlist to use for the computation

        // Model parameters
        Scalar m_ch; //!< Smoothing length to use if constant for all particles
        Scalar m_rcut; //!< Cut-off length to use if constant for all particles
        Scalar m_rcutsq; //!< Square cut-off length to use if constant for all particles
        DensityMethod m_density_method; //!< Density approach to use
        ViscosityMethod m_viscosity_method; //!< Viscosity approach to use

        // Physical variables
        Scalar m_rho0; //!< Rest density (Read from equation of state class)
        Scalar m_rhoS; //!< Rest density for solids (Read from equation of state class)
        Scalar m_c; //!< Speed of sound (Read from equation of state class)
        Scalar m_kappa; //!< Kernel scaling factor (Read from kernel class)
        Scalar m_mu; //!< Viscosity ( Must be set by user )
        Scalar m_avalpha; //!< Volumetric diffusion coefficient for artificial viscosity operator
        Scalar m_avbeta; //!< Shock diffusion coefficient for artificial viscosity operator
        Scalar m_ddiff; //!< Diffusion coefficient for Molteni type density diffusion
        unsigned int m_shepardfreq; //!< Time step frequency for Shepard reinitialization

        // Auxiliary variables
        std::vector<unsigned int> m_fluidtypes; //!< Fluid type numbers
        std::vector<unsigned int> m_solidtypes; //!< Solid wall type numbers
        std::vector<unsigned int> m_aggregatetypes; //!< Solid aggregate type numbers
        GPUArray<unsigned int> m_type_property_map; //!< to check if a particle type is soli d or fluid

        // Flags
        bool m_const_slength; //!< True if using constant smoothing length
        bool m_compute_solid_forces; //!< Set to true if forces acting on solid particle are to be computed
        bool m_artificial_viscosity; //!< Set to true if Monaghan type artificial viscosity is to be used
        bool m_density_diffusion; //!< Set to true if Molteni type density diffusion is to be used
        bool m_shepard_renormalization; //!< Set to true if Shepard type density reinitialization is to be used
        bool m_params_set; //!< True if parameters are set
        bool m_solid_removed; //!< True if solid Particles have been marked to remove 
        bool m_walls; //!< True if there are walls defined, solid body IDs start with 1 instead of 0


        bool m_bodies_changed;          //!< True if constituent particles have changed
        bool m_particles_added_removed; //!< True if particles have been added or removed

        /// The number of free particles in the simulation box.
        unsigned int m_n_free_particles_global;

        GlobalArray<unsigned int> m_body_types;  //!< Constituent particle types per type id (2D)
        GlobalArray<Scalar3> m_body_pos;         //!< Constituent particle offsets per type id (2D)
        // GlobalArray<Scalar4> m_body_orientation; //!< Constituent particle orientations per type id (2D)
        GlobalArray<unsigned int> m_body_len;    //!< Length of body per type id

        // std::vector<std::vector<Scalar>> m_body_charge;   //!< Constituent particle charges
        // std::vector<std::vector<Scalar>> m_body_diameter; //!< Constituent particle diameters
        Index2D m_body_idx;                               //!< Indexer for body parameters

        // std::vector<Scalar> m_d_max;       //!< Maximum body diameter per constituent particle type
        // std::vector<bool> m_d_max_changed; //!< True if maximum body diameter changed (per type)

        // Log parameters
        uint64_t m_log_computed_last_timestep; //!< Last time step where log quantities were computed

        // Timestep parameters
        std::vector<double> m_timestep_list = std::vector<double>(7);  //!< Cache all generated timestep quantities names

        void mark_solid_particles_toremove(uint64_t timestep);

        /*! Helper function to compute particle number density
         * \post For fluid particles, compute number density. For solid particles,
                 compute fluid normalization constant.
         */
        void compute_ndensity(uint64_t timestep);
        
        void compute_particlenumberdensity(uint64_t timestep);

        /*! Helper function to compute particle pressures
         *  \post Pressure of fluid particle computed
         */
        void compute_pressure(uint64_t timestep);

        /*! Helper function to compute fictitious solid particle properties (pressures and velocities)
        * \pre Ghost particle number densities (i.e. density array) must be up-to-date
        * \pre Solid normalization constant \sum_j w_ij must be computed and stored in density array
        * \post Fictitious particle properties are computed and stored in aux1 array
        */
        void compute_noslip(uint64_t timestep);

        /*! Helper function to apply Shepard density filter
        * \post Fluid particle densities are recomputed based on the Shepard renormalization
        */
        void renormalize_density(uint64_t timestep);

        /*! Helper function where the actual force computation takes place
         * \pre Number densities and fictitious solid particle properties must be up-to-date
         * \post h_force stores forces acting on fluid particles and .w component stores rate of change of density
         */
        void forcecomputation(uint64_t timestep);

        /*! Helper function to set communication flags and update ghosts densities
        * \param timestep The time step
        * \post Ghost particle density array is up-to-date
        */
        void update_ghost_density(uint64_t timestep);

        /*! Helper function to set communication flags and update ghosts densities and pressures
        * \param timestep The time step
        * \post Ghost particle density and pressue array is up-to-date
        */
        void update_ghost_density_pressure(uint64_t timestep);

        /*! Helper function to set communication flags and update ghosts auxiliary array 1
        * \param timestep The time step
        * \post Ghost particle auxiliary array 1 is up-to-date
        */
        void update_ghost_aux1(uint64_t timestep);

        /*! Helper function to set communication flags and update ghosts auxiliary array 1
        * \param timestep The time step
        * \post Ghost particle auxiliary array 1 is up-to-date
        */
        void update_ghost_aux2(uint64_t timestep);

        // eventually add:
        /* - equivalentradii - lubrication -centerofmasses - repulsiveforce -contactforce
        */

        //! Method to be called when particles are added or removed
        void slotPtlsAddedRemoved()
            {
            m_particles_added_removed = true;
            }

        /// Return the requested minimum ghost layer width for a body's central particle.
        virtual Scalar requestBodyGhostLayerWidth(unsigned int type, Scalar* h_r_ghost);

    private:

    };


namespace detail 
{
template<SmoothingKernelType KT_, StateEquationType SET_>
void export_SuspensionFlow(pybind11::module& m, std::string name);

} // end namespace detail
} // end namespace sph
} // end namespace hoomd

#endif // __SuspensionFlow_H__