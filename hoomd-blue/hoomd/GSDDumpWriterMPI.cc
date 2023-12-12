// Copyright (c) 2009-2023 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#include "GSDDumpWriterMPI.h"
#include "Filesystem.h"
#include "PGSD.h"
#include "HOOMDVersion.h"

#ifdef ENABLE_MPI
#include "Communicator.h"
#endif

#include <pybind11/numpy.h>
#include <pybind11/stl_bind.h>

#include <limits>
#include <list>
#include <sstream>
#include <stdexcept>
#include <string.h>
using namespace std;
using namespace hoomd::detail;

namespace hoomd
    {
std::list<std::string> GSDDumpWriterMPI::particle_chunks {"particles/position",
                                                       "particles/typeid",
                                                       "particles/mass",
                                                       "particles/slength"
                                                       // "particles/charge",
                                                       // "particles/diameter",
                                                       "particles/body",
                                                       // "particles/moment_inertia",
                                                       // "particles/orientation",
                                                       "particles/velocity",
                                                       // "particles/dpe",
                                                       "particles/density",
                                                       "particles/pressure",
                                                       "particles/energy",
                                                       "particles/auxiliary1",
                                                       "particles/auxiliary2",
                                                       "particles/auxiliary3",
                                                       "particles/auxiliary4",
                                                       // "particles/angmom",
                                                       "particles/image"};

/*! Constructs the GSDDumpWriterMPI. After construction, settings are set. No file operations are
    attempted until analyze() is called.

    \param sysdef SystemDefinition containing the ParticleData to dump
    \param fname File name to write data to
    \param group Group of particles to include in the output
    \param mode File open mode ("wb", "xb", or "ab")
    \param truncate If true, truncate the file to 0 frames every time analyze() called, then write
   out one frame

    If the group does not include all particles, then topology information cannot be written to the
   file.
*/
// GSDDumpWriterMPI::GSDDumpWriterMPI(std::shared_ptr<SystemDefinition> sysdef,
//                              std::shared_ptr<Trigger> trigger,
//                              const std::string& fname,
//                              std::shared_ptr<ParticleGroup> group,
//                              std::string mode,
//                              bool truncate)
//     : GSDDumpWriter(sysdef, trigger, fname, group, mode, truncate), m_fname(fname), m_mode(mode), m_truncate(truncate), m_group(group)
GSDDumpWriterMPI::GSDDumpWriterMPI(std::shared_ptr<SystemDefinition> sysdef,
                             std::shared_ptr<Trigger> trigger,
                             const std::string& fname,
                             std::shared_ptr<ParticleGroup> group,
                             std::string mode,
                             bool truncate)
    : Analyzer(sysdef, trigger), m_fname(fname), m_mode(mode), m_truncate(truncate), m_group(group)
    {
    m_exec_conf->msg->notice(5) << "Constructing GSDDumpWriterMPI: " << m_fname << " " << mode << " "
                                << truncate << endl;
    if (mode != "wb" && mode != "xb" && mode != "ab")
        {
        throw std::invalid_argument("Invalid PGSD file mode: " + mode);
        }
    m_log_writer = pybind11::none();

    m_dynamic.reset();
    m_dynamic[pgsd_flag::particles_position] = true;
    m_exec_conf->msg->notice(5) << "GSDDumpWriterMPI: start init File IO" << endl;
    initFileIO();
    std::cout << "GSDDumpWriterMPI: Done init File IO" << std::endl;
    }

pybind11::tuple GSDDumpWriterMPI::getDynamic()
    {
    pybind11::list result;

    if (m_dynamic[pgsd_flag::configuration_box])
        {
        result.append("configuration/box");
        }
    if (m_dynamic[pgsd_flag::particles_N])
        {
        result.append("particles/N");
        }
    if (m_dynamic[pgsd_flag::particles_position])
        {
        result.append("particles/position");
        }
    if (m_dynamic[pgsd_flag::particles_velocity])
        {
        result.append("particles/velocity");
        }
    if (m_dynamic[pgsd_flag::particles_image])
        {
        result.append("particles/image");
        }
    if (m_dynamic[pgsd_flag::particles_types])
        {
        result.append("particles/types");
        }
    if (m_dynamic[pgsd_flag::particles_type])
        {
        result.append("particles/typeid");
        }
    if (m_dynamic[pgsd_flag::particles_mass])
        {
        result.append("particles/mass");
        }
    if (m_dynamic[pgsd_flag::particles_slength])
        {
        result.append("particles/slength");
        }
    if (m_dynamic[pgsd_flag::particles_density])
        {
        result.append("particles/density");
        }
    if (m_dynamic[pgsd_flag::particles_pressure])
        {
        result.append("particles/pressure");
        }
    if (m_dynamic[pgsd_flag::particles_energy])
        {
        result.append("particles/energy");
        }
    if (m_dynamic[pgsd_flag::particles_aux1])
        {
        result.append("particles/aux1");
        }
    if (m_dynamic[pgsd_flag::particles_aux2])
        {
        result.append("particles/aux2");
        }
    if (m_dynamic[pgsd_flag::particles_aux3])
        {
        result.append("particles/aux3");
        }
    if (m_dynamic[pgsd_flag::particles_aux4])
        {
        result.append("particles/aux4");
        }

    return pybind11::tuple(result);
    }

void GSDDumpWriterMPI::setDynamic(pybind11::object dynamic)
    {
    pybind11::list dynamic_list = dynamic;
    m_dynamic.reset();
    m_write_topology = false;

    for (const auto& s_py : dynamic_list)
        {
        std::string s = s_py.cast<std::string>();
        if (s == "configuration/box" || s == "property")
            {
            m_dynamic[pgsd_flag::configuration_box] = true;
            }
        if (s == "particles/N" || s == "property")
            {
            m_dynamic[pgsd_flag::particles_N] = true;
            }
        if (s == "particles/position" || s == "property")
            {
            m_dynamic[pgsd_flag::particles_position] = true;
            }
        if (s == "particles/velocity" || s == "momentum")
            {
            m_dynamic[pgsd_flag::particles_velocity] = true;
            }
        if (s == "particles/image" || s == "momentum")
            {
            m_dynamic[pgsd_flag::particles_image] = true;
            }
        if (s == "particles/types" || s == "attribute")
            {
            m_dynamic[pgsd_flag::particles_types] = true;
            }
        if (s == "particles/typeid" || s == "attribute")
            {
            m_dynamic[pgsd_flag::particles_type] = true;
            }
        if (s == "particles/mass" || s == "attribute")
            {
            m_dynamic[pgsd_flag::particles_mass] = true;
            }
        if (s == "particles/slength" || s == "attribute")
            {
            m_dynamic[pgsd_flag::particles_slength] = true;
            }
        if (s == "particles/density" || s == "property")
            {
            m_dynamic[pgsd_flag::particles_density] = true;
            }
        if (s == "particles/pressure" || s == "property")
            {
            m_dynamic[pgsd_flag::particles_pressure] = true;
            }
        if (s == "particles/energy" || s == "property")
            {
            m_dynamic[pgsd_flag::particles_energy] = true;
            }
        if (s == "particles/aux1" || s == "momentum")
            {
            m_dynamic[pgsd_flag::particles_aux1] = true;
            }
        if (s == "particles/aux2" || s == "momentum")
            {
            m_dynamic[pgsd_flag::particles_aux2] = true;
            }
        if (s == "particles/aux3" || s == "momentum")
            {
            m_dynamic[pgsd_flag::particles_aux3] = true;
            }
        if (s == "particles/aux4" || s == "momentum")
            {
            m_dynamic[pgsd_flag::particles_aux4] = true;
            }
        }
    }

void GSDDumpWriterMPI::flush()
    {
    unsigned int rank = m_exec_conf->getRank();
    printf("Flush with rank %i\n", rank);
    // bool root = false;
    // if (m_exec_conf->isRoot())
    //     {
    m_exec_conf->msg->notice(5) << "PGSD: flush gsd file " << m_fname << endl;
        // root = true;
        // }
    int retval = pgsd_flush(&m_handle);
    PGSDUtils::checkError(retval, m_fname);
    }

void GSDDumpWriterMPI::setMaximumWriteBufferSize(uint64_t size)
    {
    // if (m_exec_conf->isRoot())
    //     {
    int retval = pgsd_set_maximum_write_buffer_size(&m_handle, size);
    PGSDUtils::checkError(retval, m_fname);

    // Scale the index buffer entires to write with the write buffer.
    retval = pgsd_set_index_entries_to_buffer(&m_handle, size / 256);
    PGSDUtils::checkError(retval, m_fname);
        // }
    }

uint64_t GSDDumpWriterMPI::getMaximumWriteBufferSize()
    {
    return pgsd_get_maximum_write_buffer_size(&m_handle);
    }

//! Initializes the output file for writing
void GSDDumpWriterMPI::initFileIO()
    {
    // create a new file or overwrite an existing one
    if (m_mode == "wb" || m_mode == "xb" || (m_mode == "ab" && !filesystem::exists(m_fname)))
        {
        ostringstream o;
        o << "HOOMD-blue " << HOOMD_VERSION;

        m_exec_conf->msg->notice(3) << "PGSD: create or overwrite gsd file " << m_fname << endl;
        m_exec_conf->msg->notice(3) << "PGSD: create or overwrite gsd file " << &m_handle << endl;
        int retval = pgsd_create_and_open(&m_handle,
                                         m_fname.c_str(),
                                         o.str().c_str(),
                                         "hoomd",
                                         pgsd_make_version(1, 4),
                                         PGSD_OPEN_APPEND,
                                         m_mode == "xb");
        printf("GSDDumpWriterMPI::initFileIO create and open retval: %i\n", retval);
        PGSDUtils::checkError(retval, m_fname);

        // in a created or overwritten file, all quantities are default
        for (auto const& chunk : particle_chunks)
            {
            m_nondefault[chunk] = false;
            }
        }
    // else if (m_mode == "ab")
    //     {
    //     // populate the non-default map
    //     populateNonDefault();

    //     // open the file in append mode
    //     m_exec_conf->msg->notice(3) << "PGSD: open gsd file " << m_fname << endl;
    //     int retval = pgsd_open(&m_handle, m_fname.c_str(), PGSD_OPEN_APPEND);
    //     PGSDUtils::checkError(retval, m_fname);

    //     // validate schema
    //     if (string(m_handle.header.schema) != string("hoomd"))
    //         {
    //         std::ostringstream s;
    //         s << "PGSD: "
    //           << "Invalid schema in " << m_fname;
    //         throw runtime_error("Error opening GSD file");
    //         }
    //     if (m_handle.header.schema_version >= pgsd_make_version(2, 0))
    //         {
    //         std::ostringstream s;
    //         s << "PGSD: "
    //           << "Invalid schema version in " << m_fname;
    //         throw runtime_error("Error opening GSD file");
    //         }
    //     }
    else
        {
        throw std::invalid_argument("Invalid PGSD file mode: " + m_mode);
        }
    MPI_Barrier(MPI_COMM_WORLD);
    m_nframes = pgsd_get_nframes(&m_handle);
    MPI_Barrier(MPI_COMM_WORLD);
    }

GSDDumpWriterMPI::~GSDDumpWriterMPI()
    {
    m_exec_conf->msg->notice(5) << "Destroying GSDDumpWriterMPI" << endl;

    // if (m_exec_conf->isRoot())
    //     {
    MPI_Barrier(MPI_COMM_WORLD);
    m_exec_conf->msg->notice(5) << "PGSD: close gsd file " << m_fname << endl;
    pgsd_close(&m_handle);
        // }
    }

//! Get the logged data for the current frame if any.
pybind11::dict GSDDumpWriterMPI::getLogData() const
    {
    if (!m_log_writer.is_none())
        {
        return m_log_writer.attr("log")().cast<pybind11::dict>();
        }
    return pybind11::dict();
    }

/*! \param timestep Current time step of the simulation

    The first call to analyze() will create or overwrite the file and write out the current system
   configuration as frame 0. Subsequent calls will append frames to the file, or keep overwriting
   frame 0 if m_truncate is true.
*/
void GSDDumpWriterMPI::analyze(uint64_t timestep)
    {
    unsigned int rank = m_exec_conf->getRank();
    printf("Analyse with rank %i\n", rank);
    m_exec_conf->msg->notice(5) << "GSDDumpWriterMPI: analyse" << endl;
    Analyzer::analyze(timestep);
    int retval;
    // bool root = true;

    // truncate the file if requested
    if (m_truncate)
        {
        // if (m_exec_conf->isRoot())
        //     {
        throw runtime_error("Error truncating GSD file: not implemented!");
        m_exec_conf->msg->notice(10) << "PGSD: truncating file" << endl;

        // retval = pgsd_truncate(&m_handle);
        // PGSDUtils::checkError(retval, m_fname);
            // }

        m_nframes = 0;
        }

    m_exec_conf->msg->notice(5) << "GSDDumpWriterMPI before Populate local frame N_global, N " << m_group->getNumMembersGlobal() << " " << m_group->getNumMembers() << endl;
    populateLocalFrame(m_local_frame, timestep);
    m_exec_conf->msg->notice(5) << "GSDDumpWriterMPI Populate local frame N_global, N " << m_group->getNumMembersGlobal() << " " << m_group->getNumMembers() << endl;
    auto log_data = getLogData();
    write(m_local_frame, log_data);
    }

void GSDDumpWriterMPI::write(GSDDumpWriterMPI::PGSDFrame& frame, pybind11::dict log_data)
    {

    printf("in GSDDumpWriterMPI:write Rank %i write header etc \n", m_exec_conf->getRank());
    writeFrameHeader(frame);
    MPI_Barrier(MPI_COMM_WORLD);
    writeAttributes(frame);
    MPI_Barrier(MPI_COMM_WORLD);
    writeProperties(frame);
    MPI_Barrier(MPI_COMM_WORLD);
    writeMomenta(frame);
    MPI_Barrier(MPI_COMM_WORLD);
    // writeLogQuantities(log_data);

    // topology is only meaningful if this is the all group
    // TODO
    // if (m_group->getNumMembersGlobal() == m_pdata->getNGlobal()
    //     && (m_write_topology || m_nframes == 0))
    //     {
    //     // if (m_exec_conf->isRoot())
    //     //     {
    //     writeTopology(frame.bond_data,
    //                   // frame.angle_data,
    //                   // frame.dihedral_data,
    //                   // frame.improper_data,
    //                   frame.constraint_data
    //                   //frame.pair_data
    //                   );
    //         // }
    //     }

    // if (m_exec_conf->isRoot())
    //     {
    m_exec_conf->msg->notice(10) << "PGSD: ending frame " << m_nframes << endl;
    MPI_Barrier(MPI_COMM_WORLD);
    int retval = pgsd_end_frame(&m_handle);
    PGSDUtils::checkError(retval, m_fname);
        // }

    m_nframes++;
    }

void GSDDumpWriterMPI::writeTypeMapping(std::string chunk, std::vector<std::string> type_mapping, const GSDDumpWriterMPI::PGSDFrame& frame)
    {
    uint32_t N_global = m_group->getNumMembersGlobal();
    unsigned int rank = m_exec_conf->getRank();

    int max_len = 0;
    for (unsigned int i = 0; i < type_mapping.size(); i++)
        {
        max_len = std::max(max_len, (int)type_mapping[i].size());
        }
    max_len += 1; // for null

    m_exec_conf->msg->notice(10) << "PGSD: writing " << chunk << endl;
    std::vector<char> types(max_len * type_mapping.size());
    for (unsigned int i = 0; i < type_mapping.size(); i++){
        strncpy(&types[max_len * i], type_mapping[i].c_str(), max_len);
    }

    int retval = pgsd_write_chunk(&m_handle,
                                 chunk.c_str(),
                                 PGSD_TYPE_UINT8,
                                 type_mapping.size(),
                                 max_len,
                                 type_mapping.size(),
                                 max_len,
                                 0,
                                 0,
                                 false,
                                 0,
                                 (void*)&types[0]);
    PGSDUtils::checkError(retval, m_fname);
    }

/*! Write the data chunks configuration/step, configuration/box, and particles/N. If this is frame
   0, also write configuration/dimensions.
*/
void GSDDumpWriterMPI::writeFrameHeader(const GSDDumpWriterMPI::PGSDFrame& frame)
    {
    int retval;
    m_exec_conf->msg->notice(10) << "PGSD: writing configuration/step" << endl;
    retval = pgsd_write_chunk(&m_handle,
                             "configuration/step",
                             PGSD_TYPE_UINT64,
                             1,
                             1,
                             1,
                             1,
                             0,
                             0,
                             false,
                             0,
                             (void*)&frame.timestep);
    PGSDUtils::checkError(retval, m_fname);

    if (m_nframes == 0)
        {
        m_exec_conf->msg->notice(10) << "PGSD: writing configuration/dimensions" << endl;
        uint8_t dimensions = (uint8_t)m_sysdef->getNDimensions();
        retval = pgsd_write_chunk(&m_handle,
                                 "configuration/dimensions",
                                 PGSD_TYPE_UINT8,
                                 1,
                                 1,
                                 1,
                                 1,
                                 0,
                                 0,
                                 false,
                                 0,
                                 (void*)&dimensions);
        PGSDUtils::checkError(retval, m_fname);
        }

    if (m_nframes == 0 || m_dynamic[pgsd_flag::configuration_box])
        {
        m_exec_conf->msg->notice(10) << "PGSD: writing configuration/box" << endl;
        float box_a[6];
        box_a[0] = (float)frame.global_box.getL().x;
        box_a[1] = (float)frame.global_box.getL().y;
        box_a[2] = (float)frame.global_box.getL().z;
        box_a[3] = (float)frame.global_box.getTiltFactorXY();
        box_a[4] = (float)frame.global_box.getTiltFactorXZ();
        box_a[5] = (float)frame.global_box.getTiltFactorYZ();
        std::cout << "box " << box_a[0] << " " << box_a[1] << " " << box_a[2] << std::endl;
        retval = pgsd_write_chunk(&m_handle,
                                 "configuration/box",
                                 PGSD_TYPE_FLOAT,
                                 6,
                                 1,
                                 6,
                                 1,
                                 0,
                                 0,
                                 false,
                                 0,
                                 (void*)box_a);
        PGSDUtils::checkError(retval, m_fname);
        }

    if (m_nframes == 0 || m_dynamic[pgsd_flag::particles_N])
        {
        m_exec_conf->msg->notice(10) << "PGSD: writing particles/N" << endl;
        uint32_t N = m_group->getNumMembersGlobal();
        retval = pgsd_write_chunk(&m_handle, 
                                  "particles/N", 
                                  PGSD_TYPE_UINT32, 
                                  1, 
                                  1,
                                  1,
                                  1,
                                  0,
                                  0,
                                  false, 
                                  0,
                                  (void*)&N);
        PGSDUtils::checkError(retval, m_fname);
        }
    }

/*! Writes the data chunks typeid, mass, body in
   particles/.
*/
void GSDDumpWriterMPI::writeAttributes(GSDDumpWriterMPI::PGSDFrame& frame)
    {
    uint32_t N = m_group->getNumMembers();
    uint32_t N_global = m_group->getNumMembersGlobal();
    // bool all_default = true;
    unsigned int rank = m_exec_conf->getRank();
    unsigned int size = m_exec_conf->getNRanks();
    int part_offset;

    std::vector<unsigned int> part_distribution(size);
    all_gather_v(N, part_distribution, MPI_COMM_WORLD);
    part_offset = std::accumulate(part_distribution.begin(), part_distribution.begin()+rank, 0);

    printf("rank %i part_distribution [rank] %i\n", rank, part_distribution[rank]);

    int i;
    for(i = 0; i < size; i++){
        printf("Rank %i, part_distribution [%i] %i\n", rank, i, part_distribution[i]);
    }

    printf("GSDDumpWriterMPI write attributes: rank %i\n", rank);
    printf("GSDDumpWriterMPI write attributes: N %i\n", N);
    printf("GSDDumpWriterMPI write attributes: N_global %i\n", N_global);

    printf("GSDDumpWriterMPI write attributes: part_offset %i\n", part_offset);

    int retval;
    MPI_Barrier(MPI_COMM_WORLD);
    if (m_dynamic[pgsd_flag::particles_types] || m_nframes == 0)
        {
        writeTypeMapping("particles/types", frame.particle_data.type_mapping, frame);
        }

    if (frame.particle_data.type.size() != 0)
        {
        assert(frame.particle_data.type.size() == N);

        m_exec_conf->msg->notice(10) << "PGSD: writing particles/typeid" << endl;
        retval = pgsd_write_chunk(&m_handle,
                                 "particles/typeid",
                                 PGSD_TYPE_UINT32,
                                 N,
                                 1,
                                 N_global,
                                 1,
                                 part_offset,
                                 N_global,
                                 true,
                                 0,
                                 (void*)frame.particle_data.type.data());
        PGSDUtils::checkError(retval, m_fname);
        if (m_nframes == 0)
            m_nondefault["particles/typeid"] = true;
        }

    if (frame.particle_data.mass.size() != 0)
        {
        assert(frame.particle_data.mass.size() == N);

        m_exec_conf->msg->notice(10) << "PGSD: writing particles/mass" << endl;
        retval = pgsd_write_chunk(&m_handle,
                                 "particles/mass",
                                 PGSD_TYPE_FLOAT,
                                 N,
                                 1,
                                 N_global,
                                 1,
                                 part_offset,
                                 N_global,
                                 true,
                                 0,
                                 (void*)frame.particle_data.mass.data());
        PGSDUtils::checkError(retval, m_fname);
        if (m_nframes == 0)
            m_nondefault["particles/mass"] = true;
        }

    if (frame.particle_data.slength.size() != 0)
        {
        assert(frame.particle_data.slength.size() == N);

        m_exec_conf->msg->notice(10) << "PGSD: writing particles/slength" << endl;
        retval = pgsd_write_chunk(&m_handle,
                                 "particles/slength",
                                 PGSD_TYPE_FLOAT,
                                 N,
                                 1,
                                 N_global,
                                 1,
                                 part_offset,
                                 N_global,
                                 true,
                                 0,
                                 (void*)frame.particle_data.slength.data());
        PGSDUtils::checkError(retval, m_fname);
        if (m_nframes == 0)
            m_nondefault["particles/slength"] = true;
        }


    if (frame.particle_data.body.size() != 0)
        {
        assert(frame.particle_data.body.size() == N);

        m_exec_conf->msg->notice(10) << "PGSD: writing particles/body" << endl;
        retval = pgsd_write_chunk(&m_handle,
                                 "particles/body",
                                 PGSD_TYPE_INT32,
                                 N,
                                 1,
                                 N_global,
                                 1,
                                 part_offset,
                                 N_global,
                                 true,
                                 0,
                                 (void*)frame.particle_data.body.data());
        PGSDUtils::checkError(retval, m_fname);
        if (m_nframes == 0)
            m_nondefault["particles/body"] = true;
        }
    }

/*! Writes the data chunks position and orientation in particles/.
 */
void GSDDumpWriterMPI::writeProperties(const GSDDumpWriterMPI::PGSDFrame& frame)
    {
    uint32_t N = m_group->getNumMembers();
    uint32_t N_global = m_group->getNumMembersGlobal();
    int retval;
    unsigned int rank = m_exec_conf->getRank();
    unsigned int size = m_exec_conf->getNRanks();
    int part_offset;
    
    std::vector<unsigned int> part_distribution(size);
    all_gather_v(N, part_distribution, MPI_COMM_WORLD);
    part_offset = std::accumulate(part_distribution.begin(), part_distribution.begin()+rank, 0);


    if (frame.particle_data.pos.size() != 0)
        {
        assert(frame.particle_data.pos.size() == N);

        m_exec_conf->msg->notice(10) << "PGSD: writing particles/position" << endl;
        retval = pgsd_write_chunk(&m_handle,
                                 "particles/position",
                                 PGSD_TYPE_FLOAT,
                                 N,
                                 3,
                                 N_global,
                                 3,
                                 3*part_offset,
                                 3*N_global,
                                 true,
                                 0,
                                 (void*)frame.particle_data.pos.data());
        PGSDUtils::checkError(retval, m_fname);
        if (m_nframes == 0)
            m_nondefault["particles/position"] = true;
        }

    if (frame.particle_data.density.size() != 0)
        {
        assert(frame.particle_data.density.size() == N);
        m_exec_conf->msg->notice(10) << "PGSD: writing particles/density" << endl;
        retval = pgsd_write_chunk(&m_handle,
                                 "particles/density",
                                 PGSD_TYPE_FLOAT,
                                 N,
                                 1,
                                 N_global,
                                 1,
                                 part_offset,
                                 N_global,
                                 true,
                                 0,
                                 (void*)frame.particle_data.density.data());
        PGSDUtils::checkError(retval, m_fname);
        if (m_nframes == 0)
            m_nondefault["particles/density"] = true;
        }

    if (frame.particle_data.pressure.size() != 0)
        {
        assert(frame.particle_data.pressure.size() == N);

        m_exec_conf->msg->notice(10) << "PGSD: writing particles/pressure" << endl;
        retval = pgsd_write_chunk(&m_handle,
                                 "particles/pressure",
                                 PGSD_TYPE_FLOAT,
                                 N,
                                 1,
                                 N_global,
                                 1,
                                 part_offset,
                                 N_global,
                                 true,
                                 0,
                                 (void*)frame.particle_data.pressure.data());
        PGSDUtils::checkError(retval, m_fname);
        if (m_nframes == 0)
            m_nondefault["particles/pressure"] = true;
        }

    if (frame.particle_data.energy.size() != 0)
        {
        assert(frame.particle_data.energy.size() == N);

        m_exec_conf->msg->notice(10) << "PGSD: writing particles/energy" << endl;
        retval = pgsd_write_chunk(&m_handle,
                                 "particles/energy",
                                 PGSD_TYPE_FLOAT,
                                 N,
                                 1,
                                 N_global,
                                 1,
                                 part_offset,
                                 N_global,
                                 true,
                                 0,
                                 (void*)frame.particle_data.energy.data());
        PGSDUtils::checkError(retval, m_fname);
        if (m_nframes == 0)
            m_nondefault["particles/energy"] = true;
        }
    }

/*! Writes the data chunks velocity, angmom, and image in particles/.
 */
void GSDDumpWriterMPI::writeMomenta(const GSDDumpWriterMPI::PGSDFrame& frame)
    {
    uint32_t N = m_group->getNumMembers();
    uint32_t N_global = m_group->getNumMembersGlobal();
    int retval;

    unsigned int rank = m_exec_conf->getRank();
    unsigned int size = m_exec_conf->getNRanks();
    int part_offset;
    std::vector<unsigned int> part_distribution(size);
    all_gather_v(N, part_distribution, MPI_COMM_WORLD);
    part_offset = std::accumulate(part_distribution.begin(), part_distribution.begin()+rank, 0);

    if (frame.particle_data.vel.size() != 0)
        {
        assert(frame.particle_data.vel.size() == N);

        m_exec_conf->msg->notice(10) << "PGSD: writing particles/velocity" << endl;
        retval = pgsd_write_chunk(&m_handle,
                                 "particles/velocity",
                                 PGSD_TYPE_FLOAT,
                                 N,
                                 3,
                                 N_global,
                                 3,
                                 3*part_offset,
                                 3*N_global,
                                 true,
                                 0,
                                 (void*)frame.particle_data.vel.data());

        PGSDUtils::checkError(retval, m_fname);
        if (m_nframes == 0)
            m_nondefault["particles/velocity"] = true;
        }

    if (frame.particle_data.aux1.size() != 0)
        {
        assert(frame.particle_data.aux1.size() == N);

        m_exec_conf->msg->notice(10) << "PGSD: writing particles/auxiliary1" << endl;
        retval = pgsd_write_chunk(&m_handle,
                                 "particles/auxiliary1",
                                 PGSD_TYPE_FLOAT,
                                 N,
                                 3,
                                 N_global,
                                 3,
                                 3*part_offset,
                                 3*N_global,
                                 true,
                                 0,
                                 (void*)frame.particle_data.aux1.data());
        PGSDUtils::checkError(retval, m_fname);
        if (m_nframes == 0)
            m_nondefault["particles/auxiliary1"] = true;
        }
    
    if (frame.particle_data.aux2.size() != 0)
        {
        assert(frame.particle_data.aux2.size() == N);

        m_exec_conf->msg->notice(10) << "PGSD: writing particles/auxiliary2" << endl;
        retval = pgsd_write_chunk(&m_handle,
                                 "particles/auxiliary2",
                                 PGSD_TYPE_FLOAT,
                                 N,
                                 3,
                                 N_global,
                                 3,
                                 3*part_offset,
                                 3*N_global,
                                 true,
                                 0,
                                 (void*)frame.particle_data.aux2.data());
        PGSDUtils::checkError(retval, m_fname);
        if (m_nframes == 0)
            m_nondefault["particles/auxiliary2"] = true;
        }

    if (frame.particle_data.aux3.size() != 0)
        {
        assert(frame.particle_data.aux3.size() == N);

        m_exec_conf->msg->notice(10) << "PGSD: writing particles/auxiliary3" << endl;
        retval = pgsd_write_chunk(&m_handle,
                                 "particles/auxiliary3",
                                 PGSD_TYPE_FLOAT,
                                 N,
                                 3,
                                 N_global,
                                 3,
                                 3*part_offset,
                                 3*N_global,
                                 true,
                                 0,
                                 (void*)frame.particle_data.aux3.data());
        PGSDUtils::checkError(retval, m_fname);
        if (m_nframes == 0)
            m_nondefault["particles/auxiliary3"] = true;
        }

    if (frame.particle_data.aux4.size() != 0)
        {
        assert(frame.particle_data.aux4.size() == N);

        m_exec_conf->msg->notice(10) << "PGSD: writing particles/auxiliary4" << endl;
        retval = pgsd_write_chunk(&m_handle,
                                 "particles/auxiliary4",
                                 PGSD_TYPE_FLOAT,
                                 N,
                                 3,
                                 N_global,
                                 3,
                                 3*part_offset,
                                 3*N_global,
                                 true,
                                 0,
                                 (void*)frame.particle_data.aux4.data());
        PGSDUtils::checkError(retval, m_fname);
        if (m_nframes == 0)
            m_nondefault["particles/auxiliary4"] = true;
        }

    if (frame.particle_data.image.size() != 0)
        {
        assert(frame.particle_data.image.size() == N);

        m_exec_conf->msg->notice(10) << "PGSD: writing particles/image" << endl;
        retval = pgsd_write_chunk(&m_handle,
                                 "particles/image",
                                 PGSD_TYPE_INT32,
                                 N,
                                 3,
                                 N_global,
                                 3,
                                 3*part_offset,
                                 3*N_global,
                                 true,
                                 0,
                                 (void*)frame.particle_data.image.data());
        PGSDUtils::checkError(retval, m_fname);
        if (m_nframes == 0)
            m_nondefault["particles/image"] = true;
        }
    }

/*! \param bond Bond data snapshot
    \param angle Angle data snapshot
    \param dihedral Dihedral data snapshot
    \param improper Improper data snapshot
    \param constraint Constraint data snapshot
    \param pair Special pair data snapshot

    Write out all the snapshot data to the GSD file
*/
void GSDDumpWriterMPI::writeTopology(BondData::Snapshot& bond,
                                  // AngleData::Snapshot& angle,
                                  // DihedralData::Snapshot& dihedral,
                                  // ImproperData::Snapshot& improper,
                                  ConstraintData::Snapshot& constraint
                                  // PairData::Snapshot& pair
                                  )
    {
    // if (bond.size > 0)
    //     {
    //     m_exec_conf->msg->notice(10) << "PGSD: writing bonds/N" << endl;
    //     uint32_t N = bond.size;
    //     int retval = gsd_write_chunk(&m_handle, "bonds/N", PGSD_TYPE_UINT32, 1, 1, 0, (void*)&N);
    //     PGSDUtils::checkError(retval, m_fname);

    //     writeTypeMapping("bonds/types", bond.type_mapping);

    //     m_exec_conf->msg->notice(10) << "PGSD: writing bonds/typeid" << endl;
    //     retval = gsd_write_chunk(&m_handle,
    //                              "bonds/typeid",
    //                              PGSD_TYPE_UINT32,
    //                              N,
    //                              1,
    //                              0,
    //                              (void*)&bond.type_id[0]);
    //     PGSDUtils::checkError(retval, m_fname);

    //     m_exec_conf->msg->notice(10) << "PGSD: writing bonds/group" << endl;
    //     retval = gsd_write_chunk(&m_handle,
    //                              "bonds/group",
    //                              PGSD_TYPE_UINT32,
    //                              N,
    //                              2,
    //                              0,
    //                              (void*)&bond.groups[0]);
    //     PGSDUtils::checkError(retval, m_fname);
    //     }
    // // if (angle.size > 0)
    // //     {
    // //     m_exec_conf->msg->notice(10) << "PGSD: writing angles/N" << endl;
    // //     uint32_t N = angle.size;
    // //     int retval = gsd_write_chunk(&m_handle, "angles/N", PGSD_TYPE_UINT32, 1, 1, 0, (void*)&N);
    // //     PGSDUtils::checkError(retval, m_fname);

    // //     writeTypeMapping("angles/types", angle.type_mapping);

    // //     m_exec_conf->msg->notice(10) << "PGSD: writing angles/typeid" << endl;
    // //     retval = gsd_write_chunk(&m_handle,
    // //                              "angles/typeid",
    // //                              PGSD_TYPE_UINT32,
    // //                              N,
    // //                              1,
    // //                              0,
    // //                              (void*)&angle.type_id[0]);
    // //     PGSDUtils::checkError(retval, m_fname);

    // //     m_exec_conf->msg->notice(10) << "PGSD: writing angles/group" << endl;
    // //     retval = gsd_write_chunk(&m_handle,
    // //                              "angles/group",
    // //                              PGSD_TYPE_UINT32,
    // //                              N,
    // //                              3,
    // //                              0,
    // //                              (void*)&angle.groups[0]);
    // //     PGSDUtils::checkError(retval, m_fname);
    // //     }
    // // if (dihedral.size > 0)
    // //     {
    // //     m_exec_conf->msg->notice(10) << "PGSD: writing dihedrals/N" << endl;
    // //     uint32_t N = dihedral.size;
    // //     int retval = gsd_write_chunk(&m_handle, "dihedrals/N", PGSD_TYPE_UINT32, 1, 1, 0, (void*)&N);
    // //     PGSDUtils::checkError(retval, m_fname);

    // //     writeTypeMapping("dihedrals/types", dihedral.type_mapping);

    // //     m_exec_conf->msg->notice(10) << "PGSD: writing dihedrals/typeid" << endl;
    // //     retval = gsd_write_chunk(&m_handle,
    // //                              "dihedrals/typeid",
    // //                              PGSD_TYPE_UINT32,
    // //                              N,
    // //                              1,
    // //                              0,
    // //                              (void*)&dihedral.type_id[0]);
    // //     PGSDUtils::checkError(retval, m_fname);

    // //     m_exec_conf->msg->notice(10) << "PGSD: writing dihedrals/group" << endl;
    // //     retval = gsd_write_chunk(&m_handle,
    // //                              "dihedrals/group",
    // //                              PGSD_TYPE_UINT32,
    // //                              N,
    // //                              4,
    // //                              0,
    // //                              (void*)&dihedral.groups[0]);
    // //     PGSDUtils::checkError(retval, m_fname);
    // //     }
    // // if (improper.size > 0)
    // //     {
    // //     m_exec_conf->msg->notice(10) << "PGSD: writing impropers/N" << endl;
    // //     uint32_t N = improper.size;
    // //     int retval = gsd_write_chunk(&m_handle, "impropers/N", PGSD_TYPE_UINT32, 1, 1, 0, (void*)&N);
    // //     PGSDUtils::checkError(retval, m_fname);

    // //     writeTypeMapping("impropers/types", improper.type_mapping);

    // //     m_exec_conf->msg->notice(10) << "PGSD: writing impropers/typeid" << endl;
    // //     retval = gsd_write_chunk(&m_handle,
    // //                              "impropers/typeid",
    // //                              PGSD_TYPE_UINT32,
    // //                              N,
    // //                              1,
    // //                              0,
    // //                              (void*)&improper.type_id[0]);
    // //     PGSDUtils::checkError(retval, m_fname);

    // //     m_exec_conf->msg->notice(10) << "PGSD: writing impropers/group" << endl;
    // //     retval = gsd_write_chunk(&m_handle,
    // //                              "impropers/group",
    // //                              PGSD_TYPE_UINT32,
    // //                              N,
    // //                              4,
    // //                              0,
    // //                              (void*)&improper.groups[0]);
    // //     PGSDUtils::checkError(retval, m_fname);
    // //     }

    // if (constraint.size > 0)
    //     {
    //     m_exec_conf->msg->notice(10) << "PGSD: writing constraints/N" << endl;
    //     uint32_t N = constraint.size;
    //     int retval
    //         = gsd_write_chunk(&m_handle, "constraints/N", PGSD_TYPE_UINT32, 1, 1, 0, (void*)&N);
    //     PGSDUtils::checkError(retval, m_fname);

    //     m_exec_conf->msg->notice(10) << "PGSD: writing constraints/value" << endl;
    //         {
    //         std::vector<float> data(N);
    //         data.reserve(1); //! make sure we allocate
    //         for (unsigned int i = 0; i < N; i++)
    //             data[i] = float(constraint.val[i]);

    //         retval = gsd_write_chunk(&m_handle,
    //                                  "constraints/value",
    //                                  PGSD_TYPE_FLOAT,
    //                                  N,
    //                                  1,
    //                                  0,
    //                                  (void*)&data[0]);
    //         PGSDUtils::checkError(retval, m_fname);
    //         }

    //     m_exec_conf->msg->notice(10) << "PGSD: writing constraints/group" << endl;
    //     retval = gsd_write_chunk(&m_handle,
    //                              "constraints/group",
    //                              PGSD_TYPE_UINT32,
    //                              N,
    //                              2,
    //                              0,
    //                              (void*)&constraint.groups[0]);
    //     PGSDUtils::checkError(retval, m_fname);
    // }

    // if (pair.size > 0)
    //     {
    //     m_exec_conf->msg->notice(10) << "PGSD: writing pairs/N" << endl;
    //     uint32_t N = pair.size;
    //     int retval = gsd_write_chunk(&m_handle, "pairs/N", PGSD_TYPE_UINT32, 1, 1, 0, (void*)&N);
    //     PGSDUtils::checkError(retval, m_fname);

    //     writeTypeMapping("pairs/types", pair.type_mapping);

    //     m_exec_conf->msg->notice(10) << "PGSD: writing pairs/typeid" << endl;
    //     retval = gsd_write_chunk(&m_handle,
    //                              "pairs/typeid",
    //                              PGSD_TYPE_UINT32,
    //                              N,
    //                              1,
    //                              0,
    //                              (void*)&pair.type_id[0]);
    //     PGSDUtils::checkError(retval, m_fname);

    //     m_exec_conf->msg->notice(10) << "PGSD: writing pairs/group" << endl;
    //     retval = gsd_write_chunk(&m_handle,
    //                              "pairs/group",
    //                              PGSD_TYPE_UINT32,
    //                              N,
    //                              2,
    //                              0,
    //                              (void*)&pair.groups[0]);
    //     PGSDUtils::checkError(retval, m_fname);
    //     }
    }

void GSDDumpWriterMPI::writeLogQuantities(pybind11::dict dict)
    {
    for (auto key_iter = dict.begin(); key_iter != dict.end(); ++key_iter)
        {
        std::string name = pybind11::cast<std::string>(key_iter->first);
        m_exec_conf->msg->notice(10) << "PGSD: writing " << name << endl;

        pybind11::array arr = pybind11::array::ensure(key_iter->second, pybind11::array::c_style);
        pgsd_type type = PGSD_TYPE_UINT8;
        auto dtype = arr.dtype();
        if (dtype.kind() == 'u' && dtype.itemsize() == 1)
            {
            type = PGSD_TYPE_UINT8;
            }
        else if (dtype.kind() == 'u' && dtype.itemsize() == 2)
            {
            type = PGSD_TYPE_UINT16;
            }
        else if (dtype.kind() == 'u' && dtype.itemsize() == 4)
            {
            type = PGSD_TYPE_UINT32;
            }
        else if (dtype.kind() == 'u' && dtype.itemsize() == 8)
            {
            type = PGSD_TYPE_UINT64;
            }
        else if (dtype.kind() == 'i' && dtype.itemsize() == 1)
            {
            type = PGSD_TYPE_INT8;
            }
        else if (dtype.kind() == 'i' && dtype.itemsize() == 2)
            {
            type = PGSD_TYPE_INT16;
            }
        else if (dtype.kind() == 'i' && dtype.itemsize() == 4)
            {
            type = PGSD_TYPE_INT32;
            }
        else if (dtype.kind() == 'i' && dtype.itemsize() == 8)
            {
            type = PGSD_TYPE_INT64;
            }
        else if (dtype.kind() == 'f' && dtype.itemsize() == 4)
            {
            type = PGSD_TYPE_FLOAT;
            }
        else if (dtype.kind() == 'f' && dtype.itemsize() == 8)
            {
            type = PGSD_TYPE_DOUBLE;
            }
        else if (dtype.kind() == 'b' && dtype.itemsize() == 1)
            {
            type = PGSD_TYPE_UINT8;
            }
        else
            {
            throw range_error("Invalid numpy array format in gsd log data [" + name
                              + "]: " + string(pybind11::str(arr.dtype())));
            }

        size_t M = 1;
        size_t N = 1;
        auto ndim = arr.ndim();
        if (ndim == 0)
            {
            // numpy converts scalars to arrays with zero dimensions
            // gsd treats them as 1x1 arrays.
            M = 1;
            N = 1;
            }
        if (ndim == 1)
            {
            N = arr.shape(0);
            M = 1;
            }
        if (ndim == 2)
            {
            N = arr.shape(0);
            M = arr.shape(1);
            if (M > std::numeric_limits<uint32_t>::max())
                throw runtime_error("Array dimension too large in gsd log data [" + name + "]");
            }
        if (ndim > 2)
            {
            throw invalid_argument("Invalid numpy dimension in gsd log data [" + name + "]");
            }

        int retval
            = pgsd_write_chunk(&m_handle, 
                               name.c_str(), 
                               type, 
                               N, 
                               (uint32_t)M,
                               N,
                               (uint32_t)M,
                               0,
                               0,
                               false,
                               0,
                               (void*)arr.data());
        PGSDUtils::checkError(retval, m_fname);
        }
    }

/*! Populate the m_nondefault map.
    Set entries to true when they exist in frame 0 of the file, otherwise, set them to false.
*/
void GSDDumpWriterMPI::populateNonDefault()
    {
    int retval;

    // open the file in read only mode
    m_exec_conf->msg->notice(3) << "PGSD: check frame 0 in gsd file " << m_fname << endl;
    retval = pgsd_open(&m_handle, m_fname.c_str(), PGSD_OPEN_READONLY);
    PGSDUtils::checkError(retval, m_fname);

    // validate schema
    if (string(m_handle.header.schema) != string("hoomd"))
        {
        std::ostringstream s;
        s << "PGSD: "
          << "Invalid schema in " << m_fname;
        throw runtime_error("Error opening GSD file");
        }
    if (m_handle.header.schema_version >= pgsd_make_version(2, 0))
        {
        std::ostringstream s;
        s << "PGSD: "
          << "Invalid schema version in " << m_fname;
        throw runtime_error("Error opening GSD file");
        }

    for (auto const& chunk : particle_chunks)
        {
        const pgsd_index_entry* entry = pgsd_find_chunk(&m_handle, 0, chunk.c_str());
        m_nondefault[chunk] = (entry != nullptr);
        }

    // close the file
    pgsd_close(&m_handle);
    }

void GSDDumpWriterMPI::populateLocalFrame(GSDDumpWriterMPI::PGSDFrame& frame, uint64_t timestep)
    {
    m_exec_conf->msg->notice(5) << "GSDDumpWriterMPI::populateLocalFrame" << endl;
    frame.timestep = timestep;
    frame.global_box = m_pdata->getGlobalBox();

    frame.particle_data.type_mapping = m_pdata->getTypeMapping();

    // Global Group size since could be zero on one rank
    uint32_t N = m_group->getNumMembersGlobal();

    // Assume values are all default to start, set flags to false when we find a non-default.
    std::bitset<n_pgsd_flags> all_default;
    all_default.set();
    frame.clear();

    ArrayHandle<unsigned int> h_rtag(m_pdata->getRTags(), access_location::host, access_mode::read);

    if (N > 0)
        {
        ArrayHandle<unsigned int> h_tag(m_pdata->getTags(),
                                        access_location::host,
                                        access_mode::read);

        m_index.resize(0);

        for (unsigned int group_tag_index = 0; group_tag_index < N; group_tag_index++)
            {
            unsigned int tag = m_group->getMemberTag(group_tag_index);
            unsigned int index = h_rtag.data[tag];
            if (index >= m_pdata->getN())
                {
                continue;
                }

            frame.particle_tags.push_back(h_tag.data[index]);
            m_index.push_back(index);
            }
        }

    std::cout << "Rank " << m_exec_conf->getRank() << " m_index size " << m_index.size() << std::endl;

    if (N > 0
        && (m_dynamic[pgsd_flag::particles_position] || m_dynamic[pgsd_flag::particles_type]
            || m_dynamic[pgsd_flag::particles_image] || m_nframes == 0))
        {
        ArrayHandle<Scalar4> h_postype(m_pdata->getPositions(),
                                       access_location::host,
                                       access_mode::read);
        ArrayHandle<int3> h_image(m_pdata->getImages(), access_location::host, access_mode::read);

        if (m_dynamic[pgsd_flag::particles_position] || m_nframes == 0)
            {
            frame.particle_data_present[pgsd_flag::particles_position] = true;
            }
        if (m_dynamic[pgsd_flag::particles_image] || m_nframes == 0)
            {
            frame.particle_data_present[pgsd_flag::particles_image] = true;
            }
        if (m_dynamic[pgsd_flag::particles_type] || m_nframes == 0)
            {
            frame.particle_data_present[pgsd_flag::particles_type] = true;
            }

        for (unsigned int index : m_index)
            {
            vec3<Scalar> position
                = vec3<Scalar>(h_postype.data[index]) - vec3<Scalar>(m_pdata->getOrigin());
            unsigned int type = __scalar_as_int(h_postype.data[index].w);
            int3 image = make_int3(0, 0, 0);

            if (m_dynamic[pgsd_flag::particles_image] || m_nframes == 0)
                {
                image = h_image.data[index];
                }

            frame.global_box.wrap(position, image);

            if (m_dynamic[pgsd_flag::particles_position] || m_nframes == 0)
                {
                if (position != vec3<Scalar>(0, 0, 0))
                    {
                    all_default[pgsd_flag::particles_position] = false;
                    }

                frame.particle_data.pos.push_back(vec3<float>(position));
                }

            if (m_dynamic[pgsd_flag::particles_image] || m_nframes == 0)
                {
                if (image != make_int3(0, 0, 0))
                    {
                    all_default[pgsd_flag::particles_image] = false;
                    }

                frame.particle_data.image.push_back(image);
                }

            if (m_dynamic[pgsd_flag::particles_type] || m_nframes == 0)
                {
                if (type != 0)
                    {
                    all_default[pgsd_flag::particles_type] = false;
                    }

                frame.particle_data.type.push_back(type);
                }
            }
        }


    if (N > 0
        && (m_dynamic[pgsd_flag::particles_velocity] || m_dynamic[pgsd_flag::particles_mass]
            || m_nframes == 0))
        {
        ArrayHandle<Scalar4> h_velocity_mass(m_pdata->getVelocities(),
                                             access_location::host,
                                             access_mode::read);

        if (m_dynamic[pgsd_flag::particles_mass] || m_nframes == 0)
            {
            frame.particle_data_present[pgsd_flag::particles_mass] = true;
            }
        if (m_dynamic[pgsd_flag::particles_velocity] || m_nframes == 0)
            {
            frame.particle_data_present[pgsd_flag::particles_velocity] = true;
            }

        for (unsigned int index : m_index)
            {
            vec3<float> velocity = vec3<float>(static_cast<float>(h_velocity_mass.data[index].x),
                                               static_cast<float>(h_velocity_mass.data[index].y),
                                               static_cast<float>(h_velocity_mass.data[index].z));
            float mass = static_cast<float>(h_velocity_mass.data[index].w);

            if (m_dynamic[pgsd_flag::particles_mass] || m_nframes == 0)
                {
                if (mass != 1.0f)
                    {
                    all_default[pgsd_flag::particles_mass] = false;
                    }

                frame.particle_data.mass.push_back(mass);
                }

            if (m_dynamic[pgsd_flag::particles_velocity] || m_nframes == 0)
                {
                if (velocity != vec3<float>(0, 0, 0))
                    {
                    all_default[pgsd_flag::particles_velocity] = false;
                    }

                frame.particle_data.vel.push_back(velocity);
                }
            }
        }

    if (N > 0 && (m_dynamic[pgsd_flag::particles_slength] || m_nframes == 0))
        {
        ArrayHandle<Scalar> h_slength(m_pdata->getSlengths(),
                                     access_location::host,
                                     access_mode::read);

        frame.particle_data_present[pgsd_flag::particles_slength] = true;

        for (unsigned int index : m_index)
            {
            float slength = static_cast<float>(h_slength.data[index]);
            if (slength != 0.0f)
                {
                all_default[pgsd_flag::particles_slength] = false;
                }

            frame.particle_data.slength.push_back(slength);
            }
        }

    if (N > 0 && (m_dynamic[pgsd_flag::particles_density] || m_nframes == 0))
        {
        ArrayHandle<Scalar> h_density(m_pdata->getDensities(),
                                     access_location::host,
                                     access_mode::read);

        frame.particle_data_present[pgsd_flag::particles_density] = true;

        for (unsigned int index : m_index)
            {
            float density = static_cast<float>(h_density.data[index]);
            if (density != 0.0f)
                {
                all_default[pgsd_flag::particles_density] = false;
                }

            frame.particle_data.density.push_back(density);
            }
        }

    if (N > 0 && (m_dynamic[pgsd_flag::particles_pressure] || m_nframes == 0))
        {
        ArrayHandle<Scalar> h_pressure(m_pdata->getPressures(),
                                     access_location::host,
                                     access_mode::read);

        frame.particle_data_present[pgsd_flag::particles_pressure] = true;

        for (unsigned int index : m_index)
            {
            float pressure = static_cast<float>(h_pressure.data[index]);
            if (pressure != 0.0f)
                {
                all_default[pgsd_flag::particles_pressure] = false;
                }

            frame.particle_data.pressure.push_back(pressure);
            }
        }

    if (N > 0 && (m_dynamic[pgsd_flag::particles_energy] || m_nframes == 0))
        {
        ArrayHandle<Scalar> h_energy(m_pdata->getEnergies(),
                                     access_location::host,
                                     access_mode::read);

        frame.particle_data_present[pgsd_flag::particles_energy] = true;

        for (unsigned int index : m_index)
            {
            float energy = static_cast<float>(h_energy.data[index]);
            if (energy != 0.0f)
                {
                all_default[pgsd_flag::particles_energy] = false;
                }

            frame.particle_data.energy.push_back(energy);
            }
        }

    if (N > 0 && (m_dynamic[pgsd_flag::particles_aux1] || m_nframes == 0))
        {
        ArrayHandle<Scalar3> h_aux1(m_pdata->getAuxiliaries1(),
                                       access_location::host,
                                       access_mode::read);

        frame.particle_data_present[pgsd_flag::particles_aux1] = true;

        for (unsigned int index : m_index)
            {
            vec3<float> aux1 = vec3<float>(h_aux1.data[index]);

            if (aux1 != vec3<float>(0, 0, 0))
                {
                all_default[pgsd_flag::particles_aux1] = false;
                }

            frame.particle_data.aux1.push_back(aux1);
            }
        }


    if (N > 0 && (m_dynamic[pgsd_flag::particles_aux2] || m_nframes == 0))
        {
        ArrayHandle<Scalar3> h_aux2(m_pdata->getAuxiliaries2(),
                                       access_location::host,
                                       access_mode::read);

        frame.particle_data_present[pgsd_flag::particles_aux2] = true;

        for (unsigned int index : m_index)
            {
            vec3<float> aux2 = vec3<float>(h_aux2.data[index]);

            if (aux2 != vec3<float>(0, 0, 0))
                {
                all_default[pgsd_flag::particles_aux2] = false;
                }

            frame.particle_data.aux2.push_back(aux2);
            }
        }

    if (N > 0 && (m_dynamic[pgsd_flag::particles_aux2] || m_nframes == 0))
        {
        ArrayHandle<Scalar3> h_aux2(m_pdata->getAuxiliaries2(),
                                       access_location::host,
                                       access_mode::read);

        frame.particle_data_present[pgsd_flag::particles_aux2] = true;

        for (unsigned int index : m_index)
            {
            vec3<float> aux2 = vec3<float>(h_aux2.data[index]);

            if (aux2 != vec3<float>(0, 0, 0))
                {
                all_default[pgsd_flag::particles_aux2] = false;
                }

            frame.particle_data.aux2.push_back(aux2);
            }
        }

    if (N > 0 && (m_dynamic[pgsd_flag::particles_aux4] || m_nframes == 0))
        {
        ArrayHandle<Scalar3> h_aux4(m_pdata->getAuxiliaries4(),
                                       access_location::host,
                                       access_mode::read);

        frame.particle_data_present[pgsd_flag::particles_aux4] = true;

        for (unsigned int index : m_index)
            {
            vec3<float> aux4 = vec3<float>(h_aux4.data[index]);

            if (aux4 != vec3<float>(0, 0, 0))
                {
                all_default[pgsd_flag::particles_aux4] = false;
                }

            frame.particle_data.aux4.push_back(aux4);
            }
        }


    if (N > 0 && (m_dynamic[pgsd_flag::particles_body] || m_nframes == 0))
        {
        ArrayHandle<unsigned int> h_body(m_pdata->getBodies(),
                                         access_location::host,
                                         access_mode::read);

        frame.particle_data_present[pgsd_flag::particles_body] = true;

        for (unsigned int index : m_index)
            {
            unsigned int body = h_body.data[index];

            if (body != NO_BODY)
                {
                all_default[pgsd_flag::particles_body] = false;
                }

            frame.particle_data.body.push_back(body);
            }
        }



    unsigned long v = all_default.to_ulong();

    // All default only when all ranks are all default.
    MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_LONG, MPI_BAND, m_exec_conf->getMPICommunicator());

    all_default = std::bitset<n_pgsd_flags>(v);

    // Present when any rank is present
    v = frame.particle_data_present.to_ulong();

    MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_LONG, MPI_BOR, m_exec_conf->getMPICommunicator());

        frame.particle_data_present = std::bitset<n_pgsd_flags>(v);

    // Keep data in arrays only when they are not all default or this is a non-zero frame
    // and the zeroth frame is non-default. To not keep, resize the arrays back to 0.
    // !(!all_default || (nframes > 0 && m_nondefault["value"])) <=>
    // (all_default && !(nframes > 0 && m_nondefault["value"])

    if (all_default[pgsd_flag::particles_position]
        && !(m_nframes > 0 && m_nondefault["particles/position"]))
        {
        frame.particle_data.pos.resize(0);
        frame.particle_data_present[pgsd_flag::particles_position] = false;
        }

    if (all_default[pgsd_flag::particles_type]
        && !(m_nframes > 0 && m_nondefault["particles/typeid"]))
        {
        frame.particle_data.type.resize(0);
        frame.particle_data_present[pgsd_flag::particles_type] = false;
        }

    if (all_default[pgsd_flag::particles_mass] && !(m_nframes > 0 && m_nondefault["particles/mass"]))
        {
        frame.particle_data.mass.resize(0);
        frame.particle_data_present[pgsd_flag::particles_mass] = false;
        }

    if (all_default[pgsd_flag::particles_slength]
        && !(m_nframes > 0 && m_nondefault["particles/slength"]))
        {
        frame.particle_data.slength.resize(0);
        frame.particle_data_present[pgsd_flag::particles_slength] = false;
        }

    if (all_default[pgsd_flag::particles_density]
        && !(m_nframes > 0 && m_nondefault["particles/density"]))
        {
        frame.particle_data.density.resize(0);
        frame.particle_data_present[pgsd_flag::particles_density] = false;
        }

    if (all_default[pgsd_flag::particles_pressure]
        && !(m_nframes > 0 && m_nondefault["particles/pressure"]))
        {
        frame.particle_data.pressure.resize(0);
        frame.particle_data_present[pgsd_flag::particles_pressure] = false;
        }

    if (all_default[pgsd_flag::particles_energy]
        && !(m_nframes > 0 && m_nondefault["particles/energy"]))
        {
        frame.particle_data.energy.resize(0);
        frame.particle_data_present[pgsd_flag::particles_energy] = false;
        }

    if (all_default[pgsd_flag::particles_body] && !(m_nframes > 0 && m_nondefault["particles/body"]))
        {
        frame.particle_data.body.resize(0);
        frame.particle_data_present[pgsd_flag::particles_body] = false;
        }

    // momenta
    if (all_default[pgsd_flag::particles_velocity]
        && !(m_nframes > 0 && m_nondefault["particles/velocity"]))
        {
        frame.particle_data.vel.resize(0);
        frame.particle_data_present[pgsd_flag::particles_velocity] = false;
        }

    if (all_default[pgsd_flag::particles_aux1]
        && !(m_nframes > 0 && m_nondefault["particles/aux1"]))
        {
        frame.particle_data.aux1.resize(0);
        frame.particle_data_present[pgsd_flag::particles_aux1] = false;
        }

    if (all_default[pgsd_flag::particles_aux2]
        && !(m_nframes > 0 && m_nondefault["particles/aux2"]))
        {
        frame.particle_data.aux2.resize(0);
        frame.particle_data_present[pgsd_flag::particles_aux2] = false;
        }

    if (all_default[pgsd_flag::particles_aux3]
        && !(m_nframes > 0 && m_nondefault["particles/aux3"]))
        {
        frame.particle_data.aux3.resize(0);
        frame.particle_data_present[pgsd_flag::particles_aux3] = false;
        }

    if (all_default[pgsd_flag::particles_aux4]
        && !(m_nframes > 0 && m_nondefault["particles/aux4"]))
        {
        frame.particle_data.aux4.resize(0);
        frame.particle_data_present[pgsd_flag::particles_aux4] = false;
        }

    if (all_default[pgsd_flag::particles_image]
        && !(m_nframes > 0 && m_nondefault["particles/image"]))
        {
        frame.particle_data.image.resize(0);
        frame.particle_data_present[pgsd_flag::particles_image] = false;
        }

    // // capture topology data
    // if (m_group->getNumMembersGlobal() != m_pdata->getNGlobal() && m_write_topology)
    //     {
    //     throw std::runtime_error("Cannot write topology for a portion of the system");
    //     }

    // if (m_group->getNumMembersGlobal() == m_pdata->getNGlobal()
    //     && (m_write_topology || m_nframes == 0))
    //     {
    //     printf("write topology");
    //     m_sysdef->getBondData()->takeSnapshotDistr(frame.bond_data);
    //     // m_sysdef->getAngleData()->takeSnapshot(frame.angle_data);
    //     // m_sysdef->getDihedralData()->takeSnapshot(frame.dihedral_data);
    //     // m_sysdef->getImproperData()->takeSnapshot(frame.improper_data);
    //     m_sysdef->getConstraintData()->takeSnapshotDistr(frame.constraint_data);
    //     // m_sysdef->getPairData()->takeSnapshot(frame.pair_data);
    //     }
    }

//#ifdef ENABLE_MPI

/*! Gather per-particle data from the local frame and sort it into ascending tag order in
    m_global_frame.
*/
// void GSDDumpWriterMPI::gatherGlobalFrame(const PGSDFrame& local_frame)
//     {
//     m_global_frame.clear();

//     m_global_frame.timestep = local_frame.timestep;
//     m_global_frame.global_box = local_frame.global_box;
//     m_global_frame.particle_data.type_mapping = local_frame.particle_data.type_mapping;
//     m_global_frame.particle_data_present = local_frame.particle_data_present;

//     m_gather_tag_order.setLocalTagsSorted(local_frame.particle_tags);

//     if (local_frame.particle_data_present[pgsd_flag::particles_position])
//         {
//         m_gather_tag_order.gatherArray(m_global_frame.particle_data.pos,
//                                        local_frame.particle_data.pos);
//         }

//     // if (local_frame.particle_data_present[pgsd_flag::particles_orientation])
//     //     {
//     //     m_gather_tag_order.gatherArray(m_global_frame.particle_data.orientation,
//     //                                    local_frame.particle_data.orientation);
//     //     }
//     if (local_frame.particle_data_present[pgsd_flag::particles_type])
//         {
//         m_gather_tag_order.gatherArray(m_global_frame.particle_data.type,
//                                        local_frame.particle_data.type);
//         }
//     if (local_frame.particle_data_present[pgsd_flag::particles_mass])
//         {
//         m_gather_tag_order.gatherArray(m_global_frame.particle_data.mass,
//                                        local_frame.particle_data.mass);
//         }
//     // if (local_frame.particle_data_present[pgsd_flag::particles_charge])
//     //     {
//     //     m_gather_tag_order.gatherArray(m_global_frame.particle_data.charge,
//     //                                    local_frame.particle_data.charge);
//     //     }
//     if (local_frame.particle_data_present[pgsd_flag::particles_slength])
//         {
//         m_gather_tag_order.gatherArray(m_global_frame.particle_data.slength,
//                                        local_frame.particle_data.slength);
//         }
//     if (local_frame.particle_data_present[pgsd_flag::particles_density])
//         {
//         m_gather_tag_order.gatherArray(m_global_frame.particle_data.density,
//                                        local_frame.particle_data.density);
//         }
//     if (local_frame.particle_data_present[pgsd_flag::particles_pressure])
//         {
//         m_gather_tag_order.gatherArray(m_global_frame.particle_data.pressure,
//                                        local_frame.particle_data.pressure);
//         }
//     if (local_frame.particle_data_present[pgsd_flag::particles_energy])
//         {
//         m_gather_tag_order.gatherArray(m_global_frame.particle_data.energy,
//                                        local_frame.particle_data.energy);
//         }
//     // if (local_frame.particle_data_present[pgsd_flag::particles_diameter])
//     //     {
//     //     m_gather_tag_order.gatherArray(m_global_frame.particle_data.diameter,
//     //                                    local_frame.particle_data.diameter);
//     //     }
//     if (local_frame.particle_data_present[pgsd_flag::particles_body])
//         {
//         m_gather_tag_order.gatherArray(m_global_frame.particle_data.body,
//                                        local_frame.particle_data.body);
//         }
//     // if (local_frame.particle_data_present[pgsd_flag::particles_inertia])
//     //     {
//     //     m_gather_tag_order.gatherArray(m_global_frame.particle_data.inertia,
//     //                                    local_frame.particle_data.inertia);
//     //     }
//     if (local_frame.particle_data_present[pgsd_flag::particles_velocity])
//         {
//         m_gather_tag_order.gatherArray(m_global_frame.particle_data.vel,
//                                        local_frame.particle_data.vel);
//         }
//     if (local_frame.particle_data_present[pgsd_flag::particles_aux1])
//         {
//         m_gather_tag_order.gatherArray(m_global_frame.particle_data.aux1,
//                                        local_frame.particle_data.aux1);
//         }
//     if (local_frame.particle_data_present[pgsd_flag::particles_aux2])
//         {
//         m_gather_tag_order.gatherArray(m_global_frame.particle_data.aux2,
//                                        local_frame.particle_data.aux2);
//         }
//     if (local_frame.particle_data_present[pgsd_flag::particles_aux3])
//         {
//         m_gather_tag_order.gatherArray(m_global_frame.particle_data.aux3,
//                                        local_frame.particle_data.aux3);
//         }
//     if (local_frame.particle_data_present[pgsd_flag::particles_aux4])
//         {
//         m_gather_tag_order.gatherArray(m_global_frame.particle_data.aux4,
//                                        local_frame.particle_data.aux4);
//         }
//     // if (local_frame.particle_data_present[pgsd_flag::particles_angmom])
//     //     {
//     //     m_gather_tag_order.gatherArray(m_global_frame.particle_data.angmom,
//     //                                    local_frame.particle_data.angmom);
//     //     }
//     if (local_frame.particle_data_present[pgsd_flag::particles_image])
//         {
//         m_gather_tag_order.gatherArray(m_global_frame.particle_data.image,
//                                        local_frame.particle_data.image);
//         }
//     }

//#endif

namespace detail
    {
void export_GSDDumpWriterMPI(pybind11::module& m)
    {
    // pybind11::bind_map<std::map<std::string, pybind11::function>>(m, "MapStringFunctionMPI");

    pybind11::class_<GSDDumpWriterMPI, Analyzer, std::shared_ptr<GSDDumpWriterMPI>>(m, "GSDDumpWriterMPI")
        .def(pybind11::init<std::shared_ptr<SystemDefinition>,
                            std::shared_ptr<Trigger>,
                            std::string,
                            std::shared_ptr<ParticleGroup>,
                            std::string,
                            bool>())
        .def_property("log_writer", &GSDDumpWriterMPI::getLogWriter, &GSDDumpWriterMPI::setLogWriter)
        .def_property_readonly("filename", &GSDDumpWriterMPI::getFilename)
        .def_property_readonly("mode", &GSDDumpWriterMPI::getMode)
        .def_property("dynamic", &GSDDumpWriterMPI::getDynamic, &GSDDumpWriterMPI::setDynamic)
        // .def_property_readonly("truncate", &GSDDumpWriterMPI::getTruncate)
        .def_property_readonly("filter",
                               [](const std::shared_ptr<GSDDumpWriterMPI> gsd)
                               { return gsd->getGroup()->getFilter(); })
        .def_property("write_diameter",
                      &GSDDumpWriterMPI::getWriteDiameter,
                      &GSDDumpWriterMPI::setWriteDiameter)
        .def("flush", &GSDDumpWriterMPI::flush)
        .def_property("maximum_write_buffer_size",
                      &GSDDumpWriterMPI::getMaximumWriteBufferSize,
                      &GSDDumpWriterMPI::setMaximumWriteBufferSize);
    }

    } // end namespace detail

    } // end namespace hoomd
