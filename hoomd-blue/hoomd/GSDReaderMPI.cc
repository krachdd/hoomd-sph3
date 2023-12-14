// Copyright (c) 2009-2023 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#include "GSDReaderMPI.h"
#include "ExecutionConfiguration.h"
#include "PGSD.h"
#include "SnapshotSystemData.h"
#include "hoomd/extern/pgsd.h"
#include <sstream>
#include <string.h>

#include <stdexcept>
using namespace std;
using namespace hoomd::detail;

namespace hoomd
    {
/*! \param exec_conf The execution configuration
    \param name File name to read
    \param frame Frame index to read from the file
    \param from_end Count frames back from the end of the file

    The GSDReaderMPI constructor opens the GSD file, initializes an empty snapshot, and reads the file
   into memory (on the root rank).
*/
GSDReaderMPI::GSDReaderMPI(std::shared_ptr<const ExecutionConfiguration> exec_conf,
                     const std::string& name,
                     const uint64_t frame,
                     bool from_end)
    : m_exec_conf(exec_conf), m_timestep(0), m_name(name), m_frame(frame)
    {
    m_snapshot = std::shared_ptr<SnapshotSystemData<float>>(new SnapshotSystemData<float>);

    // open the GSD file in read mode
    m_exec_conf->msg->notice(3) << "data.pgsd_snapshot: open gsd file " << name << endl;
    int retval = pgsd_open(&m_handle, name.c_str(), PGSD_OPEN_READONLY);
    PGSDUtils::checkError(retval, m_name);

    // validate schema
    if (string(m_handle.header.schema) != string("hoomd"))
        {
        std::ostringstream s;
        s << "Invalid schema in " << name << endl;
        throw runtime_error(s.str());
        }
    if (m_handle.header.schema_version >= pgsd_make_version(2, 1))
        {
        std::ostringstream s;
        s << "Invalid schema version in " << name << endl;
        throw runtime_error(s.str());
        }

    // set frame from the end of the file if requested
    uint64_t nframes = pgsd_get_nframes(&m_handle);
    if (from_end && frame <= nframes)
        m_frame = nframes - frame;

    // validate number of frames
    if (m_frame >= nframes)
        {
        std::ostringstream s;
        s << "Cannot read frame " << m_frame << " " << name << " only has "
          << pgsd_get_nframes(&m_handle) << " frames.";
        throw runtime_error(s.str());
        }
    readHeader();
    readParticles();
    readTopology();
    }

GSDReaderMPI::~GSDReaderMPI()
    {
    pgsd_close(&m_handle);
    }

/*! \param data Pointer to data to read into
    \param frame Frame index to read from
    \param name Name of the data chunk
    \param expected_size Expected size of the data chunk in bytes.
    \param cur_n N in the current frame.

    Attempts to read the data chunk of the given name at the given frame. If it is not present at
   this frame, attempt to read from frame 0. If it is also not present at frame 0, return false. If
   the found data chunk is not the expected size, throw an exception.

    Per the GSD spec, keep the default when the frame 0 N does not match the current N.

    Return true if data is actually read from the file.
*/
bool GSDReaderMPI::readChunk(void* data,
                          uint64_t frame,
                          const char* name,
                          size_t expected_size,
                          unsigned int cur_n,
                          uint64_t N_local, 
                          uint32_t M_local,
                          uint32_t offset, 
                          bool all)
    {

    const struct pgsd_index_entry* entry = pgsd_find_chunk(&m_handle, frame, name);
    bool empty_entry_indicator = false;

    if (entry == NULL && frame != 0 && is_root() )
        {
        empty_entry_indicator = true;
        }
    MPI_Bcast(&empty_entry_indicator, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);

    if ( empty_entry_indicator == true)
        {
        entry = pgsd_find_chunk(&m_handle, 0, name);
        }

    empty_entry_indicator = false;
    if ( entry == NULL && is_root() )
        empty_entry_indicator = true;

    MPI_Bcast(&empty_entry_indicator, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
    if ( empty_entry_indicator == true )
        {
        m_exec_conf->msg->notice(10) << "data.pgsd_snapshot: empty entry -> chunk not found " << name << endl;
        return false;
        }

    size_t m_N;
    uint8_t m_type;

    if ( m_exec_conf->getRank() == 0 )
        {
        m_N = entry->N;
        m_type = entry->type;
        }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&m_N, 1, my_MPI_SIZE_T, 0, MPI_COMM_WORLD);
    MPI_Bcast(&m_type, 1, MPI_UINT8_T, 0, MPI_COMM_WORLD);

    if ( cur_n != 0 && m_N != cur_n )
        {
        m_exec_conf->msg->notice(10) << "data.pgsd_snapshot: chunk not found " << name << endl;
        return false;
        }

    else
        {
        m_exec_conf->msg->notice(7) << "data.pgsd_snapshot: reading chunk " << name << endl;
        size_t actual_size = N_local * M_local * pgsd_sizeof_type((enum pgsd_type)m_type);
        if (actual_size != expected_size)
            {
            std::ostringstream s;
            s << "Expecting " << expected_size << " bytes in " << name << " but found "
              << actual_size << ".";
            throw runtime_error(s.str());
            }

        int retval = pgsd_read_chunk(&m_handle, data, entry, N_local, M_local, offset, all);
        PGSDUtils::checkError(retval, m_name);

        return true;
        }
    }

/*! \param frame Frame index to read from
    \param name Name of the data chunk

   Attempts to read the data chunk of the given name at the given frame. If it is not present at
   this frame, attempt to read from frame 0. If it is also not present at frame 0, return an empty
   list.

    If the data chunk is found in the file, return a vector of string type names.
*/
std::vector<std::string> GSDReaderMPI::readTypes(uint64_t frame, const char* name)
    {
    m_exec_conf->msg->notice(7) << "data.pgsd_snapshot: reading chunk " << name << endl;

    std::vector<std::string> type_mapping;

    // set the default particle type mapping per the GSD HOOMD Schema
    if (std::string(name) == "particles/types")
        type_mapping.push_back("A");

    const struct pgsd_index_entry* entry = pgsd_find_chunk(&m_handle, frame, name);

    bool empty_entry_indicator = false;
    if (entry == NULL && frame != 0 && is_root() )
        empty_entry_indicator = true;
    MPI_Bcast(&empty_entry_indicator, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);

    if ( empty_entry_indicator == true )
        entry = pgsd_find_chunk(&m_handle, 0, name);
    
    empty_entry_indicator = false;
    
    if (entry == NULL && is_root() )
        empty_entry_indicator = true;

    MPI_Bcast(&empty_entry_indicator, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);

    if ( empty_entry_indicator == true )
        return type_mapping;
    else
        {
        size_t m_N;
        uint32_t m_M;
        uint8_t m_type;

        if ( is_root() )
            {
            m_N = entry->N;
            m_M = entry->M;
            m_type = entry->type;
            }
        MPI_Bcast(&m_N, 1, my_MPI_SIZE_T, 0, MPI_COMM_WORLD);
        MPI_Bcast(&m_M, 1, my_MPI_SIZE_T, 0, MPI_COMM_WORLD);
        MPI_Bcast(&m_type, 1, MPI_UINT8_T, 0, MPI_COMM_WORLD);
        
        size_t actual_size = m_N * m_M * pgsd_sizeof_type((enum pgsd_type)m_type);
        
        std::vector<char> data(actual_size);
        int retval = pgsd_read_chunk(&m_handle, &data[0], entry, m_N, m_M, 0, true);
        PGSDUtils::checkError(retval, m_name);

        type_mapping.clear();
        for (unsigned int i = 0; i < m_N; i++)
            {
            size_t l = strnlen(&data[i * m_M], m_M);
            type_mapping.push_back(std::string(&data[i * m_M], l));
            }
        return type_mapping;
        }
    }

/*! Read the same data chunks written by GSDDumpWriter::writeFrameHeader
 */
void GSDReaderMPI::readHeader()
    {
    int nprocs, rank;
    rank = m_exec_conf->getRank();
    nprocs = m_exec_conf->getNRanksGlobal();

    uint32_t N_local = 1;
    uint32_t M_local = 1;
    unsigned int cur_n = 0;
    readChunk(&m_timestep, m_frame, "configuration/step", 8, cur_n, N_local, M_local, 0, true);
    MPI_Barrier(MPI_COMM_WORLD);
    uint8_t dim = 3;
    readChunk(&dim, m_frame, "configuration/dimensions", 1, cur_n, N_local, M_local, 0, true);
    MPI_Barrier(MPI_COMM_WORLD);

    m_snapshot->dimensions = dim;

    float box[6] = {1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f};
    readChunk(&box, m_frame, "configuration/box", 6 * 4, cur_n, 6 * N_local, M_local, 0, true);
    // Set Lz, xz, and yz to 0 for 2D boxes. Needed for working with hoomd v 2 GSD files.
    if (dim == 2)
        {
        box[2] = 0;
        box[4] = 0;
        box[5] = 0;
        }
    m_snapshot->global_box = std::make_shared<BoxDim>(BoxDim(box[0], box[1], box[2]));
    m_snapshot->global_box->setTiltFactors(box[3], box[4], box[5]);

    unsigned int N = 0;
    readChunk(&N, m_frame, "particles/N", 4, cur_n, N_local, M_local, 0, true);
    if (N == 0)
        {
        std::ostringstream s;
        s << "Cannot read a file with 0 particles.";
        throw runtime_error(s.str());
        }


    m_part_per_rank.resize( nprocs );

    unsigned int n = floor( N / nprocs );
    // std::fill(m_part_per_rank.begin(),m_part_per_rank.end(), n);
    unsigned int rem = N%nprocs;
    
    if( rank < rem ){
        n++;
    }

    all_gather_v(n, m_part_per_rank, MPI_COMM_WORLD);

    m_snapshot->particle_data.resize(n);
    }

/*! Read the same data chunks for particles
 */
void GSDReaderMPI::readParticles()
    {
    uint32_t N = m_snapshot->particle_data.size;

    unsigned int rank = m_exec_conf->getRank();
    unsigned int size = m_exec_conf->getNRanksGlobal();
    int offset;

    std::vector<unsigned int> part_distribution(size);
    all_gather_v(N, part_distribution, MPI_COMM_WORLD);
    offset = std::accumulate(part_distribution.begin(), part_distribution.begin()+rank, 0);
    uint32_t N_global = std::accumulate(part_distribution.begin(), part_distribution.end(),0);

    m_snapshot->particle_data.type_mapping = readTypes(m_frame, "particles/types");
    MPI_Barrier(MPI_COMM_WORLD);


    // the snapshot already has default values, if a chunk is not found, the value
    // is already at the default, and the failed read is not a problem
    readChunk(&m_snapshot->particle_data.type[0], m_frame, "particles/typeid", N * 4, N_global, N, 1, offset, true);
    readChunk(&m_snapshot->particle_data.mass[0], m_frame, "particles/mass", N * 4, N_global, N, 1, offset, true);
    readChunk(&m_snapshot->particle_data.slength[0], m_frame, "particles/slength", N * 4, N_global, N, 1, offset, true);
    readChunk(&m_snapshot->particle_data.body[0], m_frame, "particles/body", N * 4, N_global, N, 1, offset, true);
    readChunk(&m_snapshot->particle_data.pos[0], m_frame, "particles/position", N * 12, N_global, N, 3, offset, true);
    readChunk(&m_snapshot->particle_data.vel[0], m_frame, "particles/velocity", N * 12, N_global, N, 3, offset, true);
    readChunk(&m_snapshot->particle_data.density[0], m_frame, "particles/density", N * 4, N_global, N, 1, offset, true);
    readChunk(&m_snapshot->particle_data.pressure[0], m_frame, "particles/pressure", N * 4, N_global, N, 1, offset, true);
    readChunk(&m_snapshot->particle_data.energy[0], m_frame, "particles/energy", N * 4, N_global, N, 1, offset, true);
    readChunk(&m_snapshot->particle_data.aux1[0], m_frame, "particles/auxiliary1", N * 12, N_global, N, 3, offset, true);
    readChunk(&m_snapshot->particle_data.aux2[0], m_frame, "particles/auxiliary2", N * 12, N_global, N, 3, offset, true);
    readChunk(&m_snapshot->particle_data.aux3[0], m_frame, "particles/auxiliary3", N * 12, N_global, N, 3, offset, true);
    readChunk(&m_snapshot->particle_data.aux4[0], m_frame, "particles/auxiliary4", N * 12, N_global, N, 3, offset, true);
    readChunk(&m_snapshot->particle_data.image[0], m_frame, "particles/image", N * 12, N_global, N, 3, offset, true);
    }

/*! Read the same data chunks for topology
 */
void GSDReaderMPI::readTopology()
    {
    // TODO 
    // unsigned int N = 0;
    // readChunk(&N, m_frame, "bonds/N", 4);
    // if (N > 0)
    //     {
    //     m_snapshot->bond_data.resize(N);
    //     m_snapshot->bond_data.type_mapping = readTypes(m_frame, "bonds/types");
    //     readChunk(&m_snapshot->bond_data.type_id[0], m_frame, "bonds/typeid", N * 4, N);
    //     readChunk(&m_snapshot->bond_data.groups[0], m_frame, "bonds/group", N * 8, N);
    //     }

    // // N = 0;
    // // readChunk(&N, m_frame, "angles/N", 4);
    // // if (N > 0)
    // //     {
    // //     m_snapshot->angle_data.resize(N);
    // //     m_snapshot->angle_data.type_mapping = readTypes(m_frame, "angles/types");
    // //     readChunk(&m_snapshot->angle_data.type_id[0], m_frame, "angles/typeid", N * 4, N);
    // //     readChunk(&m_snapshot->angle_data.groups[0], m_frame, "angles/group", N * 12, N);
    // //     }

    // // N = 0;
    // // readChunk(&N, m_frame, "dihedrals/N", 4);
    // // if (N > 0)
    // //     {
    // //     m_snapshot->dihedral_data.resize(N);
    // //     m_snapshot->dihedral_data.type_mapping = readTypes(m_frame, "dihedrals/types");
    // //     readChunk(&m_snapshot->dihedral_data.type_id[0], m_frame, "dihedrals/typeid", N * 4, N);
    // //     readChunk(&m_snapshot->dihedral_data.groups[0], m_frame, "dihedrals/group", N * 16, N);
    // //     }

    // // N = 0;
    // // readChunk(&N, m_frame, "impropers/N", 4);
    // // if (N > 0)
    // //     {
    // //     m_snapshot->improper_data.resize(N);
    // //     m_snapshot->improper_data.type_mapping = readTypes(m_frame, "impropers/types");
    // //     readChunk(&m_snapshot->improper_data.type_id[0], m_frame, "impropers/typeid", N * 4, N);
    // //     readChunk(&m_snapshot->improper_data.groups[0], m_frame, "impropers/group", N * 16, N);
    // //     }

    // N = 0;
    // readChunk(&N, m_frame, "constraints/N", 4);
    // if (N > 0)
    //     {
    //     m_snapshot->constraint_data.resize(N);
    //     std::vector<float> data(N);
    //     readChunk(&data[0], m_frame, "constraints/value", N * 4, N);
    //     for (unsigned int i = 0; i < N; i++)
    //         m_snapshot->constraint_data.val[i] = Scalar(data[i]);

    //     readChunk(&m_snapshot->constraint_data.groups[0], m_frame, "constraints/group", N * 8, N);
    //     }

    // // if (m_handle.header.schema_version >= gsd_make_version(1, 1))
    // //     {
    // //     N = 0;
    // //     readChunk(&N, m_frame, "pairs/N", 4);
    // //     if (N > 0)
    // //         {
    // //         m_snapshot->pair_data.resize(N);
    // //         m_snapshot->pair_data.type_mapping = readTypes(m_frame, "pairs/types");
    // //         readChunk(&m_snapshot->pair_data.type_id[0], m_frame, "pairs/typeid", N * 4, N);
    // //         readChunk(&m_snapshot->pair_data.groups[0], m_frame, "pairs/group", N * 8, N);
    // //         }
    // //     }
    }

pybind11::list GSDReaderMPI::readTypeShapesPy(uint64_t frame)
    {
    std::vector<std::string> type_mapping = this->readTypes(frame, "particles/type_shapes");
    pybind11::list type_shapes;
    for (unsigned int i = 0; i < type_mapping.size(); i++)
        type_shapes.append(type_mapping[i]);
    return type_shapes;
    }

namespace detail
    {
void export_GSDReaderMPI(pybind11::module& m)
    {
    pybind11::class_<GSDReaderMPI, std::shared_ptr<GSDReaderMPI>>(m, "GSDReaderMPI")
        .def(pybind11::init<std::shared_ptr<const ExecutionConfiguration>,
                            const string&,
                            const uint64_t,
                            bool>())
        .def("getTimeStep", &GSDReaderMPI::getTimeStep)
        .def("getSnapshot", &GSDReaderMPI::getSnapshot)
        .def("clearSnapshot", &GSDReaderMPI::clearSnapshot)
        .def("readTypeShapesPy", &GSDReaderMPI::readTypeShapesPy);
    }

    } // end namespace detail

    } // end namespace hoomd
