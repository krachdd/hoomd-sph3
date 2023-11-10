# Copyright (c) 2016-2023 The Regents of the University of Michigan
# Part of PGSD, released under the BSD 2-Clause License.

from libc.stdint cimport uint8_t, int8_t, uint16_t, int16_t, uint32_t, int32_t,\
    uint64_t, int64_t
from mpi4py.MPI cimport MPI_Comm, MPI_Info, Comm, Info, MPI_File
from libcpp cimport bool

cdef extern from "pgsd.h" nogil:
    cdef enum pgsd_type:
        PGSD_TYPE_UINT8=1
        PGSD_TYPE_UINT16
        PGSD_TYPE_UINT32
        PGSD_TYPE_UINT64
        PGSD_TYPE_INT8
        PGSD_TYPE_INT16
        PGSD_TYPE_INT32
        PGSD_TYPE_INT64
        PGSD_TYPE_FLOAT
        PGSD_TYPE_DOUBLE

    cdef enum pgsd_open_flag:
        PGSD_OPEN_READWRITE=1
        PGSD_OPEN_READONLY
        PGSD_OPEN_APPEND

    cdef enum pgsd_error:
        PGSD_SUCCESS = 0
        PGSD_ERROR_IO = -1
        PGSD_ERROR_INVALID_ARGUMENT = -2
        PGSD_ERROR_NOT_A_PGSD_FILE = -3
        PGSD_ERROR_INVALID_PGSD_FILE_VERSION = -4
        PGSD_ERROR_FILE_CORRUPT = -5
        PGSD_ERROR_MEMORY_ALLOCATION_FAILED = -6
        PGSD_ERROR_NAMELIST_FULL = -7
        PGSD_ERROR_FILE_MUST_BE_WRITABLE = -8
        PGSD_ERROR_FILE_MUST_BE_READABLE = -9

    cdef struct pgsd_header:
        uint64_t magic
        uint32_t pgsd_version
        char application[64]
        char schema[64]
        uint32_t schema_version
        uint64_t index_location
        uint64_t index_allocated_entries
        uint64_t namelist_location
        uint64_t namelist_allocated_entries
        char reserved[80]

    cdef struct pgsd_index_entry:
        uint64_t frame
        uint64_t N
        int64_t location
        uint32_t M
        uint16_t id
        uint8_t type
        uint8_t flags

    cdef struct pgsd_namelist_entry:
        char name[64]

    cdef struct pgsd_index_buffer:
        pgsd_index_entry *data
        size_t size
        size_t reserved
        void *mapped_data
        size_t mapped_len

    cdef struct pgsd_name_id_map:
        void *v
        size_t size

    cdef struct pgsd_write_buffer:
        char *data
        size_t size
        size_t reserved

    cdef struct pgsd_handle:
        # int fd
        MPI_File fh
        pgsd_header header
        pgsd_index_buffer file_index
        pgsd_index_buffer frame_index
        pgsd_index_buffer buffer_index
        pgsd_write_buffer write_buffer
        pgsd_namelist_entry *namelist
        uint64_t namelist_num_entries
        uint64_t cur_frame
        long long int file_size
        pgsd_open_flag open_flags
        pgsd_name_id_map name_map
        uint64_t namelist_written_entries

    uint32_t pgsd_make_version(unsigned int major, unsigned int minor)
    int pgsd_create(const char *fname,
                   const char *application,
                   const char *schema,
                   uint32_t schema_version)
    int pgsd_create_and_open(pgsd_handle* handle,
                            const char *fname,
                            const char *application,
                            const char *schema,
                            uint32_t schema_version,
                            const pgsd_open_flag flags,
                            int exclusive_create)
    int pgsd_open(pgsd_handle* handle, const char *fname,
                 const pgsd_open_flag flags)
    int pgsd_truncate(pgsd_handle* handle)
    int pgsd_close(pgsd_handle* handle, bool all)
    int pgsd_end_frame(pgsd_handle* handle, bool all)
    int pgsd_flush(pgsd_handle* handle, bool all)
    int pgsd_write_chunk(pgsd_handle* handle,
                        const char *name,
                        pgsd_type type,
                        uint64_t N,
                        uint32_t M,
                        uint64_t N_global,
                        uint32_t M_global,
                        uint64_t offset,
                        uint64_t global_size,
                        bool all,
                        uint8_t flags,
                        const void *data)
    const pgsd_index_entry* pgsd_find_chunk(pgsd_handle* handle,
                                          uint64_t frame,
                                          const char *name,
                                          bool all)
    int pgsd_read_chunk(pgsd_handle* handle, void* data,
                       const pgsd_index_entry* chunk,
                       unsigned int *offset)
    uint64_t pgsd_get_nframes(pgsd_handle* handle)
    size_t pgsd_sizeof_type(pgsd_type type)
    const char *pgsd_find_matching_chunk_name(pgsd_handle* handle,
                                             const char *match,
                                             const char *prev,
                                             bool all)
    # int pgsd_upgrade(pgsd_handle *handle)
    uint64_t pgsd_get_maximum_write_buffer_size(pgsd_handle* handle)
    int pgsd_set_maximum_write_buffer_size(pgsd_handle* handle, uint64_t size)
    uint64_t pgsd_get_index_entries_to_buffer(pgsd_handle* handle)
    int pgsd_set_index_entries_to_buffer(pgsd_handle* handle, uint64_t number)
