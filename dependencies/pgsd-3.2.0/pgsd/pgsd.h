// Copyright (c) 2016-2023 The Regents of the University of Michigan
// Part of PGSD, released under the BSD 2-Clause License.

#ifndef PGSD_H
#define PGSD_H

// #ifdef ENABLE_MPI
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <mpi.h>

#ifdef __cplusplus
extern "C"
    {
#endif

    /*! \file pgsd.h
        \brief Declare PGSD data types and C API
    */

    /// Identifiers for the pgsd data chunk element types
    enum pgsd_type
        {
        /// Unsigned 8-bit integer.
        PGSD_TYPE_UINT8 = 1,

        /// Unsigned 16-bit integer.
        PGSD_TYPE_UINT16,

        /// Unsigned 32-bit integer.
        PGSD_TYPE_UINT32,

        /// Unsigned 53-bit integer.
        PGSD_TYPE_UINT64,

        /// Signed 8-bit integer.
        PGSD_TYPE_INT8,

        /// Signed 16-bit integer.
        PGSD_TYPE_INT16,

        /// Signed 32-bit integer.
        PGSD_TYPE_INT32,

        /// Signed 64-bit integer.
        PGSD_TYPE_INT64,

        /// 32-bit floating point number.
        PGSD_TYPE_FLOAT,

        /// 64-bit floating point number.
        PGSD_TYPE_DOUBLE
        };

    /// Flag for PGSD file open options
    enum pgsd_open_flag
        {
        /// Open for both reading and writing
        PGSD_OPEN_READWRITE = 1,

        /// Open only for reading
        PGSD_OPEN_READONLY,

        /// Open only for writing
        PGSD_OPEN_APPEND
        };

    /// Error return values
    enum pgsd_error
        {
        /// Success.
        PGSD_SUCCESS = 0,

        /// IO error. Check ``errno`` for details
        PGSD_ERROR_IO = -1,

        /// Invalid argument passed to function.
        PGSD_ERROR_INVALID_ARGUMENT = -2,

        /// The file is not a PGSD file.
        PGSD_ERROR_NOT_A_PGSD_FILE = -3,

        /// The PGSD file version cannot be read.
        PGSD_ERROR_INVALID_PGSD_FILE_VERSION = -4,

        /// The PGSD file is corrupt.
        PGSD_ERROR_FILE_CORRUPT = -5,

        /// PGSD failed to allocated memory.
        PGSD_ERROR_MEMORY_ALLOCATION_FAILED = -6,

        /// The PGSD file cannot store any additional unique data chunk names.
        PGSD_ERROR_NAMELIST_FULL = -7,

        /** This API call requires that the PGSD file opened in with the mode PGSD_OPEN_APPEND or
            PGSD_OPEN_READWRITE.
        */
        PGSD_ERROR_FILE_MUST_BE_WRITABLE = -8,

        /** This API call requires that the PGSD file opened the mode PGSD_OPEN_READ or
            PGSD_OPEN_READWRITE.
        */
        PGSD_ERROR_FILE_MUST_BE_READABLE = -9,
        };

    enum
        {
        /** v1 file: Size of a PGSD name in memory. v2 file: The name buffer size is a multiple of
            PGSD_NAME_SIZE.
        */
        PGSD_NAME_SIZE = 64
        };

    enum
        {
        /// Reserved bytes in the header structure
        PGSD_RESERVED_BYTES = 80
        };

    /** PGSD file header

        The in-memory and on-disk storage of the PGSD file header. Stored in the first 256 bytes of
        the file.

        @warning All members are **read-only** to the caller.
    */
    struct pgsd_header
        {
        /// Magic number marking that this is a PGSD file.
        uint64_t magic;

        /// Location of the chunk index in the file.
        uint64_t index_location;

        /// Number of index entries that will fit in the space allocated.
        uint64_t index_allocated_entries;

        /// Location of the name list in the file.
        uint64_t namelist_location;

        /// Number of bytes in the namelist divided by PGSD_NAME_SIZE.
        uint64_t namelist_allocated_entries;

        /// Schema version: from pgsd_make_version().
        uint32_t schema_version;

        /// PGSD file format version from pgsd_make_version().
        uint32_t pgsd_version;

        /// Name of the application that generated this file.
        char application[PGSD_NAME_SIZE];

        /// Name of data schema.
        char schema[PGSD_NAME_SIZE];

        /// Reserved for future use.
        char reserved[PGSD_RESERVED_BYTES];
        };

    /** Index entry

        An index entry for a single chunk of data.

        @warning All members are **read-only** to the caller.
    */
    struct pgsd_index_entry
        {
        /// Frame index of the chunk.
        uint64_t frame;

        /// Number of rows in the chunk.
        uint64_t N;

        /// Location of the chunk in the file.
        int64_t location;

        /// Number of columns in the chunk.
        uint32_t M;

        /// Index of the chunk name in the name list.
        uint16_t id;

        /// Data type of the chunk: one of pgsd_type.
        uint8_t type;

        /// Flags (for internal use).
        uint8_t flags;
        };

    /** Name/id mapping

        A string name paired with an ID. Used for storing sorted name/id mappings in a hash map.
    */
    struct pgsd_name_id_pair
        {
        /// Pointer to name (actual name storage is allocated in pgsd_handle)
        char* name;

        /// Next name/id pair with the same hash
        struct pgsd_name_id_pair* next;

        /// Entry id
        uint16_t id;
        };

    /** Name/id hash map

        A hash map of string names to integer identifiers.
    */
    struct pgsd_name_id_map
        {
        /// Name/id mappings
        struct pgsd_name_id_pair* v;

        /// Number of entries in the mapping
        size_t size;
        };

    /** Array of index entries

        May point to a mapped location of index entries in the file or an in-memory buffer.
    */
    struct pgsd_index_buffer
        {
        /// Indices in the buffer
        struct pgsd_index_entry* data;

        /// Number of entries in the buffer
        size_t size;

        /// Number of entries available in the buffer
        size_t reserved;

        /// Pointer to mapped data (NULL if not mapped)
        void* mapped_data;

        /// Number of bytes mapped
        size_t mapped_len;
        };

    /** Byte buffer

        Used to buffer of small data chunks held for a buffered write at the end of a frame. Also
        used to hold the names.
    */
    struct pgsd_byte_buffer
        {
        /// Data
        char* data;

        /// Number of bytes in the buffer
        size_t size;

        /// Number of bytes available in the buffer
        size_t reserved;
        };

    /** Name buffer

        Holds a list of string names in order separated by NULL terminators. In v1 files, each name
        is 64 bytes. In v2 files, only one NULL terminator is placed between each name.
    */
    struct pgsd_name_buffer
        {
        /// Data
        struct pgsd_byte_buffer data;

        /// Number of names in the list
        size_t n_names;
        };

    /** File handle

        A handle to an open PGSD file.

        This handle is obtained when opening a PGSD file and is passed into every method that
        operates on the file.

        @warning All members are **read-only** to the caller.
    */
    struct pgsd_handle
        {
        /// File descriptor
        // int fd;

        /// MPI File pointer
        MPI_File fh;

        /// The file header
        struct pgsd_header header;

        /// Mapped data chunk index
        struct pgsd_index_buffer file_index;

        /// Index entries to append to the current frame
        struct pgsd_index_buffer frame_index;

        /// Buffered index entries to append to the current frame
        struct pgsd_index_buffer buffer_index;

        /// Buffered write data
        struct pgsd_byte_buffer write_buffer;

        /// List of names stored in the file
        struct pgsd_name_buffer file_names;

        /// List of names added in the current frame
        struct pgsd_name_buffer frame_names;

        /// The index of the last frame in the file
        uint64_t cur_frame;

        /// Size of the file (in bytes)
        // int64_t file_size;
        long long int file_size;

        /// Flags passed to pgsd_open() when opening this handle
        enum pgsd_open_flag open_flags;

        /// Access the names in the namelist
        struct pgsd_name_id_map name_map;

        /// Number of index entries pending in the current frame.
        uint64_t pending_index_entries;

        /// Maximum write buffer size (bytes).
        uint64_t maximum_write_buffer_size;

        /// Number of index entries to buffer before flushing.
        uint64_t index_entries_to_buffer;
        };

    /** Specify a version.

        @param major major version
        @param minor minor version

        @return a packed version number aaaa.bbbb suitable for storing in a pgsd file version entry.
    */
    uint32_t pgsd_make_version(unsigned int major, unsigned int minor);

    /** Create a PGSD file.

        @param fname File name (UTF-8 encoded).
        @param application Generating application name (truncated to 63 chars).
        @param schema Schema name for data to be written in this PGSD file (truncated to 63 chars).
        @param schema_version Version of the scheme data to be written (make with
        pgsd_make_version()).

        @post Create an empty pgsd file in a file of the given name. Overwrite any existing file at
        that location.

        The generated pgsd file is not opened. Call pgsd_open() to open it for writing.

        @return
          - PGSD_SUCCESS (0) on success. Negative value on failure:
          - PGSD_ERROR_IO: IO error (check errno).
    */
    // int pgsd_create(const char* fname,
    //                const char* application,
    //                const char* schema,
    //                uint32_t schema_version);

    /** Create and open a PGSD file.

        @param handle Handle to open.
        @param fname File name (UTF-8 encoded).
        @param application Generating application name (truncated to 63 chars).
        @param schema Schema name for data to be written in this GSD file (truncated to 63 chars).
        @param schema_version Version of the scheme data to be written (make with
            pgsd_make_version()).
        @param flags Either PGSD_OPEN_READWRITE, or PGSD_OPEN_APPEND.
        @param exclusive_create Set to non-zero to force exclusive creation of the file.

        @post Create an empty pgsd file with the given name. Overwrite any existing file at that
        location.

        Open the generated pgsd file in *handle*.

        The file descriptor is closed if there when an error opening the file.

        @return
          - PGSD_SUCCESS (0) on success. Negative value on failure:
          - PGSD_ERROR_IO: IO error (check errno).
          - PGSD_ERROR_NOT_A_PGSD_FILE: Not a PGSD file.
          - PGSD_ERROR_INVALID_PGSD_FILE_VERSION: Invalid PGSD file version.
          - PGSD_ERROR_FILE_CORRUPT: Corrupt file.
          - PGSD_ERROR_MEMORY_ALLOCATION_FAILED: Unable to allocate memory.
    */
    int pgsd_create_and_open(struct pgsd_handle* handle,
                            const char* fname,
                            const char* application,
                            const char* schema,
                            uint32_t schema_version,
                            enum pgsd_open_flag flags,
                            int exclusive_create);

    /** Open a PGSD file.

        @param handle Handle to open.
        @param fname File name to open (UTF-8 encoded).
        @param flags Either PGSD_OPEN_READWRITE, PGSD_OPEN_READONLY, or PGSD_OPEN_APPEND.

        @pre The file name *fname* is a PGSD file.

        @post Open a PGSD file and populates the handle for use by API calls.

        The file descriptor is closed if there is an error opening the file.

        @return
          - PGSD_SUCCESS (0) on success. Negative value on failure:
          - PGSD_ERROR_IO: IO error (check errno).
          - PGSD_ERROR_NOT_A_PGSD_FILE: Not a PGSD file.
          - PGSD_ERROR_INVALID_PGSD_FILE_VERSION: Invalid PGSD file version.
          - PGSD_ERROR_FILE_CORRUPT: Corrupt file.
          - PGSD_ERROR_MEMORY_ALLOCATION_FAILED: Unable to allocate memory.
    */
    int pgsd_open(struct pgsd_handle* handle, const char* fname, enum pgsd_open_flag flags);

    /** Truncate a PGSD file.

        @param handle Open PGSD file to truncate.

        After truncating, a file will have no frames and no data chunks. The file size will be that
        of a newly created pgsd file. The application, schema, and schema version metadata will be
        kept. Truncate does not close and reopen the file, so it is suitable for writing restart
        files on Lustre file systems without any metadata access.

        @return
          - PGSD_SUCCESS (0) on success. Negative value on failure:
          - PGSD_ERROR_IO: IO error (check errno).
          - PGSD_ERROR_NOT_A_PGSD_FILE: Not a PGSD file.
          - PGSD_ERROR_INVALID_PGSD_FILE_VERSION: Invalid PGSD file version.
          - PGSD_ERROR_FILE_CORRUPT: Corrupt file.
          - PGSD_ERROR_MEMORY_ALLOCATION_FAILED: Unable to allocate memory.
    */
    int pgsd_truncate(struct pgsd_handle* handle);

    /** Close a PGSD file.

        @param handle PGSD file to close.

        @pre *handle* was opened by pgsd_open().
        
        @post Writable files: All data and index entries buffered before the previous call to
              pgsd_end_frame() is written to the file (see pgsd_flush()).
        @post The file is closed.
        @post *handle* is freed and can no longer be used.

        @warning Ensure that all pgsd_write_chunk() calls are completed with pgsd_end_frame() before
        closing the file.

        @return
          - PGSD_SUCCESS (0) on success. Negative value on failure:
          - PGSD_ERROR_IO: IO error (check errno).
          - PGSD_ERROR_INVALID_ARGUMENT: *handle* is NULL.
    */
    int pgsd_close(struct pgsd_handle* handle, bool all);

    /** Complete the current frame.

        @param handle Handle to an open PGSD file

        @pre *handle* was opened by pgsd_open().

        @post The current frame counter is increased by 1.
        @post Flush the write buffer if it has overflowed. See pgsd_flush().
        
        @return
          - PGSD_SUCCESS (0) on success. Negative value on failure:
          - PGSD_ERROR_IO: IO error (check errno).
          - PGSD_ERROR_INVALID_ARGUMENT: *handle* is NULL.
          - PGSD_ERROR_FILE_MUST_BE_WRITABLE: The file was opened read-only.
          - PGSD_ERROR_MEMORY_ALLOCATION_FAILED: Unable to allocate memory.
    */
    int pgsd_end_frame(struct pgsd_handle* handle, bool all);

    /** Flush the write buffer.

        @param handle Handle to an open GSD file

        @pre *handle* was opened by pgsd_open().

        @post All data buffered by pgsd_write_chunk() are present in the file.
        @post All index entries buffered by pgsd_write_chunk() prior to the last call to
              pgsd_end_frame() are present in the file.

        @return
          - PGSD_SUCCESS (0) on success. Negative value on failure:
          - PGSD_ERROR_IO: IO error (check errno).
          - PGSD_ERROR_INVALID_ARGUMENT: *handle* is NULL.
          - PGSD_ERROR_FILE_MUST_BE_WRITABLE: The file was opened read-only.
          - PGSD_ERROR_MEMORY_ALLOCATION_FAILED: Unable to allocate memory.
    */
    int pgsd_flush(struct pgsd_handle* handle, bool all);

    /** Add a data chunk to the current frame.

        @param handle Handle to an open PGSD file.
        @param name Name of the data chunk.
        @param type type ID that identifies the type of data in *data*.
        @param N Number of rows in the data.
        @param M Number of columns in the data.
        @param flags set to 0, non-zero values reserved for future use.
        @param data Data buffer.

        @pre *handle* was opened by pgsd_open().
        @pre *name* is a unique name for data chunks in the given frame.
        @pre data is allocated and contains at least `N * M * pgsd_sizeof_type(type)` bytes.

        @post When there is space in the buffer: The given data is present in the write buffer.
              Otherwise, the data is present at the end of the file.
        @post The index is present in the buffer.

        @note If the PGSD file is version 1.0, the chunk name is truncated to 63 bytes. PGSD version
        2.0 files support arbitrarily long names.

        @note *N* == 0 is allowed. When *N* is 0, *data* may be NULL.

        @return
          - PGSD_SUCCESS (0) on success. Negative value on failure:
          - PGSD_ERROR_IO: IO error (check errno).
          - PGSD_ERROR_INVALID_ARGUMENT: *handle* is NULL, *N* == 0, *M* == 0, *type* is invalid, or
            *flags* != 0.
          - PGSD_ERROR_FILE_MUST_BE_WRITABLE: The file was opened read-only.
          - PGSD_ERROR_NAMELIST_FULL: The file cannot store any additional unique chunk names.
          - PGSD_ERROR_MEMORY_ALLOCATION_FAILED: failed to allocate memory.
    */
    int pgsd_write_chunk(struct pgsd_handle* handle,
                        const char* name,
                        enum pgsd_type type,
                        uint64_t N,
                        uint32_t M,

                        uint64_t N_global, // added DK
                        uint32_t M_global,
                        uint64_t offset,
                        uint64_t global_size,
                        bool all,
                        
                        uint8_t flags,
                        const void* data);

    /** Find a chunk in the PGSD file.

        @param handle Handle to an open PGSD file
        @param frame Frame to look for chunk
        @param name Name of the chunk to find

        @pre *handle* was opened by pgsd_open() in read or readwrite mode.

        The found entry contains size and type metadata and can be passed to pgsd_read_chunk() to
        read the data.

        @return A pointer to the found chunk, or NULL if not found.

        @note pgsd_find_chunk() calls gsd_flush() when the file is writable.
    */
    const struct pgsd_index_entry*
    pgsd_find_chunk(struct pgsd_handle* handle, uint64_t frame, const char* name, bool all);

    /** Read a chunk from the PGSD file.

        @param handle Handle to an open PGSD file.
        @param data Data buffer to read into.
        @param chunk Chunk to read.

        @pre *handle* was opened in read or readwrite mode.
        @pre *chunk* was found by pgsd_find_chunk().
        @pre *data* points to an allocated buffer with at least `N * M * pgsd_sizeof_type(type)`
       bytes.

        @return
          - PGSD_SUCCESS (0) on success. Negative value on failure:
          - PGSD_ERROR_IO: IO error (check errno).
          - PGSD_ERROR_INVALID_ARGUMENT: *handle* is NULL, *data* is NULL, or *chunk* is NULL.
          - PGSD_ERROR_FILE_MUST_BE_READABLE: The file was opened in append mode.
          - PGSD_ERROR_FILE_CORRUPT: The PGSD file is corrupt.
    @note pgsd_read_chunk() calls pgsd_flush() when the file is writable.
    */
    
    int pgsd_read_chunk(struct pgsd_handle* handle, void* data, const struct pgsd_index_entry* chunk, unsigned int *offset, bool all);

    /** Get the number of frames in the PGSD file.

        @param handle Handle to an open PGSD file

        @pre *handle* was opened by pgsd_open().

        @return The number of frames in the file, or 0 on error.
    */
    uint64_t pgsd_get_nframes(struct pgsd_handle* handle);

    /** Query size of a PGSD type ID.

        @param type Type ID to query.

        @return Size of the given type in bytes, or 0 for an unknown type ID.
    */
    size_t pgsd_sizeof_type(enum pgsd_type type);

    /** Search for chunk names in a pgsd file.

        @param handle Handle to an open PGSD file.
        @param match String to match.
        @param prev Search starting point.

        @pre *handle* was opened by pgsd_open()
        @pre *prev* was returned by a previous call to pgsd_find_matching_chunk_name()

        To find the first matching chunk name, pass NULL for prev. Pass in the previous found string
        to find the next after that, and so on. Chunk names match if they begin with the string in
        *match*. Chunk names returned by this function may be present in at least one frame.

        @return Pointer to a string, NULL if no more matching chunks are found found, or NULL if
        *prev* is invalid
    
    @note  pgsd_find_matching_chunk_name() calls pgsd_flush() when the file is writable.

    */
    const char*
    pgsd_find_matching_chunk_name(struct pgsd_handle* handle, const char* match, const char* prev, bool all);

    /** Upgrade a PGSD file to the latest specification.

        @param handle Handle to an open PGSD file

        @pre *handle* was opened by pgsd_open() with a writable mode.
        @pre There are no pending data to write to the file in pgsd_end_frame()

        @return
          - PGSD_SUCCESS (0) on success. Negative value on failure:
          - PGSD_ERROR_IO: IO error (check errno).
          - PGSD_ERROR_INVALID_ARGUMENT: *handle* is NULL
          - PGSD_ERROR_FILE_MUST_BE_WRITABLE: The file was opened in read-only mode.
    */
    // int pgsd_upgrade(struct pgsd_handle* handle);


    /** Get the maximum write buffer size.

        @param handle Handle to an open GSD file

        @pre *handle* was opened by pgsd_open().

        @return The maximum write buffer size in bytes, or 0 on error.
    */
    uint64_t pgsd_get_maximum_write_buffer_size(struct pgsd_handle* handle);

    /** Set the maximum write buffer size.

        @param handle Handle to an open GSD file
        @param size Maximum number of bytes to allocate in the write buffer (must be greater than
        0).

        @pre *handle* was opened by pgsd_open().

        @return
          - PGSD_SUCCESS (0) on success. Negative value on failure:
          - PGSD_ERROR_INVALID_ARGUMENT: *handle* is NULL
          - PGSD_ERROR_INVALID_ARGUMENT: size == 0
    */
    int pgsd_set_maximum_write_buffer_size(struct pgsd_handle* handle, uint64_t size);

    /** Get the number of index entries to buffer.

        @param handle Handle to an open GSD file

        @pre *handle* was opened by pgsd_open().

        @return The number of index entries to buffer, or 0 on error.
    */
    uint64_t pgsd_get_index_entries_to_buffer(struct pgsd_handle* handle);

    /** Set the number of index entries to buffer.

        @param handle Handle to an open GSD file
        @param number Number of index entries to buffer before automatically flushing in
        `pgsd_end_frame()` (must be greater than 0).

        @pre *handle* was opened by pgsd_open().

        @note GSD may allocate more than this number of entries in the buffer, as needed to store
        all index entries for the already buffered frames and the current frame.

        @return
          - PGSD_SUCCESS (0) on success. Negative value on failure:
          - PGSD_ERROR_INVALID_ARGUMENT: *handle* is NULL
          - PGSD_ERROR_INVALID_ARGUMENT: number == 0
    */
    int pgsd_set_index_entries_to_buffer(struct pgsd_handle* handle, uint64_t number);

#ifdef __cplusplus
    }
#endif

// #endif // #ifdef ENABLE_MPI

#endif // #ifndef PGSD_H
