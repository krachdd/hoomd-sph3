// Copyright (c) 2016-2023 The Regents of the University of Michigan
// Part of PGSD, released under the BSD 2-Clause License.

// #ifdef ENABLE_MPI

#include <sys/stat.h>
#ifdef _WIN32

#pragma warning(push)
#pragma warning(disable : 4996)

#define PGSD_USE_MMAP 0
#define WIN32_LEAN_AND_MEAN
#include <io.h>
#include <windows.h>

#else // linux / mac

#define _XOPEN_SOURCE 500
#include <sys/mman.h>
#include <unistd.h>
// #define PGSD_USE_MMAP 1
#define PGSD_USE_MMAP 0

// for sys/mman.h and mmap see:
// https://pubs.opengroup.org/onlinepubs/000095399/functions/mmap.html

#endif

#ifdef __APPLE__
#include <limits.h>
#endif

#include <errno.h>
#include <fcntl.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

// Only testing
#include <inttypes.h>

#include "pgsd.h"

/** @file pgsd.c
    @brief Implements the PGSD C API
*/

/// Magic value identifying a PGSD file
const uint64_t PGSD_MAGIC_ID = 0x65DF65DF65DF65DF;

/// Initial index size
enum
    {
    PGSD_INITIAL_INDEX_SIZE = 128
    };

/// Initial namelist size
enum
    {
    PGSD_INITIAL_NAME_BUFFER_SIZE = 1024
    };

/// Size of initial frame index
enum
    {
    PGSD_INITIAL_FRAME_INDEX_SIZE = 16
    };

/// Initial size of write buffer
enum
    {
    PGSD_INITIAL_WRITE_BUFFER_SIZE = 1024
    };

/// Default maximum size of write buffer
enum
    {
    PGSD_DEFAULT_MAXIMUM_WRITE_BUFFER_SIZE = 64 * 1024 * 1024
    };

/// Default number of index entries to buffer
enum
    {
    PGSD_DEFAULT_INDEX_ENTRIES_TO_BUFFER = 256 * 1024
    };

/// Size of hash map
enum
    {
    PGSD_NAME_MAP_SIZE = 57557
    };

/// Current PGSD file specification
enum
    {
    PGSD_CURRENT_FILE_VERSION = 2
    };

// define windows wrapper functions
// #ifdef _WIN32
// #define lseek _lseeki64
// #define ftruncate _chsize
// #define fsync _commit
// typedef int64_t ssize_t;

// int S_IRUSR = _S_IREAD;
// int S_IWUSR = _S_IWRITE;
// int S_IRGRP = _S_IREAD;
// int S_IWGRP = _S_IWRITE;

// inline ssize_t pread(int fd, void* buf, size_t count, int64_t offset)
//     {
//     // Note: _read only accepts unsigned int values
//     if (count > UINT_MAX)
//         return PGSD_ERROR_IO;

//     int64_t oldpos = _telli64(fd);
//     _lseeki64(fd, offset, SEEK_SET);
//     ssize_t result = _read(fd, buf, (unsigned int)count);
//     _lseeki64(fd, oldpos, SEEK_SET);
//     return result;
//     }

// inline ssize_t pwrite(int fd, const void* buf, size_t count, int64_t offset)
//     {
//     // Note: _write only accepts unsigned int values
//     if (count > UINT_MAX)
//         return PGSD_ERROR_IO;

//     int64_t oldpos = _telli64(fd);
//     _lseeki64(fd, offset, SEEK_SET);
//     ssize_t result = _write(fd, buf, (unsigned int)count);
//     _lseeki64(fd, oldpos, SEEK_SET);
//     return result;
//     }

// #endif

/** Zero memory

    @param d pointer to memory region
    @param size_to_zero size of the area to zero in bytes
*/
inline static void pgsd_util_zero_memory(void* d, size_t size_to_zero)
    {
    memset(d, 0, size_to_zero);
    }

/** @internal
    @brief Write large data buffer to file

    The system call pwrite() fails to write very large data buffers. This method calls pwrite() as
    many times as necessary to completely write a large buffer.

    @param fd File descriptor.
    @param buf Data buffer.
    @param count Number of bytes to write.
    @param offset Location in the file to start writing.

    @returns The total number of bytes written or a negative value on error.
*/
// inline static ssize_t pgsd_io_pwrite_retry(MPI_File *f, const void* buf, size_t count, int64_t offset)
//     {
//     size_t total_bytes_written = 0;
//     const char* ptr = (char*)buf;

//     // perform multiple pwrite calls to complete a large write successfully
//     while (total_bytes_written < count)
//         {
//         size_t to_write = count - total_bytes_written;
// #if defined(_WIN32) || defined(__APPLE__)
//         // win32 and apple raise an error for writes greater than INT_MAX
//         if (to_write > INT_MAX / 2)
//             to_write = INT_MAX / 2;
// #endif

//         errno = 0;

//         ssize_t bytes_written
//             = pwrite(fd, ptr + total_bytes_written, to_write, offset + total_bytes_written);
//         if (bytes_written == -1 || (bytes_written == 0 && errno != 0))
//             {
//             return PGSD_ERROR_IO;
//             }

//         total_bytes_written += bytes_written;
//         }

//     return total_bytes_written;
//     }

/** @internal
    @brief Read large data buffer to file

    The system call pread() fails to read very large data buffers. This method calls pread() as many
    times as necessary to completely read a large buffer.

    @param fd File descriptor.
    @param buf Data buffer.
    @param count Number of bytes to read.
    @param offset Location in the file to start reading.

    @returns The total number of bytes read or a negative value on error.
*/
// inline static ssize_t pgsd_io_pread_retry(int fd, void* buf, size_t count, int64_t offset)
//     {
//     size_t total_bytes_read = 0;
//     char* ptr = (char*)buf;

//     // perform multiple pread calls to complete a large write successfully
//     while (total_bytes_read < count)
//         {
//         size_t to_read = count - total_bytes_read;
// #if defined(_WIN32) || defined(__APPLE__)
//         // win32 and apple raise errors for reads greater than INT_MAX
//         if (to_read > INT_MAX / 2)
//             to_read = INT_MAX / 2;
// #endif

//         errno = 0;
//         ssize_t bytes_read = pread(fd, ptr + total_bytes_read, to_read, offset + total_bytes_read);
//         if (bytes_read == -1 || (bytes_read == 0 && errno != 0))
//             {
//             return PGSD_ERROR_IO;
//             }

//         total_bytes_read += bytes_read;

//         // handle end of file
//         if (bytes_read == 0)
//             {
//             return total_bytes_read;
//             }
//         }

//     return total_bytes_read;
//     }

/** @internal
    @brief Allocate a name/id map

    @param map Map to allocate.
    @param size Number of entries in the map.

    @returns PGSD_SUCCESS on success, PGSD_* error codes on error.
*/
inline static int pgsd_name_id_map_allocate(struct pgsd_name_id_map* map, size_t size)
    {
    if (map == NULL || map->v || size == 0 || map->size != 0)
        {
        return PGSD_ERROR_INVALID_ARGUMENT;
        }

    map->v = calloc(size, sizeof(struct pgsd_name_id_pair));
    if (map->v == NULL)
        {
        return PGSD_ERROR_MEMORY_ALLOCATION_FAILED;
        }

    map->size = size;

    return PGSD_SUCCESS;
    }

/** @internal
    @brief Free a name/id map

    @param map Map to free.

    @returns PGSD_SUCCESS on success, PGSD_* error codes on error.
*/
inline static int pgsd_name_id_map_free(struct pgsd_name_id_map* map)
    {
    if (map == NULL || map->v == NULL || map->size == 0)
        {
        return PGSD_ERROR_INVALID_ARGUMENT;
        }

    // free all of the linked lists
    size_t i;
    for (i = 0; i < map->size; i++)
        {
        free(map->v[i].name);

        struct pgsd_name_id_pair* cur = map->v[i].next;
        while (cur != NULL)
            {
            struct pgsd_name_id_pair* prev = cur;
            cur = cur->next;
            free(prev->name);
            free(prev);
            }
        }

    // free the main map
    free(map->v);

    map->v = 0;
    map->size = 0;

    return PGSD_SUCCESS;
    }

/** @internal
    @brief Hash a string

    @param str String to hash

    @returns Hashed value of the string.
*/
inline static unsigned long pgsd_hash_str(const unsigned char* str)
    {
    unsigned long hash = 5381; // NOLINT
    int c;

    while ((c = *str++))
        {
        hash = ((hash << 5) + hash) + c; /* hash * 33 + c NOLINT */
        }

    return hash;
    }

/** @internal
    @brief Insert a string into a name/id map

    @param map Map to insert into.
    @param str String to insert.
    @param id ID to associate with the string.

    @returns PGSD_SUCCESS on success, PGSD_* error codes on error.
*/
inline static int pgsd_name_id_map_insert(struct pgsd_name_id_map* map, const char* str, uint16_t id)
    {
    if (map == NULL || map->v == NULL || map->size == 0)
        {
        return PGSD_ERROR_INVALID_ARGUMENT;
        }

    size_t hash = pgsd_hash_str((const unsigned char*)str) % map->size;

    // base case: no conflict
    if (map->v[hash].name == NULL)
        {
        map->v[hash].name = calloc(strlen(str) + 1, sizeof(char));
        if (map->v[hash].name == NULL)
            {
            return PGSD_ERROR_MEMORY_ALLOCATION_FAILED;
            }
        memcpy(map->v[hash].name, str, strlen(str) + 1);
        map->v[hash].id = id;
        map->v[hash].next = NULL;
        }
    else
        {
        // go to the end of the conflict list
        struct pgsd_name_id_pair* insert_point = map->v + hash;

        while (insert_point->next != NULL)
            {
            insert_point = insert_point->next;
            }

        // allocate and insert a new entry
        insert_point->next = malloc(sizeof(struct pgsd_name_id_pair));
        if (insert_point->next == NULL)
            {
            return PGSD_ERROR_MEMORY_ALLOCATION_FAILED;
            }

        insert_point->next->name = calloc(strlen(str) + 1, sizeof(char));
        if (insert_point->next->name == NULL)
            {
            return PGSD_ERROR_MEMORY_ALLOCATION_FAILED;
            }
        memcpy(insert_point->next->name, str, strlen(str) + 1);
        insert_point->next->id = id;
        insert_point->next->next = NULL;
        }

    return PGSD_SUCCESS;
    }

/** @internal
    @brief Find an ID in a name/id mapping

    @param map Map to search.
    @param str String to search.

    @returns The ID if found, or UINT16_MAX if not found.
*/
inline static uint16_t pgsd_name_id_map_find(struct pgsd_name_id_map* map, const char* str)
    {
    if (map == NULL || map->v == NULL || map->size == 0)
        {
        return UINT16_MAX;
        }

    size_t hash = pgsd_hash_str((const unsigned char*)str) % map->size;

    struct pgsd_name_id_pair* cur = map->v + hash;

    while (cur != NULL)
        {
        if (cur->name == NULL)
            {
            // not found
            return UINT16_MAX;
            }

        if (strcmp(str, cur->name) == 0)
            {
            // found
            return cur->id;
            }

        // keep looking
        cur = cur->next;
        }

    // not found in any conflict
    return UINT16_MAX;
    }

/** @internal
    @brief Utility function to validate index entry
    @param handle handle to the open pgsd file
    @param idx index of entry to validate

    @returns 1 if the entry is valid, 0 if it is not
*/
inline static int pgsd_is_entry_valid(struct pgsd_handle* handle, size_t idx)
    {
    const struct pgsd_index_entry entry = handle->file_index.data[idx];

    // check for valid type
    if (pgsd_sizeof_type((enum pgsd_type)entry.type) == 0)
        {
        return 0;
        }

    // validate that we don't read past the end of the file
    size_t size = entry.N * entry.M * pgsd_sizeof_type((enum pgsd_type)entry.type);
    if ((entry.location + size) > (uint64_t)handle->file_size)
        {
        return 0;
        }

    // check for valid frame (frame cannot be more than the number of index entries)
    if (entry.frame >= handle->header.index_allocated_entries)
        {
        return 0;
        }

    // check for valid id
    if (entry.id >= (handle->file_names.n_names + handle->frame_names.n_names))
        {
        return 0;
        }

    // check for valid flags
    if (entry.flags != 0)
        {
        return 0;
        }

    return 1;
    }

/** @internal
    @brief Allocate a write buffer

    @param buf Buffer to allocate.
    @param reserve Number of bytes to allocate.

    @returns PGSD_SUCCESS on success, PGSD_* error codes on error.
*/
inline static int pgsd_byte_buffer_allocate(struct pgsd_byte_buffer* buf, size_t reserve)
    {
    if (buf == NULL || buf->data || reserve == 0 || buf->reserved != 0 || buf->size != 0)
        {
        return PGSD_ERROR_INVALID_ARGUMENT;
        }

    buf->data = calloc(reserve, sizeof(char));
    if (buf->data == NULL)
        {
        return PGSD_ERROR_MEMORY_ALLOCATION_FAILED;
        }

    buf->size = 0;
    buf->reserved = reserve;

    return PGSD_SUCCESS;
    }

/** @internal
    @brief Append bytes to a byte buffer

    @param buf Buffer to append to.
    @param data Data to append.
    @param size Number of bytes in *data*.

    @returns PGSD_SUCCESS on success, PGSD_* error codes on error.

    DK : per rank
*/
inline static int pgsd_byte_buffer_append(struct pgsd_byte_buffer* buf, const char* data, size_t size)
    {
    if (buf == NULL || buf->data == NULL || size == 0 || buf->reserved == 0)
        {
        return PGSD_ERROR_INVALID_ARGUMENT;
        }

    if (buf->size + size > buf->reserved)
        {
        // reallocate by doubling
        size_t new_reserved = buf->reserved * 2;
        while (buf->size + size >= new_reserved)
            {
            new_reserved = new_reserved * 2;
            }

        char* old_data = buf->data;
        buf->data = realloc(buf->data, sizeof(char) * new_reserved);
        if (buf->data == NULL)
            {
            // this free should not be necessary, but clang-tidy disagrees
            free(old_data);
            return PGSD_ERROR_MEMORY_ALLOCATION_FAILED;
            }

        // zero the new memory, but only the portion after the end of the new section to be appended
        pgsd_util_zero_memory(buf->data + (buf->size + size),
                             sizeof(char) * (new_reserved - (buf->size + size)));
        buf->reserved = new_reserved;
        }

    memcpy(buf->data + buf->size, data, size);
    buf->size += size;

    return PGSD_SUCCESS;
    }

/** @internal
    @brief Free the memory allocated by the write buffer or unmap the mapped memory.

    @param buf Buffer to free.

    @returns PGSD_SUCCESS on success, PGSD_* error codes on error.

    DK : per rank

*/
inline static int pgsd_byte_buffer_free(struct pgsd_byte_buffer* buf)
    {
    if (buf == NULL || buf->data == NULL)
        {
        return PGSD_ERROR_INVALID_ARGUMENT;
        }

    free(buf->data);

    pgsd_util_zero_memory(buf, sizeof(struct pgsd_byte_buffer));
    return PGSD_SUCCESS;
    }

/** @internal
    @brief Allocate a buffer of index entries

    @param buf Buffer to allocate.
    @param reserve Number of entries to allocate.

    @post The buffer's data element has *reserve* elements allocated in memory.

    @returns PGSD_SUCCESS on success, PGSD_* error codes on error.

    DK : per rank
*/
inline static int pgsd_index_buffer_allocate(struct pgsd_index_buffer* buf, size_t reserve)
    {
    if (buf == NULL || buf->mapped_data || buf->data || reserve == 0 || buf->reserved != 0
        || buf->size != 0)
        {
        printf("Failed: invalid argument");
        return PGSD_ERROR_INVALID_ARGUMENT;
        }

    buf->data = calloc(reserve, sizeof(struct pgsd_index_entry));
    if (buf->data == NULL)
        {
        printf("Failed: Mem alloc error");
        return PGSD_ERROR_MEMORY_ALLOCATION_FAILED;
        }

    buf->size = 0;
    buf->reserved = reserve;
    buf->mapped_data = NULL;
    buf->mapped_len = 0;

    return PGSD_SUCCESS;
    }

/** @internal
    @brief Map index entries from the file

    @param buf Buffer to map.
    @param handle PGSD file handle to map.

    @post The buffer's data element contains the index data from the file.

    On some systems, this will use mmap to efficiently access the file. On others, it may result in
    an allocation and read of the entire index from the file.

    @returns PGSD_SUCCESS on success, PGSD_* error codes on error.
*/
inline static int pgsd_index_buffer_map(struct pgsd_index_buffer* buf, struct pgsd_handle* handle)
    {
    if (buf == NULL || buf->mapped_data || buf->data || buf->reserved != 0 || buf->size != 0)
        {
        printf("index buffer map: invalid arg\n");
        return PGSD_ERROR_INVALID_ARGUMENT;
        }

    // validate that the index block exists inside the file
    if (handle->header.index_location
            + sizeof(struct pgsd_index_entry) * handle->header.index_allocated_entries
        > (uint64_t)handle->file_size)
        {
        printf("index File corrupt\n");
        return PGSD_ERROR_FILE_CORRUPT;
        }
    printf("index buffer map: buffer correct\n");

#if PGSD_USE_MMAP
    printf("index buffer map: use mmap\n");
    // map the index in read only mode
    size_t page_size = getpagesize();
    printf("index buffer map: pagesizes %ld\n", page_size);
    size_t index_size = sizeof(struct pgsd_index_entry) * handle->header.index_allocated_entries;
    printf("index buffer map: sindexizes %ld\n", index_size);
    size_t offset = (handle->header.index_location / page_size) * page_size;
    printf("index buffer map: offset %ld\n", offset);
    printf("index buffer map: index_size + (handle->header.index_location - offset) %ld\n", index_size + (handle->header.index_location - offset));
    printf("index buffer map: handle->header.index_location %ld\n", handle->header.index_location);
    buf->mapped_data = mmap(NULL,
                            index_size + (handle->header.index_location - offset),
                            PROT_READ,
                            MAP_SHARED,
                            handle->fh,
                            offset);

    if (buf->mapped_data == MAP_FAILED)
        {
        return PGSD_ERROR_IO;
        }

    printf("index buffer map: mmap success\n");

    buf->data = (struct pgsd_index_entry*)(((char*)buf->mapped_data)
                                          + (handle->header.index_location - offset));

    buf->mapped_len = index_size + (handle->header.index_location - offset);
    buf->reserved = handle->header.index_allocated_entries;
#else
    // mmap not supported, read the data from the disk
    int retval = pgsd_index_buffer_allocate(buf, handle->header.index_allocated_entries);
    if (retval != PGSD_SUCCESS)
        {
        printf("index allocate failed\n");
        return retval;
        }
    printf("index buffer allocated");
    // size_t bytes_read = sizeof(struct pgsd_index_entry)* handle->header.index_allocated_entries;
    // MPI_File_get_size(handle->fh, &(handle->header.index_location));
    // MPI_File_seek(handle->fh, 0, MPI_SEEK_END);
    // MPI_Barrier(MPI_COMM_WORLD);
    
    MPI_File_read_at(handle->fh, handle->header.index_location, buf->data, sizeof(struct pgsd_index_entry)* handle->header.index_allocated_entries, MPI_BYTE, MPI_STATUS_IGNORE);

    // ssize_t bytes_read = pgsd_io_pread_retry(handle->fd,
    //                                         buf->data,
    //                                         sizeof(struct pgsd_index_entry)
    //                                             * handle->header.index_allocated_entries,
    //                                         handle->header.index_location);

    // if (bytes_read == -1
    //     || bytes_read != sizeof(struct pgsd_index_entry) * handle->header.index_allocated_entries)
    //     {
    //     return PGSD_ERROR_IO;
    //     }
#endif

    // determine the number of index entries in the list
    // file is corrupt if first index entry is invalid
    if (buf->data[0].location != 0 && !pgsd_is_entry_valid(handle, 0))
        {
        printf("index File corrupt after endif \n");
        return PGSD_ERROR_FILE_CORRUPT;
        }

    if (buf->data[0].location == 0)
        {
        buf->size = 0;
        }
    else
        {
        // determine the number of index entries (marked by location = 0)
        // binary search for the first index entry with location 0
        size_t L = 0;
        size_t R = buf->reserved;

        // progressively narrow the search window by halves
        do
            {
            size_t m = (L + R) / 2;
            
            // file is corrupt if any index entry is invalid or frame does not increase
            // monotonically
            if (buf->data[m].location != 0
                && (!pgsd_is_entry_valid(handle, m) || buf->data[m].frame < buf->data[L].frame))
                {
                printf("index File corrupt in loop\n");
                return PGSD_ERROR_FILE_CORRUPT;
                }

            if (buf->data[m].location != 0)
                {
                L = m;
                }
            else
                {
                R = m;
                }
            } while ((R - L) > 1);

        // this finds R = the first index entry with location = 0
        buf->size = R;
        }

    return PGSD_SUCCESS;
    }

/** @internal
    @brief Free the memory allocated by the index buffer or unmap the mapped memory.

    @param buf Buffer to free.

    @returns PGSD_SUCCESS on success, PGSD_* error codes on error.

    DK : per rank
*/
inline static int pgsd_index_buffer_free(struct pgsd_index_buffer* buf)
    {
    if (buf == NULL || buf->data == NULL)
        {
        return PGSD_ERROR_INVALID_ARGUMENT;
        }

#if PGSD_USE_MMAP
    if (buf->mapped_data)
        {
        int retval = munmap(buf->mapped_data, buf->mapped_len);

        if (retval != 0)
            {
            return PGSD_ERROR_IO;
            }
        }
    else
#endif
        {
        free(buf->data);
        }

    pgsd_util_zero_memory(buf, sizeof(struct pgsd_index_buffer));
    return PGSD_SUCCESS;
    }

/** @internal
    @brief Add a new index entry and provide a pointer to it.

    @param buf Buffer to add too.
    @param entry [out] Pointer to set to the new entry.

    Double the size of the reserved space as needed to hold the new entry. Does not accept mapped
    indices.

    @returns PGSD_SUCCESS on success, PGSD_* error codes on error.

    DK : per rank
*/
inline static int pgsd_index_buffer_add(struct pgsd_index_buffer* buf, struct pgsd_index_entry** entry)
    {
    if (buf == NULL || buf->mapped_data || entry == NULL || buf->reserved == 0)
        {
        return PGSD_ERROR_INVALID_ARGUMENT;
        }

    if (buf->size == buf->reserved)
        {
        // grow the array
        size_t new_reserved = buf->reserved * 2;
        buf->data = realloc(buf->data, sizeof(struct pgsd_index_entry) * new_reserved);
        if (buf->data == NULL)
            {
            return PGSD_ERROR_MEMORY_ALLOCATION_FAILED;
            }

        // zero the new memory
        pgsd_util_zero_memory(buf->data + buf->reserved,
                             sizeof(struct pgsd_index_entry) * (new_reserved - buf->reserved));
        buf->reserved = new_reserved;
        }

    size_t insert_pos = buf->size;
    buf->size++;
    *entry = buf->data + insert_pos;

    return PGSD_SUCCESS;
    }

inline static int pgsd_cmp_index_entry(const struct pgsd_index_entry* a,
                                      const struct pgsd_index_entry* b)
    {
    int result = 0;

    if (a->frame < b->frame)
        {
        result = -1;
        }

    if (a->frame > b->frame)
        {
        result = 1;
        }

    if (a->frame == b->frame)
        {
        if (a->id < b->id)
            {
            result = -1;
            }

        if (a->id > b->id)
            {
            result = 1;
            }

        if (a->id == b->id)
            {
            result = 0;
            }
        }

    return result;
    }

/** @internal
    @brief Compute heap parent node.
    @param i Node index.
*/
inline static size_t pgsd_heap_parent(size_t i)
    {
    return (i - 1) / 2;
    }

/** @internal
    @brief Compute heap left child.
    @param i Node index.
*/
inline static size_t pgsd_heap_left_child(size_t i)
    {
    return 2 * i + 1;
    }

/** @internal
    @brief Swap the nodes *a* and *b* in the buffer
    @param buf Buffer.
    @param a First index to swap.
    @param b Second index to swap.
*/
inline static void pgsd_heap_swap(struct pgsd_index_buffer* buf, size_t a, size_t b)
    {
    struct pgsd_index_entry tmp = buf->data[a];
    buf->data[a] = buf->data[b];
    buf->data[b] = tmp;
    }

/** @internal
    @brief Shift heap node downward
    @param buf Buffer.
    @param start First index of the valid heap in *buf*.
    @param end Last index of the valid hep in *buf*.
*/
inline static void pgsd_heap_shift_down(struct pgsd_index_buffer* buf, size_t start, size_t end)
    {
    size_t root = start;

    while (pgsd_heap_left_child(root) <= end)
        {
        size_t child = pgsd_heap_left_child(root);
        size_t swap = root;

        if (pgsd_cmp_index_entry(buf->data + swap, buf->data + child) < 0)
            {
            swap = child;
            }
        if (child + 1 <= end && pgsd_cmp_index_entry(buf->data + swap, buf->data + child + 1) < 0)
            {
            swap = child + 1;
            }

        if (swap == root)
            {
            return;
            }

        pgsd_heap_swap(buf, root, swap);
        root = swap;
        }
    }

/** @internal
    @brief Convert unordered index buffer to a heap
    @param buf Buffer.
*/
inline static void pgsd_heapify(struct pgsd_index_buffer* buf)
    {
    ssize_t start = pgsd_heap_parent(buf->size - 1);

    while (start >= 0)
        {
        pgsd_heap_shift_down(buf, start, buf->size - 1);
        start--;
        }
    }

/** @internal
    @brief Sort the index buffer.

    @param buf Buffer to sort.

    Sorts an in-memory index buffer. Does not accept mapped indices.

    @returns PGSD_SUCCESS on success, PGSD_* error codes on error.
*/
inline static int pgsd_index_buffer_sort(struct pgsd_index_buffer* buf)
    {
    if (buf == NULL || buf->mapped_data || buf->reserved == 0)
        {
        return PGSD_ERROR_INVALID_ARGUMENT;
        }

    // arrays of size 0 or 1 are already sorted
    if (buf->size <= 1)
        {
        return PGSD_SUCCESS;
        }

    pgsd_heapify(buf);

    size_t end = buf->size - 1;
    while (end > 0)
        {
        pgsd_heap_swap(buf, end, 0);
        end = end - 1;
        pgsd_heap_shift_down(buf, 0, end);
        }

    return PGSD_SUCCESS;
    }

/** @internal
    @brief Utility function to expand the memory space for the index block in the file.

    @param handle Handle to the open pgsd file.
    @param size_required The new index must be able to hold at least this many elements.

    @returns PGSD_SUCCESS on success, PGSD_* error codes on error.

    DK : only used in flush with root loop 
*/
inline static int pgsd_expand_file_index(struct pgsd_handle* handle, size_t size_required)
    {
    if (handle->open_flags == PGSD_OPEN_READONLY)
        {
        return PGSD_ERROR_FILE_MUST_BE_WRITABLE;
        }

    // multiply the index size each time it grows
    // this allows the index to grow rapidly to accommodate new frames
    const int multiplication_factor = 2;

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    bool root = false;
    if ( rank == 0 ){
        root = true;
    }
    bool all = true;

    // MPI_Barrier(MPI_COMM_WORLD);
    // save the old size and update the new size
    size_t size_old = handle->header.index_allocated_entries;
    size_t size_new = size_old * multiplication_factor;

    while (size_new <= size_required)
        {
        size_new *= multiplication_factor;
        }

    // Mac systems deadlock when writing from a mapped region into the tail end of that same region
    // unmap the index first and copy it over by chunks
    int retval = pgsd_index_buffer_free(&handle->file_index);
    if (retval != 0)
        {
        return retval;
        }

    // allocate the copy buffer
    uint64_t copy_buffer_size
        = PGSD_DEFAULT_INDEX_ENTRIES_TO_BUFFER * sizeof(struct pgsd_index_entry);
    if (copy_buffer_size > size_old * sizeof(struct pgsd_index_entry))
        {
        copy_buffer_size = size_old * sizeof(struct pgsd_index_entry);
        }
    char* buf = malloc(copy_buffer_size);


    // write the current index to the end of the file
    // int64_t new_index_location = lseek(handle->fd, 0, SEEK_END);
    // int64_t new_index_location;
    long long int new_index_location;
    int64_t old_index_location = handle->header.index_location;
    
    MPI_File_get_size(handle->fh, &new_index_location);
    MPI_File_seek(handle->fh, 0, MPI_SEEK_END);

    size_t total_bytes_written = 0;
    size_t old_index_bytes = size_old * sizeof(struct pgsd_index_entry);

    while (total_bytes_written < old_index_bytes)
        {
        size_t bytes_to_copy = copy_buffer_size;
        if (old_index_bytes - total_bytes_written < copy_buffer_size)
            {
            bytes_to_copy = old_index_bytes - total_bytes_written;
            }

        size_t bytes_read = bytes_to_copy;
        // ssize_t bytes_read = pgsd_io_pread_retry(handle->fd,
        //                                         buf,
        //                                         bytes_to_copy,
        //                                         old_index_location + total_bytes_written);
        MPI_File_read_at(handle->fh, old_index_location + total_bytes_written, buf, bytes_to_copy, MPI_BYTE, MPI_STATUS_IGNORE);


        if (bytes_read == -1 || bytes_read != bytes_to_copy)
            {
            free(buf);
            return PGSD_ERROR_IO;
            }

        size_t bytes_written = bytes_to_copy;
        // ssize_t bytes_written = pgsd_io_pwrite_retry(handle->fd,
        //                                             buf,
        //                                             bytes_to_copy,
        //                                             new_index_location + total_bytes_written);
        if( all == true  || root == true ){
            MPI_File_write_at(handle->fh, new_index_location + total_bytes_written, buf, bytes_to_copy, MPI_BYTE, MPI_STATUS_IGNORE);
        }

        if (bytes_written == -1 || bytes_written != bytes_to_copy)
            {
            free(buf);
            return PGSD_ERROR_IO;
            }

        total_bytes_written += bytes_written;
        }

    // fill the new index space with 0s
    pgsd_util_zero_memory(buf, copy_buffer_size);

    size_t new_index_bytes = size_new * sizeof(struct pgsd_index_entry);
    while (total_bytes_written < new_index_bytes)
        {
        size_t bytes_to_copy = copy_buffer_size;

        if (new_index_bytes - total_bytes_written < copy_buffer_size)
            {
            bytes_to_copy = new_index_bytes - total_bytes_written;
            }

        size_t bytes_written = bytes_to_copy;
        if( all == true  || root == true ){
            MPI_File_write_at(handle->fh, new_index_location + total_bytes_written, buf, bytes_to_copy, MPI_BYTE, MPI_STATUS_IGNORE);
        }
        // ssize_t bytes_written = pgsd_io_pwrite_retry(handle->fd,
        //                                             buf,
        //                                             bytes_to_copy,
        //                                             new_index_location + total_bytes_written);

        if (bytes_written == -1 || bytes_written != bytes_to_copy)
            {
            free(buf);
            return PGSD_ERROR_IO;
            }

        total_bytes_written += bytes_written;
        }

    // sync the expanded index
    // retval = fsync(handle->fd);
    // if (retval != 0)
    //     {
    //     free(buf);
    //     return PGSD_ERROR_IO;
    //     }

    // free the copy buffer
    free(buf);

    // update the header
    handle->header.index_location = new_index_location;
    handle->file_size = handle->header.index_location + total_bytes_written;
    handle->header.index_allocated_entries = size_new;

    // write the new header out
    size_t bytes_written = sizeof(struct pgsd_header);
    if( all == true  || root == true ){
        MPI_File_write(handle->fh, &(handle->header), sizeof(struct pgsd_header), MPI_BYTE, MPI_STATUS_IGNORE);
    }
    // ssize_t bytes_written
    //     = pgsd_io_pwrite_retry(handle->fd, &(handle->header), sizeof(struct pgsd_header), 0);
    if (bytes_written != sizeof(struct pgsd_header))
        {
        return PGSD_ERROR_IO;
        }

    // sync the updated header
    // retval = fsync(handle->fd);
    // if (retval != 0)
    //     {
    //     return PGSD_ERROR_IO;
    //     }

    // remap the file index
    retval = pgsd_index_buffer_map(&handle->file_index, handle);
    if (retval != 0)
        {
        return retval;
        }

    return PGSD_SUCCESS;
    }

/** @internal
    @brief Flush the write buffer.

    pgsd_write_frame() writes small data chunks into the write buffer. It adds index entries for
    these chunks to pgsd_handle::buffer_index with locations offset from the start of the write
    buffer. pgsd_flush_write_buffer() writes the buffer to the end of the file, moves the index
    entries to pgsd_handle::frame_index and updates the location to reference the beginning of the
    file.

    @param handle Handle to flush the write buffer.
    @returns PGSD_SUCCESS on success or PGSD_* error codes on error

    DK : collectively. Use write buffer sizes as offset. 

*/
inline static int pgsd_flush_write_buffer(struct pgsd_handle* handle)
    {
    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    bool root = false;
    if ( rank == 0 ){
        root = true;
    }

    size_t allbuffers[nprocs];
    MPI_Allgather(&handle->write_buffer.size, 1, MPI_UNSIGNED_LONG, allbuffers, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    // printf("rank %i allbuffers %lu\n", rank, allbuffers[rank]);

    
    // handle->write_buffer.size

    if (handle == NULL)
        {
        return PGSD_ERROR_INVALID_ARGUMENT;
        }

    if (handle->write_buffer.size == 0 && handle->buffer_index.size == 0)
        {
        // nothing to do
        return PGSD_SUCCESS;
        }

    if (handle->write_buffer.size > 0 && handle->buffer_index.size == 0)
        {
        // error: bytes in buffer, but no index for them
        return PGSD_ERROR_INVALID_ARGUMENT;
        }


    // write the buffer to the end of the file
    uint64_t offset = handle->file_size;
    int j;
    for( j = 0; j < rank; j++ ){
        offset += allbuffers[j];
    }

    // printf("pgsd_flush_write_buffer rank %i: write buffer with offset %lu\n", rank, offset);
    // printf("pgsd_flush_write_buffer rank %i: write biffer size %lu\n",rank, handle->write_buffer.size);
    // printf("pgsd_flush_write_buffer rank %i: handle_buffer index size %lu\n",rank, handle->buffer_index.size);

    // ssize_t bytes_written = pgsd_io_pwrite_retry(handle->fd,
    //                                             handle->write_buffer.data,
    //                                             handle->write_buffer.size,
    //                                             offset);

    size_t bytes_written = handle->write_buffer.size;
    MPI_File_write_at(handle->fh, offset, handle->write_buffer.data, handle->write_buffer.size, MPI_BYTE, MPI_STATUS_IGNORE);
    

    if (bytes_written == -1 || bytes_written != handle->write_buffer.size)
        {
        return PGSD_ERROR_IO;
        }

    handle->file_size += handle->write_buffer.size;

    // reset write_buffer for new data
    handle->write_buffer.size = 0;

    // Move buffer_index entries to frame_index.
    size_t i;
    for (i = 0; i < handle->buffer_index.size; i++)
        {
        struct pgsd_index_entry* new_index;
        int retval = pgsd_index_buffer_add(&handle->frame_index, &new_index);
        if (retval != PGSD_SUCCESS)
            {
            return retval;
            }

        *new_index = handle->buffer_index.data[i];
        new_index->location += offset;
        }

    // clear the buffer index for new entries
    handle->buffer_index.size = 0;

    return PGSD_SUCCESS;
    }

/** @internal
    @brief Flush the name buffer.

    pgsd_write_frame() adds new names to the frame_names buffer. pgsd_flush_name_buffer() flushes
    this buffer at the end of a frame write and commits the new names to the file. If necessary,
    the namelist is written to a new location in the file.

    @param handle Handle to flush the write buffer.
    @returns PGSD_SUCCESS on success or PGSD_* error codes on error


    DK :
*/
inline static int pgsd_flush_name_buffer(struct pgsd_handle* handle)
    {

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    bool root = false;
    if ( rank == 0 ){
        root = true;
    }

    if (root == true)
        { 
        if (handle == NULL)
            {
            return PGSD_ERROR_INVALID_ARGUMENT;
            }

        if (handle->frame_names.n_names == 0)
            {
            // nothing to do
            return PGSD_SUCCESS;
            }

        if (handle->frame_names.data.size == 0)
            {
            // error: bytes in buffer, but no names for them
            return PGSD_ERROR_INVALID_ARGUMENT;
            }

        size_t old_reserved = handle->file_names.data.reserved;
        size_t old_size = handle->file_names.data.size;

        int i;
        if (root){
            for (i = 0; i < handle->frame_names.data.size; i++){
                printf("handle->frame_names.data.data %s\n", handle->frame_names.data.data);
            }
        }

        // add the new names to the file name list and zero the frame list
        int retval = pgsd_byte_buffer_append(&handle->file_names.data,
                                            handle->frame_names.data.data,
                                            handle->frame_names.data.size);
        if (retval != PGSD_SUCCESS)
            {
            return retval;
            }
        handle->file_names.n_names += handle->frame_names.n_names;
        printf("rank %i : n_names1: %i\n", rank ,handle->frame_names.n_names);
        printf("rank %i : names.size: %i\n", rank ,handle->frame_names.data.size);
        // int i;
        if (root){
            for (i = 0; i < handle->frame_names.n_names; i++){
                printf("handle->frame_names.data.data %s\n", handle->frame_names.data.data);
            }
        }

        handle->frame_names.n_names = 0;
        handle->frame_names.data.size = 0;

        pgsd_util_zero_memory(handle->frame_names.data.data, handle->frame_names.data.reserved);
        printf("rank %i : n_names2: %i\n", rank ,handle->frame_names.n_names);
        // MPI_Barrier(MPI_COMM_WORLD);

        // reserved space must be a multiple of the PGSD name size
        if (handle->file_names.data.reserved % PGSD_NAME_SIZE != 0)
            {
            return PGSD_ERROR_INVALID_ARGUMENT;
            }

        if (handle->file_names.data.reserved > old_reserved)
            {
            // write the new name list to the end of the file
            uint64_t offset = handle->file_size;

            size_t bytes_written = handle->file_names.data.reserved;
            if( root == true ){
                MPI_File_write_at(handle->fh, offset, handle->file_names.data.data, handle->file_names.data.reserved, MPI_BYTE, MPI_STATUS_IGNORE);
            }
            // ssize_t bytes_written = pgsd_io_pwrite_retry(handle->fd,
            //                                             handle->file_names.data.data,
            //                                             handle->file_names.data.reserved,
            //                                             offset);

            if (bytes_written == -1 || bytes_written != handle->file_names.data.reserved)
                {
                return PGSD_ERROR_IO;
                }

            // sync the updated name list
            // retval = fsync(handle->fd);
            // if (retval != 0)
            //     {
            //     return PGSD_ERROR_IO;
            //     }

            handle->file_size += handle->file_names.data.reserved;
            handle->header.namelist_location = offset;
            handle->header.namelist_allocated_entries
                = handle->file_names.data.reserved / PGSD_NAME_SIZE;

            // write the new header out
            bytes_written = sizeof(struct pgsd_header);
            if( root == true ){
                MPI_File_write(handle->fh, &(handle->header), sizeof(struct pgsd_header), MPI_BYTE, MPI_STATUS_IGNORE);
            }
            
            // bytes_written
            //     = pgsd_io_pwrite_retry(handle->fd, &(handle->header), sizeof(struct pgsd_header), 0);


            if (bytes_written != sizeof(struct pgsd_header))
                {
                return PGSD_ERROR_IO;
                }
            }
        else
            {
            // write the new name list to the old index location
            uint64_t offset = handle->header.namelist_location;

            size_t bytes_written = (handle->file_names.data.reserved - old_size);
            if( root == true ){
                MPI_File_write_at(handle->fh, offset + old_size, handle->file_names.data.data + old_size, handle->file_names.data.reserved + old_size, MPI_BYTE, MPI_STATUS_IGNORE);
            }

            // ssize_t bytes_written = pgsd_io_pwrite_retry(handle->fd,
            //                                             handle->file_names.data.data + old_size,
            //                                             handle->file_names.data.reserved - old_size,
            //                                             offset + old_size);
            if (bytes_written != (handle->file_names.data.reserved - old_size))
                {
                return PGSD_ERROR_IO;
                }
            }

        // sync the updated name list or header
        // retval = fsync(handle->fd);
        // if (retval != 0)
        //     {
        //     return PGSD_ERROR_IO;
        //     }
        }
    return PGSD_SUCCESS;
    }

/** @internal
    @brief utility function to append a name to the namelist

    @param id [out] ID of the new name
    @param handle handle to the open pgsd file
    @param name string name

    Append a name to the names in the current frame. pgsd_end_frame() will add this list to the
    file names.

    @return
      - PGSD_SUCCESS (0) on success. Negative value on failure:
      - PGSD_ERROR_IO: IO error (check errno).
      - PGSD_ERROR_MEMORY_ALLOCATION_FAILED: Unable to allocate memory.
      - PGSD_ERROR_FILE_MUST_BE_WRITABLE: File must not be read only.
*/
inline static int pgsd_append_name(uint16_t* id, struct pgsd_handle* handle, const char* name)
    {
    if (handle->open_flags == PGSD_OPEN_READONLY)
        {
        return PGSD_ERROR_FILE_MUST_BE_WRITABLE;
        }

    if (handle->file_names.n_names + handle->frame_names.n_names == UINT16_MAX)
        {
        // no more names may be added
        return PGSD_ERROR_NAMELIST_FULL;
        }

    // Provide the ID of the new name
    *id = (uint16_t)(handle->file_names.n_names + handle->frame_names.n_names);

    if (handle->header.pgsd_version < pgsd_make_version(2, 0))
        {
        // v1 files always allocate PGSD_NAME_SIZE bytes for each name and put a NULL terminator
        // at address 63
        char name_v1[PGSD_NAME_SIZE];
        strncpy(name_v1, name, PGSD_NAME_SIZE - 1);
        name_v1[PGSD_NAME_SIZE - 1] = 0;
        pgsd_byte_buffer_append(&handle->frame_names.data, name_v1, PGSD_NAME_SIZE);
        handle->frame_names.n_names++;

        // update the name/id mapping with the truncated name
        int retval = pgsd_name_id_map_insert(&handle->name_map, name_v1, *id);
        if (retval != PGSD_SUCCESS)
            {
            return retval;
            }
        }
    else
        {
        pgsd_byte_buffer_append(&handle->frame_names.data, name, strlen(name) + 1);
        handle->frame_names.n_names++;

        // update the name/id mapping
        int retval = pgsd_name_id_map_insert(&handle->name_map, name, *id);
        if (retval != PGSD_SUCCESS)
            {
            return retval;
            }
        }

    return PGSD_SUCCESS;
    }

/** @internal
    @brief Cross-platform wrapper for the POSIX open() system function.
    @param pathname file path using UTF-8 encoding on all platforms
    @return file descriptor
*/
// inline static int pgsd_open_file(const char* pathname, int flags, int mode)
//     {
// #ifndef _WIN32
//     return open(pathname, flags, mode);
// #else
//     // On Windows, we call the _wopen() function, which requires converting the UTF-8 input path to
//     // UTF-16 wide-character encoding.
//     int count_wchars;
//     wchar_t* wpathname;
//     int fd;

//     // First, determine the number of wide characters needed to represent the input string.
//     count_wchars = MultiByteToWideChar(CP_UTF8, 0, pathname, -1, NULL, 0);
//     // Then allocate temporary wchar_t buffer and perform the string conversion.
//     wpathname = malloc(sizeof(wchar_t) * count_wchars);
//     MultiByteToWideChar(CP_UTF8, 0, pathname, -1, wpathname, count_wchars);
//     fd = _wopen(wpathname, flags, mode);
//     free(wpathname);
//     return fd;
// #endif
//     }

/** @internal
    @brief Truncate the file and write a new pgsd header.

    @param fd file descriptor to initialize
    @param application Generating application name (truncated to 63 chars)
    @param schema Schema name for data to be written in this PGSD file (truncated to 63 chars)
    @param schema_version Version of the scheme data to be written (make with pgsd_make_version())
*/
inline static int
pgsd_initialize_file(MPI_File fh, const char* application, const char* schema, uint32_t schema_version)
    {
    printf("Initialize File\n");
    // check if the file was created
    // if (fd == -1)
    //     {
    //     return PGSD_ERROR_IO;
    //     }

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
    bool root = false;
    if ( rank == 0 ){
        root = true;
    }

    printf("MPI_File %p\n", fh);

    // check if the file was created
    MPI_Offset file_size = 0;
    MPI_File_set_size(fh, file_size);
    MPI_File_get_size(fh, &file_size);
    int retval = file_size;
    // MPI_File_seek_shared(fh, 0, MPI_SEEK_SET);
    MPI_File_seek(fh, 0, MPI_SEEK_SET);
    printf("In initil File: handle fh : %p\n", fh);

    // retval = ftruncate(fh, 0);
    if (retval != 0)
        {
        return PGSD_ERROR_IO;
        }
    if( root == true ){
    // MPI_Barrier(MPI_COMM_WORLD);
    // populate header fields
        printf("Populate Header on rank %i\n", rank);
        struct pgsd_header header;
        pgsd_util_zero_memory(&header, sizeof(header));

        header.magic = PGSD_MAGIC_ID;
        // printf("%lu \n", PGSD_MAGIC_ID);
        // printf("%lu \n", header.magic);
        header.pgsd_version = pgsd_make_version(PGSD_CURRENT_FILE_VERSION, 0);
        strncpy(header.application, application, sizeof(header.application) - 1);
        header.application[sizeof(header.application) - 1] = 0;
        strncpy(header.schema, schema, sizeof(header.schema) - 1);
        header.schema[sizeof(header.schema) - 1] = 0;
        header.schema_version = schema_version;
        header.index_location = sizeof(header);
        header.index_allocated_entries = PGSD_INITIAL_INDEX_SIZE;
        header.namelist_location
            = header.index_location + sizeof(struct pgsd_index_entry) * header.index_allocated_entries;
        header.namelist_allocated_entries = PGSD_INITIAL_NAME_BUFFER_SIZE / PGSD_NAME_SIZE;
        pgsd_util_zero_memory(header.reserved, sizeof(header.reserved));

        MPI_File_write(fh, &header, sizeof(header), MPI_BYTE, MPI_STATUS_IGNORE);
    // write the header out
    // if( rank == 0 ){

        // MPI_File_close(&fh);
        printf("Write File header on rank %i\n", rank);
        printf("print header.magic %lu\n", header.magic);
    // }
    // ssize_t bytes_written = pgsd_io_pwrite_retry(fd, &header, sizeof(header), 0);
    // if (bytes_written != sizeof(header))
    //     {
    //     return PGSD_ERROR_IO;
    //     }

    // allocate and zero default index memory
        struct pgsd_index_entry index[PGSD_INITIAL_INDEX_SIZE];
        pgsd_util_zero_memory(index, sizeof(index));
    
    // if( rank == 0 ){
        MPI_File_write(fh, &index, sizeof(index), MPI_BYTE, MPI_STATUS_IGNORE);
        printf("Write File index on rank %i\n", rank);
    // }

    // write the empty index out
    // bytes_written = pgsd_io_pwrite_retry(fd, index, sizeof(index), sizeof(header));
    // if (bytes_written != sizeof(index))
    //     {
    //     return PGSD_ERROR_IO;
    //     }

    // allocate and zero the namelist memory
        char names[PGSD_INITIAL_NAME_BUFFER_SIZE];
        pgsd_util_zero_memory(names, sizeof(char) * PGSD_INITIAL_NAME_BUFFER_SIZE);


    // if( rank == 0 ){
        MPI_File_write(fh, &names, sizeof(names), MPI_BYTE, MPI_STATUS_IGNORE);
    }

    // // write the namelist out
    // bytes_written = pgsd_io_pwrite_retry(fd, names, sizeof(names), sizeof(header) + sizeof(index));
    // if (bytes_written != sizeof(names))
    //     {
    //     return PGSD_ERROR_IO;
    //     }

    // sync file
    // retval = fsync(fh);

    // if (retval != 0)
    //     {
    //     return PGSD_ERROR_IO;
    //     }
    MPI_Barrier(MPI_COMM_WORLD);
    return PGSD_SUCCESS;
    }

/** @internal
    @brief Read in the file index and initialize the handle.

    @param handle Handle to read the header

    @pre handle->fd is an open file.
    @pre handle->open_flags is set.
*/
inline static int 
pgsd_initialize_handle(struct pgsd_handle* handle)
    {
    printf("Initialize File Handle\n");

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // check if the file was created
    if (handle->fh == NULL)
        {
        return PGSD_ERROR_IO;
        // printf("Initialize handle file doesnt exist.");
        }

    // printf("First handle->fh %p\n", handle->fh);

    // read the header
    MPI_File_seek(handle->fh, 0, MPI_SEEK_SET);
    // MPI_File_seek_shared(handle->fh, 0, MPI_SEEK_SET);
    MPI_File_read(handle->fh, &(handle->header), sizeof(struct pgsd_header), MPI_BYTE, MPI_STATUS_IGNORE);
    
    // printf("Second handle->fh %p\n", handle->fh);
    // printf("handle->header.magic %lu\n", handle->header.magic);
    // printf("handle->header.version %lu\n", handle->header.pgsd_version);
    // ssize_t bytes_read
    //     = pgsd_io_pread_retry(handle->fd, &handle->header, sizeof(struct pgsd_header), 0);
    // if (bytes_read == -1)
    //     {
    //     return PGSD_ERROR_IO;
    //     }
    // if (bytes_read != sizeof(struct pgsd_header))
    //     {
    //     return PGSD_ERROR_NOT_A_PGSD_FILE;
    //     }

    // validate the header
    if (handle->header.magic != PGSD_MAGIC_ID)
        {
        // printf("handle->header.magic %lu\n", handle->header.magic);
        return PGSD_ERROR_NOT_A_PGSD_FILE;
        }
    if (handle->header.pgsd_version < pgsd_make_version(1, 0)
        && handle->header.pgsd_version != pgsd_make_version(0, 3))
        {
        return PGSD_ERROR_INVALID_PGSD_FILE_VERSION;
        }
    // printf("PGSD File Version correct!\n");

    if (handle->header.pgsd_version >= pgsd_make_version(3, 0))
        {
        return PGSD_ERROR_INVALID_PGSD_FILE_VERSION;
        }
    // printf("PGSD Version correct!\n");

    // determine the file size
    // handle->file_size = lseek(handle->fd, 0, SEEK_END);
    // printf("Get Size!\n");
    MPI_File_get_size(handle->fh, &(handle->file_size));
    // printf("File see end !\n");
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_File_seek(handle->fh, 0, MPI_SEEK_END);
    // printf("rank %i: File see end correct!\n", rank);
    
    // validate that the namelist block exists inside the file
    if (handle->header.namelist_location
            + (PGSD_NAME_SIZE * handle->header.namelist_allocated_entries)
        > (uint64_t)handle->file_size)
        {
        printf("Here is the return value for corrupt file 1\n");
        return PGSD_ERROR_FILE_CORRUPT;
        }

    // printf("rank %i: Namelist block validated!\n", rank);
    // allocate the hash map
    int retval = pgsd_name_id_map_allocate(&handle->name_map, PGSD_NAME_MAP_SIZE);
    if (retval != PGSD_SUCCESS)
        {
        printf("error in allocate has map\n");

        return retval;
        }

    printf("allocate hashmap done!\n");
    // read the namelist block
    size_t namelist_n_bytes = PGSD_NAME_SIZE * handle->header.namelist_allocated_entries;
    retval = pgsd_byte_buffer_allocate(&handle->file_names.data, namelist_n_bytes);
    if (retval != PGSD_SUCCESS)
        {
        printf("error in allocate has byte buffer\n");
        return retval;
        }
    printf("read namelist block done!\n");

    // ssize_t bytes_read = namelist_n_bytes;
    MPI_File_read_at(handle->fh, handle->header.namelist_location, handle->file_names.data.data, namelist_n_bytes, MPI_BYTE, MPI_STATUS_IGNORE);
    
    // bytes_read = pgsd_io_pread_retry(handle->fd,
    //                                 handle->file_names.data.data,
    //                                 namelist_n_bytes,
    //                                 handle->header.namelist_location);


    // if (bytes_read == -1 || bytes_read != namelist_n_bytes)
    //     {
    //     return PGSD_ERROR_IO;
    //     }
    printf("Test namelist bytes read!\n");

    // The name buffer must end in a NULL terminator or else the file is corrupt
    if (handle->file_names.data.data[handle->file_names.data.reserved - 1] != 0)
        {
        printf("Here is the return value for corrupt file 1\n");
        return PGSD_ERROR_FILE_CORRUPT;
        }
    printf("Test name buffer corruoptnes!\n");

    // Add the names to the hash map. Also determine the number of used bytes in the namelist.
    size_t name_start = 0;
    handle->file_names.n_names = 0;
    while (name_start < handle->file_names.data.reserved)
        {
        char* name = handle->file_names.data.data + name_start;

        // an empty name notes the end of the list
        if (name[0] == 0)
            {
            break;
            }

        retval
            = pgsd_name_id_map_insert(&handle->name_map, name, (uint16_t)handle->file_names.n_names);
        if (retval != PGSD_SUCCESS)
            {
            printf("error in id map insert\n");
            return retval;
            }
        handle->file_names.n_names++;

        if (handle->header.pgsd_version < pgsd_make_version(2, 0))
            {
            // pgsd v1 stores names in fixed 64 byte segments
            name_start += PGSD_NAME_SIZE;
            }
        else
            {
            size_t len = strnlen(name, handle->file_names.data.reserved - name_start);
            name_start += len + 1;
            }
        }

    handle->file_names.data.size = name_start;
    printf("ADD names to hashmap done successfully!\n");

    // read in the file index
    retval = pgsd_index_buffer_map(&handle->file_index, handle);
    if (retval != PGSD_SUCCESS)
        {
        printf("error in index buffer map\n");
        return retval;
        }

    printf("read in the file index done successfully!\n");
    // determine the current frame counter
    if (handle->file_index.size == 0)
        {
        handle->cur_frame = 0;
        }
    else
        {
        handle->cur_frame = handle->file_index.data[handle->file_index.size - 1].frame + 1;
        }
    printf("determine current file counter successfully!\n");

    // if this is a write mode, allocate the initial frame index and the name buffer
    if (handle->open_flags != PGSD_OPEN_READONLY)
        {
        retval = pgsd_index_buffer_allocate(&handle->frame_index, PGSD_INITIAL_FRAME_INDEX_SIZE);
        if (retval != PGSD_SUCCESS)
            {
            // printf("error in buffer allocate\n");
            return retval;
            }

        retval = pgsd_index_buffer_allocate(&handle->buffer_index, PGSD_INITIAL_FRAME_INDEX_SIZE);
        if (retval != PGSD_SUCCESS)
            {
            // printf("error in buffer allocate2\n");
            return retval;
            }

        retval = pgsd_byte_buffer_allocate(&handle->write_buffer, PGSD_INITIAL_WRITE_BUFFER_SIZE);
        if (retval != PGSD_SUCCESS)
            {
            // printf("error in buffer allocate3\n");
            return retval;
            }

        handle->frame_names.n_names = 0;
        retval = pgsd_byte_buffer_allocate(&handle->frame_names.data, PGSD_NAME_SIZE);
        if (retval != PGSD_SUCCESS)
            {
            // printf("error in buffer allocate4\n");
            return retval;
            }
        }

    handle->pending_index_entries = 0;
    handle->maximum_write_buffer_size = PGSD_DEFAULT_MAXIMUM_WRITE_BUFFER_SIZE;
    handle->index_entries_to_buffer = PGSD_DEFAULT_INDEX_ENTRIES_TO_BUFFER;

    return PGSD_SUCCESS;
    }

uint32_t pgsd_make_version(unsigned int major, unsigned int minor)
    {
    return major << (sizeof(uint32_t) * 4) | minor;
    }

// int pgsd_create(const char* fname,
//                const char* application,
//                const char* schema,
//                uint32_t schema_version)
//     {
//     printf("Create pgsd\n");
//     int extra_flags = 0;
// #ifdef _WIN32
//     extra_flags = _O_BINARY;
// #endif

//     // create the file
//     // int fd = pgsd_open_file(fname,
//     //                        O_RDWR | O_CREAT | O_TRUNC | extra_flags,
//     //                        S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP);
//     MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &(handle->fh));
    
//     int retval = pgsd_initialize_file(&(handle->fh), application, schema, schema_version);
//     if (fh != -1)
//         {
//         close(fh);
//         }
//     return retval;
//     }

int pgsd_create_and_open(struct pgsd_handle* handle,
                        const char* fname,
                        const char* application,
                        const char* schema,
                        uint32_t schema_version,
                        const enum pgsd_open_flag flags,
                        int exclusive_create)
    {
    // zero the handle
    printf("Start pgsd until zero memory in create anbd open\n");
    pgsd_util_zero_memory(handle, sizeof(struct pgsd_handle));

    int extra_flags = 0;
#ifdef _WIN32
    extra_flags = _O_BINARY;
#endif

    // set the open flags in the handle
    if (flags == PGSD_OPEN_READWRITE)
        {
        handle->open_flags = PGSD_OPEN_READWRITE;
        }
    else if (flags == PGSD_OPEN_READONLY)
        {
        return PGSD_ERROR_FILE_MUST_BE_WRITABLE;
        }
    else if (flags == PGSD_OPEN_APPEND)
        {
        handle->open_flags = PGSD_OPEN_APPEND;
        }

    // set the exclusive create bit
    if (exclusive_create)
        {
        extra_flags |= O_EXCL;
        }

    // create the file
    // handle->fd = pgsd_open_file(fname,
    //                            O_RDWR | O_CREAT | O_TRUNC | extra_flags,
    //                            S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP);

    printf("Initialize MPI FILE\n");
    MPI_File fh;
    
    handle->fh = fh;
    printf("Initialize\n");
    // printf("File handle %p\n", fh);
    MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &(handle->fh));

    // printf("File handle->fh %p\n", handle->fh);
    // MPI_Barrier(MPI_COMM_WORLD);
    
    printf("Opened MPI FILE\n");
    int retval = pgsd_initialize_file(handle->fh, application, schema, schema_version);
    MPI_Barrier(MPI_COMM_WORLD);
    // printf("After initialized file retval = %i\n", retval);
    if (retval != 0)
        {
        if (retval != -1)
            {
            MPI_File_close(&(handle->fh));
            }
        return retval;
        }
    printf("retval before init handle: %i\n", retval);
    retval = pgsd_initialize_handle(handle);
    if (retval != 0)
        {
        if (handle->fh != NULL)
            {
            MPI_File_close(&(handle->fh));
            }
        }
    printf("retval after init handle: %i\n", retval);
    return retval;
    }

int pgsd_open(struct pgsd_handle* handle, const char* fname, const enum pgsd_open_flag flags)
    {
    // zero the handle
    pgsd_util_zero_memory(handle, sizeof(struct pgsd_handle));
    printf("Filename ind pgsd_open %s\n", fname);
//     int extra_flags = 0;
// #ifdef _WIN32
//     extra_flags = _O_BINARY;
// #endif

    // open the file
    if (flags == PGSD_OPEN_READWRITE)
        {
        MPI_File_open(MPI_COMM_WORLD ,fname, MPI_MODE_RDWR, MPI_INFO_NULL, &(handle->fh));
        // handle->fd = pgsd_open_file(fname, O_RDWR | extra_flags, 0);
        handle->open_flags = PGSD_OPEN_READWRITE;
        }
    else if (flags == PGSD_OPEN_READONLY)
        {
        MPI_File_open(MPI_COMM_WORLD ,fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &(handle->fh));
        // handle->fd = pgsd_open_file(fname, O_RDONLY | extra_flags, 0);
        handle->open_flags = PGSD_OPEN_READONLY;
        }
    else if (flags == PGSD_OPEN_APPEND)
        {
        MPI_File_open(MPI_COMM_WORLD ,fname, MPI_MODE_RDWR, MPI_INFO_NULL, &(handle->fh));
        // handle->fd = pgsd_open_file(fname, O_RDWR | extra_flags, 0);
        handle->open_flags = PGSD_OPEN_APPEND;
        }

    int retval = pgsd_initialize_handle(handle);
    printf("retval after open file write only %i\n", retval);
    if (retval != 0)
        {
        if (handle->fh != NULL)
            {
            MPI_File_close(&(handle->fh));
            }
        }

    return retval;
    }

// int pgsd_truncate(struct pgsd_handle* handle)
//     {
//     if (handle == NULL)
//         {
//         return PGSD_ERROR_INVALID_ARGUMENT;
//         }
//     if (handle->open_flags == PGSD_OPEN_READONLY)
//         {
//         return PGSD_ERROR_FILE_MUST_BE_WRITABLE;
//         }

//     int retval = 0;

//     // deallocate indices
//     if (handle->frame_names.data.reserved > 0)
//         {
//         retval = pgsd_byte_buffer_free(&handle->frame_names.data);
//         if (retval != PGSD_SUCCESS)
//             {
//             return retval;
//             }
//         }

//     if (handle->file_names.data.reserved > 0)
//         {
//         retval = pgsd_byte_buffer_free(&handle->file_names.data);
//         if (retval != PGSD_SUCCESS)
//             {
//             return retval;
//             }
//         }

//     retval = pgsd_name_id_map_free(&handle->name_map);
//     if (retval != PGSD_SUCCESS)
//         {
//         return retval;
//         }

//     retval = pgsd_index_buffer_free(&handle->file_index);
//     if (retval != PGSD_SUCCESS)
//         {
//         return retval;
//         }

//     if (handle->frame_index.reserved > 0)
//         {
//         retval = pgsd_index_buffer_free(&handle->frame_index);
//         if (retval != PGSD_SUCCESS)
//             {
//             return retval;
//             }
//         }

//     if (handle->buffer_index.reserved > 0)
//         {
//         retval = pgsd_index_buffer_free(&handle->buffer_index);
//         if (retval != PGSD_SUCCESS)
//             {
//             return retval;
//             }
//         }

//     if (handle->write_buffer.reserved > 0)
//         {
//         retval = pgsd_byte_buffer_free(&handle->write_buffer);
//         if (retval != PGSD_SUCCESS)
//             {
//             return retval;
//             }
//         }

//     // keep a copy of the old header
//     struct pgsd_header old_header = handle->header;
//     retval = pgsd_initialize_file(handle->fh,
//                                  old_header.application,
//                                  old_header.schema,
//                                  old_header.schema_version);

//     if (retval != PGSD_SUCCESS)
//         {
//         return retval;
//         }

//     return pgsd_initialize_handle(handle);
//     }

int pgsd_close(struct pgsd_handle* handle)
    {
    if (handle == NULL)
        {
        return PGSD_ERROR_INVALID_ARGUMENT;
        }

    int retval;

    if (handle->open_flags != PGSD_OPEN_READONLY)
        {
        retval = pgsd_flush(handle);
        if (retval != PGSD_SUCCESS)
            {
            return retval;
            }
        }
    MPI_Barrier(MPI_COMM_WORLD);
    // save the fd so we can use it after freeing the handle
    // int fd = handle->fd;
    // MPI_File fh = handle->fh;

    retval = pgsd_index_buffer_free(&handle->file_index);
    if (retval != PGSD_SUCCESS)
        {
        return retval;
        }

    if (handle->frame_index.reserved > 0)
        {
        retval = pgsd_index_buffer_free(&handle->frame_index);
        if (retval != PGSD_SUCCESS)
            {
            return retval;
            }
        }

    if (handle->buffer_index.reserved > 0)
        {
        retval = pgsd_index_buffer_free(&handle->buffer_index);
        if (retval != PGSD_SUCCESS)
            {
            return retval;
            }
        }

    if (handle->write_buffer.reserved > 0)
        {
        retval = pgsd_byte_buffer_free(&handle->write_buffer);
        if (retval != PGSD_SUCCESS)
            {
            return retval;
            }
        }

    retval = pgsd_name_id_map_free(&handle->name_map);
    if (retval != PGSD_SUCCESS)
        {
        return retval;
        }

    if (handle->frame_names.data.reserved > 0)
        {
        handle->frame_names.n_names = 0;
        retval = pgsd_byte_buffer_free(&handle->frame_names.data);
        if (retval != PGSD_SUCCESS)
            {
            return retval;
            }
        }

    if (handle->file_names.data.reserved > 0)
        {
        handle->file_names.n_names = 0;
        retval = pgsd_byte_buffer_free(&handle->file_names.data);
        if (retval != PGSD_SUCCESS)
            {
            return retval;
            }
        }

    // close the file
    // retval = close(fh);
    MPI_Barrier(MPI_COMM_WORLD);
    retval = MPI_File_close(&(handle->fh));

    if (retval != 0)
        {
        return PGSD_ERROR_IO;
        }

    return PGSD_SUCCESS;
    }

int pgsd_end_frame(struct pgsd_handle* handle)
    {
    if (handle == NULL)
        {
        return PGSD_ERROR_INVALID_ARGUMENT;
        }
    if (handle->open_flags == PGSD_OPEN_READONLY)
        {
        return PGSD_ERROR_FILE_MUST_BE_WRITABLE;
        }

    // increment the frame counter
    handle->cur_frame++;
    handle->pending_index_entries = 0;

    if (handle->frame_index.size > 0 || handle->buffer_index.size > handle->index_entries_to_buffer)
        {
        return pgsd_flush(handle);
        }

    return PGSD_SUCCESS;
    }

int pgsd_flush(struct pgsd_handle* handle)
    {
    if (handle == NULL)
        {
        return PGSD_ERROR_INVALID_ARGUMENT;
        }
    if (handle->open_flags == PGSD_OPEN_READONLY)
        {
        return PGSD_ERROR_FILE_MUST_BE_WRITABLE;
        }

    int rank;
    bool root = false;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if( rank == 0 ){
        root = true;
    }

    // flush the namelist buffer
    int retval = pgsd_flush_name_buffer(handle);
    if (retval != PGSD_SUCCESS)
        {
        return retval;
        }

    // flush the write buffer
    retval = pgsd_flush_write_buffer(handle);
    if (retval != PGSD_SUCCESS)
        {
        return retval;
        }

    // sync the data before writing the index
    // retval = fsync(handle->fh);
    // if (retval != 0)
    //     {
    //     return PGSD_ERROR_IO;
    //     }
    if ( root == true ){
        // Write the frame index to the file, excluding the index entries that are part of the current
        // frame.
        if (handle->pending_index_entries > handle->frame_index.size)
            {
            return PGSD_ERROR_INVALID_ARGUMENT;
            }
        uint64_t index_entries_to_write = handle->frame_index.size - handle->pending_index_entries;

        if (index_entries_to_write > 0)
            {
            // ensure there is enough space in the index
            if ((handle->file_index.size + index_entries_to_write) > handle->file_index.reserved)
                {
                pgsd_expand_file_index(handle, handle->file_index.size + index_entries_to_write);
                }

            // sort the index before writing
            retval = pgsd_index_buffer_sort(&handle->frame_index);
            if (retval != 0)
                {
                return retval;
                }

            // write the frame index entries to the file
            int64_t write_pos = handle->header.index_location
                                + sizeof(struct pgsd_index_entry) * handle->file_index.size;

            size_t bytes_to_write = sizeof(struct pgsd_index_entry) * index_entries_to_write;
            // ssize_t bytes_written
            //     = pgsd_io_pwrite_retry(handle->fd, handle->frame_index.data, bytes_to_write, write_pos);

            size_t bytes_written = sizeof(struct pgsd_index_entry) * index_entries_to_write;
            MPI_File_write_at(handle->fh, write_pos, handle->frame_index.data, sizeof(struct pgsd_index_entry) * handle->frame_index.size, MPI_BYTE, MPI_STATUS_IGNORE);

            if (bytes_written == -1 || bytes_written != bytes_to_write)
                {
                return PGSD_ERROR_IO;
                }

    #if !PGSD_USE_MMAP
            // add the entries to the file index
            memcpy(handle->file_index.data + handle->file_index.size,
                   handle->frame_index.data,
                   sizeof(struct pgsd_index_entry) * handle->frame_index.size);
    #endif

            // update size of file index
            handle->file_index.size += index_entries_to_write;

            // Clear the frame index, keeping those in the current unfinished frame.
            if (handle->pending_index_entries > 0)
                {
                for (uint64_t i = 0; i < handle->pending_index_entries; i++)
                    {
                    handle->frame_index.data[i]
                        = handle->frame_index
                              .data[handle->frame_index.size - handle->pending_index_entries];
                    }
                }

            handle->frame_index.size = handle->pending_index_entries;
            }
        }

    return PGSD_SUCCESS;
    }

int pgsd_write_chunk(struct pgsd_handle* handle,
                    const char* name,
                    enum pgsd_type type,
                    uint64_t N,
                    uint32_t M,
                    uint64_t N_global,
                    uint32_t M_global,
                    uint64_t offset,
                    uint64_t global_size,
                    bool all,
                    uint8_t flags,
                    const void* data)
    {
    // validate input
    if (N > 0 && data == NULL)
        {
        return PGSD_ERROR_INVALID_ARGUMENT;
        }
    if (M == 0)
        {
        return PGSD_ERROR_INVALID_ARGUMENT;
        }
    if (handle->open_flags == PGSD_OPEN_READONLY)
        {
        return PGSD_ERROR_FILE_MUST_BE_WRITABLE;
        }
    if (flags != 0)
        {
        return PGSD_ERROR_INVALID_ARGUMENT;
        }

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    bool root = false;
    if ( rank == 0 ){
        root = true;
    }

    uint16_t id = pgsd_name_id_map_find(&handle->name_map, name);
    printf("Rank %i write chunk name:%s\n",rank, name);
    if (id == UINT16_MAX)
        {
        // not found, append to the index
        if ( root == true ){ 
            int retval = pgsd_append_name(&id, handle, name);
            if (retval != PGSD_SUCCESS)
                {
                return retval;
                }
            if (id == UINT16_MAX)
                {
                // this should never happen
                return PGSD_ERROR_NAMELIST_FULL;
                }
            }
        }

    struct pgsd_index_entry entry;
    // populate fields in the entry's data
    pgsd_util_zero_memory(&entry, sizeof(struct pgsd_index_entry));
    entry.frame = handle->cur_frame;
    entry.id = id;
    entry.type = (uint8_t)type;
    entry.N = N_global;
    entry.M = M_global;
    size_t size = N * M * pgsd_sizeof_type(type);

    global_size *= pgsd_sizeof_type(type);
    offset *= pgsd_sizeof_type(type);
    if(global_size == 0 && offset ==0){
        global_size = size;
    }

    // decide whether to write this chunk to the buffer or straight to disk
    if (size < handle->maximum_write_buffer_size)
        {
        // flush the buffer if this entry won't fit
        if (size > (handle->maximum_write_buffer_size - handle->write_buffer.size))
            {
            pgsd_flush_write_buffer(handle);
            }

        entry.location = handle->write_buffer.size + offset;

        int retval = 0;

        // add an entry to the buffer index
        struct pgsd_index_entry* index_entry;

        retval = pgsd_index_buffer_add(&handle->buffer_index, &index_entry);
        if (retval != PGSD_SUCCESS)
            {
            return retval;
            }
        *index_entry = entry;

        // add the data to the write buffer
        if (size > 0)
            {
            retval = pgsd_byte_buffer_append(&handle->write_buffer, data, size);
            if (retval != PGSD_SUCCESS)
                {
                return retval;
                }
            }
        }
    else
        {
        // add an entry to the frame index
        struct pgsd_index_entry* index_entry;

        int retval = pgsd_index_buffer_add(&handle->frame_index, &index_entry);
        if (retval != PGSD_SUCCESS)
            {
            return retval;
            }
        *index_entry = entry;

        // find the location at the end of the file for the chunk
        index_entry->location = handle->file_size + offset;


        // write the data
        // ssize_t bytes_written = pgsd_io_pwrite_retry(handle->fd, data, size, index_entry->location);
        
        size_t bytes_written = size;
        if( all == true || root == true ){
            printf("rank %i: Write %s data at index entry location: %lu\n", rank, name, index_entry->location);
            printf("rank %i: Write %s data at index entry file_size: %llu\n", rank, name, handle->file_size);
            printf("rank %i: Write %s data at index entry offset: %lu\n", rank, name, offset);
            MPI_File_write_at(handle->fh, index_entry->location, data, size, MPI_BYTE, MPI_STATUS_IGNORE);
        }
        if (bytes_written == -1 || bytes_written != size)
            {
            return PGSD_ERROR_IO;
            }

        // update the file_size in the handle
        handle->file_size += bytes_written;
        }
    handle->pending_index_entries++;
    return PGSD_SUCCESS;
    }

uint64_t pgsd_get_nframes(struct pgsd_handle* handle)
    {
    if (handle == NULL)
        {
        return 0;
        }
    return handle->cur_frame;
    }

const struct pgsd_index_entry*
pgsd_find_chunk(struct pgsd_handle* handle, uint64_t frame, const char* name)
    {
    if (handle == NULL)
        {
        return NULL;
        }
    if (name == NULL)
        {
        return NULL;
        }
    if (frame >= pgsd_get_nframes(handle))
        {
        return NULL;
        }
    if (handle->open_flags != PGSD_OPEN_READONLY)
        {
        int retval = pgsd_flush(handle);
        if (retval != PGSD_SUCCESS)
            {
            return NULL;
            }
        }

    // find the id for the given name
    uint16_t match_id = pgsd_name_id_map_find(&handle->name_map, name);
    if (match_id == UINT16_MAX)
        {
        return NULL;
        }

    if (handle->header.pgsd_version >= pgsd_make_version(2, 0))
        {
        // pgsd 2.0 files sort the entire index
        // binary search for the index entry
        ssize_t L = 0;
        ssize_t R = handle->file_index.size - 1;
        struct pgsd_index_entry T;
        T.frame = frame;
        T.id = match_id;

        while (L <= R)
            {
            size_t m = (L + R) / 2;
            int cmp = pgsd_cmp_index_entry(handle->file_index.data + m, &T);
            if (cmp == -1)
                {
                L = m + 1;
                }
            else if (cmp == 1)
                {
                R = m - 1;
                }
            else
                {
                return &(handle->file_index.data[m]);
                }
            }
        }
    else
        {
        // pgsd 1.0 file: use binary search to find the frame and linear search to find the entry
        size_t L = 0;
        size_t R = handle->file_index.size;

        // progressively narrow the search window by halves
        do
            {
            size_t m = (L + R) / 2;

            if (frame < handle->file_index.data[m].frame)
                {
                R = m;
                }
            else
                {
                L = m;
                }
            } while ((R - L) > 1);

        // this finds L = the rightmost index with the desired frame
        int64_t cur_index;

        // search all index entries with the matching frame
        for (cur_index = L; (cur_index >= 0) && (handle->file_index.data[cur_index].frame == frame);
             cur_index--)
            {
            // if the frame matches, check the id
            if (match_id == handle->file_index.data[cur_index].id)
                {
                return &(handle->file_index.data[cur_index]);
                }
            }
        }

    // if we got here, we didn't find the specified chunk
    return NULL;
    }

int pgsd_read_chunk(struct pgsd_handle* handle, void* data, const struct pgsd_index_entry* chunk, uint32_t * offset, bool all)
    {
    if (handle == NULL)
        {
        return PGSD_ERROR_INVALID_ARGUMENT;
        }
    if (data == NULL)
        {
        return PGSD_ERROR_INVALID_ARGUMENT;
        }
    if (chunk == NULL)
        {
        return PGSD_ERROR_INVALID_ARGUMENT;
        }
    if (handle->open_flags != PGSD_OPEN_READONLY)
        {
        int retval = pgsd_flush(handle);
        if (retval != PGSD_SUCCESS)
            {
            return retval;
            }
        }

    uint64_t stride = 0;
    int rank;
    int nprocs;
    unsigned int loc_N = chunk->N;
    size_t size;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    // printf("pgsd read chuck: rank %i, size %i\n", rank, nprocs);
    // printf("pgsd read chuck: offset[0] %i, stride %i\n", offset[0], stride);

    if (all == false){
        stride = 0;
        size = chunk->N * chunk->M * pgsd_sizeof_type((enum pgsd_type)chunk->type);
    }
    else if (all == true){
        if(offset != NULL) {
            unsigned int i;
            for( i = 0; i < rank; i++){
                stride += offset[i];
            }
            loc_N = offset[rank];
        }
        // printf("pgsd_sizeof_type(chunk->type) %i\n", pgsd_sizeof_type(chunk->type));
        stride *= chunk->M *pgsd_sizeof_type(chunk->type);
        size = loc_N * chunk->M * pgsd_sizeof_type(chunk->type);
        // printf("pgsd read chuck: stride %i, size %i\n", stride, size);
        // printf("pgsd read chuck: chunk->M %i, loc_N %i\n", chunk->M, loc_N);

    }

    // size_t size = chunk->N * chunk->M * pgsd_sizeof_type((enum pgsd_type)chunk->type);
    if (size == 0)
        {
        printf("Here is the return value for corrupt file 3\n");
        return PGSD_ERROR_FILE_CORRUPT;
        }
    if (chunk->location == 0)
        {
        printf("Here is the return value for corrupt file 4\n");
        return PGSD_ERROR_FILE_CORRUPT;
        }

    // validate that we don't read past the end of the file
    if ((chunk->location + size + stride) > (uint64_t)handle->file_size)
        {
        printf("Here is the return value for corrupt file 5\n");
        return PGSD_ERROR_FILE_CORRUPT;
        }
    // printf("handle->fh: %lu\n", handle->fh);
    // printf("size: %i\n", size);
    // printf("chunk->location+stride: %lu\n", chunk->location+stride);
    // printf("chunk->location: %lu\n", chunk->location);
    // printf("before read file at\n");
    MPI_File_read_at(handle->fh, chunk->location+stride, data, size, MPI_BYTE, MPI_STATUS_IGNORE);
    // ssize_t bytes_read = pgsd_io_pread_retry(handle->fd, data, size, chunk->location);
    // if (bytes_read == -1 || bytes_read != size)
    //     {
    //     return PGSD_ERROR_IO;
    //     }
    printf("after read file at\n");

    return PGSD_SUCCESS;
    }

size_t pgsd_sizeof_type(enum pgsd_type type)
    {
    size_t val = 0;
    if (type == PGSD_TYPE_UINT8)
        {
        val = sizeof(uint8_t);
        }
    else if (type == PGSD_TYPE_UINT16)
        {
        val = sizeof(uint16_t);
        }
    else if (type == PGSD_TYPE_UINT32)
        {
        val = sizeof(uint32_t);
        }
    else if (type == PGSD_TYPE_UINT64)
        {
        val = sizeof(uint64_t);
        }
    else if (type == PGSD_TYPE_INT8)
        {
        val = sizeof(int8_t);
        }
    else if (type == PGSD_TYPE_INT16)
        {
        val = sizeof(int16_t);
        }
    else if (type == PGSD_TYPE_INT32)
        {
        val = sizeof(int32_t);
        }
    else if (type == PGSD_TYPE_INT64)
        {
        val = sizeof(int64_t);
        }
    else if (type == PGSD_TYPE_FLOAT)
        {
        val = sizeof(float);
        }
    else if (type == PGSD_TYPE_DOUBLE)
        {
        val = sizeof(double);
        }
    else
        {
        return 0;
        }
    return val;
    }

const char*
pgsd_find_matching_chunk_name(struct pgsd_handle* handle, const char* match, const char* prev)
    {
    if (handle == NULL)
        {
        return NULL;
        }
    if (match == NULL)
        {
        return NULL;
        }
    if (handle->file_names.n_names == 0)
        {
        return NULL;
        }
    if (handle->open_flags != PGSD_OPEN_READONLY)
        {
        int retval = pgsd_flush(handle);
        if (retval != PGSD_SUCCESS)
            {
            return NULL;
            }
        }

    // return nothing found if the name buffer is corrupt
    if (handle->file_names.data.data[handle->file_names.data.reserved - 1] != 0)
        {
        return NULL;
        }

    // determine search start index
    const char* search_str;
    if (prev == NULL)
        {
        search_str = handle->file_names.data.data;
        }
    else
        {
        // return not found if prev is not in range
        if (prev < handle->file_names.data.data)
            {
            return NULL;
            }
        if (prev >= (handle->file_names.data.data + handle->file_names.data.reserved))
            {
            return NULL;
            }

        if (handle->header.pgsd_version < pgsd_make_version(2, 0))
            {
            search_str = prev + PGSD_NAME_SIZE;
            }
        else
            {
            search_str = prev + strlen(prev) + 1;
            }
        }

    size_t match_len = strlen(match);

    while (search_str < (handle->file_names.data.data + handle->file_names.data.reserved))
        {
        if (search_str[0] != 0 && 0 == strncmp(match, search_str, match_len))
            {
            return search_str;
            }

        if (handle->header.pgsd_version < pgsd_make_version(2, 0))
            {
            search_str += PGSD_NAME_SIZE;
            }
        else
            {
            search_str += strlen(search_str) + 1;
            }
        }

    // searched past the end of the list, return NULL
    return NULL;
    }

// int pgsd_upgrade(struct pgsd_handle* handle)
//     {
//     if (handle == NULL)
//         {
//         return PGSD_ERROR_INVALID_ARGUMENT;
//         }
//     if (handle->open_flags == PGSD_OPEN_READONLY)
//         {
//         return PGSD_ERROR_INVALID_ARGUMENT;
//         }
//     if (handle->frame_index.size > 0 || handle->frame_names.n_names > 0)
//         {
//         return PGSD_ERROR_INVALID_ARGUMENT;
//         }

//     if (handle->header.pgsd_version < pgsd_make_version(2, 0))
//         {
//         if (handle->file_index.size > 0)
//             {
//             // make a copy of the file index
//             struct pgsd_index_buffer buf;
//             pgsd_util_zero_memory(&buf, sizeof(struct pgsd_index_buffer));
//             int retval = pgsd_index_buffer_allocate(&buf, handle->file_index.size);
//             if (retval != PGSD_SUCCESS)
//                 {
//                 return retval;
//                 }
//             memcpy(buf.data,
//                    handle->file_index.data,
//                    sizeof(struct pgsd_index_entry) * handle->file_index.size);
//             buf.size = handle->file_index.size;

//             // sort the copy and write it back out to the file
//             retval = pgsd_index_buffer_sort(&buf);
//             if (retval != PGSD_SUCCESS)
//                 {
//                 pgsd_index_buffer_free(&buf);
//                 return retval;
//                 }

//             ssize_t bytes_written = pgsd_io_pwrite_retry(handle->fd,
//                                                         buf.data,
//                                                         sizeof(struct pgsd_index_entry) * buf.size,
//                                                         handle->header.index_location);

//             if (bytes_written == -1 || bytes_written != sizeof(struct pgsd_index_entry) * buf.size)
//                 {
//                 pgsd_index_buffer_free(&buf);
//                 return PGSD_ERROR_IO;
//                 }

//             retval = pgsd_index_buffer_free(&buf);
//             if (retval != PGSD_SUCCESS)
//                 {
//                 return retval;
//                 }

//             // sync the updated index
//             retval = fsync(handle->fd);
//             if (retval != 0)
//                 {
//                 return PGSD_ERROR_IO;
//                 }
//             }

//         if (handle->file_names.n_names > 0)
//             {
//             // compact the name list without changing its size or position on the disk
//             struct pgsd_byte_buffer new_name_buf;
//             pgsd_util_zero_memory(&new_name_buf, sizeof(struct pgsd_byte_buffer));
//             int retval = pgsd_byte_buffer_allocate(&new_name_buf, handle->file_names.data.reserved);
//             if (retval != PGSD_SUCCESS)
//                 {
//                 return retval;
//                 }

//             const char* name = pgsd_find_matching_chunk_name(handle, "", NULL);
//             while (name != NULL)
//                 {
//                 retval = pgsd_byte_buffer_append(&new_name_buf, name, strlen(name) + 1);
//                 if (retval != PGSD_SUCCESS)
//                     {
//                     pgsd_byte_buffer_free(&new_name_buf);
//                     return retval;
//                     }
//                 name = pgsd_find_matching_chunk_name(handle, "", name);
//                 }

//             if (new_name_buf.reserved != handle->file_names.data.reserved)
//                 {
//                 pgsd_byte_buffer_free(&new_name_buf);
//                 return PGSD_ERROR_FILE_CORRUPT;
//                 }

//             // write the new names out to disk
//             ssize_t bytes_written = pgsd_io_pwrite_retry(handle->fd,
//                                                         new_name_buf.data,
//                                                         new_name_buf.reserved,
//                                                         handle->header.namelist_location);

//             if (bytes_written == -1 || bytes_written != new_name_buf.reserved)
//                 {
//                 pgsd_byte_buffer_free(&new_name_buf);
//                 return PGSD_ERROR_IO;
//                 }

//             // swap in the re-organized name buffer
//             retval = pgsd_byte_buffer_free(&handle->file_names.data);
//             if (retval != PGSD_SUCCESS)
//                 {
//                 pgsd_byte_buffer_free(&new_name_buf);
//                 return retval;
//                 }
//             handle->file_names.data = new_name_buf;

//             // sync the updated name list
//             retval = fsync(handle->fd);
//             if (retval != 0)
//                 {
//                 pgsd_byte_buffer_free(&new_name_buf);
//                 return PGSD_ERROR_IO;
//                 }
//             }

//         // label the file as a v2.0 file
//         handle->header.pgsd_version = pgsd_make_version(PGSD_CURRENT_FILE_VERSION, 0);

//         // write the new header out
//         ssize_t bytes_written
//             = pgsd_io_pwrite_retry(handle->fd, &(handle->header), sizeof(struct pgsd_header), 0);
//         if (bytes_written != sizeof(struct pgsd_header))
//             {
//             return PGSD_ERROR_IO;
//             }

//         // sync the updated header
//         int retval = fsync(handle->fd);
//         if (retval != 0)
//             {
//             return PGSD_ERROR_IO;
//             }

//         // remap the file index
//         retval = pgsd_index_buffer_free(&handle->file_index);
//         if (retval != 0)
//             {
//             return retval;
//             }

//         retval = pgsd_index_buffer_map(&handle->file_index, handle);
//         if (retval != 0)
//             {
//             return retval;
//             }
//         }

//     return PGSD_SUCCESS;
//     }

uint64_t pgsd_get_maximum_write_buffer_size(struct pgsd_handle* handle)
    {
    if (handle == NULL)
        {
        return 0;
        }
    return handle->maximum_write_buffer_size;
    }

int pgsd_set_maximum_write_buffer_size(struct pgsd_handle* handle, uint64_t size)
    {
    if (handle == NULL || size == 0)
        {
        return PGSD_ERROR_INVALID_ARGUMENT;
        }

    handle->maximum_write_buffer_size = size;

    return PGSD_SUCCESS;
    }

uint64_t pgsd_get_index_entries_to_buffer(struct pgsd_handle* handle)
    {
    if (handle == NULL)
        {
        return 0;
        }
    return handle->index_entries_to_buffer;
    }

int pgsd_set_index_entries_to_buffer(struct pgsd_handle* handle, uint64_t number)
    {
    if (handle == NULL || number == 0)
        {
        return PGSD_ERROR_INVALID_ARGUMENT;
        }

    handle->index_entries_to_buffer = number;

    return PGSD_SUCCESS;
    }

// undefine windows wrapper macros
// #ifdef _WIN32
// #undef lseek
// #undef write
// #undef read
// #undef open
// #undef ftruncate
// #pragma warning(pop)

// #endif

// #endif // End ENABLE_MPI
