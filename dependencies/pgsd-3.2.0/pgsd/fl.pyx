# Copyright (c) 2016-2023 The Regents of the University of Michigan
# Part of PGSD, released under the BSD 2-Clause License.

"""PGSD file layer API.

Low level access to pgsd files. :py:mod:`pgsd.fl` allows direct access to create,
read, and write ``pgsd`` files. The module is implemented in C and is optimized.
See :ref:`fl-examples` for detailed example code.

* :py:class:`PGSDFile` - Class interface to read and write pgsd files.
* :py:func:`open` - Open a pgsd file.

"""

import logging
import numpy
import os
from pickle import PickleError
import warnings
from libc.stdint cimport uint8_t, int8_t, uint16_t, int16_t, uint32_t, int32_t,uint64_t, int64_t
from libc.errno cimport errno
cimport pgsd.libpgsd as libpgsd
cimport numpy
from libcpp cimport bool

logger = logging.getLogger('pgsd.fl')

####################
# Helper functions #

cdef __format_errno(fname):
    """Return a tuple for constructing an IOError."""
    return (errno, os.strerror(errno), fname)

cdef __raise_on_error(retval, extra):
    """Raise the appropriate error type.

    Args:
        retval: Return value from a pgsd C API call
        extra: Extra string to pass along with the exception
    """
    if retval == libpgsd.PGSD_ERROR_IO:
        raise IOError(*__format_errno(extra))
    elif retval == libpgsd.PGSD_ERROR_NOT_A_PGSD_FILE:
        raise RuntimeError("Not a PGSD file: " + extra)
    elif retval == libpgsd.PGSD_ERROR_INVALID_PGSD_FILE_VERSION:
        raise RuntimeError("Unsupported PGSD file version: " + extra)
    elif retval == libpgsd.PGSD_ERROR_FILE_CORRUPT:
        raise RuntimeError("Corrupt PGSD file: " + extra)
    elif retval == libpgsd.PGSD_ERROR_MEMORY_ALLOCATION_FAILED:
        raise MemoryError("Memory allocation failed: " + extra)
    elif retval == libpgsd.PGSD_ERROR_NAMELIST_FULL:
        raise RuntimeError("PGSD namelist is full: " + extra)
    elif retval == libpgsd.PGSD_ERROR_FILE_MUST_BE_WRITABLE:
        raise RuntimeError("File must be writable: " + extra)
    elif retval == libpgsd.PGSD_ERROR_FILE_MUST_BE_READABLE:
        raise RuntimeError("File must be readable: " + extra)
    elif retval == libpgsd.PGSD_ERROR_INVALID_ARGUMENT:
        raise RuntimeError("Invalid pgsd argument: " + extra)
    elif retval != 0:
        raise RuntimeError("Unknown error: " + extra)

# Getter methods for 2D numpy arrays of all supported types
# cython needs strongly typed numpy arrays to get a void *
# to the data, so we implement each by hand here and dispacth
# from the read and write routines

cdef void * __get_ptr_uint8(data):
    cdef numpy.ndarray[uint8_t, ndim=2, mode="c"] data_array_uint8
    data_array_uint8 = data
    if (data.size == 0):
        return NULL
    else:
        return <void*>&data_array_uint8[0, 0]

cdef void * __get_ptr_uint16(data):
    cdef numpy.ndarray[uint16_t, ndim=2, mode="c"] data_array_uint16
    data_array_uint16 = data
    if (data.size == 0):
        return NULL
    else:
        return <void*>&data_array_uint16[0, 0]

cdef void * __get_ptr_uint32(data):
    cdef numpy.ndarray[uint32_t, ndim=2, mode="c"] data_array_uint32
    data_array_uint32 = data
    if (data.size == 0):
        return NULL
    else:
        return <void*>&data_array_uint32[0, 0]

cdef void * __get_ptr_uint64(data):
    cdef numpy.ndarray[uint64_t, ndim=2, mode="c"] data_array_uint64
    data_array_uint64 = data
    if (data.size == 0):
        return NULL
    else:
        return <void*>&data_array_uint64[0, 0]

cdef void * __get_ptr_int8(data):
    cdef numpy.ndarray[int8_t, ndim=2, mode="c"] data_array_int8
    data_array_int8 = data
    if (data.size == 0):
        return NULL
    else:
        return <void*>&data_array_int8[0, 0]

cdef void * __get_ptr_int16(data):
    cdef numpy.ndarray[int16_t, ndim=2, mode="c"] data_array_int16
    data_array_int16 = data
    if (data.size == 0):
        return NULL
    else:
        return <void*>&data_array_int16[0, 0]

cdef void * __get_ptr_int32(data):
    cdef numpy.ndarray[int32_t, ndim=2, mode="c"] data_array_int32
    data_array_int32 = data
    if (data.size == 0):
        return NULL
    else:
        return <void*>&data_array_int32[0, 0]

cdef void * __get_ptr_int64(data):
    cdef numpy.ndarray[int64_t, ndim=2, mode="c"] data_array_int64
    data_array_int64 = data
    if (data.size == 0):
        return NULL
    else:
        return <void*>&data_array_int64[0, 0]

cdef void * __get_ptr_float32(data):
    cdef numpy.ndarray[float, ndim=2, mode="c"] data_array_float32
    data_array_float32 = data
    if (data.size == 0):
        return NULL
    else:
        return <void*>&data_array_float32[0, 0]

cdef void * __get_ptr_float64(data):
    cdef numpy.ndarray[double, ndim=2, mode="c"] data_array_float64
    data_array_float64 = data
    if (data.size == 0):
        return NULL
    else:
        return <void*>&data_array_float64[0, 0]


def open(name, mode, application=None, schema=None, schema_version=None):
    """open(name, mode, application=None, schema=None, schema_version=None)

    :py:func:`open` opens a PGSD file and returns a :py:class:`PGSDFile` instance.
    The return value of :py:func:`open` can be used as a context manager.

    Args:
        name (str): File name to open.

        mode (str): File access mode.

        application (str): Name of the application creating the file.

        schema (str): Name of the data schema.

        schema_version (tuple[int, int]): Schema version number
            (major, minor).

    Valid values for mode:

Valid values for ``mode``:

    +------------------+---------------------------------------------+
    | mode             | description                                 |
    +==================+=============================================+
    | ``'r'``          | Open an existing file for reading.          |
    +------------------+---------------------------------------------+
    | ``'r+'``         | Open an existing file for reading and       |
    |                  | writing.                                    |
    +------------------+---------------------------------------------+
    | ``'w'``          | Open a file for reading and writing.        |
    |                  | Creates the file if needed, or overwrites   |
    |                  | an existing file.                           |
    +------------------+---------------------------------------------+
    | ``'x'``          | Create a gsd file exclusively and opens it  |
    |                  | for reading and writing.                    |
    |                  | Raise :py:exc:`FileExistsError`             |
    |                  | if it already exists.                       |
    +------------------+---------------------------------------------+
    | ``'a'``          | Open a file for reading and writing.        |
    |                  | Creates the file if it doesn't exist.       |
    +------------------+---------------------------------------------+

    When opening a file for reading (``'r'`` and ``'r+'`` modes):
    ``application`` and ``schema_version`` are ignored and may be ``None``.
    When ``schema`` is not ``None``, :py:func:`open` throws an exception if the
    file's schema does not match ``schema``.

    When opening a file for writing (``'w'``, ``'x'``, or ``'a'`` modes): The
    given ``application``, ``schema``, and ``schema_version`` must not be None.

    Example:

        .. ipython:: python

            with pgsd.fl.open(name='file.pgsd', mode='w',
                             application="My application", schema="My Schema",
                             schema_version=[1,0]) as f:
                f.write_chunk(name='chunk1',
                              data=numpy.array([1,2,3,4], dtype=numpy.float32))
                f.write_chunk(name='chunk2',
                              data=numpy.array([[5,6],[7,8]],
                                               dtype=numpy.float32))
                f.end_frame()
                f.write_chunk(name='chunk1',
                              data=numpy.array([9,10,11,12],
                                               dtype=numpy.float32))
                f.write_chunk(name='chunk2',
                              data=numpy.array([[13,14],[15,16]],
                                               dtype=numpy.float32))
                f.end_frame()

            f = pgsd.fl.open(name='file.pgsd', mode='r')
            if f.chunk_exists(frame=0, name='chunk1'):
                data = f.read_chunk(frame=0, name='chunk1')
            data
            f.close()
    """

    return PGSDFile(str(name), mode, application, schema, schema_version)


cdef class PGSDFile:
    """PGSDFile

    PGSD file access interface.

    Args:

        name (str): Name of the open file.

        mode (str): Mode of the open file.

        pgsd_version (tuple[int, int]): PGSD file layer version number
            (major, minor).

        application (str): Name of the generating application.

        schema (str): Name of the data schema.

        schema_version (tuple[int, int]): Schema version number
            (major, minor).

        nframes (int): Number of frames.

    :py:class:`PGSDFile` implements an object oriented class interface to the PGSD
    file layer. Use :py:func:`open` to open a PGSD file and obtain a
    :py:class:`PGSDFile` instance. :py:class:`PGSDFile` can be used as a context
    manager.

    Attributes:

        name (str): Name of the open file.

        mode (str): Mode of the open file.

        pgsd_version (tuple[int, int]): PGSD file layer version number
            (major, minor).

        application (str): Name of the generating application.

        schema (str): Name of the data schema.

        schema_version (tuple[int, int]): Schema version number
            (major, minor).

        nframes (int): Number of frames.


        maximum_write_buffer_size (int): The Maximum write buffer size (bytes).

        index_entries_to_buffer (int): Number of index entries to buffer before
            flushing.
    """

    cdef libpgsd.pgsd_handle __handle
    cdef bint __is_open
    cdef str mode
    cdef str name

    def __init__(self,
                 name,
                 mode,
                 application,
                 schema,
                 schema_version):
        cdef libpgsd.pgsd_open_flag c_flags
        cdef int exclusive_create = 0
        cdef int overwrite = 0

        self.mode = mode

        if mode == 'w':
            c_flags = libpgsd.PGSD_OPEN_READWRITE
            overwrite = 1
        elif mode == 'r':
            c_flags = libpgsd.PGSD_OPEN_READONLY
        elif mode == 'r+':
            c_flags = libpgsd.PGSD_OPEN_READWRITE
        elif mode == 'x':
            c_flags = libpgsd.PGSD_OPEN_READWRITE
            overwrite = 1
            exclusive_create = 1
        elif mode == 'a':
            c_flags = libpgsd.PGSD_OPEN_READWRITE
            if not os.path.exists(name):
                overwrite = 1
        else:
            raise ValueError("Invalid mode: " + mode)

        self.name = name

        cdef char * c_name
        cdef char * c_application
        cdef char * c_schema
        cdef int _c_schema_version
        cdef str schema_truncated

        if overwrite:
            if application is None:
                raise ValueError("Provide application when creating a file")
            if schema is None:
                raise ValueError("Provide schema when creating a file")
            if schema_version is None:
                raise ValueError("Provide schema_version when creating a file")

            # create a new file or overwrite an existing one
            logger.info('overwriting file: ' + name + ' with mode: ' + mode
                        + ', application: ' + application
                        + ', schema: ' + schema
                        + ', and schema_version: ' + str(schema_version))
            name_e = name.encode('utf-8')
            c_name = name_e

            application_e = application.encode('utf-8')
            c_application = application_e

            schema_e = schema.encode('utf-8')
            c_schema = schema_e

            _c_schema_version = libpgsd.pgsd_make_version(schema_version[0],
                                                        schema_version[1])

            with nogil:
                retval = libpgsd.pgsd_create_and_open(&self.__handle,
                                                    c_name,
                                                    c_application,
                                                    c_schema,
                                                    _c_schema_version,
                                                    c_flags,
                                                    exclusive_create)
        else:
            # open an existing file
            logger.info('opening file: ' + name + ' with mode: ' + mode)
            name_e = name.encode('utf-8')
            c_name = name_e

            with nogil:
                retval = libpgsd.pgsd_open(&self.__handle, c_name, c_flags)

        __raise_on_error(retval, name)

        # validate schema
        if schema is not None:
            schema_truncated = schema
            if len(schema_truncated) > 64:
                schema_truncated = schema_truncated[0:63]
            if self.schema != schema_truncated:
                raise RuntimeError('file ' + name + ' has incorrect schema: '
                                   + self.schema)

        self.__is_open = True

    def close(self, write_all):
        """close()

        Close the file.

        Once closed, any other operation on the file object will result in a
        `ValueError`. :py:meth:`close()` may be called more than once.
        The file is automatically closed when garbage collected or when
        the context manager exits.

        Example:
            .. ipython:: python
               :okexcept:

                f = pgsd.fl.open(name='file.pgsd', mode='w',
                                application="My application",
                                schema="My Schema", schema_version=[1,0])
                f.write_chunk(name='chunk1',
                              data=numpy.array([1,2,3,4], dtype=numpy.float32))
                f.end_frame()
                data = f.read_chunk(frame=0, name='chunk1')

                f.close()
                # Read fails because the file is closed
                data = f.read_chunk(frame=0, name='chunk1')

        """

        cdef bool wr_all;
        wr_all = write_all

        if self.__is_open:
            logger.info('closing file: ' + self.name)
            with nogil:
                retval = libpgsd.pgsd_close(&self.__handle)
            self.__is_open = False

            __raise_on_error(retval, self.name)

    # def truncate(self):
    #     """truncate()

    #     Truncate all data from the file. After truncation, the file has no
    #     frames and no data chunks. The application, schema, and schema version
    #     remain the same.

    #     Example:
    #         .. ipython:: python

    #             with pgsd.fl.open(name='file.pgsd', mode='w',
    #                              application="My application",
    #                              schema="My Schema", schema_version=[1,0]) as f:
    #                 for i in range(10):
    #                     f.write_chunk(name='chunk1',
    #                                   data=numpy.array([1,2,3,4],
    #                                                    dtype=numpy.float32))
    #                     f.end_frame()

    #             f = pgsd.fl.open(name='file.pgsd', mode='r+',
    #                             application="My application",
    #                             schema="My Schema", schema_version=[1,0])
    #             f.nframes
    #             f.schema, f.schema_version, f.application
    #             f.truncate()
    #             f.nframes
    #             f.schema, f.schema_version, f.application
    #             f.close()
    #     """

    #     if not self.__is_open:
    #         raise ValueError("File is not open")

    #     logger.info('truncating file: ' + self.name)
    #     with nogil:
    #         retval = libpgsd.pgsd_truncate(&self.__handle)

    #     __raise_on_error(retval, self.name)

    def end_frame(self, write_all):
        """end_frame()

        Complete writing the current frame. After calling :py:meth:`end_frame()`
        future calls to :py:meth:`write_chunk()` will write to the **next**
        frame in the file.

        .. danger::
            Call :py:meth:`end_frame()` to complete the current frame
            **before** closing the file. If you fail to call
            :py:meth:`end_frame()`, the last frame will not be written
            to disk.

        Example:
            .. ipython:: python

                f = pgsd.fl.open(name='file.pgsd', mode='w',
                                application="My application",
                                schema="My Schema", schema_version=[1,0])

                f.write_chunk(name='chunk1',
                              data=numpy.array([1,2,3,4], dtype=numpy.float32))
                f.end_frame()
                f.write_chunk(name='chunk1',
                              data=numpy.array([9,10,11,12],
                                               dtype=numpy.float32))
                f.end_frame()
                f.write_chunk(name='chunk1',
                              data=numpy.array([13,14],
                                               dtype=numpy.float32))
                f.end_frame()
                f.nframes
                f.close()

        """
        cdef bool wr_all;
        wr_all = write_all;

        if not self.__is_open:
            raise ValueError("File is not open")

        logger.debug('end frame: ' + self.name)

        with nogil:
            retval = libpgsd.pgsd_end_frame(&self.__handle)

        __raise_on_error(retval, self.name)

    def flush(self, write_all):
        """flush()

        Flush all buffered frames to the file.
        """
        cdef bool wr_all;
        wr_all = write_all;

        if not self.__is_open:
            raise ValueError("File is not open")

        logger.debug('flush: ' + self.name)

        with nogil:
            retval = libpgsd.pgsd_flush(&self.__handle)

        __raise_on_error(retval, self.name)

    def write_chunk(self, name, offset, rank, write_all, data):
        """write_chunk(name, data)

        Write a data chunk to the file. After writing all chunks in the
        current frame, call :py:meth:`end_frame()`.

        Args:
            name (str): Name of the chunk
            data: Data to write into the chunk. Must be a numpy
                  array, or array-like, with 2 or fewer
                  dimensions.

        Warning:
            :py:meth:`write_chunk()` will implicitly converts array-like and
            non-contiguous numpy arrays to contiguous numpy arrays with
            ``numpy.ascontiguousarray(data)``. This may or may not produce
            desired data types in the output file and incurs overhead.

        Example:
            .. ipython:: python

                f = pgsd.fl.open(name='file.pgsd', mode='w',
                                application="My application",
                                schema="My Schema", schema_version=[1,0])

                f.write_chunk(name='float1d',
                              data=numpy.array([1,2,3,4],
                                               dtype=numpy.float32))
                f.write_chunk(name='float2d',
                              data=numpy.array([[13,14],[15,16],[17,19]],
                                               dtype=numpy.float32))
                f.write_chunk(name='double2d',
                              data=numpy.array([[1,4],[5,6],[7,9]],
                                               dtype=numpy.float64))
                f.write_chunk(name='int1d',
                              data=numpy.array([70,80,90],
                                               dtype=numpy.int64))
                f.end_frame()
                f.nframes
                f.close()
        """

        if not self.__is_open:
            raise ValueError("File is not open")

        data_array = numpy.ascontiguousarray(data)
        if data_array is not data:
            logger.warning('implicit data copy when writing chunk: ' + name)
        data_array = data_array.view()

        cdef uint64_t N
        cdef uint32_t M

        cdef uint64_t N_global; # added DK 
        cdef uint64_t stride;
        cdef bool wr_all;

        if len(data_array.shape) > 2:
            raise ValueError("PGSD can only write 1 or 2 dimensional arrays: "
                             + name)

        if len(data_array.shape) == 1:
            data_array = data_array.reshape([data_array.shape[0], 1])

        N = data_array.shape[0]
        M = data_array.shape[1]

        wr_all = write_all; # added DK 
        N_global = N;
        stride = 0;
        if offset is not None:
            N_global = offset.sum();
            stride = M * offset[0:rank].sum();

        cdef libpgsd.pgsd_type pgsd_type
        cdef void *data_ptr
        if data_array.dtype == numpy.uint8:
            pgsd_type = libpgsd.PGSD_TYPE_UINT8
            data_ptr = __get_ptr_uint8(data_array)
        elif data_array.dtype == numpy.uint16:
            pgsd_type = libpgsd.PGSD_TYPE_UINT16
            data_ptr = __get_ptr_uint16(data_array)
        elif data_array.dtype == numpy.uint32:
            pgsd_type = libpgsd.PGSD_TYPE_UINT32
            data_ptr = __get_ptr_uint32(data_array)
        elif data_array.dtype == numpy.uint64:
            pgsd_type = libpgsd.PGSD_TYPE_UINT64
            data_ptr = __get_ptr_uint64(data_array)
        elif data_array.dtype == numpy.int8:
            pgsd_type = libpgsd.PGSD_TYPE_INT8
            data_ptr = __get_ptr_int8(data_array)
        elif data_array.dtype == numpy.int16:
            pgsd_type = libpgsd.PGSD_TYPE_INT16
            data_ptr = __get_ptr_int16(data_array)
        elif data_array.dtype == numpy.int32:
            pgsd_type = libpgsd.PGSD_TYPE_INT32
            data_ptr = __get_ptr_int32(data_array)
        elif data_array.dtype == numpy.int64:
            pgsd_type = libpgsd.PGSD_TYPE_INT64
            data_ptr = __get_ptr_int64(data_array)
        elif data_array.dtype == numpy.float32:
            pgsd_type = libpgsd.PGSD_TYPE_FLOAT
            data_ptr = __get_ptr_float32(data_array)
        elif data_array.dtype == numpy.float64:
            pgsd_type = libpgsd.PGSD_TYPE_DOUBLE
            data_ptr = __get_ptr_float64(data_array)
        else:
            raise ValueError("invalid type for chunk: " + name)

        logger.debug('write chunk: ' + self.name + ' - ' + name)

        cdef char * c_name
        name_e = name.encode('utf-8')
        c_name = name_e
        with nogil:
            retval = libpgsd.pgsd_write_chunk(&self.__handle,
                                            c_name,
                                            pgsd_type,
                                            N,
                                            M,
                                            N_global,
                                            M,
                                            stride,
                                            N_global*M,
                                            wr_all,
                                            0,
                                            data_ptr)

        __raise_on_error(retval, self.name)

    def chunk_exists(self, frame, name, write_all):
        """chunk_exists(frame, name)

        Test if a chunk exists.

        Args:
            frame (int): Index of the frame to check
            name (str): Name of the chunk

        Returns:
            bool: ``True`` if the chunk exists in the file at the given frame.\
              ``False`` if it does not.

        Example:
            .. ipython:: python

                with pgsd.fl.open(name='file.pgsd', mode='w',
                                 application="My application",
                                 schema="My Schema", schema_version=[1,0]) as f:
                    f.write_chunk(name='chunk1',
                                  data=numpy.array([1,2,3,4],
                                                   dtype=numpy.float32))
                    f.write_chunk(name='chunk2',
                                  data=numpy.array([[5,6],[7,8]],
                                                   dtype=numpy.float32))
                    f.end_frame()
                    f.write_chunk(name='chunk1',
                                  data=numpy.array([9,10,11,12],
                                                   dtype=numpy.float32))
                    f.write_chunk(name='chunk2',
                                  data=numpy.array([[13,14],[15,16]],
                                                   dtype=numpy.float32))
                    f.end_frame()

                f = pgsd.fl.open(name='file.pgsd', mode='r',
                                application="My application",
                                schema="My Schema", schema_version=[1,0])

                f.chunk_exists(frame=0, name='chunk1')
                f.chunk_exists(frame=0, name='chunk2')
                f.chunk_exists(frame=0, name='chunk3')
                f.chunk_exists(frame=10, name='chunk1')
                f.close()
        """
        cdef bool wr_all;
        wr_all = write_all;

        cdef const libpgsd.pgsd_index_entry* index_entry
        cdef char * c_name
        name_e = name.encode('utf-8')
        c_name = name_e
        cdef int64_t c_frame
        c_frame = frame

        logger.debug('chunk exists: ' + self.name + ' - ' + name)

        with nogil:
            index_entry = libpgsd.pgsd_find_chunk(&self.__handle, c_frame, c_name)

        return index_entry != NULL

    def read_chunk(self, frame, name, uint32_t[:] offset, r_all):
        """read_chunk(frame, name)

        Read a data chunk from the file and return it as a numpy array.

        Args:
            frame (int): Index of the frame to read
            name (str): Name of the chunk

        Returns:
            ``(N,M)`` or ``(N,)`` `numpy.ndarray` of ``type``: Data read from
            file. ``N``, ``M``, and ``type`` are determined by the chunk
            metadata. If the data is NxM in the file and M > 1, return a 2D
            array. If the data is Nx1, return a 1D array.

        .. tip::
            Each call invokes a disk read and allocation of a
            new numpy array for storage. To avoid overhead, call
            :py:meth:`read_chunk()` on the same chunk only once.

        Example:
            .. ipython:: python
               :okexcept:

                with pgsd.fl.open(name='file.pgsd', mode='w',
                                 application="My application",
                                 schema="My Schema", schema_version=[1,0]) as f:
                    f.write_chunk(name='chunk1',
                                  data=numpy.array([1,2,3,4],
                                                   dtype=numpy.float32))
                    f.write_chunk(name='chunk2',
                                  data=numpy.array([[5,6],[7,8]],
                                                   dtype=numpy.float32))
                    f.end_frame()
                    f.write_chunk(name='chunk1',
                                  data=numpy.array([9,10,11,12],
                                                   dtype=numpy.float32))
                    f.write_chunk(name='chunk2',
                                  data=numpy.array([[13,14],[15,16]],
                                                   dtype=numpy.float32))
                    f.end_frame()

                f = pgsd.fl.open(name='file.pgsd', mode='r',
                                application="My application",
                                schema="My Schema", schema_version=[1,0])
                f.read_chunk(frame=0, name='chunk1')
                f.read_chunk(frame=1, name='chunk1')
                f.read_chunk(frame=2, name='chunk1')
                f.close()
        """
        # raise NotImplementedError;

        if not self.__is_open:
            raise ValueError("File is not open")

        cdef const libpgsd.pgsd_index_entry* index_entry
        cdef char * c_name
        name_e = name.encode('utf-8')
        print("in file fl.pyx name: " +str(name_e))
        c_name = name_e
        cdef int64_t c_frame
        c_frame = frame
        # cdef uint32_t * c_offset;
        # c_offset  = offset;
        cdef bool wr_all
        wr_all = r_all

        with nogil:
            index_entry = libpgsd.pgsd_find_chunk(&self.__handle, c_frame, c_name)

        if index_entry == NULL:
            raise KeyError("frame " + str(frame) + " / chunk " + name
                           + " not found in: " + self.name)

        cdef libpgsd.pgsd_type pgsd_type
        pgsd_type = <libpgsd.pgsd_type>index_entry.type
        print("in file fl.pyx pgsd_type: " +str(pgsd_type))


        cdef void *data_ptr
        if pgsd_type == libpgsd.PGSD_TYPE_UINT8:
            data_array = numpy.empty(dtype=numpy.uint8,
                                     shape=[index_entry.N, index_entry.M])
        elif pgsd_type == libpgsd.PGSD_TYPE_UINT16:
            data_array = numpy.empty(dtype=numpy.uint16,
                                     shape=[index_entry.N, index_entry.M])
        elif pgsd_type == libpgsd.PGSD_TYPE_UINT32:
            data_array = numpy.empty(dtype=numpy.uint32,
                                     shape=[index_entry.N, index_entry.M])
        elif pgsd_type == libpgsd.PGSD_TYPE_UINT64:
            data_array = numpy.empty(dtype=numpy.uint64,
                                     shape=[index_entry.N, index_entry.M])
        elif pgsd_type == libpgsd.PGSD_TYPE_INT8:
            data_array = numpy.empty(dtype=numpy.int8,
                                     shape=[index_entry.N, index_entry.M])
        elif pgsd_type == libpgsd.PGSD_TYPE_INT16:
            data_array = numpy.empty(dtype=numpy.int16,
                                     shape=[index_entry.N, index_entry.M])
        elif pgsd_type == libpgsd.PGSD_TYPE_INT32:
            data_array = numpy.empty(dtype=numpy.int32,
                                     shape=[index_entry.N, index_entry.M])
        elif pgsd_type == libpgsd.PGSD_TYPE_INT64:
            data_array = numpy.empty(dtype=numpy.int64,
                                     shape=[index_entry.N, index_entry.M])
        elif pgsd_type == libpgsd.PGSD_TYPE_FLOAT:
            data_array = numpy.empty(dtype=numpy.float32,
                                     shape=[index_entry.N, index_entry.M])
        elif pgsd_type == libpgsd.PGSD_TYPE_DOUBLE:
            data_array = numpy.empty(dtype=numpy.float64,
                                     shape=[index_entry.N, index_entry.M])
        else:
            raise ValueError("invalid type for chunk: " + name)

        logger.debug('read chunk: ' + self.name + ' - '
                     + str(frame) + ' - ' + name)

        # only read chunk if we have data
        if index_entry.N != 0 and index_entry.M != 0:
            if pgsd_type == libpgsd.PGSD_TYPE_UINT8:
                data_ptr = __get_ptr_uint8(data_array)
            elif pgsd_type == libpgsd.PGSD_TYPE_UINT16:
                data_ptr = __get_ptr_uint16(data_array)
            elif pgsd_type == libpgsd.PGSD_TYPE_UINT32:
                data_ptr = __get_ptr_uint32(data_array)
            elif pgsd_type == libpgsd.PGSD_TYPE_UINT64:
                data_ptr = __get_ptr_uint64(data_array)
            elif pgsd_type == libpgsd.PGSD_TYPE_INT8:
                data_ptr = __get_ptr_int8(data_array)
            elif pgsd_type == libpgsd.PGSD_TYPE_INT16:
                data_ptr = __get_ptr_int16(data_array)
            elif pgsd_type == libpgsd.PGSD_TYPE_INT32:
                data_ptr = __get_ptr_int32(data_array)
            elif pgsd_type == libpgsd.PGSD_TYPE_INT64:
                data_ptr = __get_ptr_int64(data_array)
            elif pgsd_type == libpgsd.PGSD_TYPE_FLOAT:
                data_ptr = __get_ptr_float32(data_array)
            elif pgsd_type == libpgsd.PGSD_TYPE_DOUBLE:
                data_ptr = __get_ptr_float64(data_array)
            else:
                raise ValueError("invalid type for chunk: " + name)

            with nogil:
                retval = libpgsd.pgsd_read_chunk(&self.__handle,
                                               data_ptr,
                                               index_entry,
                                               &offset[0],
                                               wr_all)

            __raise_on_error(retval, self.name)

        if index_entry.M == 1:
            return data_array.reshape([index_entry.N])
        else:
            return data_array

    def find_matching_chunk_names(self, match, write_all):
        """find_matching_chunk_names(match)

        Find all the chunk names in the file that start with the string *match*.

        Args:
            match (str): Start of the chunk name to match

        Returns:
            list[str]: Matching chunk names

        Example:
            .. ipython:: python

                with pgsd.fl.open(name='file.pgsd', mode='w',
                                 application="My application",
                                 schema="My Schema", schema_version=[1,0]) as f:
                    f.write_chunk(name='data/chunk1',
                                  data=numpy.array([1,2,3,4],
                                                   dtype=numpy.float32))
                    f.write_chunk(name='data/chunk2',
                                  data=numpy.array([[5,6],[7,8]],
                                                   dtype=numpy.float32))
                    f.write_chunk(name='input/chunk3',
                                  data=numpy.array([9, 10],
                                                   dtype=numpy.float32))
                    f.end_frame()
                    f.write_chunk(name='input/chunk4',
                                  data=numpy.array([11, 12, 13, 14],
                                                   dtype=numpy.float32))
                    f.end_frame()

                f = pgsd.fl.open(name='file.pgsd', mode='r',
                                application="My application",
                                schema="My Schema", schema_version=[1,0])

                f.find_matching_chunk_names('')
                f.find_matching_chunk_names('data')
                f.find_matching_chunk_names('input')
                f.find_matching_chunk_names('other')
                f.close()
        """

        cdef bool wr_all;
        wr_all = write_all;

        if not self.__is_open:
            raise ValueError("File is not open")

        cdef const char * c_found
        cdef char * c_match
        match_e = match.encode('utf-8')
        c_match = match_e

        retval = []

        with nogil:
            c_found = libpgsd.pgsd_find_matching_chunk_name(&self.__handle,
                                                          c_match,
                                                          NULL)

        while c_found != NULL:
            retval.append(c_found.decode('utf-8'))

            with nogil:
                c_found = libpgsd.pgsd_find_matching_chunk_name(&self.__handle,
                                                              c_match,
                                                              c_found)

        return retval

    # def upgrade(self):
    #     """upgrade()

    #     Upgrade a GSD file to the v2 specification in place. The file must be
    #     open in a writable mode.

    #     """

    #     if not self.__is_open:
    #         raise ValueError("File is not open")

    #     logger.info('upgrading file: ' + self.name)

    #     with nogil:
    #         retval = libpgsd.pgsd_upgrade(&self.__handle)

    #     __raise_on_error(retval, self.name)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def __reduce__(self):
        """Allows filehandles to be pickled when in read only mode."""
        if self.mode not in ['rb', 'r']:
            raise PickleError("Only read only GSDFiles can be pickled.")
        return (PGSDFile,
                (self.name, self.mode, self.application,
                    self.schema, self.schema_version),
                )

    property name:
        def __get__(self):
            return self.name

    property mode:
        def __get__(self):
            return self.mode

    property pgsd_version:
        def __get__(self):
            cdef uint32_t v = self.__handle.header.pgsd_version
            return (v >> 16, v & 0xffff)

    property schema_version:
        def __get__(self):
            cdef uint32_t v = self.__handle.header.schema_version
            return (v >> 16, v & 0xffff)

    property schema:
        def __get__(self):
            return self.__handle.header.schema.decode('utf-8')

    property application:
        def __get__(self):
            return self.__handle.header.application.decode('utf-8')

    property nframes:
        def __get__(self):
            if not self.__is_open:
                raise ValueError("File is not open")

            return libpgsd.pgsd_get_nframes(&self.__handle)

    property maximum_write_buffer_size:
        def __get__(self):
            if not self.__is_open:
                raise ValueError("File is not open")

            return libpgsd.pgsd_get_maximum_write_buffer_size(&self.__handle)

        def __set__(self, size):
            if not self.__is_open:
                raise ValueError("File is not open")

            retval = libpgsd.pgsd_set_maximum_write_buffer_size(&self.__handle, size)
            __raise_on_error(retval, self.name)

    property index_entries_to_buffer:
        def __get__(self):
            if not self.__is_open:
                raise ValueError("File is not open")

            return libpgsd.pgsd_get_index_entries_to_buffer(&self.__handle)

        def __set__(self, number):
            if not self.__is_open:
                raise ValueError("File is not open")

            retval = libpgsd.pgsd_set_index_entries_to_buffer(&self.__handle, number)
            __raise_on_error(retval, self.name)

    def __dealloc__(self):
        if self.__is_open:
            logger.info('closing file: ' + self.name)
            libpgsd.pgsd_close(&self.__handle)
            self.__is_open = False
