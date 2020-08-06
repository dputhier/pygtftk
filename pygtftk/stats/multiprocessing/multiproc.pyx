
"""
Using Cython and OpenMP to multiprocess any arbitrary function over NumPy arrays.
As of January 2020, the functions here are not yet used in a meaningful way in the code.

OLOGRAM implements some examples of use, but in parts of the code where it does not bring much meaningful improvement.
Consider these as examples of what multiprocessing can do, it you need it later.

    Author : Quentin Ferré <quentin.ferre@gmail.com>

Based on : https://cython.readthedocs.io/en/latest/src/userguide/parallelism.html


IMPLEMENTATION NOTES :

The main function of this module is apply_func_multiproc_cython() ; it applies a
user-defined nogil function over a list of NumPy arrays. The list will first be
converted to a list of C arrays.
The function to be applied must have a very specific signature (see comments of
apply_func_multiproc_cython() for details) ; in broad strokes, it must return a
result object containing a single 2D result array.

These functions are designed for 2D arrays because this is what I use here.

I use long long all over the place for compatibility reasons. The arrays will
be long longs and that will not change, but the shapes could have been int.
However there will be few shapes so that's okay.
"""

cimport cython
cimport openmp
import numpy as np
cimport numpy as np
cimport libc.stdio as stdio
from libc.stdlib cimport abort, malloc, free
from cython.parallel import parallel, prange, threadid

import time

# Import structures
from pygtftk.stats.multiprocessing.multiproc_structs cimport *


cdef extern from "sched.h":
    int sched_getcpu() nogil

np.import_array() # Init numpy ?


# To set memory ownership flags
cdef extern from "numpy/arrayobject.h":
    void PyArray_ENABLEFLAGS(np.ndarray arr, int flags)




# ------------------- Converting NumPy arrays to C objects ------------------- #
# Takes a list of numpy arrays, convert it into two objects : a list of C
# pointers to the same arrays, and a list of array shapes because C pointers do
# not know the shape.

# This "list_of_arrays" we created in the previous step is a C object, meaning
# it can be used without the GIL.

cdef long long  ** list_of_numpy_arrays_to_c_object(list list_of_arrays_Py):

    # Allocate a look-up-table of pointers.
    numPointers = len(list_of_arrays_Py)
    cdef long long ** list_of_arrays = <long long **> malloc(sizeof(long long *) * numPointers)

    # Rq : this is for 2D arrays. To make it suitable for n-dimensional arrays,
    # write ndim = n AND write &temp[0,0,...] with n zeros.

    # Now replace the pointers with the points of the NumPy arrays
    cdef np.ndarray[long long, ndim=2, mode="c"] temp
    for i in range(numPointers):
        temp = list_of_arrays_Py[i]
        list_of_arrays[i] = &temp[0,0]

    return list_of_arrays




# ------------------------------ Multiprocessing ----------------------------- #

cdef list apply_func_multiproc_cython(list python_list_of_numpy_arrays, list python_list_of_parameters, MULTIPROC_CUSTOM_FUNCTYPE cython_func, int num_threads):
    """
    Uses OpenMP multiprocessing to apply a Cython function a a list of NumPy arrays

    The function must have the following signature :
        cdef FUNC_2DARRAY_RESULT func(long long * array, long long * array_shape, long long parameter) nogil

        It takes as input a numpy array, and a long long parameter. You can encode anything you want in the parameter.

        where :
        cdef struct FUNC_2DARRAY_RESULT:
            long long * result_array                  # Pointer to a C array as flattned 1D array (!)
            long long[2] result_shape                 # x, y respectively
        See comments in code for more info
        The result_array is 1D, so to get result_array[x,y] write result_array[x*result_shape[1]+y]

    And it must only contain Cython code, no Python, so that the GIL may be released.
    """


    # Convert the Pyton list of numpy array into a C list of C arrays (ie. pointers to arrays)
    cdef long long ** list_of_arrays = list_of_numpy_arrays_to_c_object(python_list_of_numpy_arrays)
    # We must also remember the shapes ! Artificially make them 2D arrays for compatibility
    shapes = [np.array(x.shape).reshape((2,1)) for x in python_list_of_numpy_arrays]
    cdef long long ** list_of_shapes = list_of_numpy_arrays_to_c_object(shapes)

    # Further convert the list of parameters to a 1D C array
    cdef long long  * list_of_parameters
    list_of_parameters = <long long *>malloc(len(python_list_of_parameters)*sizeof(long long))
    if list_of_parameters is NULL: raise MemoryError()
    for pk in xrange(len(python_list_of_parameters)):
        list_of_parameters[pk] = python_list_of_parameters[pk]


    # Calling code sequences (function) in parallel over a "list" of objects.
    # Each object we work on will be stocked as one element of a buffer.
    cdef int n_threads_cy = num_threads

    # Remember the number of elements in the lists, as the 'n' variable
    cdef Py_ssize_t i, n = len(python_list_of_numpy_arrays)


    # Allocate a list of structures for the results
    cdef FUNC_2DARRAY_RESULT* results = <FUNC_2DARRAY_RESULT *> malloc(sizeof(FUNC_2DARRAY_RESULT) * n)


    with nogil, parallel(num_threads = n_threads_cy):



        # The "local buffer" are the lists we created above in a sequential loop
        # It is designed to be used as a list of thread-local buffers to share
        # the work using prange

        # TODO : which schedule to use ?
        # Likely dynamic as the tasks do not always have the same size (I often multiproc by chromosome)

        for i in prange(n, schedule='dynamic'):

            # Debug print
            # TODO RE-ADD THIS !!
            #stdio.printf("tid: %d   cpuid: %d\n", openmp.omp_get_thread_num(), sched_getcpu())



            results[i] = cython_func(list_of_arrays[i], list_of_shapes[i], list_of_parameters[i])
            # NOTE : Must use `ptr[0]` instead of `*ptr` in Cython ! Also as we see here above, 'i' can be used to iterate over
            # multiple lists of arrays, as created above.




    # Convert results from a list of pointers to a list of NumPyarrays
    cdef results_return = list()
    cdef np.npy_intp[2] myshape
    cdef long long * arr_point
    cdef np.ndarray[np.int64_t, ndim=2] arr_adding

    for r in range(n):
        myshape = <np.npy_intp*> results[r].result_shape
        # Read as 1D array, then reshape
        arr_point = results[r].result_array
        arr_adding = np.PyArray_SimpleNewFromData(2, myshape, np.NPY_LONGLONG, arr_point)

        # Tell Python that it can deallocate the memory when the object gets garbage collected
        
        # np.PyArray_UpdateFlags(arr_adding, arr_adding.flags.num | np.NPY_OWNDATA) # NO ! THIS DOES NOT WORK, AND DOES NOT SEND AN ERROR MESSAGE
        
        PyArray_ENABLEFLAGS(arr_adding, np.NPY_OWNDATA)
    
        results_return += [arr_adding]

    # Finally, return and free memory
    try:
        return results_return

    finally:
        free(list_of_arrays)
        free(list_of_shapes)
        free(list_of_parameters)
        free(results)
        # NOTE Do not free the pointers ! They have been recast as NumPy arrays ! 









# -------------------------------- Testing ----------------------------------- #

# Example of function that operates on a 2D array. It doubles the values.
# Remember to declare it as nogil and specify its return type!
# Also be careful to not use a single Python command in it.
cdef FUNC_2DARRAY_RESULT example_func(long long* arr, long long* shape, long long parameter) nogil:

    # Create a result array of the same shape
    cdef size_t k
    cdef size_t xmax = shape[0]
    cdef size_t ymax = shape[1]

    # Create a result object
    cdef FUNC_2DARRAY_RESULT result
    result.result_shape[0] = xmax
    result.result_shape[1] = ymax
    result.result_array = <long long *> malloc(xmax*ymax*sizeof(long long));

    # Operation

    # WARNING : 2d arrays are flattened as 1d when put into C arrays.
    # When making a query, stride must be taken into account !
    # And also, unlike in 1d arrays, you must declare the iterators as size_t
    # to prevent it from falling back into Python function
    cdef size_t i
    cdef size_t j


    for i in range(xmax):
        for j in range(ymax):
            result.result_array[(i*ymax)+j] = parameter * arr[(i*ymax)+j]

    # Control : very CPU intensive but useless operation, just to check that the
    # threads are being properly used.
    cdef int n = 1
    cdef double tmp = 0.0
    while n > 0:
        tmp = (n ** 0.5) ** 0.5
        n -= 1


    return result
    # TODO and what about the shape of the result ?


cdef printcustom(FUNC_2DARRAY_RESULT result):
    cdef long long[2] shape = result.result_shape
    cdef long long * arr = result.result_array
    for i in range(shape[0]):
        for j in range(shape[1]):
            print(arr[i*shape[1]+j])










"""
# TODO Make this a proper unitary test
data = list()
for d in range(10):
    if d % 2 == 0 : data += [np.array([[d,d,d,d],[d,d,d,d]])]
    else: data += [np.array([[d,d,d],[d,d,d]])]
parameters = [2] * len(data)

start = time.time()
r = apply_func_multiproc_cython(data, parameters, example_func, num_threads = 2)
stop = time.time()
print('Took '+str(stop-start)+' s.')
print(r)
"""






"""
# ------- Multiprocessing type 2

# This is OpenCL-based multiprocessing.
# This answers a different question : above we compute purely independant operations
# on several threads (MIMD parallelism). This part here is dedicated to SIMD parallelism.

# Implementation note : there are a few matrix calculations steps in OLOGRAM that could benefit from this.
# I will leave it here as commented code in case such an opportunity comes.

import pyopencl as cl
import numpy

A = numpy.random.rand(1000).astype(numpy.float32)
B = numpy.random.rand(1000).astype(numpy.float32)
C = numpy.empty_like(A)

cl.get_plaftorms() # Check available OpenCL devices

# Create a context and a job queue. The context is a device. Pased when creating buffers.
ctx = cl.Context()
queue = cl.CommandQueue(ctx)

## Read kernel
# The Kernel is an OpenCL program that will be run. It must be compiled first.
# __global means we draw from the global context memory
# OpenCL syntax is similar to C
kernel_as_string = '''
  __kernel
  void sum(__global const float* a, __global const float* b, __global float* c){
    int i = get_global_id(0);
    c[i] = a[i] + b[i];
  }
  '''
# 'i' here will give the id of the workgroup the kernel is in (equivalent in CUDA would be id of Block)
# get_local_id() would give the id within the workgroup. You can use __shared in kernel code to declare memory that will be shared by all kernels of a workgroup

## Input buffers
# Requires NumPy : the arrays are transfered into the device memory
# Specify memory flags and host buffer to draw from
A_buf = cl.Buffer(ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=A)
B_buf = cl.Buffer(ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=B)

# Output buffer
C_buf = cl.Buffer(ctx, cl.mem_flags.WRITE_ONLY, A.nbytes)

# Build program
prg = cl.Program(ctx, kernel_as_string).build()

# Apply program
# Remember that OpenCL computing is SIMD in essence.
# The following line of code will call the `sum` kernel for indices up to A.shape (indices are the i
# in the kernel above)
prg.sum(queue, A.shape, None, # queue, global_size, local_size
    A_buf, B_buf, C_buf) # argumenrs

    # global_size gives the number of workgroups, or Blocks in CUDA terms.
    # local_size is the size of each of them in terms of threads

    # TODO : just need to find that out. Maybe assign one Block per array for multiprocessing !!
    # It's not so simple. Or is it ? Yeah it is ! Use the individual array size as local_size, and the number of arrays as global_size ?
    # NO ! not recommended ! Although it would work, that is not the intended design. it would try to allocate too many resources at once I think.


cl.enqueue_copy(queue, C_buf, C) # Transpose C_buf into C
"""
