"""
Functions used to perform multiprocessing operations with the pygtftk.stats.multiprocessing.multiproc module
should return this structure each time, which consists of a pointer to a 2D C array, and its shape.
"""
cdef struct FUNC_2DARRAY_RESULT:
    # The 2D array will be represented flattened, as a 1D array
    long long * result_array
    long long[2] result_shape

# Function type to be used in the multiprocessing function, to be passed to to apply_func_multiproc_cython()
ctypedef FUNC_2DARRAY_RESULT (*MULTIPROC_CUSTOM_FUNCTYPE) (long long *, long long*, long long) nogil
