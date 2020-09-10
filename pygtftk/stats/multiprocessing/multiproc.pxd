# Import structures
from pygtftk.stats.multiprocessing.multiproc_structs cimport *

# --------- Multiproceessing header
cdef list apply_func_multiproc_cython(list, list, MULTIPROC_CUSTOM_FUNCTYPE, int)
