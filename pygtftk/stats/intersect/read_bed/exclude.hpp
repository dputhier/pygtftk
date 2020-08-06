// Header for exclude.cpp ; see the read_bed_as_list.pyx file for its
// integration in the main Cython code.
namespace exclusion {
  void cpp_excludeConcatenateForThisChrom(long long* bedfile_starts, long long* bedfile_ends,
                                          long long* exclusion_starts, long long* exclusion_ends,
                                          long long* result_starts, long long* result_ends,
                                          long bed_size, long excl_size);
}
