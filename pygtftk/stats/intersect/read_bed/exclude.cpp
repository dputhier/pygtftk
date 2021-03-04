/* This function will perform exclusion for a single chromosome.
 *
 * For more details about role of this function, please see the documentation in
 * read_bed_as_list.pyx
 *
 * Author : Quentin FERRE <quentin.ferre@gmail.com>
 */

 /* IMPLEMENTATION NOTES
  *
  * We use hand-crafted C++ here because of the computational cost. It is meant
  * to alleviate the bottleneck observed in the old pandas code. For the other
  * parts of the code, manually integrating C++ code was not worth the hassle
  * compared to the advantages of easy writing with Cython.
  *
  * Some helper variables are declared to make the code more legible (eg. excl_length)
  *
  * The reason why I use a result array and do not perform the operations directly
  * on the bedfile_starts and bedfile_ends is because this means that after the
  * first exclusion has been performed we would be comparing positions in two
  * different coordinates systems (between 'exclusion' with the original ones
  * and the shifted 'bedfile'). Same reason why the conditions are computed from
  * the bedfile but the results expressed as a function of the previous value in
  * result : because the algorithm is iterative with one exclusion after the other.
  *
  * Code is based on the Cython C++ tutorial and NumPy binding tutorial.
  *
  * // Vectorisation
  *
  * In summary, I compute several intermediate arrays holding logical conditions
  * beause such computations can be SIMD-vectorized easily by the compiler.
  * Then, the results are used as masks for the processing.
  *
  * This also motivates the use of C++ because, although possible, it is harder
  * to write SIMD code in Cython.
  *
  * Masks and operations are written as for loops for better vectorization : we
  * compile with `-O3`, so loop unrolling will be used by gcc. As such, 'else if'
  * constructions are avoided. Same reason I use array indices instead of pointers.
  * Somewhat counterintuitively, to ensure better vectorisation, it is better to
  * write five for loops.
  *
  * Source : https://software.intel.com/en-us/cpp-compiler-developer-guide-and-reference-programming-guidelines-for-vectorization
  */



////////////////////////////////////////////////////////////////////////////////

#include "exclude.hpp"

namespace exclusion {

/* Template of functions to write in a boolean logical array (result_array)
 * whether each value of a given other array is under a threshold. */
template <class T>
void arrayWiseIsLower(const T * array, T value, long array_size, bool * result_array){

  const T * const end = array + array_size;
  bool * p_result = result_array;

  for(const T * p = array; p != end; ++p){
    *p_result = (*p < value);
    p_result++;
  }

  // TODO: If the array is sorted, this can be done more efficiently with
  // std::lower_bound ; but the current code is already fast enough in most cases.

  return;
}







/* This function is passed numpy arrays. In other words, it is passsed pointers
 * to C arrays :) Remember that the operations are performed in-place
 * (we are passing references) so this function returns void. */
void cpp_excludeConcatenateForThisChrom(long long* bedfile_starts, long long* bedfile_ends,
                                        long long* exclusion_starts, long long* exclusion_ends,
                                        long long* result_starts, long long* result_ends,
                                        long bed_size, long excl_size)
{
  // Iterating over all exclusion regions will definitely use a pointer.
  const long long* p_start_excl = exclusion_starts;   // pointer onto exclusion_starts
  const long long* p_end_excl = exclusion_ends;       // pointer onto exclusion_ends

  const long long* const p_excl_END = p_start_excl + excl_size; // pointer onto end of all_excl_starts

  long long start_excl_value;
  long long end_excl_value;



  // For each region in 'exclusion' :
  while (p_start_excl != p_excl_END){


    // Query the coordinates through the pointer
    start_excl_value = *p_start_excl;
    end_excl_value = *p_end_excl;

    long long excl_length = end_excl_value - start_excl_value;

    // Is region start under exclusion start ?
    bool* rs_lt_es = new bool[bed_size]();
    arrayWiseIsLower(bedfile_starts, start_excl_value, bed_size, rs_lt_es);

    // Is region end under exclusion start ?
    bool* re_lt_es = new bool[bed_size]();
    arrayWiseIsLower(bedfile_ends, start_excl_value, bed_size, re_lt_es);

    // Is region start under exclusion end ?
    bool* rs_lt_ee = new bool[bed_size]();
    arrayWiseIsLower(bedfile_starts, end_excl_value, bed_size, rs_lt_ee);

    // Is region end under exclusion end ?
    bool* re_lt_ee = new bool[bed_size]();
    arrayWiseIsLower(bedfile_ends, end_excl_value, bed_size, re_lt_ee);




    /*************************** Determining cases ****************************/
    /* The cases are arranged in priority order ! So be careful if you reorder
     * them. In the original code, they were "else ifs". Here this means we must
     * specify that a case can only be considered if all previous are false. */

    bool* case_0 = new bool[bed_size]();
    for(long i = 0; i < bed_size; ++i){
      case_0[i] = re_lt_es[i];
    }

    bool* case_1 = new bool[bed_size]();
    for(long i = 0; i < bed_size; ++i){
      if (!(case_0[i])){
        case_1[i] = (rs_lt_es[i] && (!(re_lt_es[i])) && re_lt_ee[i]);
      }
    }

    bool* case_2 = new bool[bed_size]();
    for(long i = 0; i < bed_size; ++i){
      if (!(case_0[i] || case_1[i])){
        case_2[i] = (rs_lt_es[i] && (!(re_lt_ee[i])));
      }
    }

    bool* case_3 = new bool[bed_size]();
    for(long i = 0; i < bed_size; ++i){
      if (!(case_0[i] || case_1[i] || case_2[i])){
        case_3[i] = ((!(rs_lt_es[i])) && re_lt_ee[i]);
      }
    }

    bool* case_4 = new bool[bed_size]();
    for(long i = 0; i < bed_size; ++i){
      if (!(case_0[i] || case_1[i] || case_2[i] || case_3[i])){
        case_4[i] = ((!(rs_lt_es[i])) && rs_lt_ee[i] && (!(re_lt_ee[i])));
      }
    }

    bool* case_5 = new bool[bed_size]();
    for(long i = 0; i < bed_size; ++i){
      if (!(case_0[i] || case_1[i] || case_2[i] || case_3[i] || case_4[i])){
        case_5[i] = ((!(rs_lt_ee[i])) && (!(re_lt_ee[i])));
      }
    }

    // Cleanup
    // TODO Rewrite this code with smart pointers and no manual deletion
    delete[] rs_lt_es;
    delete[] re_lt_es;
    delete[] rs_lt_ee;
    delete[] re_lt_ee;

    /***************************** Main processing ****************************/

    // TODO Try putting all `if`s in the same `for` loop. Without elses, they
    // are still treatable as masked assignments so the compiler should like it.



    // Case 0 : regions before exclusion. Do nothing.

    // Case 1 :
    /* All regions where region_start is under exclu_start but region_end is
     * higher than exclu_start BUT lower than excl_end: truncate by setting
     * region_end to exclu_start */
    for(long i = 0; i < bed_size; i++){
      if (case_1[i]){
        long long truncate_by = bedfile_ends[i] - start_excl_value;
        result_ends[i] -= truncate_by;
      }
    }

    // Case 2 :
    /* All which contain the excluded region (start before and end after) :
     * shorten the end by the region length */
    for(long i = 0; i < bed_size; i++){
      if (case_2[i]){
        result_ends[i] -= excl_length;
      }
    }

    // Case 3 :
    /* All regions where region_start > excl_start but region_end < excl_end
     * (so are included) : eliminate those by setting 0 0 */
    for(long i = 0; i < bed_size; i++){
      if (case_3[i]){
        result_starts[i] = 0;
        result_ends[i] = 0;
      }
    }

    // Case 4 :
    /* All regions where region_start is higher than excl_start but lower than
     * excl_end and region_end is higher than excl_end : truncate by setting
     * region_start to excl_end and also
     * region_end = region_end - nb_of_nt_of_region_that_are_in_excl */
    for(long i = 0; i < bed_size; i++){
      if (case_4[i]){

        // Compute some utils
        long long region_length_before_truncating = result_ends[i] - result_starts[i];
        long long nb_of_bp_of_region_that_are_in_excl = end_excl_value - bedfile_starts[i];

        // Move start point
        long long forward_by = bedfile_starts[i] - start_excl_value;
        result_starts[i] -= forward_by;

        // Move end point to 'new start point + new length'
        long long new_length = region_length_before_truncating - nb_of_bp_of_region_that_are_in_excl;
        result_ends[i] = result_starts[i] + new_length;
      }
    }

    // Case 5 :
    /* All regions where region_start and region_end are both higher than
     * excl_end : move by setting region_start = region_start - excl_length
     * and region_end = region_end - excl_length */
    for(long i = 0; i < bed_size; i++){
      if (case_5[i]){
        result_starts[i] -= excl_length;
        result_ends[i] -= excl_length;
      }
    }

    // Cleanup
    // TODO Rewrite this code with smart pointers and no manual deletion
    delete[] case_0;
    delete[] case_1;
    delete[] case_2;
    delete[] case_3;
    delete[] case_4;
    delete[] case_5;

    // Finally, move to the next exclusion region
    p_start_excl++;
    p_end_excl++;
  }

  // Return once we have processed all exclusion regions
  return;
}





}
