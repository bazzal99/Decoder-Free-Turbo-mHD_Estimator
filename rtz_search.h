/**
 * @file rtz_search.h
 * @brief Graph-girth RTZ search kernels for minimum distance estimation.
 *
 * The minimum distance estimator works by searching for short cycles in the
 * correlation graph of the QPP interleaver. A cycle of girth g corresponds
 * to an input sequence of weight proportional to g that is simultaneously
 * RTZ for both component encoders — i.e., a low-weight turbo codeword.
 *
 * Two kernel families are provided, each in ZT and TB variants:
 *
 *  RTZ-2 kernels (Min_Girth_One_Layer_RTZ_2 / _Circular):
 *    Exploit the periodic structure of weight-2 RTZ sequences.
 *    Neighbor steps are multiples of LTE_ENCODER_PERIOD (= 7).
 *    Used for input weights iw = 2, 4, 6, ...
 *
 *  RTZ-Multiple kernels (Min_Girth_One_Layer_RTZ_Multiple_with_End / _Circular):
 *    General search with unit steps.
 *    Used for input weight iw = 3.
 *
 * Each "Layer" function is the public entry point for one starting position.
 * The "Ref" functions are the recursive implementation.
 *
 * Global operation counters (total_operations_multiple_iw2,
 * total_operations_not_multiple) are defined in rtz_search.c and
 * declared extern here for reporting by the main driver.
 */

#ifndef RTZ_SEARCH_H
#define RTZ_SEARCH_H

/* =========================================================================
 * Operation counters (for complexity reporting)
 * ========================================================================= */

extern long long total_operations_multiple_iw2;
extern long long total_operations_not_multiple;


/* =========================================================================
 * RTZ-2 kernel — Zero-Termination (ZT) mode
 * ========================================================================= */

/**
 * Search for weight-2 RTZ cycles starting from interleaver position `layer`.
 *
 * @param layer        Starting position in the interleaver domain.
 * @param pi           QPP forward map:  pi[i]  = π(i).
 * @param depi         QPP inverse map:  depi[i] = π⁻¹(i).
 * @param threshold    Current best distance (updated in place if improved).
 * @param mult         Multiplicity counter (updated in place).
 * @param max_RTZ_size Maximum half-period budget (controls search depth).
 * @param mingirth     Girth bound (controls maximum cycle length).
 * @param framesize    Frame length K.
 * @return             0 if a new minimum was found, 1 otherwise.
 */
int Min_Girth_One_Layer_RTZ_2(int layer, int pi[], int depi[],
                               int *threshold, int *mult,
                               int max_RTZ_size, int mingirth, int framesize);

int Min_Girth_One_Ref_With_Layer_RTZ_2(int layer, int pi[], int depi[],
                                        int girth, int ref, int the_ref, int starting,
                                        int path[], int pathLen,
                                        int *threshold, int *mult,
                                        int max_RTZ_size, int mingirth, int framesize);


/* =========================================================================
 * RTZ-Multiple kernel — Zero-Termination (ZT) mode
 * ========================================================================= */

/**
 * Search for weight-3+ RTZ cycles starting from interleaver position `layer`.
 * (Same parameter semantics as Min_Girth_One_Layer_RTZ_2.)
 */
int Min_Girth_One_Layer_RTZ_Multiple_with_End(int layer, int pi[], int depi[],
                                               int *threshold, int *mult,
                                               int max_RTZ_size, int mingirth, int framesize);

int Min_Girth_One_Ref_With_Layer_RTZ_Multiple_with_End(int layer, int pi[], int depi[],
                                                        int girth, int ref, int the_ref,
                                                        int starting, int path[], int pathLen,
                                                        int *threshold, int *mult,
                                                        int max_RTZ_size, int mingirth, int framesize);


/* =========================================================================
 * RTZ-2 kernel — Circular / Tail-Biting (TB) mode
 * ========================================================================= */

int Min_Girth_One_Layer_RTZ_2_Circular(int layer, int pi[], int depi[],
                                        int *threshold, int *mult,
                                        int max_RTZ_size, int mingirth, int framesize);

int Min_Girth_One_Ref_With_Layer_RTZ_2_Circular(int layer, int pi[], int depi[],
                                                  int girth, int ref, int the_ref, int starting,
                                                  int path[], int pathLen,
                                                  int *threshold, int *mult,
                                                  int max_RTZ_size, int mingirth, int framesize);


/* =========================================================================
 * RTZ-Multiple kernel — Circular / Tail-Biting (TB) mode
 * ========================================================================= */

int Min_Girth_One_Layer_RTZ_Multiple_Circular(int layer, int pi[], int depi[],
                                               int *threshold, int *mult,
                                               int max_RTZ_size, int mingirth, int framesize);

int Min_Girth_One_Ref_With_Layer_RTZ_Multiple_Circular(int layer, int pi[], int depi[],
                                                        int girth, int ref, int the_ref,
                                                        int starting, int path[], int pathLen,
                                                        int *threshold, int *mult,
                                                        int max_RTZ_size, int mingirth, int framesize);

#endif /* RTZ_SEARCH_H */
