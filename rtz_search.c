/**
 * @file rtz_search.c
 * @brief Graph-girth RTZ search kernel implementations.
 *
 * Each kernel recursively explores the correlation graph of the QPP
 * interleaver, alternating between the original domain (even depths) and
 * the interleaved domain (odd depths). When a closed cycle is found,
 * the corresponding input pattern is passed to the distance computation
 * functions in distance.c.
 *
 * The four kernels differ along two axes:
 *   - Step policy: multiples of LTE_ENCODER_PERIOD (RTZ-2) vs unit steps (RTZ-Multiple)
 *   - Distance oracle: ZT impulse-response (Without_circular) vs TB trellis (circular)
 */

#include <stdio.h>
#include "trellis.h"
#include "distance.h"
#include "rtz_search.h"
#include "utils.h"

/* =========================================================================
 * Operation counters
 * ========================================================================= */

long long total_operations_multiple_iw2 = 0;
long long total_operations_not_multiple  = 0;


/* =========================================================================
 * Shared helper: write a candidate codeword to the output file.
 *
 * If d < *threshold: overwrite the file with this new best result.
 * If d == *threshold: append this codeword (record the multiplicity).
 * ========================================================================= */

static void record_candidate(int d, int *threshold, int *mult,
                              int path[], int pathLen, int framesize)
{
    char filename[32];
    FILE *file;
    int   hid[20];

    sprintf(filename, "output_%d.txt", framesize);

    if (d == *threshold) {
        (*mult)++;
        duplicateArray(path, hid, pathLen);
        bubbleSort(hid, pathLen);
        file = fopen(filename, "a");
        fprintf(file, "d=%d, ", d);
        for (int i = 0; i < pathLen; i++) fprintf(file, "%d ", hid[i]);
        fprintf(file, "\n");
        fclose(file);
    } else if (d < *threshold) {
        *mult      = 0;
        *threshold = d;
        duplicateArray(path, hid, pathLen);
        bubbleSort(hid, pathLen);
        file = fopen(filename, "w");
        fprintf(file, "d=%d, ", d);
        for (int i = 0; i < pathLen; i++) fprintf(file, "%d ", hid[i]);
        fprintf(file, "\n");
        fclose(file);
        printf("\nNew best: d=%d  positions:", d);
        for (int i = 0; i < pathLen; i++) printf(" %d", hid[i]);
        printf("\n");
    }
}


/* =========================================================================
 * RTZ-2 kernel — Zero-Termination (ZT) mode
 *
 * Neighbors are spaced by multiples of LTE_ENCODER_PERIOD (=7), exploiting
 * the periodic structure of weight-2 RTZ sequences for the LTE encoder.
 * ========================================================================= */

int Min_Girth_One_Layer_RTZ_2(int layer, int pi[], int depi[],
                               int *threshold, int *mult,
                               int max_RTZ_size, int mingirth, int framesize)
{
    int path[20];
    return Min_Girth_One_Ref_With_Layer_RTZ_2(
               layer, pi, depi, 0, layer, layer, 0,
               path, 0, threshold, mult, max_RTZ_size, mingirth, framesize);
}

int Min_Girth_One_Ref_With_Layer_RTZ_2(int layer, int pi[], int depi[],
                                        int girth, int ref, int the_ref, int starting,
                                        int path[], int pathLen,
                                        int *threshold, int *mult,
                                        int max_RTZ_size, int mingirth, int framesize)
{
    if (ref < 0 || ref >= framesize) return 1;
    /* Symmetry pruning: only explore ref >= the_ref in the original domain */
    if (starting == 0 && ref < the_ref) return 1;

    path[pathLen] = (starting == 1) ? depi[ref] : ref;
    pathLen++;

    if (girth > mingirth) return 1;

    /* Closed cycle detected: evaluate distance */
    if (ref == the_ref && starting == 0 && girth != 0 && girth <= mingirth) {
        int d = Turbo_parity_distance_Without_circular_With_Puncturing(
                    path, pathLen - 1, depi, framesize);
        total_operations_multiple_iw2++;
        record_candidate(d, threshold, mult, path, pathLen - 1, framesize);
        return (d < *threshold) ? 0 : 1;
    }

    int x = 1;
    for (int i = LTE_ENCODER_PERIOD; i <= max_RTZ_size * LTE_ENCODER_PERIOD; i += LTE_ENCODER_PERIOD) {
        if (starting == 0) {
            if (ref - i >= 0)
                x &= Min_Girth_One_Ref_With_Layer_RTZ_2(layer, pi, depi, girth + 1,
                        pi[ref - i], the_ref, 1, path, pathLen, threshold, mult, max_RTZ_size, mingirth, framesize);
            if (ref + i < framesize)
                x &= Min_Girth_One_Ref_With_Layer_RTZ_2(layer, pi, depi, girth + 1,
                        pi[ref + i], the_ref, 1, path, pathLen, threshold, mult, max_RTZ_size, mingirth, framesize);
        } else {
            if (ref - i >= 0)
                x &= Min_Girth_One_Ref_With_Layer_RTZ_2(layer, pi, depi, girth + 1,
                        depi[ref - i], the_ref, 0, path, pathLen, threshold, mult, max_RTZ_size, mingirth, framesize);
            if (ref + i < framesize)
                x &= Min_Girth_One_Ref_With_Layer_RTZ_2(layer, pi, depi, girth + 1,
                        depi[ref + i], the_ref, 0, path, pathLen, threshold, mult, max_RTZ_size, mingirth, framesize);
        }
    }
    return x;
}


/* =========================================================================
 * RTZ-Multiple kernel — Zero-Termination (ZT) mode
 *
 * Neighbors are explored with unit steps, allowing detection of RTZ cycles
 * of arbitrary input weight (iw = 3, 5, ...).
 * ========================================================================= */

int Min_Girth_One_Layer_RTZ_Multiple_with_End(int layer, int pi[], int depi[],
                                               int *threshold, int *mult,
                                               int max_RTZ_size, int mingirth, int framesize)
{
    int path[20];
    return Min_Girth_One_Ref_With_Layer_RTZ_Multiple_with_End(
               layer, pi, depi, 0, layer, layer, 0,
               path, 0, threshold, mult, max_RTZ_size, mingirth, framesize);
}

int Min_Girth_One_Ref_With_Layer_RTZ_Multiple_with_End(int layer, int pi[], int depi[],
                                                        int girth, int ref, int the_ref,
                                                        int starting, int path[], int pathLen,
                                                        int *threshold, int *mult,
                                                        int max_RTZ_size, int mingirth, int framesize)
{
    if (ref < 0 || ref >= framesize) return 1;
    /* Ordering pruning: only explore ref <= the_ref in the original domain */
    if (starting == 0 && ref > the_ref) return 1;

    path[pathLen] = (starting == 1) ? depi[ref] : ref;
    /* Skip duplicate consecutive positions */
    if (path[pathLen] != path[pathLen - 1]) pathLen++;

    if (girth > mingirth) return 1;
    /* Cycle self-intersection check (excluding the target) */
    if (path[pathLen - 1] != the_ref && search(path[pathLen - 1], path, pathLen - 1) != -1)
        return 1;

    if (path[pathLen - 1] == the_ref && girth != 0 && girth <= mingirth) {
        int d = Turbo_parity_distance_Without_circular_With_Puncturing(
                    path, pathLen - 1, pi, framesize);
        total_operations_not_multiple++;
        record_candidate(d, threshold, mult, path, pathLen - 1, framesize);
        return (d < *threshold) ? 0 : 1;
    }

    int x = 1;
    for (int i = 1; i <= max_RTZ_size; i++) {
        if (starting == 0) {
            if (ref - i >= 0)
                x &= Min_Girth_One_Ref_With_Layer_RTZ_Multiple_with_End(layer, pi, depi, girth + 1,
                        pi[ref - i], the_ref, 1, path, pathLen, threshold, mult, max_RTZ_size, mingirth, framesize);
            if (ref + i < framesize)
                x &= Min_Girth_One_Ref_With_Layer_RTZ_Multiple_with_End(layer, pi, depi, girth + 1,
                        pi[ref + i], the_ref, 1, path, pathLen, threshold, mult, max_RTZ_size, mingirth, framesize);
        } else {
            if (ref - i >= 0)
                x &= Min_Girth_One_Ref_With_Layer_RTZ_Multiple_with_End(layer, pi, depi, girth + 1,
                        depi[ref - i], the_ref, 0, path, pathLen, threshold, mult, max_RTZ_size, mingirth, framesize);
            if (ref + i < framesize)
                x &= Min_Girth_One_Ref_With_Layer_RTZ_Multiple_with_End(layer, pi, depi, girth + 1,
                        depi[ref + i], the_ref, 0, path, pathLen, threshold, mult, max_RTZ_size, mingirth, framesize);
        }
    }
    return x;
}


/* =========================================================================
 * RTZ-2 kernel — Circular / Tail-Biting (TB) mode
 * ========================================================================= */

int Min_Girth_One_Layer_RTZ_2_Circular(int layer, int pi[], int depi[],
                                        int *threshold, int *mult,
                                        int max_RTZ_size, int mingirth, int framesize)
{
    int path[20];
    return Min_Girth_One_Ref_With_Layer_RTZ_2_Circular(
               layer, pi, depi, 0, layer, layer, 0,
               path, 0, threshold, mult, max_RTZ_size, mingirth, framesize);
}

int Min_Girth_One_Ref_With_Layer_RTZ_2_Circular(int layer, int pi[], int depi[],
                                                  int girth, int ref, int the_ref, int starting,
                                                  int path[], int pathLen,
                                                  int *threshold, int *mult,
                                                  int max_RTZ_size, int mingirth, int framesize)
{
    if (ref < 0 || ref >= framesize) return 1;
    if (starting == 0 && ref < the_ref) return 1;

    path[pathLen] = (starting == 1) ? depi[ref] : ref;
    pathLen++;

    if (girth > mingirth) return 1;

    total_operations_multiple_iw2++;

    if (ref == the_ref && starting == 0 && girth != 0 && girth <= mingirth) {
        int d = Turbo_parity_distance_circular_With_Puncturing(
                    path, pathLen - 1, pi, framesize);
        record_candidate(d, threshold, mult, path, pathLen - 1, framesize);
        return (d < *threshold) ? 0 : 1;
    }

    int x = 1;
    for (int i = LTE_ENCODER_PERIOD; i <= max_RTZ_size * LTE_ENCODER_PERIOD; i += LTE_ENCODER_PERIOD) {
        int cm;
        if (starting == 0) {
            cm = circular_modulus(ref - i, framesize);
            if (cm != -1)
                x &= Min_Girth_One_Ref_With_Layer_RTZ_2_Circular(layer, pi, depi, girth + 1,
                        pi[cm], the_ref, 1, path, pathLen, threshold, mult, max_RTZ_size, mingirth, framesize);
            cm = circular_modulus(ref + i, framesize);
            if (cm != -1)
                x &= Min_Girth_One_Ref_With_Layer_RTZ_2_Circular(layer, pi, depi, girth + 1,
                        pi[cm], the_ref, 1, path, pathLen, threshold, mult, max_RTZ_size, mingirth, framesize);
        } else {
            cm = circular_modulus(ref - i, framesize);
            if (cm != -1)
                x &= Min_Girth_One_Ref_With_Layer_RTZ_2_Circular(layer, pi, depi, girth + 1,
                        depi[cm], the_ref, 0, path, pathLen, threshold, mult, max_RTZ_size, mingirth, framesize);
            cm = circular_modulus(ref + i, framesize);
            if (cm != -1)
                x &= Min_Girth_One_Ref_With_Layer_RTZ_2_Circular(layer, pi, depi, girth + 1,
                        depi[cm], the_ref, 0, path, pathLen, threshold, mult, max_RTZ_size, mingirth, framesize);
        }
    }
    return x;
}


/* =========================================================================
 * RTZ-Multiple kernel — Circular / Tail-Biting (TB) mode
 * ========================================================================= */

int Min_Girth_One_Layer_RTZ_Multiple_Circular(int layer, int pi[], int depi[],
                                               int *threshold, int *mult,
                                               int max_RTZ_size, int mingirth, int framesize)
{
    int path[20];
    return Min_Girth_One_Ref_With_Layer_RTZ_Multiple_Circular(
               layer, pi, depi, 0, layer, layer, 0,
               path, 0, threshold, mult, max_RTZ_size, mingirth, framesize);
}

int Min_Girth_One_Ref_With_Layer_RTZ_Multiple_Circular(int layer, int pi[], int depi[],
                                                        int girth, int ref, int the_ref,
                                                        int starting, int path[], int pathLen,
                                                        int *threshold, int *mult,
                                                        int max_RTZ_size, int mingirth, int framesize)
{
    if (ref < 0 || ref >= framesize) return 1;

    path[pathLen] = (starting == 1) ? depi[ref] : ref;
    pathLen++;

    if (girth > mingirth) return 1;
    if (path[pathLen - 1] != the_ref && search(path[pathLen - 1], path, pathLen - 1) != -1)
        return 1;

    total_operations_not_multiple++;

    if (path[pathLen - 1] == the_ref && girth != 0 && girth <= mingirth) {
        int d = Turbo_parity_distance_circular_With_Puncturing(
                    path, pathLen - 1, pi, framesize);
        record_candidate(d, threshold, mult, path, pathLen - 1, framesize);
        return (d < *threshold) ? 0 : 1;
    }

    int x = 1;
    for (int i = 1; i <= max_RTZ_size; i++) {
        int cm;
        if (starting == 0) {
            cm = circular_modulus(ref - i, framesize);
            if (cm != -1)
                x &= Min_Girth_One_Ref_With_Layer_RTZ_Multiple_Circular(layer, pi, depi, girth + 1,
                        pi[cm], the_ref, 1, path, pathLen, threshold, mult, max_RTZ_size, mingirth, framesize);
            cm = circular_modulus(ref + i, framesize);
            if (cm != -1)
                x &= Min_Girth_One_Ref_With_Layer_RTZ_Multiple_Circular(layer, pi, depi, girth + 1,
                        pi[cm], the_ref, 1, path, pathLen, threshold, mult, max_RTZ_size, mingirth, framesize);
        } else {
            cm = circular_modulus(ref - i, framesize);
            if (cm != -1)
                x &= Min_Girth_One_Ref_With_Layer_RTZ_Multiple_Circular(layer, pi, depi, girth + 1,
                        depi[cm], the_ref, 0, path, pathLen, threshold, mult, max_RTZ_size, mingirth, framesize);
            cm = circular_modulus(ref + i, framesize);
            if (cm != -1)
                x &= Min_Girth_One_Ref_With_Layer_RTZ_Multiple_Circular(layer, pi, depi, girth + 1,
                        depi[cm], the_ref, 0, path, pathLen, threshold, mult, max_RTZ_size, mingirth, framesize);
        }
    }
    return x;
}
