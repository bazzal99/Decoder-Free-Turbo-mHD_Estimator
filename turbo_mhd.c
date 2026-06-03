/**
 * @file turbo_mhd.c
 * @brief Decoder-independent minimum Hamming distance estimator for LTE Turbo Codes.
 *
 * Entry point and top-level estimation driver.
 *
 * Supports all 188 LTE block sizes (K = 40 … 6144, 3GPP TS 36.212) and
 * both termination modes:
 *   - Zero-Termination (ZT)  — circular = 0
 *   - Tail-Biting (TB)       — circular = 1
 *
 * Usage:
 *   Edit the configuration block in main(), compile, and run:
 *     make
 *     ./turbo_mhd
 *
 * Output files (written to the current directory):
 *   output_<K>.txt          — minimum-distance codewords: "d=<d>, <pos0> <pos1> ..."
 *   results.txt             — summary: K, elapsed time, multiplicity, best distance
 *   Number_of_Operations.txt — operation count breakdown (iw=2 vs iw=3+)
 *
 * Reference:
 *   M. Bazzal, J. Nadal, S. Weithoffer, C. Abdel Nour, C. Douillard,
 *   "Efficient decoder-free minimum distance estimation for concatenated
 *   convolutional codes," IEEE VTC, Oslo, July 2025.
 *
 * Author:  Mohammad Bazzal, IMT Atlantique, Lab-STICC, Brest, France
 * License: MIT
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "trellis.h"
#include "distance.h"
#include "rtz_search.h"
#include "utils.h"

/* =========================================================================
 * LTE QPP interleaver parameter tables (3GPP TS 36.212 Table 5.1.3-3)
 * ========================================================================= */

#define NUM_LTE_BLOCK_SIZES 188

static const int LTE_CODE_BLOCK_SIZES[NUM_LTE_BLOCK_SIZES] = {
    40,  48,  56,  64,  72,  80,  88,  96, 104, 112, 120, 128, 136, 144, 152,
   160, 168, 176, 184, 192, 200, 208, 216, 224, 232, 240, 248, 256, 264, 272,
   280, 288, 296, 304, 312, 320, 328, 336, 344, 352, 360, 368, 376, 384, 392,
   400, 408, 416, 424, 432, 440, 448, 456, 464, 472, 480, 488, 496, 504, 512,
   528, 544, 560, 576, 592, 608, 624, 640, 656, 672, 688, 704, 720, 736, 752,
   768, 784, 800, 816, 832, 848, 864, 880, 896, 912, 928, 944, 960, 976, 992,
  1008,1024,1056,1088,1120,1152,1184,1216,1248,1280,1312,1344,1376,1408,1440,
  1472,1504,1536,1568,1600,1632,1664,1696,1728,1760,1792,1824,1856,1888,1920,
  1952,1984,2016,2048,2112,2176,2240,2304,2368,2432,2496,2560,2624,2688,2752,
  2816,2880,2944,3008,3072,3136,3200,3264,3328,3392,3456,3520,3584,3648,3712,
  3776,3840,3904,3968,4032,4096,4160,4224,4288,4352,4416,4480,4544,4608,4672,
  4736,4800,4864,4928,4992,5056,5120,5184,5248,5312,5376,5440,5504,5568,5632,
  5696,5760,5824,5888,5952,6016,6080,6144
};

static const int LTE_F1[NUM_LTE_BLOCK_SIZES] = {
    3,  7, 19,  7,  7, 11,  5, 11,  7, 41,103, 15,  9, 17,  9, 21,101, 21,
   57, 23, 13, 27, 11, 27, 85, 29, 33, 15, 17, 33,103, 19, 19, 37, 19, 21,
   21,115,193, 21,133, 81, 45, 23,243,151,155, 25, 51, 47, 91, 29, 29,247,
   29, 89, 91,157, 55, 31, 17, 35,227, 65, 19, 37, 41, 39,185, 43, 21,155,
   79,139, 23,217, 25, 17,127, 25,239, 17,137,215, 29, 15,147, 29, 59, 65,
   55, 31, 17,171, 67, 35, 19, 39, 19,199, 21,211, 21, 43,149, 45, 49, 71,
   13, 17, 25,183, 55,127, 27, 29, 29, 57, 45, 31, 59,185,113, 31, 17,171,
  209,253,367,265,181, 39, 27,127,143, 43, 29, 45,157, 47, 13,111,443, 51,
   51,451,257, 57,313,271,179,331,363,375,127, 31, 33, 43, 33,477, 35,233,
  357,337, 37, 71, 71, 37, 39,127, 39, 39, 31,113, 41,251, 43, 21, 43, 45,
   45,161, 89,323, 47, 23, 47,263
};

static const int LTE_F2[NUM_LTE_BLOCK_SIZES] = {
   10, 12, 42, 16, 18, 20, 22, 24, 26, 84, 90, 32, 34,108, 38,120, 84, 44,
   46, 48, 50, 52, 36, 56, 58, 60, 62, 32,198, 68,210, 36, 74, 76, 78,120,
   82, 84, 86, 44, 90, 46, 94, 48, 98, 40,102, 52,106, 72,110,168,114, 58,
  118,180,122, 62, 84, 64, 66, 68,420, 96, 74, 76,234, 80, 82,252, 86, 44,
  120, 92, 94, 48, 98, 80,102, 52,106, 48,110,112,114, 58,118, 60,122,124,
   84, 64, 66,204,140, 72, 74, 76, 78,240, 82,252, 86, 88, 60, 92,846, 48,
   28, 80,102,104,954, 96,110,112,114,116,354,120,610,124,420, 64, 66,136,
  420,216,444,456,468, 80,164,504,172, 88,300, 92,188, 96, 28,240,204,104,
  212,192,220,336,228,232,236,120,244,248,168, 64,130,264,134,408,138,280,
  142,480,146,444,120,152,462,234,158, 80, 96,902,166,336,170, 86,174,176,
  178,120,182,184,186, 94,190,480
};


/* =========================================================================
 * min_distance_CHOOSE
 * ========================================================================= */
void min_distance_CHOOSE(int f1, int f2, int threshold, int framesize,
                         int max_input_weight_test, int circular)
{
    static const int g_prime_2[] = {2, 4, 6, 8};
    static const int p_2[]       = {14, 10, 7, 4};
    static const int g_prime_M[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    static const int p_M[]       = {-1, -1, -1, 21, 7, 10, 14, 7, 5, 4};

    int OriginalAddress[framesize];
    int InterleavedAddress[framesize];
    for (long long idx = 0; idx < framesize; idx++) {
        long long pos            = (f1 * idx + f2 * idx * idx) % framesize;
        OriginalAddress[idx]     = (int)pos;
        InterleavedAddress[(int)pos] = (int)idx;
    }

    char filename[32];
    sprintf(filename, "output_%d.txt", framesize);
    FILE *file = fopen(filename, "w");
    if (file) fclose(file);

    int mult = 0;

    /* ZT: quick single-impulse threshold initialiser */
    if (circular == 0) {
        int x[1];
        for (int i = framesize - 1; i > framesize / 2; i--) {
            x[0] = i;
            int d = Turbo_parity_distance_Without_circular_With_Puncturing(
                        x, 1, OriginalAddress, framesize);
            if (d < threshold) {
                file = fopen(filename, "w");
                fprintf(file, "d=%d, %d \n", d, i);
                fclose(file);
                threshold = d;
            }
        }
    }

    /* Weight-2 RTZ search */
    for (int iw = 0; iw < 3; iw++) {
        if (g_prime_2[iw] > max_input_weight_test) break;
        printf("K=%d  iw=%d  [RTZ-2 kernel]\n", framesize, g_prime_2[iw]);
        for (int i = framesize - 1; i >= 0; i--) {
            if (circular)
                Min_Girth_One_Layer_RTZ_2_Circular(i, InterleavedAddress, OriginalAddress,
                                                   &threshold, &mult, p_2[iw], g_prime_2[iw], framesize);
            else
                Min_Girth_One_Layer_RTZ_2(i, InterleavedAddress, OriginalAddress,
                                          &threshold, &mult, p_2[iw], g_prime_2[iw], framesize);
        }
    }

    /* Weight-3 RTZ search */
    for (int iw = 3; iw <= 3; iw++) {
        printf("K=%d  iw=%d  [RTZ-Multiple kernel]\n", framesize, iw);
        for (int i = framesize - 1; i >= 0; i--) {
            if (circular)
                Min_Girth_One_Layer_RTZ_Multiple_Circular(i, InterleavedAddress, OriginalAddress,
                                                          &threshold, &mult, p_M[iw], g_prime_M[iw], framesize);
            else
                Min_Girth_One_Layer_RTZ_Multiple_with_End(i, InterleavedAddress, OriginalAddress,
                                                          &threshold, &mult, p_M[iw], g_prime_M[iw], framesize);
        }
    }

    (void)max_input_weight_test;

    printf("Total operations = %.2e\n",
           (double)(total_operations_not_multiple + total_operations_multiple_iw2));
    FILE *f_ops = fopen("Number_of_Operations.txt", "w");
    if (f_ops) {
        fprintf(f_ops, "number of iw_2 operations  = %.2e\n", (double)total_operations_multiple_iw2);
        fprintf(f_ops, "number of other operations = %.2e\n\n", (double)total_operations_not_multiple);
        fprintf(f_ops, "Total number of operations = %.2e\n",
                (double)(total_operations_not_multiple + total_operations_multiple_iw2));
        fclose(f_ops);
    }

    removeDuplicateLines(filename);
}


/* =========================================================================
 * main — edit the four constants below, then: make && ./turbo_mhd
 * ========================================================================= */
int main(void)
{
    const int K         = 6144; /* LTE frame size in bits                              */
    const int iw        = 7;    /* Maximum input weight to search                      */
    const int threshold = 70;   /* Initial distance upper bound (tightened at runtime) */
    const int circular  = 0;    /* 0 = zero-termination (ZT),  1 = tail-biting (TB)   */

    int index = -1;
    for (int i = 0; i < NUM_LTE_BLOCK_SIZES; i++) {
        if (LTE_CODE_BLOCK_SIZES[i] == K) { index = i; break; }
    }
    if (index < 0) {
        fprintf(stderr, "Error: K=%d is not a valid LTE block size.\n", K);
        return 1;
    }

    Initialize(CircularStates, TrellisNextState, TrellisOutput,
               BitValue_ForTheDistance, K);

    time_t t_start = time(NULL);
    min_distance_CHOOSE(LTE_F1[index], LTE_F2[index], threshold, K, iw, circular);
    double elapsed = difftime(time(NULL), t_start);

    char filename[32];
    sprintf(filename, "output_%d.txt", K);
    FILE *file = fopen(filename, "r");
    if (!file) { perror("Cannot open output file"); return 1; }

    char line[64] = {0};
    int  multiplicity = 0;
    while (fgets(line, sizeof(line), file)) multiplicity++;
    fclose(file);

    printf("\nK = %d\n", K);
    printf("Elapsed time:  %.2f seconds\n", elapsed);
    printf("Multiplicity:  %d\n", multiplicity);
    printf("Best result:   %s", line);

    file = fopen("results.txt", "w");
    if (file) {
        fprintf(file, "K = %d\n", K);
        fprintf(file, "Elapsed time: %.2f seconds\n", elapsed);
        fprintf(file, "multiplicity = %d\n", multiplicity);
        fprintf(file, "%s", line);
        fclose(file);
    }

    return 0;
}
