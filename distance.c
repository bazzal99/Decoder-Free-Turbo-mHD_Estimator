/**
 * @file distance.c
 * @brief Turbo codeword distance computation and RTZ detection.
 */

#include <stdlib.h>
#include "trellis.h"
#include "distance.h"
#include "utils.h"

/**
 * Error pattern (syndrome) lookup table for weight-2 RTZ detection.
 *
 * These are the 16 valid syndrome values for the LTE RSC encoder
 * G(1, 15/13)₈. A set of bit positions is RTZ iff its XOR-accumulated
 * residues (mod LTE_ENCODER_PERIOD) produce one of these syndromes after
 * XOR-ing with a tail-bit correction pattern.
 */
static const int ERROR_LIST_2[] = {
     0, 11, 22, 29, 39, 44, 49,  58,
    69, 78, 83, 88, 98, 105, 116, 127
};
#define ERROR_LIST_2_SIZE 16


/* =========================================================================
 * RTZ detection
 * ========================================================================= */

int is_RTZ(int x[], int size, int framesize)
{
    int y[LTE_ENCODER_PERIOD];
    int z[3];
    int xframesize[3];
    int news[LTE_ENCODER_PERIOD];

    for (int j = 0; j < LTE_ENCODER_PERIOD; j++) y[j] = 0;
    for (int i = 0; i < 3; i++) xframesize[i] = (framesize + i) % LTE_ENCODER_PERIOD;

    /* Accumulate XOR of position residues modulo LTE_ENCODER_PERIOD */
    for (int i = 0; i < size; i++)
        y[x[i] % LTE_ENCODER_PERIOD] ^= 1;

    int temp = binary_to_int(y, LTE_ENCODER_PERIOD);

    /* Test all 8 possible tail-bit patterns */
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < LTE_ENCODER_PERIOD; j++) news[j] = 0;
        int_to_binary(i, 3, z);
        for (int j = 0; j < 3; j++)
            if (z[j] == 1) news[xframesize[j]] = 1;
        int temp2 = binary_to_int(news, LTE_ENCODER_PERIOD);
        if (search(temp ^ temp2, ERROR_LIST_2, ERROR_LIST_2_SIZE) != -1)
            return i;
    }
    return 0;
}


/* =========================================================================
 * Impulse response convolution (ZT mode)
 * ========================================================================= */

void ImpulseResponse(int impulse_response[], int size, int position)
{
    /*
     * Periodic impulse response of the LTE RSC encoder G(1, 15/13)₈.
     * The encoder has period 7; this is the first period of its IIR output
     * when driven by a single 1 at time 0.
     */
    static const int iir[LTE_ENCODER_PERIOD] = {0, 1, 1, 1, 0, 0, 1};
    int j = 1;
    for (int i = 0; i < size; i++) {
        if      (i < position)  impulse_response[i] = 0;
        else if (i == position) impulse_response[i] = 1;
        else                  { impulse_response[i] = iir[j % LTE_ENCODER_PERIOD]; j++; }
    }
}

void Bit_X_Many_Impulse_Response(int InformationFrame[], int size,
                                  int impulse_response[])
{
    int a[size];
    for (int i = 0; i < size; i++) impulse_response[i] = 0;
    for (int i = 0; i < size; i++) {
        if (InformationFrame[i] == 1) {
            ImpulseResponse(a, size, i);
            Array_XOR(impulse_response, a, impulse_response, size);
        }
    }
}


/* =========================================================================
 * Distance computation
 * ========================================================================= */

int Turbo_parity_distance_Without_circular_With_Puncturing(int x[], int size,
                                                            int Original_address[],
                                                            int framesize)
{
    int interleaver[framesize + 3];
    int x_prime[size];
    int impulse_response[framesize + 3];
    int positions[3];
    int d = size;  /* systematic weight = input weight */

    /* --- First RSC encoder --- */
    for (int i = 0; i < framesize + 3; i++) interleaver[i] = 0;
    for (int i = 0; i < size; i++)          interleaver[x[i]] = 1;

    int pos = is_RTZ(x, size, framesize);
    int_to_binary(pos, 3, positions);
    for (int i = 0; i < 3; i++) {
        interleaver[framesize + i] = positions[i];
        d += positions[i];
    }
    Bit_X_Many_Impulse_Response(interleaver, framesize + 3, impulse_response);
    d += Number_Of_Ones(impulse_response, framesize + 3);

    /* --- Second RSC encoder (interleaved input) --- */
    for (int i = 0; i < framesize + 3; i++) interleaver[i] = 0;
    for (int i = 0; i < size; i++) {
        x_prime[i] = Original_address[x[i]];
        interleaver[x_prime[i]] = 1;
    }
    pos = is_RTZ(x_prime, size, framesize);
    int_to_binary(pos, 3, positions);
    for (int i = 0; i < 3; i++) {
        interleaver[framesize + i] = positions[i];
        d += positions[i];
    }
    Bit_X_Many_Impulse_Response(interleaver, framesize + 3, impulse_response);
    d += Number_Of_Ones(impulse_response, framesize + 3);

    return d;
}

int Turbo_parity_distance_circular_With_Puncturing(int x[], int size,
                                                    int Original_address[],
                                                    int framesize)
{
    int interleaver[framesize];
    int impulse_response[framesize];
    int d = size;

    /* --- First RSC encoder --- */
    for (int i = 0; i < framesize; i++) interleaver[i] = 0;
    for (int i = 0; i < size; i++)      interleaver[x[i]] = 1;
    Encode(interleaver, impulse_response, framesize);
    d += Number_Of_Ones(impulse_response, framesize);

    /* --- Second RSC encoder (interleaved input) --- */
    for (int i = 0; i < framesize; i++) interleaver[i] = 0;
    for (int i = 0; i < size; i++)      interleaver[Original_address[x[i]]] = 1;
    Encode(interleaver, impulse_response, framesize);
    d += Number_Of_Ones(impulse_response, framesize);

    return d;
}
