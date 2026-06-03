/**
 * @file distance.h
 * @brief Turbo codeword distance computation and RTZ detection.
 *
 * Two distance functions are provided, one for each termination mode:
 *
 *  - Turbo_parity_distance_Without_circular_With_Puncturing()
 *      Zero-termination (ZT): computes distance analytically via impulse
 *      response convolution. Tail bits are handled through is_RTZ().
 *
 *  - Turbo_parity_distance_circular_With_Puncturing()
 *      Tail-biting (TB): computes distance by running the circular trellis
 *      encoder directly (calls Encode() from trellis.h).
 *
 * Both functions count systematic bits + first RSC parity + second RSC parity
 * (after QPP interleaving), matching the rate-1/3 LTE turbo code structure.
 */

#ifndef DISTANCE_H
#define DISTANCE_H

/* =========================================================================
 * RTZ detection
 * ========================================================================= */

/**
 * Detect whether the input pattern x[0..size-1] is a Return-To-Zero (RTZ)
 * sequence for the LTE RSC encoder, and determine which tail bits are needed.
 *
 * The detection is based on XOR accumulation of position residues modulo
 * LTE_ENCODER_PERIOD, compared against a precomputed syndrome table.
 *
 * @param x         Array of bit positions where the input equals 1.
 * @param size      Number of 1s in the input (input weight).
 * @param framesize Frame length K.
 * @return          Tail-bit index (0–7) for the matching syndrome, or 0
 *                  if the sequence is not RTZ.
 */
int is_RTZ(int x[], int size, int framesize);


/* =========================================================================
 * Impulse response convolution (ZT mode)
 * ========================================================================= */

/**
 * Compute the parity output produced by a single impulse at position
 * `position` in a zero-terminated RSC frame of length `size`.
 *
 * Uses the periodic impulse response of the LTE RSC encoder:
 *   {0, 1, 1, 1, 0, 0, 1} (period = LTE_ENCODER_PERIOD = 7).
 */
void ImpulseResponse(int impulse_response[], int size, int position);

/**
 * Compute the parity output for an arbitrary binary input frame
 * using superposition of individual impulse responses.
 *
 * @param InformationFrame  Binary input (0s and 1s), length `size`.
 * @param size              Frame length (may include tail-bit positions).
 * @param impulse_response  Output parity bits (same length `size`).
 */
void Bit_X_Many_Impulse_Response(int InformationFrame[], int size,
                                  int impulse_response[]);


/* =========================================================================
 * Distance computation
 * ========================================================================= */

/**
 * Compute the Hamming weight of the turbo codeword produced by the sparse
 * input x[] under the zero-terminated LTE turbo encoder.
 *
 * Contribution = input weight (systematic) +
 *                first RSC parity weight   +
 *                tail bits (both encoders) +
 *                second RSC parity weight (interleaved input)
 *
 * @param x                Binary-1 positions in the input frame.
 * @param size             Input weight (number of 1s).
 * @param Original_address QPP interleaver mapping: Original_address[i] = π(i).
 * @param framesize        Frame length K.
 * @return                 Total codeword Hamming weight.
 */
int Turbo_parity_distance_Without_circular_With_Puncturing(int x[], int size,
                                                            int Original_address[],
                                                            int framesize);

/**
 * Compute the Hamming weight of the turbo codeword produced by the sparse
 * input x[] under the circular (tail-biting) LTE turbo encoder.
 *
 * Uses Encode() (from trellis.h) directly; no tail-bit correction needed.
 *
 * @param x                Binary-1 positions in the input frame.
 * @param size             Input weight.
 * @param Original_address QPP interleaver mapping.
 * @param framesize        Frame length K.
 * @return                 Total codeword Hamming weight.
 */
int Turbo_parity_distance_circular_With_Puncturing(int x[], int size,
                                                    int Original_address[],
                                                    int framesize);

#endif /* DISTANCE_H */
