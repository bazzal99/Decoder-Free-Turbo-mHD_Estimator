/**
 * @file utils.h
 * @brief General-purpose utility functions used across the estimator.
 *
 * Covers integer/binary conversion, array operations, sorting,
 * circular modular arithmetic, and output deduplication.
 */

#ifndef UTILS_H
#define UTILS_H

/* =========================================================================
 * Integer / binary conversion
 * ========================================================================= */

/** Convert an integer to its binary representation (MSB first). */
void int_to_binary(int num, int size, int binary_arr[]);

/** Convert a binary array (MSB first) back to an integer. */
int  binary_to_int(int binary_arr[], int size);


/* =========================================================================
 * Array utilities
 * ========================================================================= */

/**
 * Search for `target` in `arr[0..n-1]`.
 * Returns the index of the first match, or -1 if not found.
 */
int  search(int target, const int arr[], int n);

/** Count the number of 1s in an integer array. */
int  Number_Of_Ones(int a[], int size);

/** Element-wise XOR: c[i] = a[i] ^ b[i]. */
void Array_XOR(int a[], int b[], int c[], int size);

/** Copy `size` elements from `source` into `destination`. */
void duplicateArray(int source[], int destination[], int size);

/** In-place bubble sort (ascending). */
void bubbleSort(int arr[], int size);


/* =========================================================================
 * Circular / modular arithmetic
 * ========================================================================= */

/**
 * Non-negative modulus: returns (value % range + range) % range.
 * Always in [0, range-1], unlike C's `%` operator for negative values.
 */
int  circular_modulus(int value, int range);


/* =========================================================================
 * Output deduplication
 * ========================================================================= */

/**
 * Remove duplicate lines from `filename` in-place, preserving order.
 * Uses a chained hash table for O(n) average-case performance.
 */
void removeDuplicateLines(const char *filename);

#endif /* UTILS_H */
