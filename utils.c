/**
 * @file utils.c
 * @brief General-purpose utility function implementations.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "utils.h"

/* Hash table size for duplicate-line removal */
#define TABLE_SIZE    1009
#define MAX_LINE_LEN  60

/* =========================================================================
 * Integer / binary conversion
 * ========================================================================= */

void int_to_binary(int num, int size, int binary_arr[])
{
    for (int i = size - 1; i >= 0; i--) {
        binary_arr[i] = num % 2;
        num /= 2;
    }
}

int binary_to_int(int binary_arr[], int size)
{
    int decimal = 0;
    for (int i = 0; i < size; i++)
        decimal += binary_arr[i] * (int)pow(2, size - 1 - i);
    return decimal;
}


/* =========================================================================
 * Array utilities
 * ========================================================================= */

int search(int target, const int arr[], int n)
{
    for (int i = 0; i < n; i++)
        if (arr[i] == target) return i;
    return -1;
}

int Number_Of_Ones(int a[], int size)
{
    int s = 0;
    for (int i = 0; i < size; i++) s += a[i];
    return s;
}

void Array_XOR(int a[], int b[], int c[], int size)
{
    for (int i = 0; i < size; i++) c[i] = a[i] ^ b[i];
}

void duplicateArray(int source[], int destination[], int size)
{
    for (int i = 0; i < size; i++) destination[i] = source[i];
}

void bubbleSort(int arr[], int size)
{
    for (int i = 0; i < size - 1; i++)
        for (int j = 0; j < size - i - 1; j++)
            if (arr[j] > arr[j + 1]) {
                int tmp = arr[j];
                arr[j]  = arr[j + 1];
                arr[j + 1] = tmp;
            }
}


/* =========================================================================
 * Circular / modular arithmetic
 * ========================================================================= */

int circular_modulus(int value, int range)
{
    return (value % range + range) % range;
}


/* =========================================================================
 * Output deduplication (chained hash table)
 * ========================================================================= */

static unsigned long hash_line(const char *str)
{
    unsigned long h = 5381;
    int c;
    while ((c = *str++)) h = ((h << 5) + h) + c;
    return h % TABLE_SIZE;
}

struct Node {
    char line[MAX_LINE_LEN];
    struct Node *next;
};

static void ht_insert(struct Node *table[], const char *line)
{
    unsigned long idx = hash_line(line);
    struct Node *node = malloc(sizeof(struct Node));
    if (!node) { perror("malloc"); exit(1); }
    strncpy(node->line, line, MAX_LINE_LEN - 1);
    node->line[MAX_LINE_LEN - 1] = '\0';
    node->next  = table[idx];
    table[idx]  = node;
}

static int ht_find(struct Node *table[], const char *line)
{
    unsigned long idx = hash_line(line);
    struct Node *cur  = table[idx];
    while (cur) {
        if (strcmp(cur->line, line) == 0) return 1;
        cur = cur->next;
    }
    return 0;
}

void removeDuplicateLines(const char *filename)
{
    FILE *file = fopen(filename, "r");
    if (!file) { perror("fopen"); return; }

    struct Node *table[TABLE_SIZE];
    for (int i = 0; i < TABLE_SIZE; i++) table[i] = NULL;

    FILE *tmp = fopen("temp_dedup.txt", "w");
    if (!tmp) { perror("fopen temp"); fclose(file); return; }

    char line[MAX_LINE_LEN];
    while (fgets(line, sizeof(line), file)) {
        if (!ht_find(table, line)) {
            ht_insert(table, line);
            fputs(line, tmp);
        }
    }
    fclose(file);
    fclose(tmp);

    if (remove(filename) != 0)               { perror("remove"); return; }
    if (rename("temp_dedup.txt", filename) != 0) perror("rename");
}
