# Decoder-Free-Turbo-mHD_Estimator

A decoder-independent, computationally efficient tool for estimating the **minimum Hamming distance (mHD)** of LTE Turbo Codes using the Quadratic Permutation Polynomial (QPP) interleaver.

**Code:** https://github.com/bazzal99/Decoder-Free-Turbo-mHD_Estimator

Supports both termination modes defined in **3GPP TS 36.212**:
- **Zero-Termination (ZT)** — default LTE mode
- **Tail-Biting / Circular (TB)** — used in some LTE variants

---

## Background

The minimum Hamming distance governs the error floor at high SNR. Computing it exactly for large turbo codes is intractable by brute force, and classical impulse-based methods (SIM, DIM, EVIM) require a full decoder and are computationally expensive.

This implementation uses a **decoder-free approach** based on identifying **Return-To-Zero (RTZ)** input sequences directly from the encoder structure. By searching for short cycles in the **correlation graph** of the QPP interleaver, the algorithm finds the minimum-distance codewords orders of magnitude faster than decoder-based methods.

**Reference:**

> M. Bazzal, J. Nadal, S. Weithoffer, C. Abdel Nour, C. Douillard,  
> *"Efficient decoder-free minimum distance estimation for concatenated convolutional codes,"*  
> IEEE Vehicular Technology Conference (VTC), Oslo, July 2025.

---

## Repository structure

```
Decoder-Free-Turbo-mHD_Estimator/
├── turbo_mhd.c     — main() driver + top-level estimation loop (min_distance_CHOOSE)
├── rtz_search.c/h  — graph-girth RTZ search kernels (ZT and TB variants)
├── distance.c/h    — codeword distance computation + RTZ detection (is_RTZ)
├── trellis.c/h     — trellis construction, encoder, circular state computation
├── utils.c/h       — sorting, binary conversion, array ops, output deduplication
├── Makefile
├── LICENSE
└── README.md
```

### Module responsibilities

| File | Responsibility |
|------|---------------|
| `turbo_mhd.c` | Configuration, QPP table lookup, search orchestration, result reporting |
| `rtz_search.c` | Four recursive graph-girth kernels (RTZ-2 and RTZ-Multiple × ZT and TB) |
| `distance.c` | `Turbo_parity_distance_*` functions + `is_RTZ()` + impulse response convolution |
| `trellis.c` | LTE RSC trellis tables, `Encode()`, circular state discovery |
| `utils.c` | `bubbleSort`, `search`, `circular_modulus`, `duplicateArray`, `removeDuplicateLines` |

---

## Build

Requires a C99 compiler and the math library:

```bash
make
```

Or manually:

```bash
gcc -O2 -std=c99 -o turbo_mhd turbo_mhd.c trellis.c distance.c rtz_search.c utils.c -lm
```

---

## Usage

Edit the four constants at the top of `main()` in `turbo_mhd.c`:

```c
const int K         = 6144;  /* LTE frame size in bits (must be a valid LTE block size) */
const int iw        = 7;     /* Maximum input weight to search                          */
const int threshold = 70;    /* Initial distance upper bound (tightened at runtime)     */
const int circular  = 0;     /* 0 = zero-termination (ZT),  1 = tail-biting (TB)       */
```

Then:

```bash
make run       # or: make && ./turbo_mhd
```

### Output files

| File | Content |
|------|---------|
| `output_<K>.txt` | Minimum-distance codewords: `d=<dist>, <pos1> <pos2> ...` |
| `results.txt` | Summary: K, elapsed time, multiplicity, best distance |
| `Number_of_Operations.txt` | Operation count breakdown (iw=2 vs iw=3+) |

### Example (K = 6144, ZT mode, iw = 7)

```
K = 6144
Elapsed time: 4583.00 seconds
Multiplicity: 1
Best result:  d=26, 6124 6141
```

---

## Algorithm overview

The LTE turbo code concatenates two identical RSC encoders (G(1, 15/13)₈, constraint length 3) through a QPP interleaver π(i) = (f₁·i + f₂·i²) mod K.

A codeword has low Hamming weight when its input is simultaneously RTZ for both component encoders. The algorithm searches for such inputs by exploring **girth-limited cycles** in the QPP correlation graph:

- **RTZ-2 kernel** — exploits the period-7 structure of weight-2 RTZ sequences; steps by multiples of 7
- **RTZ-Multiple kernel** — general search with unit steps; catches weight-3+ sequences

For ZT codes, tail bits are handled analytically via `is_RTZ()` and a syndrome lookup table.  
For TB codes, the correct circular starting state is found and used in `Encode()`.

---

## Complexity

For K = 6144 (ZT mode, iw ≤ 7):
- ~2.57 × 10⁶ iw=2 operations, ~1.76 × 10⁴ other operations
- Runtime: ~76 minutes on a standard PC
- Reduction of **several orders of magnitude** vs decoder-based methods

---

## License

MIT License. See [LICENSE](LICENSE).

---

## Citation

If you use this code in your research, please cite:

```bibtex
@inproceedings{bazzal2025turbo,
  author    = {Bazzal, Mohammad and Nadal, Jeremy and Weithoffer, Stefan
               and {Abdel Nour}, Charbel and Douillard, Catherine},
  title     = {Efficient decoder-free minimum distance estimation for
               concatenated convolutional codes},
  booktitle = {IEEE Vehicular Technology Conference (VTC)},
  address   = {Oslo, Norway},
  month     = jul,
  year      = {2025}
}
```
