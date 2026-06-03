/**
 * @file trellis.h
 * @brief Trellis construction, encoder, and circular state computation
 *        for the LTE RSC component encoder G(1, 15/13)₈.
 *
 * The LTE turbo code uses two identical recursive systematic convolutional
 * (RSC) encoders with:
 *   - Generator polynomial: G(1, 15/13) in octal notation
 *   - Constraint length:    3  (4 trellis states)
 *   - Encoder period:       7
 *
 * Global trellis tables (TrellisNextState, TrellisOutput, CircularStates,
 * BitValue_ForTheDistance) are declared here and defined in trellis.c.
 * They are populated once by Initialize() before any distance computation.
 */

#ifndef TRELLIS_H
#define TRELLIS_H

/* =========================================================================
 * Dimension constants
 * ========================================================================= */

#define MAX_NUMBER_OF_INPUTS       2
#define MAX_NUMBER_OF_OUTPUTS      7
#define MAX_CONSTRAINED_LENGTH     5
#define MAX_BPS                    8
#define MAX_DISTANCE               256
#define MAX_INPUT_STATES           4
#define MAX_NUMBER_OF_CODE_STATES  32

/**
 * Period of the LTE RSC encoder impulse response.
 * G(1, 15/13)₈ has a periodic free response of period 7.
 */
#define LTE_ENCODER_PERIOD  7


/* =========================================================================
 * Global trellis tables (populated by Initialize)
 * ========================================================================= */

extern int BitValue_ForTheDistance[MAX_BPS][MAX_DISTANCE];
extern int TrellisNextState[MAX_NUMBER_OF_CODE_STATES][MAX_INPUT_STATES];
extern int TrellisOutput[MAX_NUMBER_OF_CODE_STATES][MAX_INPUT_STATES];
extern int CircularStates[MAX_NUMBER_OF_CODE_STATES];


/* =========================================================================
 * Public API
 * ========================================================================= */

/**
 * Initialize all global trellis tables for the LTE RSC encoder.
 *
 * Must be called once before any call to Encode() or
 * Turbo_parity_distance_circular_With_Puncturing().
 * For ZT mode it is harmless but still recommended for completeness.
 *
 * @param size_frame  Frame length K (needed to find the circular states).
 */
void Initialize(int CircularStates[MAX_NUMBER_OF_CODE_STATES],
                int TrellisNextState[MAX_NUMBER_OF_CODE_STATES][MAX_INPUT_STATES],
                int TrellisOutput[MAX_NUMBER_OF_CODE_STATES][MAX_INPUT_STATES],
                int BitValue_ForTheDistance[MAX_BPS][MAX_DISTANCE],
                int size_frame);

/**
 * Encode one frame using the global trellis tables.
 *
 * In circular (tail-biting) mode the encoder is run twice: once to find
 * the correct starting state (via CircularStates[]), and once to produce
 * the actual parity output.
 *
 * @param InformationData     Binary input frame of length `framesize`.
 * @param DataToBeTransmitted Output parity bits (length `framesize`).
 * @param framesize           Frame length K.
 * @return                    Parity weight (number of 1s in the output).
 */
int Encode(int InformationData[], int DataToBeTransmitted[], int framesize);

/* --- Internal helpers (used by Initialize; not needed externally) --- */
void TrellisCreationAndFindCircular(int NumberOfInputs, int NumberOfOutputs,
                                    int NumberOfRedundancyOutputs,
                                    int InputPolynomial[MAX_NUMBER_OF_INPUTS][MAX_CONSTRAINED_LENGTH],
                                    int ConstrainedLength,
                                    int RecursionPolynomial[MAX_CONSTRAINED_LENGTH][MAX_CONSTRAINED_LENGTH],
                                    int FeedForwardPolynomial[MAX_NUMBER_OF_OUTPUTS][MAX_NUMBER_OF_INPUTS+MAX_CONSTRAINED_LENGTH+1],
                                    int OctalInputPolynomial[MAX_NUMBER_OF_INPUTS],
                                    int OctalRecursionPolynomial,
                                    int OctalOutputPolynomial[MAX_NUMBER_OF_OUTPUTS],
                                    int BitValue_ForTheDistance[MAX_BPS][MAX_DISTANCE],
                                    int DistanceValue_ForTheBit[MAX_BPS][2][MAX_DISTANCE/2],
                                    int TrellisNextState[MAX_NUMBER_OF_CODE_STATES][MAX_INPUT_STATES],
                                    int TrellisOutput[MAX_NUMBER_OF_CODE_STATES][MAX_INPUT_STATES],
                                    int NumberOfCodeStates,
                                    int CircularStates[MAX_NUMBER_OF_CODE_STATES],
                                    int Information_FrameSize_K,
                                    int Damier, int size_frame);

void PrepareConvolutionalInputOneRecursion(int NumberOfInputs, int NumberOfOutputs,
                                           int NumberOfRedundancyOutputs,
                                           int InputPolynomial[MAX_NUMBER_OF_INPUTS][MAX_CONSTRAINED_LENGTH],
                                           int ConstrainedLength, int SourceAlphabetSize,
                                           int RecursionPolynomial[MAX_CONSTRAINED_LENGTH][MAX_CONSTRAINED_LENGTH],
                                           int FeedForwardPolynomial[MAX_NUMBER_OF_OUTPUTS][MAX_NUMBER_OF_INPUTS+MAX_CONSTRAINED_LENGTH+1],
                                           int OctalInputPolynomial[MAX_NUMBER_OF_INPUTS],
                                           int OctalRecursionPolynomial,
                                           int OctalOutputPolynomial[MAX_NUMBER_OF_OUTPUTS],
                                           int BitValue_ForTheDistance[MAX_BPS][MAX_DISTANCE],
                                           int DistanceValue_ForTheBit[MAX_BPS][2][MAX_DISTANCE/2]);

void ConvolutionalTrellisConstruction(int NumberOfInputs, int NumberOfOutputs,
                                      int NumberOfRedundancyOutputs,
                                      int InputPolynomial[MAX_NUMBER_OF_INPUTS][MAX_CONSTRAINED_LENGTH],
                                      int ConstrainedLength, int SourceAlphabetSize,
                                      int RecursionPolynomial[MAX_CONSTRAINED_LENGTH][MAX_CONSTRAINED_LENGTH],
                                      int FeedForwardPolynomial[MAX_NUMBER_OF_OUTPUTS][MAX_NUMBER_OF_INPUTS+MAX_CONSTRAINED_LENGTH+1],
                                      int TrellisNextState[MAX_NUMBER_OF_CODE_STATES][MAX_INPUT_STATES],
                                      int TrellisOutput[MAX_NUMBER_OF_CODE_STATES][MAX_INPUT_STATES]);

void Marginalisation(int Bit_Value_For_The_Distance[MAX_BPS][MAX_DISTANCE],
                     int Distance_Value_For_The_Bit[MAX_BPS][2][MAX_DISTANCE/2],
                     int BitPerSymbol);

void EncodingInputData(int Damierr, int NumInputs, int SourceAlphabetSize,
                       int StartingState, int *EndingState, int Sizee,
                       int InformationData[], int DataToBeTransmitted[],
                       int FindTrellisNextState[MAX_NUMBER_OF_CODE_STATES][MAX_INPUT_STATES],
                       int FindTrellisOutput[MAX_NUMBER_OF_CODE_STATES][MAX_INPUT_STATES],
                       int framesize);

void GenerateRandomInput(int InformationData[], int TailleTrame, int AlphabetSize);

#endif /* TRELLIS_H */
