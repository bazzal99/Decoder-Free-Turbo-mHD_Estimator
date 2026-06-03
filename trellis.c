/**
 * @file trellis.c
 * @brief Trellis construction, encoder, and circular state computation.
 *
 * Implements the LTE RSC component encoder G(1, 15/13)₈ (octal notation),
 * constraint length 3, using a generic state-machine framework that supports
 * arbitrary polynomial specifications in octal.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "trellis.h"

/* =========================================================================
 * Global trellis tables (defined here, declared extern in trellis.h)
 * ========================================================================= */

int BitValue_ForTheDistance[MAX_BPS][MAX_DISTANCE];
int TrellisNextState[MAX_NUMBER_OF_CODE_STATES][MAX_INPUT_STATES];
int TrellisOutput[MAX_NUMBER_OF_CODE_STATES][MAX_INPUT_STATES];
int CircularStates[MAX_NUMBER_OF_CODE_STATES];


/* =========================================================================
 * Initialize
 * ========================================================================= */

void Initialize(int CircularStates[MAX_NUMBER_OF_CODE_STATES],
                int TrellisNextState[MAX_NUMBER_OF_CODE_STATES][MAX_INPUT_STATES],
                int TrellisOutput[MAX_NUMBER_OF_CODE_STATES][MAX_INPUT_STATES],
                int BitValue_ForTheDistance[MAX_BPS][MAX_DISTANCE],
                int size_frame)
{
    int InputPolynomial[MAX_NUMBER_OF_INPUTS][MAX_CONSTRAINED_LENGTH];
    int RecursionPolynomial[MAX_CONSTRAINED_LENGTH][MAX_CONSTRAINED_LENGTH];
    int FeedForwardPolynomial[MAX_NUMBER_OF_OUTPUTS][MAX_NUMBER_OF_INPUTS + MAX_CONSTRAINED_LENGTH + 1];
    int OctalInputPolynomial[MAX_NUMBER_OF_INPUTS];
    int OctalOutputPolynomial[MAX_NUMBER_OF_OUTPUTS];
    int DistanceValue_ForTheBit[MAX_BPS][2][MAX_DISTANCE / 2];

    /* LTE RSC encoder: G(1, 15/13)₈ */
    const int NumberOfInputs            = 1;
    const int NumberOfOutputs           = 2;
    const int NumberOfRedundancyOutputs = 1;
    const int ConstrainedLength         = 3;
    const int OctalRecursionPolynomial  = 13;  /* denominator 1 + D² + D³  (octal 13) */

    OctalInputPolynomial[0]  = 4;
    OctalInputPolynomial[1]  = 0;
    OctalOutputPolynomial[0] = 15;  /* numerator  1 + D + D² + D³  (octal 15) */
    for (int i = 1; i < MAX_NUMBER_OF_OUTPUTS; i++) OctalOutputPolynomial[i] = 0;

    Marginalisation(BitValue_ForTheDistance, DistanceValue_ForTheBit, MAX_BPS);
    TrellisCreationAndFindCircular(
        NumberOfInputs, NumberOfOutputs, NumberOfRedundancyOutputs,
        InputPolynomial, ConstrainedLength, RecursionPolynomial,
        FeedForwardPolynomial,
        OctalInputPolynomial, OctalRecursionPolynomial, OctalOutputPolynomial,
        BitValue_ForTheDistance, DistanceValue_ForTheBit,
        TrellisNextState, TrellisOutput,
        (int)pow(2, ConstrainedLength), CircularStates,
        size_frame, 0, size_frame);
}


/* =========================================================================
 * TrellisCreationAndFindCircular
 * ========================================================================= */

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
                                    int Damier, int size_frame)
{
    int InformationData[size_frame];
    int EncoderOutput[size_frame];
    int StartingState, StartingStateMiddle, EndingState;
    int BoolCircular;
    int CodeNumberOfStates;

    PrepareConvolutionalInputOneRecursion(
        NumberOfInputs, NumberOfOutputs, NumberOfRedundancyOutputs,
        InputPolynomial, ConstrainedLength, 2,
        RecursionPolynomial, FeedForwardPolynomial,
        OctalInputPolynomial, OctalRecursionPolynomial, OctalOutputPolynomial,
        BitValue_ForTheDistance, DistanceValue_ForTheBit);

    ConvolutionalTrellisConstruction(
        NumberOfInputs, NumberOfOutputs, NumberOfRedundancyOutputs,
        InputPolynomial, ConstrainedLength, 2,
        RecursionPolynomial, FeedForwardPolynomial,
        TrellisNextState, TrellisOutput);

    /* Initialise circular state table */
    CodeNumberOfStates  = NumberOfCodeStates;
    CircularStates[0]   = 0;
    for (int i = 1; i < CodeNumberOfStates; i++) CircularStates[i] = 9999;

    /*
     * Find all circular (tail-biting) starting states by random probing.
     * For each distinct ending state encountered, record the corresponding
     * starting state. Repeat until all states are covered.
     */
    do {
        BoolCircular = 1;
        GenerateRandomInput(InformationData, Information_FrameSize_K, 2);
        StartingState = 0;

        EncodingInputData(Damier, NumberOfInputs, 2, StartingState,
                          &StartingStateMiddle, Information_FrameSize_K,
                          InformationData, EncoderOutput,
                          TrellisNextState, TrellisOutput, size_frame);

        EncodingInputData(Damier, NumberOfInputs, 2, StartingStateMiddle,
                          &EndingState, Information_FrameSize_K,
                          InformationData, EncoderOutput,
                          TrellisNextState, TrellisOutput, size_frame);

        CircularStates[EndingState] = StartingStateMiddle;

        for (int j = 1; j < CodeNumberOfStates; j++)
            if (CircularStates[j] == 9999) { BoolCircular = 0; break; }

    } while (BoolCircular == 0);
}


/* =========================================================================
 * PrepareConvolutionalInputOneRecursion
 *
 * Converts octal polynomial specifications into binary coefficient arrays
 * used by ConvolutionalTrellisConstruction().
 * ========================================================================= */

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
                                           int DistanceValue_ForTheBit[MAX_BPS][2][MAX_DISTANCE/2])
{
    (void)SourceAlphabetSize;          /* unused in this specialisation */
    (void)DistanceValue_ForTheBit;

    int l, o, i, k;
    int MemoryOctal;
    int SSValues[100];
    int SS[100];
    int Statet;

    /* Helper macro: determine octal digit count and expand into SS[] */
#define EXPAND_OCTAL(poly)                                          \
    do {                                                            \
        MemoryOctal = 0;                                            \
        if ((poly) >     7) MemoryOctal = 2;                       \
        if ((poly) >    77) MemoryOctal = 3;                       \
        if ((poly) >   777) MemoryOctal = 4;                       \
        if ((poly) >  7777) MemoryOctal = 5;                       \
        if ((poly) > 77777) MemoryOctal = 6;                       \
        for (k = 0; k < 100; k++) SS[k] = 0;                      \
        for (k = 0; k < MemoryOctal; k++)                          \
            SSValues[k] = (int)pow(10., k);                         \
        Statet = (poly);                                            \
        for (k = MemoryOctal - 1; k > 0; k--) {                   \
            SS[k] = Statet / SSValues[k];                           \
            Statet %= SSValues[k];                                  \
        }                                                           \
        SS[0] = Statet;                                             \
    } while (0)

    /* --- Recursion (feedback) polynomial --- */
    EXPAND_OCTAL(OctalRecursionPolynomial);
    for (o = 0; o < ConstrainedLength; o++)
        RecursionPolynomial[0][ConstrainedLength - 1 - o] =
            BitValue_ForTheDistance[o % 3][SS[o / 3]];
    for (l = 1; l < ConstrainedLength; l++)
        for (o = 0; o < ConstrainedLength; o++)
            RecursionPolynomial[l][o] = (o == l - 1) ? 1 : 0;

    /* --- Input polynomials --- */
    for (i = 0; i < NumberOfInputs; i++) {
        EXPAND_OCTAL(OctalInputPolynomial[i]);
        for (o = 0; o < ConstrainedLength; o++)
            InputPolynomial[i][ConstrainedLength - 1 - o] =
                BitValue_ForTheDistance[o % 3][SS[o / 3]];
    }

    /* --- Systematic (pass-through) output connections --- */
    for (i = 0; i < NumberOfOutputs - NumberOfRedundancyOutputs; i++) {
        for (o = 0; o < NumberOfInputs; o++)
            FeedForwardPolynomial[i + NumberOfRedundancyOutputs][ConstrainedLength + 1 + o] =
                (i == o) ? 1 : 0;
        for (o = 0; o < ConstrainedLength + 1; o++)
            FeedForwardPolynomial[i + NumberOfRedundancyOutputs][o] = 0;
    }

    /* --- Parity (feed-forward) output polynomials --- */
    for (i = 0; i < NumberOfRedundancyOutputs; i++) {
        EXPAND_OCTAL(OctalOutputPolynomial[i]);
        for (o = 0; o < ConstrainedLength + 1; o++)
            FeedForwardPolynomial[i][ConstrainedLength - o] =
                BitValue_ForTheDistance[o % 3][SS[o / 3]];
        for (o = 0; o < NumberOfOutputs - NumberOfRedundancyOutputs; o++)
            FeedForwardPolynomial[i][ConstrainedLength + 1 + o] = 0;
    }

#undef EXPAND_OCTAL
}


/* =========================================================================
 * ConvolutionalTrellisConstruction
 *
 * Fills TrellisNextState[state][input] and TrellisOutput[state][input]
 * for all states and inputs.
 * ========================================================================= */

void ConvolutionalTrellisConstruction(int NumberOfInputs, int NumberOfOutputs,
                                      int NumberOfRedundancyOutputs,
                                      int InputPolynomial[MAX_NUMBER_OF_INPUTS][MAX_CONSTRAINED_LENGTH],
                                      int ConstrainedLength, int SourceAlphabetSize,
                                      int RecursionPolynomial[MAX_CONSTRAINED_LENGTH][MAX_CONSTRAINED_LENGTH],
                                      int FeedForwardPolynomial[MAX_NUMBER_OF_OUTPUTS][MAX_NUMBER_OF_INPUTS+MAX_CONSTRAINED_LENGTH+1],
                                      int TrellisNextState[MAX_NUMBER_OF_CODE_STATES][MAX_INPUT_STATES],
                                      int TrellisOutput[MAX_NUMBER_OF_CODE_STATES][MAX_INPUT_STATES])
{
    int i, j, k, l, o;
    int SS[MAX_CONSTRAINED_LENGTH], SSValues[MAX_CONSTRAINED_LENGTH];
    int UpdatedSS[MAX_CONSTRAINED_LENGTH];
    int InputValue[MAX_NUMBER_OF_INPUTS], IIValues[MAX_NUMBER_OF_INPUTS];
    int OOValues[MAX_NUMBER_OF_INPUTS + MAX_CONSTRAINED_LENGTH + 1];
    int Statet, IIStatet, Intermediate, Intermediate2;

    int ActualNumberOfStates = (int)pow(SourceAlphabetSize, ConstrainedLength);
    int InputAtTime_t_Max    = (int)pow(SourceAlphabetSize, NumberOfInputs);

    for (k = 0; k < NumberOfOutputs;   k++) OOValues[k] = (int)pow(SourceAlphabetSize, k);
    for (k = 0; k < ConstrainedLength; k++) SSValues[k] = (int)pow(SourceAlphabetSize, k);
    for (k = 0; k < NumberOfInputs;    k++) IIValues[k] = (int)pow(SourceAlphabetSize, k);

    for (i = 0; i < ActualNumberOfStates; i++) {
        Statet = i;
        for (k = ConstrainedLength - 1; k > 0; k--) {
            SS[ConstrainedLength - 1 - k] = Statet / SSValues[k];
            Statet %= SSValues[k];
        }
        SS[ConstrainedLength - 1] = Statet;

        for (j = 0; j < InputAtTime_t_Max; j++) {
            IIStatet = j;
            for (o = 0; o < ConstrainedLength; o++) UpdatedSS[o] = 0;
            for (k = NumberOfInputs - 1; k > 0; k--) {
                InputValue[k] = IIStatet / IIValues[k];
                IIStatet %= IIValues[k];
            }
            InputValue[0] = IIStatet;

            /* Compute next state */
            for (o = 0; o < ConstrainedLength; o++) {
                for (l = 0; l < ConstrainedLength; l++)
                    UpdatedSS[o] = (UpdatedSS[o] + RecursionPolynomial[o][l] * SS[l]) % SourceAlphabetSize;
                for (l = 0; l < NumberOfInputs; l++)
                    UpdatedSS[o] = (UpdatedSS[o] + InputPolynomial[l][o] * InputValue[l]) % SourceAlphabetSize;
            }
            TrellisNextState[i][j] = 0;
            for (k = 0; k < ConstrainedLength; k++)
                TrellisNextState[i][j] += UpdatedSS[k] * SSValues[ConstrainedLength - 1 - k];

            /* Compute output symbol */
            TrellisOutput[i][j] = 0;
            for (k = 0; k < NumberOfOutputs - NumberOfRedundancyOutputs; k++) {
                Intermediate = 0;
                for (l = 0; l < ConstrainedLength; l++)
                    Intermediate = (Intermediate + FeedForwardPolynomial[k + NumberOfRedundancyOutputs][l] * SS[l]) % SourceAlphabetSize;
                for (l = 0; l < NumberOfInputs; l++)
                    Intermediate = (Intermediate + FeedForwardPolynomial[k + NumberOfRedundancyOutputs][ConstrainedLength + 1 + l] * InputValue[l]) % SourceAlphabetSize;
                TrellisOutput[i][j] += Intermediate * OOValues[k + NumberOfRedundancyOutputs];
            }
            for (k = 0; k < NumberOfRedundancyOutputs; k++) {
                Intermediate2 = 0;
                for (l = 0; l < ConstrainedLength + 1; l++) {
                    int SSVal = (l == 0) ? InputValue[0] : SS[l - 1];
                    Intermediate2 = (Intermediate2 + FeedForwardPolynomial[k][ConstrainedLength - l] * SSVal) % SourceAlphabetSize;
                }
                TrellisOutput[i][j] += Intermediate2 * OOValues[k];
            }
        }
    }
}


/* =========================================================================
 * Encoder functions
 * ========================================================================= */

void EncodingInputData(int Damierr, int NumInputs, int SourceAlphabetSize,
                       int StartingState, int *EndingState, int Sizee,
                       int InformationData[], int DataToBeTransmitted[],
                       int FindTrellisNextState[MAX_NUMBER_OF_CODE_STATES][MAX_INPUT_STATES],
                       int FindTrellisOutput[MAX_NUMBER_OF_CODE_STATES][MAX_INPUT_STATES],
                       int framesize)
{
    (void)Damierr; /* Damier path not used for the LTE encoder */
    (void)framesize;

    int i, k;
    int NextStates[Sizee + 1];
    int OValues[MAX_NUMBER_OF_INPUTS];
    int CodedCurrentSymbol;
    int CodingSize = Sizee / NumInputs;

    for (k = 0; k < NumInputs; k++) OValues[k] = (int)pow(SourceAlphabetSize, k);
    NextStates[0] = StartingState;

    for (i = 1; i < CodingSize + 1; i++) {
        CodedCurrentSymbol = 0;
        for (k = 0; k < NumInputs; k++)
            CodedCurrentSymbol += OValues[k] * InformationData[NumInputs * (i - 1) + k];
        NextStates[i] = FindTrellisNextState[NextStates[i - 1]][CodedCurrentSymbol];
        DataToBeTransmitted[i - 1] = FindTrellisOutput[NextStates[i - 1]][CodedCurrentSymbol];
    }
    *EndingState = NextStates[CodingSize];
}

/**
 * Encode() — circular encoding for the tail-biting distance computation.
 *
 * Pass 1: run the trellis from state 0 to find the ending state.
 * Pass 2: look up the circular starting state and re-encode to get the
 *         actual parity output.
 */
int Encode(int InformationData[], int DataToBeTransmitted[], int framesize)
{
    int i, k, Count, distance;
    int NextStates[framesize + 1];
    int OValues[MAX_NUMBER_OF_INPUTS];
    int CodedCurrentSymbol;

    const int NumInputs              = 1;
    const int NumberOfRedundancyBits = 1;
    const int SourceAlphabetSize     = 2;

    for (k = 0; k < NumInputs; k++) OValues[k] = (int)pow(SourceAlphabetSize, k);

    /* Pass 1: determine ending state */
    NextStates[0] = 0;
    for (i = 1; i < framesize + 1; i++) {
        CodedCurrentSymbol = 0;
        for (k = 0; k < NumInputs; k++)
            CodedCurrentSymbol += OValues[k] * InformationData[NumInputs * (i - 1) + k];
        NextStates[i] = TrellisNextState[NextStates[i - 1]][CodedCurrentSymbol];
    }

    /* Pass 2: encode from the circular starting state */
    NextStates[0] = CircularStates[NextStates[framesize]];
    Count = distance = 0;
    for (i = 1; i < framesize + 1; i++) {
        CodedCurrentSymbol = 0;
        for (k = 0; k < NumInputs; k++)
            CodedCurrentSymbol += OValues[k] * InformationData[NumInputs * (i - 1) + k];
        NextStates[i] = TrellisNextState[NextStates[i - 1]][CodedCurrentSymbol];
        for (k = 0; k < NumberOfRedundancyBits; k++) {
            DataToBeTransmitted[Count] =
                BitValue_ForTheDistance[k][TrellisOutput[NextStates[i - 1]][CodedCurrentSymbol]];
            distance += DataToBeTransmitted[Count];
            Count++;
        }
    }
    return distance;
}

void GenerateRandomInput(int InformationData[], int TailleTrame, int AlphabetSize)
{
    for (int j = 0; j < TailleTrame; j++)
        InformationData[j] = rand() % AlphabetSize;
}

/**
 * Marginalisation — builds the bit-value and distance-value lookup tables.
 *
 * BitValue_ForTheDistance[b][s] = value of bit b in symbol s.
 * Distance_Value_For_The_Bit[b][v][k] = k-th symbol index where bit b = v.
 */
void Marginalisation(int Bit_Value_For_The_Distance[MAX_BPS][MAX_DISTANCE],
                     int Distance_Value_For_The_Bit[MAX_BPS][2][MAX_DISTANCE/2],
                     int BitPerSymbol)
{
    int TotalSlots = (int)pow(2, BitPerSymbol);
    for (int j = 0; j < BitPerSymbol; j++) {
        int stride = (int)pow(2, j);
        int n = TotalSlots / (2 * stride);
        int o = 0, p = 0;
        for (int m = 0; m < n; m++) {
            for (int l = 0; l < stride; l++) {
                Bit_Value_For_The_Distance[j][2 * m * stride + l] = 0;
                Distance_Value_For_The_Bit[j][0][o++] = 2 * m * stride + l;
            }
            for (int l = stride; l < 2 * stride; l++) {
                Bit_Value_For_The_Distance[j][2 * m * stride + l] = 1;
                Distance_Value_For_The_Bit[j][1][p++] = m * 2 * stride + l;
            }
        }
    }
}
