// ****************************************************************************
// This file is part of VocalTractLab.
// Copyright (C) 2020, Peter Birkholz, Dresden, Germany
// www.vocaltractlab.de
// author: Peter Birkholz
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.
//
// ****************************************************************************

// ****************************************************************************
// This file defines the entry point for the DLL application, and the functions
// defined here are C-compatible so that they can be used with the 
// MATLAB and Python shared library interface.
// ****************************************************************************

// Make an extern "C" section so that the functions can be accessed from Matlab

#ifdef __cplusplus
extern "C"{ /* start extern "C" */
#endif

// Definition for function export, if the file is compiled as part of a dll.

#ifdef WIN32
  #ifdef _USRDLL
    #define C_EXPORT __declspec(dllexport)
  #else
    #define C_EXPORT
  #endif        // DLL
#else
  #define C_EXPORT
#endif  // WIN32

#include <stdbool.h>

// ****************************************************************************
// The exported C-compatible functions.
// IMPORTANT: 
// All the functions defined below must be named in the VocalTractLabApi.def 
// file in the project folder, so that they are usable from MATLAB !!!
// ****************************************************************************

// ****************************************************************************
// Init. the synthesis with the given speaker file name, e.g. "JD3.speaker".
// This function should be called before any other function of this API.
// Return values:
// 0: success.
// 1: Loading the speaker file failed.
// ****************************************************************************

C_EXPORT int vtlInitialize(const char *speakerFileName);

// ****************************************************************************
// Save the current speaker information (vocal tract and glottis shape) in
// a speaker file (e.g., "JD3.speaker")
// Return values:
// 0: success.
// 1: Saving the speaker file failed.
// ****************************************************************************

C_EXPORT int vtlSaveSpeaker(const char* speakerFileName);

// ****************************************************************************
// Clean up the memory and shut down the synthesizer.
// Return values:
// 0: success.
// 1: The API was not initialized.
// ****************************************************************************

C_EXPORT int vtlClose();


// ****************************************************************************
// Switch to turn off/on the automatic calculation of the tongue root 
// parameters TRX and TRY.
//
// Return values:
// 0: success.
// 1: The API was not initialized.
// ****************************************************************************

C_EXPORT int vtlCalcTongueRootAutomatically(bool automaticCalculation);


// ****************************************************************************
// Returns the version of this API as a string that contains the compile data.
// Reserve at least 32 chars for the string.
// ****************************************************************************

C_EXPORT void vtlGetVersion(char *version);


// ****************************************************************************
// Returns a couple of constants:
// o The audio sampling rate of the synthesized signal.
// o The number of supraglottal tube sections.
// o The number of vocal tract model parameters.
// o The number of glottis model parameters.
// o The number of audio samples that correspond to a single tract state sample
// o The number of evaluated tract state samples per second (internal sampling rate).
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// ****************************************************************************

C_EXPORT int vtlGetConstants(int *audioSamplingRate, int *numTubeSections,
  int *numVocalTractParams, int *numGlottisParams, int *numAudioSamplesPerTractState, double *internalSamplingRate);


// ****************************************************************************
// Returns for each supra glottal parameter the minimum value, the maximum value,
// and the standard (default) value. Each array passed to this function must have at 
// least as many elements as the number of supra glottal parameters.
// The "names" string receives the names of the parameters separated
// by tabs. This string should have at least 10*numParams elements.
// The "descriptions" string receives the descriptions of the parameters separated
// by tabs. This string should have at least 100*numParams elements.
// The "units" string receives the names of the parameter units separated
// by tabs. This string should have at least 10*numParams elements.
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// ****************************************************************************

C_EXPORT int vtlGetTractParamInfo(char* names, char* descriptions, char* units,
    double* paramMin, double* paramMax, double* paramStandard);


// ****************************************************************************
// Returns for each glottis model parameter the minimum value, the maximum value,
// and the standard (default) value. Each array passed to this function must have at 
// least as many elements as the number of glottis model parameters.
// The "names" string receives the names of the parameters separated
// by tabs. This string should have at least 10*numParams elements.
// The "descriptions" string receives the descriptions of the parameters separated
// by tabs. This string should have at least 100*numParams elements.
// The "units" string receives the names of the parameter units separated
// by tabs. This string should have at least 10*numParams elements.
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// ****************************************************************************

C_EXPORT int vtlGetGlottisParamInfo(char* names, char* descriptions, char* units,
    double* paramMin, double* paramMax, double* paramStandard);


// ****************************************************************************
// Returns the sub-glottal parameters for the given shape as defined in the
// speaker file.
// The array passed to this function must have at least as many elements as 
// the number of glottis model parameters.
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// 2: A shape with the given name does not exist.
// ****************************************************************************

C_EXPORT int vtlGetGlottisParams(const char* shapeName, double *glottisParams);


// ****************************************************************************
// Returns the supra-glottal parameters for the given shape as defined in the
// speaker file.
// The array passed to this function must have at least as many elements as 
// the number of vocal tract model parameters.
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// 2: A shape with the given name does not exist.
// ****************************************************************************

C_EXPORT int vtlGetTractParams(const char *shapeName, double *tractParams);


// ****************************************************************************
// Exports the vocal tract contours for the given vector of vocal tract
// parameters as a SVG file (scalable vector graphics).
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// 2: Writing the SVG file failed.
// ****************************************************************************

C_EXPORT int vtlExportTractSvg(double *tractParams, const char *fileName);


// ****************************************************************************
// Provides the tube data (especially the area function) for the given vector
// of tractParams. The vectors tubeLength_cm, tubeArea_cm2, and tubeArticulator, 
// must each have as many elements as tube sections.
// The values incisorPos_cm, tongueTipSideElevation, and velumOpening_cm2 are 
// one double value each.
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// ****************************************************************************

C_EXPORT int vtlTractToTube(double* tractParams,
  double* tubeLength_cm, double* tubeArea_cm2, int* tubeArticulator,
  double* incisorPos_cm, double* tongueTipSideElevation, double* velumOpening_cm2);


// ****************************************************************************
// Provides the tube data (especially the area function) for the given vector
// of tractParams. This version does NOT store and restore the previous vocal
// tract state. That means it is at least twice as fast as vtlTractToTube.
// However, this breaks incremental synthesis. That means you should not call
// this method during synthesis via the methods vtlSynthesisAddTube or
// vtlSynthesisAddTract (otherwise the current tract state will be changed).
// It has no negative impact on all other synthesis methods.
// 
// The vectors tubeLength_cm, tubeArea_cm2, and tubeArticulator, 
// must each have as many elements as tube sections.
// The values incisorPos_cm, tongueTipSideElevation, and velumOpening_cm2 are 
// one double value each.
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// ****************************************************************************

C_EXPORT int vtlFastTractToTube(double* tractParams,
  double* tubeLength_cm, double* tubeArea_cm2, int* tubeArticulator,
  double* incisorPos_cm, double* tongueTipSideElevation, double* velumOpening_cm2);


// ****************************************************************************
// Enumerates the different options to model radiation of the sound wave
// from the mouth
// ****************************************************************************

typedef enum
{
    NO_RADIATION,
    PISTONINSPHERE_RADIATION,
    PISTONINWALL_RADIATION,
    PARALLEL_RADIATION,
    NUM_RADIATION_OPTIONS
} RadiationType;


// ****************************************************************************
// Enumerates the different types of spectra (or transfer functions)
// ****************************************************************************

typedef enum
{
    SPECTRUM_UU,  // Output flow over input flow
    SPECTRUM_PU   // Output output pressure over input flow
} SpectrumType;


// ****************************************************************************
// A struct containing the various options available for the calculation 
// of the vocal tract transfer function
// ****************************************************************************

typedef struct
{
    SpectrumType spectrumType;      // What kind of transfer function to calculate
    RadiationType radiationType;        // Radiation model
    bool boundaryLayer;             // Consider boundary layer resistance
    bool heatConduction;            // Consider heat conduction losses
    bool softWalls;                 // Consider soft walls
    bool hagenResistance;           // Consider Hagen-Poiseuille resistance
    bool innerLengthCorrections;    // Make inner (tube) length corrections
    bool lumpedElements;            // Use lumped elements in T-sections
    bool paranasalSinuses;          // Include the paranasal sinuses
    bool piriformFossa;             // Include the piriform fossa
    bool staticPressureDrops;       // Consider static pressure drops
} TransferFunctionOptions;


// ****************************************************************************
// Returns the default options for the transfer function calculation. 
// 
// Parameters out:
// o opts: A struct containing the default values for the options available for
// the transfer function calculation.
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// ****************************************************************************
C_EXPORT int vtlGetDefaultTransferFunctionOptions(TransferFunctionOptions* opts);


// ****************************************************************************
// Calculates the volume velocity transfer function of the vocal tract between 
// the glottis and the lips for the given vector of vocal tract parameters and
// returns the spectrum in terms of magnitude and phase.
//
// Parameters in:
// o tractParams: Is a vector of vocal tract parameters with 
//     numVocalTractParams elements.
// o numSpectrumSamples: The number of samples (points) in the requested 
//     spectrum. This number of samples includes the negative frequencies and
//     also determines the frequency spacing of the returned magnitude and
//     phase vectors. The frequency spacing is 
//     deltaFreq = SAMPLING_RATE / numSpectrumSamples.
//     For example, with the sampling rate of 44100 Hz and 
//     numSpectrumSamples = 512, the returned magnitude and phase values are 
//     at the frequencies 0.0, 86.13, 172.3, ... Hz.
//     The value of numSpectrumSamples should not be greater than 16384,
//     otherwise the returned spectrum will be bandlimited to below 10 kHz.
// o opts: The options to use for the transfer function calculation. If NULL 
//     is passed, the default options will be used (see 
//     vtlGetDefaultTransferFunctionOptions()).
// 
// Parameters out:
// o magnitude: Vector of spectral magnitudes at equally spaced discrete 
//     frequencies. This vector mus have at least numSpectrumSamples elements.
// o phase_rad: Vector of the spectral phase in radians at equally 
//     spaced discrete frequencies. This vector must have at least 
//     numSpectrumSamples elements.
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// ****************************************************************************

C_EXPORT int vtlGetTransferFunction(double* tractParams, int numSpectrumSamples,
    TransferFunctionOptions* opts, double* magnitude, double* phase_rad);


// ****************************************************************************
// Calculates the real limited tract params (the ones that are actually used
// in the synthesis) from a given arbitrary set of tract parameters
//
// Parameters:
// o inTractParams (in): Is a vector of vocal tract parameters with 
//     numVocalTractParams elements.
// o outTractParams (out): Is a vector of vocal tract parameters with 
//     numVocalTractParams elements.
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// ****************************************************************************

C_EXPORT int vtlInputTractToLimitedTract(double* inTractParams, double* outTractParams);


// ****************************************************************************
// Resets the time-domain synthesis of continuous speech (using the functions
// vtlSynthesisAddTube() or vtlSynthesisAddTract()). This function must be 
// called every time you start a new synthesis.
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// ****************************************************************************

C_EXPORT int vtlSynthesisReset();


// ****************************************************************************
// Synthesize a part of a speech signal with numNewSamples samples, during 
// which the vocal tract tube changes linearly from the tube shape passed to
// the previous call of this function to the tube shape passed to this call.
// To synthesize parts of 5 ms duration, call this function with 
// numNewSamples = 220. The synthesized signal part is written to the array 
// audio (the caller must allocate the memory for the array).
// During the *first* call of this function after vtlSynthesisReset(), no audio
// is synthesized, and numNewSamples should be 0. During the first call, only 
// the initial tube state is set.
//
// The new tube state is given in terms of the following parameters:
// o tubeLength_cm: Vector of tube sections lengths from the glottis (index 0)
//     to the mouth (index numTubeSections; see vtlGetConstants()).
// o tubeArea_cm2: According vector of tube section areas in cm^2.
// o tubeArticulator: Vector of characters (letters) that denote the articulator 
//     that confines the vocal tract at the position of the tube. We discriminate
//     1 (tongue), 2 (lower incisors), 3 (lower lip), 4 (other articulator).
// o incisorPos_cm: Position of the incisors from the glottis.
// o velumOpening_cm2: Opening of the velo-pharyngeal port in cm^2.
// o tongueTipSideElevation: Corresponds to the TS3 parameter of the vocal tract.
// o newGlottisParams: vector with parameters of the glottis model.
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// 2: Number of generated audio samples is wrong (may happen when 
//    numNewSamples != 0 during the first call of this function after reset).
// ****************************************************************************

C_EXPORT int vtlSynthesisAddTube(int numNewSamples, double *audio,
  double *tubeLength_cm, double *tubeArea_cm2, int *tubeArticulator,
  double incisorPos_cm, double velumOpening_cm2, double tongueTipSideElevation,
  double *newGlottisParams);


// ****************************************************************************
// Synthesize a part of a speech signal with numNewSamples samples, during 
// which the vocal tract changes linearly from the tract shape passed to
// the previous call of this function to the tract shape passed to this call.
// To synthesize parts of 5 ms duration, call this function with 
// numNewSamples = 220. The synthesized signal part is written to the array 
// audio (the caller must allocate the memory for the array).
// During the *first* call of this function after vtlSynthesisReset(), no audio
// is synthesized, and numNewSamples should be 0. During the first call, only 
// the initial tube state is set.
//
// The new vocal tract state is given in terms of the following parameters:
// o tractParams: Vector of vocal tract parameters.
// o glottisParams: Vector of vocal fold model parameters.
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// 2: Number of generated audio samples is wrong (may happen when 
//    numNewSamples != 0 during the first call of this function after reset).
// ****************************************************************************

C_EXPORT int vtlSynthesisAddTract(int numNewSamples, double *audio,
  double *tractParams, double *glottisParams);


// ****************************************************************************
// Synthesize speech with a given sequence of vocal tract model states and 
// glottis model states, and return the corresponding audio signal.
// This function makes successive calls to the function vtlSynthesisAddTract().
//
// Parameters (in/out):
// o tractParams (in): Is a concatenation of vocal tract parameter vectors
//     with the total length of (numVocalTractParams*numFrames) elements.
// o glottisParams (in): Is a concatenation of glottis parameter vectors
//     with the total length of (numGlottisParams*numFrames) elements.
// o numFrames (in): Number of successive states of the glottis and vocal tract
//     that are going to be concatenated.
// o frameStep_samples (in): The number of audio samples between adjacent 
//     frames (states). A typical value is 220, which corresponds to 5 ms.
// o audio (out): The resulting audio signal with sample values in the range 
//     [-1, +1] and with the sampling rate audioSamplingRate. The signal
//     will have (numFrames-1) * frameStep_samples samples, so the array must
//     be at least of this size.
// o enableConsoleOutput (in): Set to 1, if you want to allow output about the
//   synthesis progress in the console window. Otherwise, set it to 0.
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// ****************************************************************************

C_EXPORT int vtlSynthBlock(double *tractParams, double *glottisParams,
  int numFrames, int frameStep_samples, double *audio, bool enableConsoleOutput);


// ****************************************************************************
// Test function for this API.
// Audio should contain at least 44100 double values.
// Run this WITHOUT calling vtlInitialize() !
// ****************************************************************************

C_EXPORT int vtlApiTest(const char *speakerFileName, double *audio, int *numSamples);


// ****************************************************************************
// This function converts a segment sequence file (a TXT file containing the 
// sequence of speech segments in SAMPA and the associated durations) with the 
// name segFileName into a gestural score file (gesFileName).
// The f0 tier in the gestural score is set to a "standard" f0.
//
// o enableConsoleOutput (in): Set to 1, if you want to allow output about the
//   synthesis progress in the console window. Otherwise, set it to 0.
//
// Function return value:
// 0: success.
// 1: The API was not initialized.
// 2: Loading the segment sequence file failed.
// 3: Saving the gestural score file failed.
// ****************************************************************************

C_EXPORT int vtlSegmentSequenceToGesturalScore(const char *segFileName, const char *gesFileName, bool enableConsoleOutput);


// ****************************************************************************
// This function directly converts a gestural score to an audio signal or file.
//
// Parameters:
// o gesFileName (in): Name of the gestural score file to synthesize.
// o wavFileName (in): Name of the audio file with the resulting speech signal.
//     This can be the empty string "" if you do not want to save a WAV file.
// o audio (out): The resulting audio signal with sample values in the range 
//     [-1, +1] and with the sampling rate audioSamplingRate. Make sure that
//     this buffer is big enough for the synthesized signal. If you are not 
//     interested in the audio signal, set this pointer to NULL.
// o numSamples (out): The number of audio samples in the synthesized signal.
//     If you are not interested in this value, set this pointer to NULL.
// o enableConsoleOutput (in): Set to 1, if you want to allow output about the
//   synthesis progress in the console window. Otherwise, set it to 0.
//
// Function return value:
// 0: success.
// 1: The API was not initialized.
// 2: Loading the gestural score file failed.
// 3: Values in the gestural score file are out of range.
// 4: The WAV file could not be saved.
// ****************************************************************************

C_EXPORT int vtlGesturalScoreToAudio(const char* gesFileName, const char* wavFileName,
  double* audio, int* numSamples, bool enableConsoleOutput);


// ****************************************************************************
// This function directly converts a gestural score to a tract sequence file.
// The latter is a text file containing the parameters of the vocal fold and 
// vocal tract models in steps of 5 ms.
//
// Parameters:
// o gesFileName (in): Name of the gestural score file to convert.
// o tractSequenceFileName (in): Name of the tract sequence file.
//
// Function return value:
// 0: success.
// 1: The API was not initialized.
// 2: Loading the gestural score file failed.
// 3: Values in the gestural score file are out of range.
// 4: The tract sequence file could not be saved.
// ****************************************************************************

C_EXPORT int vtlGesturalScoreToTractSequence(const char* gesFileName, 
  const char* tractSequenceFileName);


// ****************************************************************************
// This function gets the duration from a gestural score.
//
// Parameters:
// o gesFileName (in): Name of the gestural score file.
// o audioFileDuration (out): The number of audio samples, the audio file would
//   have, if the gestural score was synthesized. This number can be slightly 
//   larger than the length of the gestural score because the audio is 
//   synthesized in chunks of a constant size. If not wanted, set to NULL.
// o gesFileDuration (out): The duration of the gestural score (in samples).
//   If not wanted, set to NULL.
//
// Function return value:
// 0: success.
// 1: The API was not initialized.
// 2: Loading the gestural score file failed.
// 3: Values in the gestural score file are out of range.
// ****************************************************************************

C_EXPORT int vtlGetGesturalScoreDuration(const char* gesFileName, int* numAudioSamples,
    int* numGestureSamples);


// ****************************************************************************
// This function converts a tract sequence file into an audio signal or file.
//
// Parameters:
// o tractSequenceFileName (in): Name of the tract sequence file to synthesize.
// o wavFileName (in): Name of the audio file with the resulting speech signal.
//     This can be the empty string "" if you do not want to save a WAV file.
// o audio (out): The resulting audio signal with sample values in the range 
//     [-1, +1] and with the sampling rate audioSamplingRate. Make sure that
//     this buffer is big enough for the synthesized signal. If you are not 
//     interested in the audio signal, set this pointer to NULL.
// o numSamples (out): The number of audio samples in the synthesized signal.
//     If you are not interested in this value, set this pointer to NULL.
//
// Function return value:
// 0: success.
// 1: The API was not initialized.
// 2: Synthesis of the tract sequence file failed.
// 3: The WAV file could not be saved.
// ****************************************************************************

C_EXPORT int vtlTractSequenceToAudio(const char* tractSequenceFileName,
  const char* wavFileName, double* audio, int* numSamples);


// ****************************************************************************
// This function calculates and exports EMA points.
//
// Parameters:
// o gestureFileName: Name of the gestural score file to synthesize.
// o emaFileName: Name of a text file to which the sequence of EMA point Coordinates
//      and other feedback data will be written.
//
// The return value is 0 if successful, and otherwise an error code >= 1.
// Error codes:
// 0: success.
// 1: The API was not initialized.
// 2: Loading the gestural score file failed.
// 3: Values in the gestural score files are out of range.
// ****************************************************************************

C_EXPORT int vtlGesturalScoreToEma(const char *gestureFileName, const char *emaFileName);


// ****************************************************************************
// Calculate and export selected EMA points and mesh data with a given sequence
// of vocal tract model states and glottis model states. For each frame in
// the incoming model, a 3D-mesh of the vocal tract is calculated and exported
// in an .obj-file and a corresponding .mtl file. The files' names consist of
// the name handed in "fileName" and the number of the current frame. These
// files are stored in a subfolder "fileName-meshes" of the given directory
// "filePath". The EMA points are exported into a .txt file named
// "fileName-ema" in the directory "filePath". It contains the time of the
// frame and the 3D-coordinates of all selected EMA points per frame in a row.
//
// Parameters in:
// o tractParams: Is a concatenation of vocal tract parameter vectors
//     with the total length of (numVocalTractParams*numFrames) elements.
// o glottisParams: Is a concatenation of glottis parameter vectors
//     with the total length of (numGlottisParams*numFrames) elements.
// o numTractParams: length of tractParams
// o numGlottisParams: length of glottisParams
// o numFrames: number of frames the model data contain
// o numEmaPoints: number of selected EMA points
// o surf: Array with indices of surfaces of selected EMA points
//            (UPPER_TEETH = 0, LOWER_TEETH = 1, UPPER_COVER = 2,
//             LOWER_COVER = 3, UPPER_LIP = 4, LOWER_LIP = 5,
//             PALATE = 6, MANDIBLE = 7, LOWER_TEETH_ORIGINAL = 8,
//             LOW_VELUM = 9, MID_VELUM = 10, HIGH_VELUM = 11,
//             NARROW_LARYNX_FRONT = 12, NARROW_LARYNX_BACK = 13,
//             WIDE_LARYNX_FRONT = 14, WIDE_LARYNX_BACK = 15,
//             TONGUE = 16, UPPER_COVER_TWOSIDE = 17,
//             LOWER_COVER_TWOSIDE = 18, UPPER_TEETH_TWOSIDE = 19,
//             LOWER_TEETH_TWOSIDE = 20, UPPER_LIP_TWOSIDE = 21,
//             LOWER_LIP_TWOSIDE = 22, LEFT_COVER = 23,
//             RIGHT_COVER = 24, UVULA_ORIGINAL = 25, UVULA = 26,
//             UVULA_TWOSIDE = 27, EPIGLOTTIS_ORIGINAL = 28,
//             EPIGLOTTIS = 29, EPIGLOTTIS_TWOSIDE = 30)
// o vert: Array with index of vertices of selected EMA points
//            (Predefined EMA point vertex indices:
//             Tongue Back (TB) = 115
//             Tongue Middle (TM) = 225
//             Tongue Tip (TT) = 335
//             Upper Lip (UL) = 89
//             Lower Lip (LL) = 89
//             Lower Cover (JAW) = 148)
// o filePath: path leading to the directory where EMA and mesh files shall be stored.
// o fileName: name of all exported datasets
//
// The return value is 0 if successful, and otherwise an error code >= 1.
// Error codes:
// 0: success.
// 1: numEmaPoints <= 0
// 2: surf == NULL
// 3: vert == NULL
// 4: surface index of a selected EmaPoints exceeds 30
// 5: vertex index of a selected EmaPoint exceeds possible range
// 6: filePath is not valid
// 7: mesh folder already exist: prevents overwriting
// 8: EMA file already exists: prevents overwriting
// 9: EMA file could not be opened
// 10: API has not been initialized
// ****************************************************************************

C_EXPORT int vtlTractSequenceToEmaAndMesh(double *tractParams, double *glottisParams, int numTractParams, int numGlottisParams, int numFrames, int numEmaPoints, int *surf, int *vert, const char *filePath, const char *fileName);


// ****************************************************************************
// This function calculates a sequence of vocal tract model states and
// glottis model states from a gestural score and hands them over to vtlTractandGlottisToEma().
//
// Parameters in:
// o gestureFileName: Name of the gestural score file to synthesize.
// o filePath: path leading to the directory where EMA and mesh files shall be stored.
// o fileName: name of all exported datasets
//
// The return value is 0 if successful, and otherwise an error code >= 1.
// Error codes:
// 0: success.
// 1: The API was not initialized.
// 2: Loading the gestural score file failed.
// 3: Values in the gestural score files are out of range.
// ****************************************************************************

C_EXPORT int vtlGesturalScoreToEmaAndMesh(const char *gestureFileName, const char *filePath, const char *fileName);

// ****************************************************************************

// ****************************************************************************
// Set control parameters on the current glottis model, call calcGeometry(),
// and return derived parameters and tube data.
//
// Parameters:
// o controlParams (in): Control parameter values (numGlottisParams elements).
// o derivedParams (out): Derived parameter values after calcGeometry().
//     Must have at least 8 elements for the geometric glottis.
// o numDerivedParams (out): Number of derived parameters written.
// o tubeLength_cm (out): Tube section lengths (2 elements).
// o tubeArea_cm2 (out): Tube section areas (2 elements).
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// ****************************************************************************

C_EXPORT int vtlGlottisCalcGeometry(double *controlParams,
                                    double *derivedParams,
                                    int *numDerivedParams,
                                    double *tubeLength_cm,
                                    double *tubeArea_cm2);

// ****************************************************************************
// Perform one incTime step on the current glottis model, then call
// calcGeometry() and return the updated state.
//
// Parameters:
// o timeIncrement_s (in): Time step in seconds.
// o pressure_dPa (in): Pressure values (4 elements):
//     [subglottal, lower_glottis, upper_glottis, supraglottal].
// o controlParams (in): Control parameter values for this step.
// o derivedParams (out): Derived parameters after the step.
// o numDerivedParams (out): Number of derived parameters written.
// o tubeLength_cm (out): Tube section lengths (2 elements).
// o tubeArea_cm2 (out): Tube section areas (2 elements).
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// ****************************************************************************

C_EXPORT int vtlGlottisIncTime(double timeIncrement_s, double *pressure_dPa,
                               double *controlParams, double *derivedParams,
                               int *numDerivedParams, double *tubeLength_cm,
                               double *tubeArea_cm2);

// ****************************************************************************
// Reset the motion state of the current glottis model (phase, time, filters).
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// ****************************************************************************

C_EXPORT int vtlGlottisResetMotion();

// ****************************************************************************
// Returns static parameter info for the current glottis model.
// The arrays must have at least numStaticParams elements.
//
// Parameters out:
// o names: Tab-separated parameter names (at least 10*numStaticParams chars).
// o paramMin, paramMax, paramStandard: Min, max, default values.
// o numStaticParams: Number of static parameters.
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// ****************************************************************************

C_EXPORT int vtlGetGlottisStaticParamInfo(char *names, double *paramMin,
                                          double *paramMax,
                                          double *paramStandard,
                                          int *numStaticParams);

// ****************************************************************************
// TDS (Time-Domain Simulation) Component Testing API
// ****************************************************************************

// ****************************************************************************
// Set TDS options (e.g., disable noise sources for deterministic testing).
//
// Parameters:
// o generateNoiseSources (in): Enable/disable noise source generation.
// o turbulenceLosses (in): Consider fluid dynamic losses due to turbulence.
// o softWalls (in): Consider losses due to soft walls.
// o radiationFromSkin (in): Allow sound radiation from the skin.
// o piriformFossa (in): Include the piriform fossa.
// o innerLengthCorrections (in): Additional inductivities between sections.
// o transvelarCoupling (in): Sound transmission through the velum tissue.
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// ****************************************************************************

C_EXPORT int vtlTdsSetOptions(
    bool generateNoiseSources,
    bool turbulenceLosses,
    bool softWalls,
    bool radiationFromSkin,
    bool piriformFossa,
    bool innerLengthCorrections,
    bool transvelarCoupling
);

// ****************************************************************************
// Reset the TDS model motion state.
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// ****************************************************************************

C_EXPORT int vtlTdsResetMotion();

// ****************************************************************************
// Set the piriform fossa dimensions on the global tube object.
// This allows testing the TDS/Synthesizer with non-default fossa values
// through synthesis_add_tube (which normally uses the global tube's defaults).
//
// Parameters (in):
// o length_cm: Fossa length in cm.
// o volume_cm3: Fossa volume in cm^3.
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// ****************************************************************************

C_EXPORT int vtlSetFossaDims(double length_cm, double volume_cm3);

// ****************************************************************************
// Set tube geometry on the TDS model and run one time step.
// Returns all internal state for component-level testing.
//
// Parameters (in):
// o tubeLength_cm: Pharynx+mouth section lengths (40 elements).
// o tubeArea_cm2: Pharynx+mouth section areas (40 elements).
// o tubeArticulator: Pharynx+mouth articulators (40 elements).
// o incisorPos_cm: Position of the incisors.
// o velumOpening_cm2: Naso-pharyngeal port area.
// o tongueTipSideElevation: TS3 parameter.
// o filtering: Apply area filtering (true/false).
// o pressureSourceSection: Section index for pressure source (-1 = none).
// o pressureSourceAmp: Pressure source amplitude in dPa.
//
// Parameters (out):
// o secArea: Area per section (93 elements).
// o secLength: Length per section (93 elements).
// o secR0: Left resistance per section (93 elements).
// o secR1: Right resistance per section (93 elements).
// o secL: Inductance per section (93 elements).
// o secC: Capacitance per section (93 elements).
// o secD: D value per section (93 elements).
// o secE: E value per section (93 elements).
// o secAlpha: Alpha (wall vibration) per section (93 elements).
// o secBeta: Beta (wall vibration) per section (93 elements).
// o secPressure: Pressure per section (93 elements).
// o bcMagnitude: Branch current magnitude (97 elements).
// o mouthFlow: Radiated flow from mouth (scalar).
// o nostrilFlow: Radiated flow from nostrils (scalar).
// o skinFlow: Radiated flow from skin (scalar).
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// ****************************************************************************

C_EXPORT int vtlTdsSetTubeAndRun(
    double *tubeLength_cm, double *tubeArea_cm2, int *tubeArticulator,
    double incisorPos_cm, double velumOpening_cm2,
    double tongueTipSideElevation, bool filtering, int pressureSourceSection,
    double pressureSourceAmp, double *secArea, double *secLength, double *secR0,
    double *secR1, double *secL, double *secC, double *secD, double *secE,
    double *secAlpha, double *secBeta, double *secPressure, double *bcMagnitude,
    double *mouthFlow, double *nostrilFlow, double *skinFlow);

// ****************************************************************************

// ****************************************************************************
// Returns all 93 tube sections (area, length, volume, wall properties, etc.)
// for the given vocal tract parameters.
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// ****************************************************************************

C_EXPORT int vtlTractToFullTube(double *tractParams,
                                double *tubeLength_cm,
                                double *tubeArea_cm2,
                                double *tubeVolume_cm3,
                                double *tubeWallMass_cgs,
                                double *tubeWallStiffness_cgs,
                                double *tubeWallResistance_cgs,
                                int *tubeArticulator,
                                double *incisorPos_cm,
                                double *tongueTipSideElevation,
                                double *velumOpening_cm2,
                                double *piriformFossaLength_cm,
                                double *piriformFossaVolume_cm3);

// ****************************************************************************
// Returns intermediate TL model values (matrixProduct for all 93 sections,
// fossa input impedance, radiation impedances) at a given frequency index.
// Used for analysis and debugging of the TL model computation.
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// 2: freqIndex out of range.
// ****************************************************************************

C_EXPORT int vtlGetTLIntermediateValues(
    double *tractParams,
    int numSpectrumSamples,
    TransferFunctionOptions *opts,
    int freqIndex,
    double *matrix_A_re, double *matrix_A_im,
    double *matrix_B_re, double *matrix_B_im,
    double *matrix_C_re, double *matrix_C_im,
    double *matrix_D_re, double *matrix_D_im,
    double *fossa_input_imp_re, double *fossa_input_imp_im,
    double *nose_rad_imp_re, double *nose_rad_imp_im,
    double *mouth_rad_imp_re, double *mouth_rad_imp_im);

// ****************************************************************************
// Returns the 129 cross-section areas, positions, and articulators from the
// VocalTract model after calling calculateAll on the given tract parameters.
// NUM_CENTERLINE_POINTS = (1 << 7) + 1 = 129.
//
// Parameters:
// o tractParams (in): Is a vector of vocal tract parameters with
//     numVocalTractParams elements.
// o crossSectionAreas (out): Is a vector of 129 cross-section areas in cm^2.
// o crossSectionPositions (out): Is a vector of 129 cross-section positions
//     along the center line in cm.
// o crossSectionArticulators (out): Is a vector of 129 articulator indices
//     (int cast of Tube::Articulator enum).
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// ****************************************************************************

C_EXPORT int vtlGetCrossSections(double *tractParams,
                                 double *crossSectionAreas,
                                 double *crossSectionPositions,
                                 int *crossSectionArticulators);

// ****************************************************************************
// Returns the upper and lower cross-sectional profiles at a specific
// centerline index for the given vocal tract parameters.
//
// Parameters:
// o tractParams (in): Is a vector of vocal tract parameters with
//     numVocalTractParams elements.
// o centerlineIndex (in): Index along the centerline (0..128).
// o upperProfile (out): Is a vector of 96 (NUM_PROFILE_SAMPLES) doubles
//     representing the upper profile.
// o lowerProfile (out): Is a vector of 96 (NUM_PROFILE_SAMPLES) doubles
//     representing the lower profile.
// o centerlineInfo (out): Is a vector of 6 doubles:
//     [point.x, point.y, normal.x, normal.y, area, pos].
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// 2: The centerline index is out of range.
// ****************************************************************************

C_EXPORT int vtlGetProfiles(double *tractParams, int centerlineIndex,
                             double *upperProfile, double *lowerProfile,
                             double *centerlineInfo);

// ****************************************************************************
// Returns all 129 centerline points (x, y, normal_x, normal_y, pos) after
// computing the vocal tract geometry for the given parameters.
//
// Parameters:
// o tractParams (in): Is a vector of vocal tract parameters with
//     numVocalTractParams elements.
// o centerlineData (out): Is a vector of 129*5 = 645 doubles.
//     For each point i (0..128): [x, y, normal_x, normal_y, pos].
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// ****************************************************************************

C_EXPORT int vtlGetCenterline(double *tractParams, double *centerlineData);

// ****************************************************************************
// Returns the 4 outlines (upper, lower, tongue, epiglottis) used for
// centerline computation, after calculating vocal tract geometry.
//
// Parameters:
// o tractParams (in): Vocal tract parameters (numVocalTractParams elements).
// o outlineData (out): Flat array for 4 outlines, each point = (x, y).
//     Layout: [upper_n, upper_x0, upper_y0, ..., lower_n, lower_x0, ...]
//     Max total size: 4 * (1 + 2*200) = 1604 doubles.
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// ****************************************************************************

C_EXPORT int vtlGetOutlines(double *tractParams, double *outlineData, int *outlineSizes);

C_EXPORT int vtlGetTongueRibData(double *tractParams, double *ribData, int *numRibs);

C_EXPORT int vtlGetTongueWidthBounds(double *tractParams, double *boundsData, int *numRibs);

C_EXPORT int vtlGetSurfaceVertices(double *tractParams, int surfaceIndex,
    double *vertexData, int *numRibs, int *numRibPoints);

C_EXPORT int vtlGetCuts(double *tractParams, int centerlineIndex,
    double *cutData, int *numCuts);

// ****************************************************************************
// Apply anatomy parameters derived from age and gender to the loaded vocal
// tract. This calls AnatomyParams::calcFromAge(), restrictParams(), and
// setFor() on the currently loaded VocalTract.
//
// Parameters:
// o ageMonths (in): Age in months (minimum 12).
// o isMale (in): true for male, false for female.
//
// Return values:
// 0: success.
// 1: The API has not been initialized.
// ****************************************************************************

C_EXPORT int vtlSetAnatomyFromAge(int ageMonths, bool isMale);

// ****************************************************************************
// Get the 13 anatomy parameters from the currently loaded vocal tract.
//
// Parameters:
// o anatomyParams (out): Array of 13 doubles to receive the values.
//
// Return values:
// 0: success.
// 1: The API has not been initialized.
// ****************************************************************************

C_EXPORT int vtlGetAnatomyParams(double *anatomyParams);

// ****************************************************************************
// Set the 13 anatomy parameters on the currently loaded vocal tract.
// This calls AnatomyParams::restrictParams() and setFor().
//
// Parameters:
// o anatomyParams (in): Array of 13 doubles with the anatomy values.
//
// Return values:
// 0: success.
// 1: The API has not been initialized.
// ****************************************************************************

C_EXPORT int vtlSetAnatomyParams(double *anatomyParams);

// ****************************************************************************

#ifdef __cplusplus
} /* end extern "C" */
#endif

