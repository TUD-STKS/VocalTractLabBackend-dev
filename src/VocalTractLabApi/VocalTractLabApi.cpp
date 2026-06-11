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

#include "VocalTractLabApi.h"
#include "VocalTractLabBackend/Dsp.h"
#include "VocalTractLabBackend/AudioFile.h"
#include "VocalTractLabBackend/Synthesizer.h"
#include "VocalTractLabBackend/SegmentSequence.h"

#include "VocalTractLabBackend/GeometricGlottis.h"
#include "VocalTractLabBackend/TwoMassModel.h"
#include "VocalTractLabBackend/TriangularGlottis.h"

#include "VocalTractLabBackend/VocalTract.h"
#include "VocalTractLabBackend/TdsModel.h"
#include "VocalTractLabBackend/GesturalScore.h"

#include "VocalTractLabBackend/Speaker.h"
#include "VocalTractLabBackend/AnatomyParams.h"

#include "VocalTractLabBackend/TlModel.h"

#include <iostream>
#include <fstream>
#include <iomanip>  // for 'setprecision'
#include <filesystem>
namespace fs = std::filesystem;


enum GlottisModel
{
  GEOMETRIC_GLOTTIS,
  TWO_MASS_MODEL,
  TRIANGULAR_GLOTTIS,
  NUM_GLOTTIS_MODELS
};

static Glottis *glottis[NUM_GLOTTIS_MODELS];
static int selectedGlottis;

static VocalTract *vocalTract = NULL;
static TdsModel *tdsModel = NULL;
static Synthesizer *synthesizer = NULL;
static Tube *tube = NULL;

static bool vtlApiInitialized = false;


#if defined(WIN32) && defined(_USRDLL) 

// ****************************************************************************
/// Windows entry point for the DLL.
// ****************************************************************************

// Windows Header Files
#include <windows.h>

BOOL APIENTRY DllMain( HANDLE hModule, DWORD  ul_reason_for_call, LPVOID lpReserved)
{
  switch (ul_reason_for_call)
	{
		case DLL_PROCESS_ATTACH:
		case DLL_THREAD_ATTACH:
		case DLL_THREAD_DETACH:
		case DLL_PROCESS_DETACH:
	  break;
  }
  return TRUE;
}

#endif  // WIN32 && _USRDLL


// ****************************************************************************
// Loads the VT anatomy and the configurations for the different glottis 
// models from a speaker file.
// This function is not visible in the interface.
// ****************************************************************************

bool vtlLoadSpeaker(const char *speakerFileName, VocalTract *vocalTract, 
  Glottis *glottis[], int &selectedGlottis)
{
    //// ****************************************************************
	//// Load the speaker from XML file
	//// ****************************************************************
    Speaker speaker(speakerFileName);

    *vocalTract = *speaker.getVocalTract();
    const auto glottisInfo = speaker.getGlottisModels();
    int i = 0;
	for (const auto& g : glottisInfo)
    {
        glottis[i++] = g;
	}
    selectedGlottis = speaker.getSelectedGlottis();


  return true;
}



// ****************************************************************************
// Init. the synthesis with the given speaker file name, e.g. "JD3.speaker".
// This function should be called before any other function of this API.
// Return values:
// 0: success.
// 1: Loading the speaker file failed.
// ****************************************************************************

int vtlInitialize(const char *speakerFileName)
{
  if (vtlApiInitialized)
  {
    vtlClose();
  }

  // ****************************************************************
  // Init the vocal tract.
  // ****************************************************************

  vocalTract = new VocalTract();
  vocalTract->calculateAll();

  // ****************************************************************
  // Init the list with glottis models
  // ****************************************************************

  glottis[GEOMETRIC_GLOTTIS] = new GeometricGlottis();
  glottis[TWO_MASS_MODEL] = new TwoMassModel();
  glottis[TRIANGULAR_GLOTTIS] = new TriangularGlottis();
  
  selectedGlottis = GEOMETRIC_GLOTTIS;

  bool ok = vtlLoadSpeaker(speakerFileName, vocalTract, glottis, selectedGlottis);

  if (ok == false)
  {
    int i;
    for (i = 0; i < NUM_GLOTTIS_MODELS; i++)
    {
      delete glottis[i];
    }
    delete vocalTract;

    printf("Error in vtlInitialize(): vtlLoadSpeaker() failed.\n");
    return 1;
  }

  // ****************************************************************
  // Init the object for the time domain simulation.
  // ****************************************************************

  tdsModel = new TdsModel();

  // ****************************************************************
  // Init the Synthesizer object.
  // ****************************************************************

  synthesizer = new Synthesizer();
  synthesizer->init(glottis[selectedGlottis], vocalTract, tdsModel);

  tube = new Tube();

  // We are now initialized!
  vtlApiInitialized = true;

  return 0;
}


// ****************************************************************************
// Save the current speaker information (vocal tract and glottis shape) in
// a speaker file (e.g., "JD3.speaker")
// Return values:
// 0: success.
// 1: Saving the speaker file failed.
// ****************************************************************************

int vtlSaveSpeaker(const char* speakerFileName)
{
    Speaker speaker(vocalTract, { std::begin(glottis), std::end(glottis) }, selectedGlottis);
    try
    {
        speaker.save(speakerFileName);
        return 0;
    }
    catch (std::exception&)
    {
        return 1;
    }
}


// ****************************************************************************
// Clean up the memory and shut down the synthesizer.
// Return values:
// 0: success.
// 1: The API was not initialized.
// ****************************************************************************

int vtlClose()
{
  if (!vtlApiInitialized)
  {
    printf("Error: The API was not initialized.\n");
    return 1;
  }

  delete tube;
  delete synthesizer;
  delete tdsModel;

  int i;
  for (i = 0; i < NUM_GLOTTIS_MODELS; i++)
  {
    delete glottis[i];
  }

  delete vocalTract;

  vtlApiInitialized = false;
  
  return 0;
}


// ****************************************************************************
// Switch to turn off/on the automatic calculation of the tongue root 
// parameters TRX and TRY.
//
// Return values:
// 0: success.
// 1: The API was not initialized.
// ****************************************************************************

int vtlCalcTongueRootAutomatically(bool automaticCalculation)
{
    if (!vtlApiInitialized)
    {
        printf("Error: The API was not initialized.\n");
        return 1;
    }

    vocalTract->anatomy.automaticTongueRootCalc = automaticCalculation;
    vocalTract->calculateAll();

    return 0;
}


// ****************************************************************************
// Returns the version of this API as a string that contains the compile data.
// Reserve at least 64 chars for the string.
// ****************************************************************************

void vtlGetVersion(char *version)
{
  strcpy(version, "API 2.4.2 ");
  strcat(version, __DATE__);
}


// ****************************************************************************
// Returns a couple of constants:
// o The audio sampling rate of the synthesized signal.
// o The number of supraglottal tube sections.
// o The number of vocal tract model parameters.
// o The number of glottis model parameters.
// o The number of audio samples per tract state.
// o The internal sampling rate.
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// ****************************************************************************

int vtlGetConstants(int *audioSamplingRate, int *numTubeSections,
                    int *numVocalTractParams, int *numGlottisParams,
                    int *numAudioSamplesPerTractState, double *internalSamplingRate)
{
  if (!vtlApiInitialized)
  {
    printf("Error: The API has not been initialized.\n");
    return 1;
  }

  *audioSamplingRate = SAMPLING_RATE;
  *numTubeSections = Tube::NUM_PHARYNX_MOUTH_SECTIONS;
  *numVocalTractParams = VocalTract::NUM_PARAMS;
  *numGlottisParams = (int)glottis[selectedGlottis]->controlParam.size();
  *numAudioSamplesPerTractState = Synthesizer::NUM_CHUNCK_SAMPLES;
  *internalSamplingRate = (double)SAMPLING_RATE / (double)Synthesizer::NUM_CHUNCK_SAMPLES;


  return 0;
}


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

int vtlGetTractParamInfo(char* names, char* descriptions, char* units,
    double* paramMin, double* paramMax, double* paramStandard)
{
    if (!vtlApiInitialized)
    {
        printf("Error: The API has not been initialized.\n");
        return 1;
    }

    int i;

    strcpy(names, "");
    strcpy(descriptions, "");
    strcpy(units, "");

    for (i = 0; i < VocalTract::NUM_PARAMS; i++)
    {
        strcat(names, vocalTract->param[i].name.c_str());
        strcat(descriptions, vocalTract->param[i].description.c_str());
        strcat(units, vocalTract->param[i].unit.c_str());
        if (i != VocalTract::NUM_PARAMS - 1)
        {
            strcat(names, "\t");
            strcat(descriptions, "\t");
            strcat(units, "\t");
        }

        paramMin[i] = vocalTract->param[i].min;
        paramMax[i] = vocalTract->param[i].max;
        paramStandard[i] = vocalTract->param[i].neutral;
    }

    return 0;
}


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

int vtlGetGlottisParamInfo(char* names, char* descriptions, char* units,
    double* paramMin, double* paramMax, double* paramStandard)
{
    if (!vtlApiInitialized)
    {
        printf("Error: The API has not been initialized.\n");
        return 1;
    }

    int i;
    int numGlottisParams = (int)glottis[selectedGlottis]->controlParam.size();

    strcpy(names, "");
    strcpy(descriptions, "");
    strcpy(units, "");

    for (i = 0; i < numGlottisParams; i++)
    {
        strcat(names, glottis[selectedGlottis]->controlParam[i].name.c_str());
        strcat(descriptions, glottis[selectedGlottis]->controlParam[i].description.c_str());
        strcat(units, glottis[selectedGlottis]->controlParam[i].cgsUnit.c_str());
        if (i != numGlottisParams - 1)
        {
            strcat(names, "\t");
            strcat(descriptions, "\t");
            strcat(units, "\t");
        }

        paramMin[i] = glottis[selectedGlottis]->controlParam[i].min;
        paramMax[i] = glottis[selectedGlottis]->controlParam[i].max;
        paramStandard[i] = glottis[selectedGlottis]->controlParam[i].neutral;
    }

    return 0;
}


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

int vtlGetGlottisParams(const char *shapeName, double *glottisParams)
{
    if (!vtlApiInitialized)
    {
        printf("Error: The API has not been initialized.\n");
        return 1;
    }

    int index = glottis[selectedGlottis]->getShapeIndex(string(shapeName));
    if (index == -1)
    {
        return 2;
    }

    int numGlottisParams = (int)glottis[selectedGlottis]->controlParam.size();
    int i;
    for (i = 0; i < numGlottisParams; i++)
    {
        glottisParams[i] = glottis[selectedGlottis]->shape[index].controlParam[i];
    }

    return 0;
}



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

int vtlGetTractParams(const char *shapeName, double *tractParams)
{
    if (!vtlApiInitialized)
    {
        printf("Error: The API has not been initialized.\n");
        return 1;
    }

    int index = vocalTract->getShapeIndex(string(shapeName));
    if (index == -1)
    {
        return 2;
    }

    int i;
    for (i = 0; i < VocalTract::NUM_PARAMS; i++)
    {
        tractParams[i] = vocalTract->shapes[index].param[i];
    }

    return 0;
}


// ****************************************************************************
// Exports the vocal tract contours for the given vector of vocal tract
// parameters as a SVG file (scalable vector graphics).
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// 2: Writing the SVG file failed.
// ****************************************************************************

int vtlExportTractSvg(double *tractParams, const char *fileName)
{
  if (!vtlApiInitialized)
  {
    printf("Error: The API has not been initialized.\n");
    return 1;
  }

  // Store the current control parameter values.
  vocalTract->storeControlParams();

  // Set the given vocal tract parameters.
  int i;
  for (i = 0; i < VocalTract::NUM_PARAMS; i++)
  {
    vocalTract->param[i].x = tractParams[i];
  }
  vocalTract->calculateAll();
  
  // Save the contour as SVG file.
  bool ok = vocalTract->exportTractContourSvg(string(fileName), false, false);
  
  // Restore the previous control parameter values and 
  // recalculate the vocal tract shape.

  vocalTract->restoreControlParams();
  vocalTract->calculateAll();

  if (ok)
  {
    return 0;
  }
  else
  {
    return 2;
  }
}


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

int vtlTractToTube(double *tractParams,
  double *tubeLength_cm, double *tubeArea_cm2, int *tubeArticulator,
  double *incisorPos_cm, double *tongueTipSideElevation, double *velumOpening_cm2)
{
  if (!vtlApiInitialized)
  {
    printf("Error: The API has not been initialized.\n");
    return 1;
  }

  // ****************************************************************
  // Store the current control parameter values.
  // ****************************************************************

  vocalTract->storeControlParams();

  // ****************************************************************
  // Set the given vocal tract parameters.
  // ****************************************************************

  int i;
  for (i = 0; i < VocalTract::NUM_PARAMS; i++)
  {
    vocalTract->param[i].x = tractParams[i];
  }

  // ****************************************************************
  // Get the tube for the new vocal tract shape.
  // ****************************************************************

  Tube tube;
  vocalTract->calculateAll();
  vocalTract->getTube(&tube);

  // ****************************************************************
  // Copy the tube parameters to the user arrays.
  // ****************************************************************

  Tube::Section *ts = NULL;
  for (i = 0; i < Tube::NUM_PHARYNX_MOUTH_SECTIONS; i++)
  {
    ts = &tube.pharynxMouthSection[i];

    tubeLength_cm[i] = ts->length_cm;
    tubeArea_cm2[i] = ts->area_cm2;
    tubeArticulator[i] = ts->articulator;
  }

  *incisorPos_cm = tube.teethPosition_cm;
  *tongueTipSideElevation = tube.tongueTipSideElevation;
  *velumOpening_cm2 = tube.getVelumOpening_cm2();

  // ****************************************************************
  // Restore the previous control parameter values and 
  // recalculate the vocal tract shape.
  // ****************************************************************

  vocalTract->restoreControlParams();

  return 0;
}


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

int vtlFastTractToTube(double *tractParams,
  double *tubeLength_cm, double *tubeArea_cm2, int *tubeArticulator,
  double *incisorPos_cm, double *tongueTipSideElevation, double *velumOpening_cm2)
{
  if (!vtlApiInitialized)
  {
    printf("Error: The API has not been initialized.\n");
    return 1;
  }

  // ****************************************************************
  // Do not store the current control parameter values.
  // ****************************************************************



  // ****************************************************************
  // Set the given vocal tract parameters.
  // ****************************************************************

  int i;
  for (i = 0; i < VocalTract::NUM_PARAMS; i++)
  {
    vocalTract->param[i].x = tractParams[i];
  }

  // ****************************************************************
  // Get the tube for the new vocal tract shape.
  // ****************************************************************

  Tube tube;
  vocalTract->calculateAll();
  vocalTract->getTube(&tube);

  // ****************************************************************
  // Copy the tube parameters to the user arrays.
  // ****************************************************************

  Tube::Section *ts = NULL;
  for (i = 0; i < Tube::NUM_PHARYNX_MOUTH_SECTIONS; i++)
  {
    ts = &tube.pharynxMouthSection[i];

    tubeLength_cm[i] = ts->length_cm;
    tubeArea_cm2[i] = ts->area_cm2;
    tubeArticulator[i] = ts->articulator;
  }

  *incisorPos_cm = tube.teethPosition_cm;
  *tongueTipSideElevation = tube.tongueTipSideElevation;
  *velumOpening_cm2 = tube.getVelumOpening_cm2();

  // ****************************************************************
  // Do not restore the previous control parameter values. This way
  // calculateAll is only called once and the processing time is 
  // reduced by 50%.
  // ****************************************************************

  return 0;
}


// ****************************************************************************
// Returns all 93 tube sections (area, length, volume, wall properties, etc.)
// for the given vocal tract parameters. Unlike vtlTractToTube which only
// returns pharynx+mouth sections, this returns the full tube including
// trachea, glottis, nose, fossa, and sinus sections.
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// ****************************************************************************

int vtlTractToFullTube(double *tractParams,
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
                       double *piriformFossaVolume_cm3)
{
  if (!vtlApiInitialized)
  {
    printf("Error: The API has not been initialized.\n");
    return 1;
  }

  int i;
  for (i = 0; i < VocalTract::NUM_PARAMS; i++)
  {
    vocalTract->param[i].x = tractParams[i];
  }

  vocalTract->calculateAll();

  Tube tube;
  vocalTract->getTube(&tube);
  tube.setGlottisArea(0.0);

  // Copy all 93 sections
  for (i = 0; i < Tube::NUM_SECTIONS; i++)
  {
    Tube::Section *ts = tube.section[i];
    tubeLength_cm[i] = ts->length_cm;
    tubeArea_cm2[i] = ts->area_cm2;
    tubeVolume_cm3[i] = ts->volume_cm3;
    tubeWallMass_cgs[i] = ts->wallMass_cgs;
    tubeWallStiffness_cgs[i] = ts->wallStiffness_cgs;
    tubeWallResistance_cgs[i] = ts->wallResistance_cgs;
    tubeArticulator[i] = ts->articulator;
  }

  *incisorPos_cm = tube.teethPosition_cm;
  *tongueTipSideElevation = tube.tongueTipSideElevation;
  *velumOpening_cm2 = tube.getVelumOpening_cm2();

  // Get static tube dimensions (fossa length and volume)
  double subglottalLength, nasalLength, fossaLength, fossaVolume;
  tube.getStaticTubeDimensions(subglottalLength, nasalLength, fossaLength, fossaVolume);

  *piriformFossaLength_cm = fossaLength;
  *piriformFossaVolume_cm3 = fossaVolume;

  return 0;
}

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

int vtlGetDefaultTransferFunctionOptions(TransferFunctionOptions* opts)
{
    opts->radiationType = PARALLEL_RADIATION;
    opts->boundaryLayer = true;
    opts->heatConduction = false;
    opts->softWalls = true;
    opts->hagenResistance = false;
    opts->lumpedElements = true;
    opts->innerLengthCorrections = false;
    opts->paranasalSinuses = true;
    opts->piriformFossa = true;
    opts->staticPressureDrops = true;
    opts->spectrumType = SPECTRUM_UU;
    return 0;
}

// ****************************************************************************
// Calculates the transfer function of the vocal tract between 
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

int vtlGetTransferFunction(double* tractParams, int numSpectrumSamples, TransferFunctionOptions* opts, double* magnitude, double* phase_rad)
{
    if (!vtlApiInitialized)
    {
        printf("Error: The API has not been initialized.\n");
        return 1;
    }

    int i;
    ComplexSignal s;

    if (numSpectrumSamples < 16)
    {
        numSpectrumSamples = 16;
    }

    // Calculate the vocal tract shape from the vocal tract parameters.

    for (i = 0; i < VocalTract::NUM_PARAMS; i++)
    {
        vocalTract->param[i].x = tractParams[i];
    }
    vocalTract->calculateAll();

    // Calculate the transfer function.

    TlModel* tlModel = new TlModel();

    // Set the options.
    // tfOpts is declared outside the if-block below, because it must
    // still be alive when opts is dereferenced further down.
    TransferFunctionOptions tfOpts;
    TlModel::Options tlOpts;
    if (opts == NULL)
    {
      vtlGetDefaultTransferFunctionOptions(&tfOpts);
      opts = &tfOpts;
    }
    tlOpts.boundaryLayer = opts->boundaryLayer;
    tlOpts.hagenResistance = opts->hagenResistance;
    tlOpts.heatConduction = opts->heatConduction;
    tlOpts.innerLengthCorrections = opts->innerLengthCorrections;
    tlOpts.lumpedElements = opts->lumpedElements;
    tlOpts.paranasalSinuses = opts->paranasalSinuses;
    tlOpts.piriformFossa = opts->piriformFossa;
    tlOpts.radiation = (TlModel::RadiationType)opts->radiationType;
    tlOpts.softWalls = opts->softWalls;
    tlOpts.staticPressureDrops = opts->staticPressureDrops;

    tlModel->options = tlOpts;
    vocalTract->getTube(&tlModel->tube);
    tlModel->tube.setGlottisArea(0.0);

    tlModel->getSpectrum(TlModel::FLOW_SOURCE_TF, &s, numSpectrumSamples, Tube::FIRST_PHARYNX_SECTION);

    if (opts->spectrumType == SPECTRUM_PU)
    {
        ComplexSignal radiationSpectrum(0);
        tlModel->getSpectrum(TlModel::RADIATION, &radiationSpectrum, numSpectrumSamples, 0);
        s *= radiationSpectrum;
    }

    // Separate the transfer function into magnitude and phase.
    for (i = 0; i < numSpectrumSamples; i++)
    {
        magnitude[i] = s.getMagnitude(i);
        phase_rad[i] = s.getPhase(i);
    }

    delete tlModel;

    return 0;
}


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

int vtlInputTractToLimitedTract(double* inTractParams, double* outTractParams)
{
    if (!vtlApiInitialized)
    {
        printf("Error: The API has not been initialized.\n");
        return 1;
    }

    // Calculate the vocal tract shape from the vocal tract parameters.
    int i;
    for (i = 0; i < VocalTract::NUM_PARAMS; i++)
    {
        vocalTract->param[i].x = inTractParams[i];
    }
    vocalTract->calculateAll();

    for (i = 0; i < VocalTract::NUM_PARAMS; i++)
    {
        outTractParams[i] = vocalTract->param[i].limitedX;
    }

    return 0;
}


// ****************************************************************************
// Resets the time-domain synthesis of continuous speech (using the functions
// vtlSynthesisAddTube() or vtlSynthesisAddTract()). This function must be 
// called every time you start a new synthesis.
//
// Function return value:
// 0: success.
// 1: The API has not been initialized.
// ****************************************************************************

int vtlSynthesisReset()
{
  if (!vtlApiInitialized)
  {
    printf("Error: The API has not been initialized.\n");
    return 1;
  }

  synthesizer->reset();
  tube->resetDynamicPart();

  return 0;
}


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

int vtlSynthesisAddTube(int numNewSamples, double* audio,
  double* tubeLength_cm, double* tubeArea_cm2, int* tubeArticulator,
  double incisorPos_cm, double velumOpening_cm2, double tongueTipSideElevation,
  double* newGlottisParams)
{
  if (!vtlApiInitialized)
  {
    printf("Error: The API has not been initialized.\n");
    return 1;
  }

  Tube::Articulator articulator[Tube::NUM_PHARYNX_MOUTH_SECTIONS];
  int i;
  for (i = 0; i < Tube::NUM_PHARYNX_MOUTH_SECTIONS; i++)
  {
    articulator[i] = (Tube::Articulator)tubeArticulator[i];
  }

  // Set the properties of the target tube.

  tube->setPharynxMouthGeometry(tubeLength_cm, tubeArea_cm2, articulator, 
    incisorPos_cm, tongueTipSideElevation);
  tube->setVelumOpening(velumOpening_cm2);
  // The aspiration strength will be set based on the glottis parameters
  // in synthesizer->add(...) below.
  tube->setAspirationStrength(0.0);

  // Synthesize the speech signal part.

  vector<double> audioVector;
  synthesizer->add(newGlottisParams, tube, numNewSamples, audioVector);

  if ((int)audioVector.size() != numNewSamples)
  {
    printf("Error in vtlSynthesisAddTube(): Number of audio samples is wrong.\n");
    return 2;
  }

  // Copy the audio samples in the given buffer.

  for (i = 0; i < numNewSamples; i++)
  {
    audio[i] = audioVector[i];
  }

  return 0;
}


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

int vtlSynthesisAddTract(int numNewSamples, double *audio,
  double *tractParams, double *glottisParams)
{
  if (!vtlApiInitialized)
  {
    printf("Error: The API has not been initialized.\n");
    return 1;
  }

  vector<double> audioVector;
  synthesizer->add(glottisParams, tractParams, numNewSamples, audioVector);

  if ((int)audioVector.size() != numNewSamples)
  {
    printf("Error in vtlSynthesisAddTube(): Number of audio samples is wrong.\n");
    return 2;
  }

  // Copy the audio samples in the given buffer.

  int i;
  for (i = 0; i < numNewSamples; i++)
  {
    audio[i] = audioVector[i];
  }

  return 0;
}


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

int vtlSynthBlock(double *tractParams, double *glottisParams,
  int numFrames, int frameStep_samples, double *audio, bool enableConsoleOutput)
{
  if (!vtlApiInitialized)
  {
    printf("Error: The API has not been initialized.\n");
    return 1;
  }

  int i;
  int samplePos = 0;
  int numGlottisParams = (int)glottis[selectedGlottis]->controlParam.size();

  if (enableConsoleOutput)
  {
    printf("Block synthesis in progress ...");
  }

  vtlSynthesisReset();

  for (i = 0; i < numFrames; i++)
  {
    if (i == 0)
    {
      // Only set the initial state of the vocal tract and glottis without generating audio.
      vtlSynthesisAddTract(0, &audio[0],
        &tractParams[i*VocalTract::NUM_PARAMS], &glottisParams[i*numGlottisParams]);
    }
    else
    {
      vtlSynthesisAddTract(frameStep_samples, &audio[samplePos],
        &tractParams[i*VocalTract::NUM_PARAMS], &glottisParams[i*numGlottisParams]);
      samplePos += frameStep_samples;
    }

    if ((enableConsoleOutput != 0) && ((i % 20) == 0))
    {
      printf(".");
    }
  }

  if (enableConsoleOutput != 0)
  {
    printf(" finished\n");
  }

  return 0;
}


// ****************************************************************************
// Test function for this API.
// Audio should contain at least 44100 double values.
// Run this WITHOUT calling vtlInitialize() !
// ****************************************************************************

int vtlApiTest(const char *speakerFileName, double *audio, int *numSamples)
{
  int failed = vtlInitialize(speakerFileName);
  if (failed != 0)
  {
    printf("Error in  in vtlApiTest(): vtlInitialize() failed.\n");
    return 1;
  }

  char version[100];
  vtlGetVersion(version);
  printf("Compile date of the library: %s\n", version);

  int audioSamplingRate = -1;
  int numTubeSections = -1;
  int numVocalTractParams = -1;
  int numGlottisParams = -1;
  int numAudioSamplesPerTractState = -1;
  double internalSamplingRate = -1.0;

  vtlGetConstants(&audioSamplingRate, &numTubeSections, &numVocalTractParams, &numGlottisParams, &numAudioSamplesPerTractState, &internalSamplingRate);

  printf("Audio sampling rate = %d\n", audioSamplingRate);
  printf("Num. of tube sections = %d\n", numTubeSections);
  printf("Num. of vocal tract parameters = %d\n", numVocalTractParams);
  printf("Num. of glottis parameters = %d\n", numGlottisParams);

  char tractParamNames[50 * 32];
  char tractParamDescriptions[500 * 32];
  char tractParamUnits[50 * 32];
  double tractParamMin[50];
  double tractParamMax[50];
  double tractParamStandard[50];

  vtlGetTractParamInfo(tractParamNames, tractParamDescriptions, tractParamUnits, tractParamMin, tractParamMax, tractParamStandard);

  char glottisParamNames[50 * 32];
  char glottisParamDescriptions[500 * 32];
  char glottisParamUnits[50 * 32];
  double glottisParamMin[50];
  double glottisParamMax[50];
  double glottisParamStandard[50];

  vtlGetGlottisParamInfo(glottisParamNames, glottisParamDescriptions, glottisParamUnits, glottisParamMin, glottisParamMax, glottisParamStandard);

  // ****************************************************************
  // Define two target tube shapes: one for /a/ and one for /i/.
  // ****************************************************************

  const int MAX_TUBES = 100;

  // These parameters are the same for both /i/ and /a/:
  double incisorPos_cm = 15.0;
  double velumOpening_cm2 = 0.0;
  double tongueTipSideElevation = 0.0;

  int i;

  // ****************************************************************
  // Define the tube for /i/.
  // ****************************************************************

  double tubeLength_cm_i[MAX_TUBES];
  double tubeArea_cm2_i[MAX_TUBES];
  int tubeArticulator_i[MAX_TUBES];

  for (i = 0; i < numTubeSections; i++)
  {
    // Full tube length is 16 cm.
    tubeLength_cm_i[i] = 16.0 / (double)numTubeSections;
    
    // Articulator is always the tongue (although not fully correct here)
    tubeArticulator_i[i] = 1;   // = tongue
    
    // Narrow mouth sections and wide pharynx sections
    if (i < numTubeSections / 2)
    {
      tubeArea_cm2_i[i] = 8.0;
    }
    else
    {
      tubeArea_cm2_i[i] = 2.0;
    }
  }

  // ****************************************************************
  // Define the tube for /a/.
  // ****************************************************************

  double tubeLength_cm_a[MAX_TUBES];
  double tubeArea_cm2_a[MAX_TUBES];
  int tubeArticulator_a[MAX_TUBES];

  for (i = 0; i < numTubeSections; i++)
  {
    // Full tube length is 16 cm.
    tubeLength_cm_a[i] = 16.0 / (double)numTubeSections;

    // Articulator is always the tongue (although not fully correct here)
    tubeArticulator_a[i] = 1;   // = tongue

    // Narrow mouth sections and wide pharynx sections
    if (i < numTubeSections / 2)
    {
      tubeArea_cm2_a[i] = 0.3;
    }
    else
    {
      tubeArea_cm2_a[i] = 8.0;
    }
  }

  // ****************************************************************
  // Set glottis parameters to default (neutral) values, which are
  // suitable for phonation.
  // ****************************************************************

  double glottisParams[Glottis::MAX_CONTROL_PARAMS];

  for (i = 0; i < numGlottisParams; i++)
  {
    glottisParams[i] = glottisParamStandard[i];
  }

  // **************************************************************************
  // Synthesize a transition from /a/ to /i/ to /a/.
  // **************************************************************************

  int numTotalSamples = 0;
  int numNewSamples = 0;

  vtlSynthesisReset();

  // Initialize with /a/ at 120 Hz.

  glottisParams[0] = 120.0;   // 120 Hz F0
  glottisParams[1] = 0.0;     // P_sub = 0 dPa.
  vtlSynthesisAddTube(0, audio, tubeLength_cm_a, tubeArea_cm2_a, tubeArticulator_a,
    incisorPos_cm, velumOpening_cm2, tongueTipSideElevation, 
    glottisParams);

  // Make 0.2 s transition to /i/ at 100 Hz.

  glottisParams[0] = 100.0;   // 100 Hz F0
  glottisParams[1] = 8000.0;  // P_sub = 8000 dPa.
  numNewSamples = (int)(0.2*audioSamplingRate);
  printf("Adding %d samples...\n", numNewSamples);

  vtlSynthesisAddTube(numNewSamples, &audio[numTotalSamples], tubeLength_cm_i, tubeArea_cm2_i, tubeArticulator_i,
    incisorPos_cm, velumOpening_cm2, tongueTipSideElevation,
    glottisParams);
  numTotalSamples += numNewSamples;

  // Make 0.2 s transition to /a/ at 80 Hz.

  glottisParams[0] = 80.0;   // 80 Hz F0
  numNewSamples = (int)(0.2*audioSamplingRate);
  printf("Adding %d samples...\n", numNewSamples);

  vtlSynthesisAddTube(numNewSamples, &audio[numTotalSamples], tubeLength_cm_a, tubeArea_cm2_a, tubeArticulator_a,
    incisorPos_cm, velumOpening_cm2, tongueTipSideElevation, glottisParams);
  numTotalSamples += numNewSamples;

  printf("Done.\n");

  *numSamples = numTotalSamples;

  // **************************************************************************
  // Clean up and close the VTL synthesis.
  // **************************************************************************

  vtlClose();

  return 0;
}


// ****************************************************************************
// This function converts a segment sequence file (a TXT file containing the 
// sequence of speech segments in SAMPA and the associated durations) with the 
// name segFileName into a gestural score file (gesFileName).
// The f0 tier in the gestural score is set to a "standard" f0.
//
// Function return value:
// 0: success.
// 1: The API was not initialized.
// 2: Loading the segment sequence file failed.
// 3: Saving the gestural score file failed.
// ****************************************************************************

int vtlSegmentSequenceToGesturalScore(const char *segFileName, const char *gesFileName, bool enableConsoleOutput)
{
  if (!vtlApiInitialized)
  {
    printf("Error: The API has not been initialized.\n");
    return 1;
  }

  // Create and load the segment sequence file.
  
  SegmentSequence *segmentSequence = new SegmentSequence();
  if (segmentSequence->readFromFile(string(segFileName)) == false)
  {
    delete segmentSequence;
    printf("Error in vtlSegmentSequenceToGesturalScore(): Segment sequence file could not be loaded.\n");
    return 2;
  }

  // Create and save the gestural score.

  GesturalScore *gesturalScore = new GesturalScore(vocalTract, glottis[selectedGlottis]);
  gesturalScore->createFromSegmentSequence(segmentSequence, enableConsoleOutput);
  if (gesturalScore->saveGesturesXml(string(gesFileName)) == false)
  {
    delete segmentSequence;
    delete gesturalScore;
    printf("Error in vtlSegmentSequenceToGesturalScore(): Gestural score file could not be saved.\n");
    return 3;
  }

  delete segmentSequence;
  delete gesturalScore;

  return 0;
}


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

int vtlGesturalScoreToAudio(const char *gesFileName, const char *wavFileName,
  double *audio, int *numSamples, bool enableConsoleOutput)
{
  if (!vtlApiInitialized)
  {
    printf("Error: The API has not been initialized.\n");
    return 1;
  }

  int i;

  // ****************************************************************
  // Init and load the gestural score.
  // ****************************************************************

  GesturalScore *gesturalScore = new GesturalScore(vocalTract, glottis[selectedGlottis]);

  bool allValuesInRange = true;
  if (gesturalScore->loadGesturesXml(string(gesFileName), allValuesInRange) == false)
  {
    printf("Error in vtlGesturalScoreToAudio(): Loading the gestural score file failed!\n");
    delete gesturalScore;
    return 2;
  }

  if (allValuesInRange == false)
  {
    printf("Error in vtlGesturalScoreToAudio(): Some values in the gestural score are out of range!\n");
    delete gesturalScore;
    return 3;
  }

  // Important !!!
  gesturalScore->calcCurves();

  // ****************************************************************
  // Do the actual synthesis.
  // ****************************************************************

  vector<double> audioVector;
  Synthesizer::synthesizeGesturalScore(gesturalScore, tdsModel, audioVector, enableConsoleOutput);
  int numVectorSamples = (int)audioVector.size();

  // ****************************************************************
  // Copy the number of audio samples to the return value numSamples.
  // ****************************************************************

  if (numSamples != NULL)
  {
    *numSamples = numVectorSamples;
  }

  // ****************************************************************
  // Copy the synthesized signal into the return buffer audio.
  // ****************************************************************

  if (audio != NULL)
  {
    for (i = 0; i < numVectorSamples; i++)
    {
      audio[i] = audioVector[i];
    }
  }
   
  // ****************************************************************
  // Save the result as WAV file (if the name is not an empty string).
  // ****************************************************************

  if (wavFileName[0] != '\0')
  {
    AudioFile<double> audioFile;
    audioFile.setAudioBufferSize(1, numVectorSamples);
    audioFile.setBitDepth(16);
    audioFile.setSampleRate(SAMPLING_RATE);

    for (i = 0; i < numVectorSamples; i++)
    {
      audioFile.samples[0][i] = audioVector[i];
    }

    if (audioFile.save(string(wavFileName)) == false)
    {
      printf("Error in vtlGesturalScoreToAudio(): The WAV file could not be saved!\n");
      delete gesturalScore;
      return 4;
    }
  }

  // ****************************************************************
  // Free the memory and return.
  // ****************************************************************

  delete gesturalScore;
  return 0;
}


// ****************************************************************************
// This function directly converts a gestural score to a tract sequence file.
// The latter is a text file containing the parameters of the vocal fold and 
// vocal tract models in steps of about 2.5 ms.
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

int vtlGesturalScoreToTractSequence(const char* gesFileName, const char* tractSequenceFileName)
{
  if (!vtlApiInitialized)
  {
    printf("Error: The API has not been initialized.\n");
    return 1;
  }

  // ****************************************************************
  // Init and load the gestural score.
  // ****************************************************************

  GesturalScore* gesturalScore = new GesturalScore(vocalTract, glottis[selectedGlottis]);

  bool allValuesInRange = true;
  if (gesturalScore->loadGesturesXml(string(gesFileName), allValuesInRange) == false)
  {
    printf("Error in vtlGesturalScoreToTractSequence(): Loading the gestural score file failed!\n");
    delete gesturalScore;
    return 2;
  }

  if (allValuesInRange == false)
  {
    printf("Error in vtlGesturalScoreToTractSequence(): Some values in the gestural score are out of range!\n");
    delete gesturalScore;
    return 3;
  }

  // Important !!!
  gesturalScore->calcCurves();

  // ****************************************************************
  // Do the actual conversion.
  // ****************************************************************

  bool ok = Synthesizer::gesturalScoreToTractSequenceFile(gesturalScore, string(tractSequenceFileName));

  if (ok == false)
  {
    printf("Error in vtlGesturalScoreToTractSequence(): Saving the tract sequence file failed!\n");
    delete gesturalScore;
    return 4;
  }
   
  // ****************************************************************
  // Free the memory and return.
  // ****************************************************************

  delete gesturalScore;
  return 0;
}



// ****************************************************************************
// This function gets the duration from a gestural score.
//
// Parameters:
// o gesFileName (in): Name of the gestural score file.
// o numAudioSamples (out): The number of audio samples, the audio file would
//   have, if the gestural score was synthesized. This number can be slightly 
//   larger than the length of the gestural score because the audio is 
//   synthesized in chunks of a constant size. If not wanted, set to NULL.
// o numGestureSamples (out): The duration of the gestural score (in samples).
//   If not wanted, set to NULL.
//
// Function return value:
// 0: success.
// 1: The API was not initialized.
// 2: Loading the gestural score file failed.
// 3: Values in the gestural score file are out of range.
// ****************************************************************************

int vtlGetGesturalScoreDuration(const char* gesFileName, int* numAudioSamples, int* numGestureSamples)
{
    if (!vtlApiInitialized)
    {
        printf("Error: The API has not been initialized.\n");
        return 1;
    }



    // ****************************************************************
    // Init and load the gestural score.
    // ****************************************************************

    GesturalScore* gesturalScore = new GesturalScore(vocalTract, glottis[selectedGlottis]);
    static const int NUM_CHUNK_SAMPLES = 110;

    bool allValuesInRange = true;
    if (gesturalScore->loadGesturesXml(string(gesFileName), allValuesInRange) == false)
    {
        printf("Error in vtlGesturalScoreToGlottisSignals: Loading the gestural score file failed!\n");
        delete gesturalScore;
        return 2;
    }

    if (allValuesInRange == false)
    {
        printf("Error in vtlGesturalScoreToGlottisSignals: Some values in the gestural score are out of range!\n");
        delete gesturalScore;
        return 3;
    }

    // Important !!!
    gesturalScore->calcCurves();

    if (numGestureSamples != NULL)
    {
        *numGestureSamples = gesturalScore->getDuration_pt();
    }

    if (numAudioSamples != NULL)
    {
        *numAudioSamples = ( (int)( ( gesturalScore->getDuration_pt() ) / NUM_CHUNK_SAMPLES ) + 1 )  * NUM_CHUNK_SAMPLES;
    }

    // ****************************************************************
    // Free the memory and return.
    // ****************************************************************

    delete gesturalScore;
    return 0;
}


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

int vtlTractSequenceToAudio(const char* tractSequenceFileName, const char* wavFileName,
  double* audio, int* numSamples)
{
  if (!vtlApiInitialized)
  {
    printf("Error: The API has not been initialized.\n");
    return 1;
  }

  int i;
  vector<double> audioVector;

  bool ok = Synthesizer::synthesizeTractSequence(string(tractSequenceFileName),
    glottis[selectedGlottis], vocalTract, tdsModel, audioVector);

  if (ok == false)
  {
    printf("Error in vtlTractSequenceToAudio(): Synthesis of the tract sequence file failed.\n");
    return 2;
  }

  int numVectorSamples = (int)audioVector.size();

  // ****************************************************************
  // Copy the number of audio samples to the return value numSamples.
  // ****************************************************************

  if (numSamples != NULL)
  {
    *numSamples = numVectorSamples;
  }

  // ****************************************************************
  // Copy the synthesized signal into the return buffer audio.
  // ****************************************************************

  if (audio != NULL)
  {
    for (i = 0; i < numVectorSamples; i++)
    {
      audio[i] = audioVector[i];
    }
  }

  // ****************************************************************
  // Save the result as WAV file (if the file name is not empty).
  // ****************************************************************

  if (wavFileName[0] != '\0')
  {
    AudioFile<double> audioFile;
    audioFile.setAudioBufferSize(1, numVectorSamples);
    audioFile.setBitDepth(16);
    audioFile.setSampleRate(SAMPLING_RATE);

    for (i = 0; i < numVectorSamples; i++)
    {
      audioFile.samples[0][i] = audioVector[i];
    }

    if (audioFile.save(string(wavFileName)) == false)
    {
      printf("Error in vtlTractSequenceToAudio(): The WAV file could not be saved!\n");
      return 3;
    }
  }

  // ****************************************************************

  return 0;
}


// ****************************************************************************
// This function calculates and exports EMA points.
//
// Parameters:
// o gestureFileName: Name of the gestural score file to synthesize.
// o emaFileName: Name of a text file to which the sequence of EMA point coordinates
//      and other feedback data will be written.
//
// The return value is 0 if successful, and otherwise an error code >= 1.
// Error codes:
// 0: success.
// 1: The API was not initialized.
// 2: Loading the gestural score file failed.
// 3: Values in the gestural score files are out of range.
// ****************************************************************************

int vtlGesturalScoreToEma(const char *gestureFileName, const char *emaFileName)
{
  if (!vtlApiInitialized)
  {
    printf("Error: The API has not been initialized.\n");
    return 1;
  }

  int i, k;

  int numEmaPoints = (int)vocalTract->emaPoints.size();

  double t_s; // time in sec
  double tractParams[VocalTract::NUM_PARAMS];
  double glottisParams[256];

  Point3D P;
  VocalTract::EmaPoint *e;

  // ****************************************************************
  // Init and load the gestural score.
  // ****************************************************************

  bool allValuesInRange = true;

  // Create the gestural score not before we know what selectedGlottis
  // is, which is determined in loadSpeakerFile().
  GesturalScore *gesturalScore = new GesturalScore(vocalTract, glottis[selectedGlottis]);

  if (gesturalScore->loadGesturesXml(string(gestureFileName), allValuesInRange) == false)
  {
    printf("Error in vtlGesturalScoreToEma(): Loading the gestural score file failed!\n");
    return 2;
  }

  if (allValuesInRange == false)
  {
    printf("Error in vtlGesturalScoreToEma(): Some values in the gestural score are out of range!\n");
    return 3;
  }

  // Important !!!
  gesturalScore->calcCurves();

  // ****************************************************************
  // Open the text file for the EMA points.
  // ****************************************************************

  ofstream emaFile;
  if (emaFileName != NULL)
  {
    emaFile.open(emaFileName);
    if (emaFile.is_open() == false)
    {
      printf("Warning in vtlGesturalScoreToEma(): The EMA file could not be opened!\n");
    }
    else
    {
      // Write some header data into the file.

      emaFile << "#comment" << endl;

    }


    // From Frontend -> Data.cpp:

    // ****************************************************************
    // Write the Header.
    // ****************************************************************

    emaFile << "time[s] ";
    for (i=0; i < numEmaPoints; i++)
    {
      e = &vocalTract->emaPoints[i];
      emaFile << e->name << "-x[cm] " << e->name << "-y[cm] " << e->name << "-z[cm]";

      if (i != numEmaPoints-1)  // No space behind last value
      {
        emaFile << " ";
      }

    }
    emaFile << endl;

    // ****************************************************************
    // Keep in mind the current vocal tract state.
    // ****************************************************************

    double oldTractParams[VocalTract::NUM_PARAMS];
    for (i=0; i < VocalTract::NUM_PARAMS; i++)
    {
      oldTractParams[i] = vocalTract->param[i].x;
    }

    // ****************************************************************
    // Write the EMA points in the file (every 0.005 sec).
    // ****************************************************************

    const double EMA_SAMPLING_RATE_HZ = 200.0;
    int numFrames = (int)(gesturalScore->getScoreDuration_s() * EMA_SAMPLING_RATE_HZ);

    emaFile << std::setprecision(8);

    for (k=0; k < numFrames; k++)
    {
      t_s = (double)k / (double)EMA_SAMPLING_RATE_HZ;
      gesturalScore->getParams(t_s, tractParams, glottisParams);

      for (i=0; i < VocalTract::NUM_PARAMS; i++)
      {
        vocalTract->param[i].x = tractParams[i];
      }
      vocalTract->calculateAll();

      // Write data to file.

      emaFile << t_s << " ";
      for (i=0; i < numEmaPoints; i++)
      {
        P = vocalTract->getEmaPointCoord(i);
        emaFile << P.x << " " << P.y << " " << P.z;

        if (i != numEmaPoints-1)  // No space behind last value
               {
                 emaFile << " ";
               }

      }
      emaFile << endl;
    }

    // ****************************************************************
    // Set back the old vocal tract state.
    // ****************************************************************

    for (i=0; i < VocalTract::NUM_PARAMS; i++)
    {
      vocalTract->param[i].x = oldTractParams[i];
    }
    vocalTract->calculateAll();

    // ****************************************************************
    // Close the file with the EMA points.
    // ****************************************************************

    if (emaFile.is_open())
    {
      emaFile.close();
    }
  }

  // ****************************************************************
  // Free the memory.
  // ****************************************************************

  delete gesturalScore;

  return 0;
}


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

int vtlTractSequenceToEmaAndMesh(double *tractParams, double *glottisParams, int numTractParams, int numGlottisParams, int numFrames, int numEmaPoints, int *surf, int *vert, const char *filePath, const char *fileName)
{
  // Return if no EMA point is selected
  if (numEmaPoints <= 0)
  {
    printf("Error in vtlTractSequenceToEmaAndMesh(): numEmaPoints <= 0");
    return 1;
  }

  if (surf == NULL)
  {
    printf("Error in vtlTractSequenceToEmaAndMesh(): surf == NULL");
    return 2;
  }

  if (vert == NULL)
  {
    printf("Error in vtlTractSequenceToEmaAndMesh(): vert == NULL");
    return 3;
  }

  if (!vtlApiInitialized)
  {
    printf("Error: The API has not been initialized.\n");
    return 10;
  }

  // ****************************************************************
  // Create mesh and EMA file paths
  // ****************************************************************

  fs::path path = filePath;

  if (!fs::exists(path))
  {
    printf("Error in vtlTractSequenceToEmaAndMesh(): filePath doesn't exist");
    return 6;
  }

  std::string fileNameString;
  fileNameString += fileName;

  std::string emaFileName;
  emaFileName += fileName;
  emaFileName += "-ema.txt";
  fs::path emaFilePath = path;
  emaFilePath /= emaFileName;  // append file name (with ending) to path

  fs::path objFilePath = path;
  std::string meshFolderName;
  meshFolderName += fileNameString;
  meshFolderName += "-meshes";
  objFilePath /= meshFolderName;

  if (fs::exists(objFilePath))
  {
    printf("Error in vtlTractSequenceToEmaAndMesh(): Mesh-Folder already exists");
    return 7;
  }
  // create mesh-directory
  fs::create_directory(objFilePath);

  // append file name (without ending) to path
  objFilePath /= fileName;

  // stores current mesh path with complete file ending
  std::string objFileName;

  string surfaceNames[] = {"UPPER_TEETH", "LOWER_TEETH",
        "UPPER_COVER", "LOWER COVER",
        "UPPER LIP", "LOWER LIP",
        "PALATE", "MANDIBLE", "LOWER_TEETH_ORIGINAL",
        "LOW_VELUM", "MID_VELUM", "HIGH_VELUM",
        "NARROW_LARYNX_FRONT", "NARROW_LARYNX_BACK",
        "WIDE_LARYNX_FRONT", "WIDE_LARYNX_BACK",
        "TONGUE",
        "UPPER_COVER_TWOSIDE",
        "LOWER_COVER_TWOSIDE",
        "UPPER_TEETH_TWOSIDE",
        "LOWER_TEETH_TWOSIDE",
        "UPPER_LIP_TWOSIDE",
        "LOWER_LIP_TWOSIDE",
        "LEFT_COVER",
        "RIGHT_COVER",
        "UVULA_ORIGINAL",
        "UVULA",
        "UVULA_TWOSIDE",
        "EPIGLOTTIS_ORIGINAL",
        "EPIGLOTTIS",
        "EPIGLOTTIS_TWOSIDE"};

  std::vector<double> glotParams(numGlottisParams, 0.0);

  Point3D P(0.0, 0.0, 0.0);

  Surface *s = NULL;

  // Store maximal index of Vertices for each surface in surf[]
  std::vector<int> maxVertex(numEmaPoints, 0);
  int  i,e;

  //Calculate max Vertex for each surface in surf (Compare with VocalTract::getEmaSurfaceVertexRange())
  for (e = 0; e < numEmaPoints; e++)
  {
    s = &vocalTract->surface[surf[e]];
    maxVertex[e] = s->numVertices - 1;
  }

  s = NULL;

  // Check if all Vertices are in range for corresponding surface
  for (e = 0; e < numEmaPoints; e++)
  {
    if (surf[e] > 30)
    {
      printf("Error in vtlTractSequenceToEmaAndMesh(): Surface at surf[%d] out of range, maximal Surface is %d\n", e, 30);
      return 4;
    }
    if (vert[e] > maxVertex[e])
    {
      printf("Error in vtlTractSequenceToEmaAndMesh(): Vertex at vert[%d] out of range, maximal Vertex of chosen Surface %d is %d\n", e, surf[e],maxVertex[e]);
      return 5;
    }
  }

  ofstream emaFile;
  if (fs::exists(emaFilePath))
  {
    printf("Error in vtlTractSequenceToEmaAndMesh(): EMA file already exists");
    return 8;
  }

  else
  {
    emaFile.open(emaFilePath.string());
    if (emaFile.is_open() == false)
    {
      printf("Warning in vtlTractSequenceToEmaAndMesh(): The Ema file could not be opened!\n");
      return 9;

    }
    else
    {
      // ****************************************************************
      // Write the Header in EMA file
      // ****************************************************************
      emaFile << "time(s) ";

      for (e = 0; e < numEmaPoints; e++)
      {
        // Take names from string surfaceNames[], append index of vertex
        emaFile << surfaceNames[surf[e]] << "_" << vert[e] << "-x[cm] " << surfaceNames[surf[e]] << "_" << vert[e] << "-y[cm] " << surfaceNames[surf[e]] << "_" << vert[e] << "-z[cm]";

        if (e != numEmaPoints-1)  // No space behind last value
        {
          emaFile << " ";
        }
      }
      emaFile << endl;
    }
    emaFile << std::setprecision(8);

    for (int f = 0; f < numFrames; f++)
    {
      fs::path objFileName = objFilePath.string() + to_string(f) + ".obj" ;  // name of current .obj file
      emaFile << f*0.005 << " "; // Write current time (frame*0.005 currently) into file

      // Set the properties of the target tube.
      for (i = 0; i < numTractParams; i++)
      {
        vocalTract->param[i].x = tractParams[i + (numTractParams*f)];
      }

      // Extract Glottis Parameters in each Frame
      for (i = 0; i < numGlottisParams; i++)
      {
        glotParams[i] = glottisParams[i + (numGlottisParams*f)];
      }
      vocalTract->calculateAll();

      // Calculate coordinates for each selected EMA points
      for (e = 0; e < numEmaPoints; e++)
      {
        s = &vocalTract->surface[surf[e]];

        P = s->vertex[vert[e]].coord;
        emaFile << P.x << " " << P.y << " " << P.z;

        if (e != numEmaPoints-1)  // No space behind last value
        {
          emaFile << " ";
        }
      }

      vocalTract->saveAsObjFile(objFileName.string(), true); // write .obj files of current vocal tract
      emaFile << endl;

    }

    // ****************************************************************
    // Close the file with the EMA points.
    // ****************************************************************
    if (emaFile.is_open())
    {
      emaFile.close();
    }
  }
  return 0;
}


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

int vtlGesturalScoreToEmaAndMesh(const char *gestureFileName, const char *filePath, const char *fileName)
{
  if (!vtlApiInitialized)
  {
    printf("Error: The API has not been initialized.\n");
    return 1;
  }

  int i, f;
  double t_s; // time in sec

  // ****************************************************************
  // Init and load the gestural score.
  // ****************************************************************

  bool allValuesInRange = true;

  // Create the gestural score not before we know what selectedGlottis
  // is, which is determined in loadSpeakerFile() (this can be overwritten later).
  GesturalScore *gesturalScore = new GesturalScore(vocalTract, glottis[selectedGlottis]);

  if (gesturalScore->loadGesturesXml(string(gestureFileName), allValuesInRange) == false)
  {
    printf("Error in vtlGesToTractAndGlottisModel(): Loading the gestural score file failed!\n");
    return 2;
  }

  if (allValuesInRange == false)
  {
    printf("Error in vtlGesToTractAndGlottisModel(): Some values in the gestural score are out of range!\n");
    return 3;
  }

  // Important !!!
  gesturalScore->calcCurves();

  // *****************************************************************************************
  // Store Tract and Glottis Parameter in tractParams[] and glottisParams[]
  // *****************************************************************************************

  const double EMA_SAMPLING_RATE_HZ = 200.0;
  int numFrames = (int)(gesturalScore->getScoreDuration_s() * EMA_SAMPLING_RATE_HZ);

  int numTractParams = VocalTract::NUM_PARAMS;
  int numGlottisParams = (int)glottis[selectedGlottis]->controlParam.size();

  std::vector<double> storeTractParams(numTractParams, 0.0);
  std::vector<double> storeGlottisParams(numGlottisParams, 0.0);


  std::vector<double> tractParams(numFrames * numTractParams, 0.0);
  std::vector<double> glottisParams(numFrames * numGlottisParams, 0.0);

  // ************************************************************************
  // Extract glottis and tract params here
  // ************************************************************************

  for (f=0; f < numFrames; f++)
  {
    t_s = (double)f / (double)EMA_SAMPLING_RATE_HZ;
    gesturalScore->getParams(t_s, storeTractParams.data(), storeGlottisParams.data());

    for (i = 0; i < numTractParams; i++)
    {
      tractParams[i + (f*numTractParams)] = storeTractParams[i];
    }

    for (i = 0; i < numGlottisParams; i++)
    {
      glottisParams[i + (f*numGlottisParams)] = storeGlottisParams[i];
    }
  }

  // *****************************************************************************************
  // Free the memory.
  // *****************************************************************************************

  delete gesturalScore;

  int surf[] = {16,16,16};   // Surface of selected EMA point
  int vertex[] = {115,225,335};  // Vertex Index of selected EMA point

  // Call of vtlTractSequenceToEmaAndMesh could also be done from outside!
  return vtlTractSequenceToEmaAndMesh(tractParams.data(), glottisParams.data(), numTractParams, numGlottisParams, numFrames, 3, surf, vertex, filePath, fileName);
}

// ****************************************************************************

// ****************************************************************************
// Set control parameters on the current glottis model, call calcGeometry(),
// and return derived parameters and tube data.
// ****************************************************************************

int vtlGlottisCalcGeometry(double *controlParams, double *derivedParams,
                           int *numDerivedParams, double *tubeLength_cm,
                           double *tubeArea_cm2)
{
  if (!vtlApiInitialized)
  {
    printf("Error: The API has not been initialized.\n");
    return 1;
  }

  Glottis *g = glottis[selectedGlottis];
  int numCP = (int)g->controlParam.size();

  // Set control parameters
  for (int i = 0; i < numCP; i++)
  {
    g->controlParam[i].x = controlParams[i];
  }

  // Calculate geometry
  g->calcGeometry();

  // Return derived parameters
  int numDP = (int)g->derivedParam.size();
  *numDerivedParams = numDP;
  for (int i = 0; i < numDP; i++)
  {
    derivedParams[i] = g->derivedParam[i].x;
  }

  // Return tube data
  g->getTubeData(tubeLength_cm, tubeArea_cm2);

  return 0;
}

// ****************************************************************************
// Perform one incTime step, then calcGeometry, and return updated state.
// ****************************************************************************

int vtlGlottisIncTime(double timeIncrement_s, double *pressure_dPa,
                      double *controlParams, double *derivedParams,
                      int *numDerivedParams, double *tubeLength_cm,
                      double *tubeArea_cm2)
{
  if (!vtlApiInitialized)
  {
    printf("Error: The API has not been initialized.\n");
    return 1;
  }

  Glottis *g = glottis[selectedGlottis];
  int numCP = (int)g->controlParam.size();

  // Set control parameters
  for (int i = 0; i < numCP; i++)
  {
    g->controlParam[i].x = controlParams[i];
  }

  // Advance time
  g->incTime(timeIncrement_s, pressure_dPa);

  // Calculate geometry with updated state
  g->calcGeometry();

  // Return derived parameters
  int numDP = (int)g->derivedParam.size();
  *numDerivedParams = numDP;
  for (int i = 0; i < numDP; i++)
  {
    derivedParams[i] = g->derivedParam[i].x;
  }

  // Return tube data
  g->getTubeData(tubeLength_cm, tubeArea_cm2);

  return 0;
}

// ****************************************************************************
// Reset the motion state of the current glottis model.
// ****************************************************************************

int vtlGlottisResetMotion()
{
  if (!vtlApiInitialized)
  {
    printf("Error: The API has not been initialized.\n");
    return 1;
  }

  glottis[selectedGlottis]->resetMotion();
  return 0;
}

// ****************************************************************************
// Returns static parameter info for the current glottis model.
// ****************************************************************************

int vtlGetGlottisStaticParamInfo(char *names, double *paramMin,
                                 double *paramMax, double *paramStandard,
                                 int *numStaticParams)
{
  if (!vtlApiInitialized)
  {
    printf("Error: The API has not been initialized.\n");
    return 1;
  }

  Glottis *g = glottis[selectedGlottis];
  int numSP = (int)g->staticParam.size();
  *numStaticParams = numSP;

  strcpy(names, "");
  for (int i = 0; i < numSP; i++)
  {
    strcat(names, g->staticParam[i].name.c_str());
    if (i != numSP - 1)
    {
      strcat(names, "\t");
    }
    paramMin[i] = g->staticParam[i].min;
    paramMax[i] = g->staticParam[i].max;
    paramStandard[i] = g->staticParam[i].neutral;
  }

  return 0;
}

// ****************************************************************************
// TDS API implementations
// ****************************************************************************

int vtlTdsSetOptions(bool generateNoiseSources, bool turbulenceLosses,
                     bool softWalls, bool radiationFromSkin, bool piriformFossa,
                     bool innerLengthCorrections, bool transvelarCoupling)
{
  if (!vtlApiInitialized)
  {
    printf("Error: The API has not been initialized.\n");
    return 1;
  }

  tdsModel->options.generateNoiseSources = generateNoiseSources;
  tdsModel->options.turbulenceLosses = turbulenceLosses;
  tdsModel->options.softWalls = softWalls;
  tdsModel->options.radiationFromSkin = radiationFromSkin;
  tdsModel->options.piriformFossa = piriformFossa;
  tdsModel->options.innerLengthCorrections = innerLengthCorrections;
  tdsModel->options.transvelarCoupling = transvelarCoupling;
  return 0;
}

int vtlTdsResetMotion()
{
  if (!vtlApiInitialized)
  {
    printf("Error: The API has not been initialized.\n");
    return 1;
  }

  tdsModel->resetMotion();
  return 0;
}

int vtlSetFossaDims(double length_cm, double volume_cm3)
{
  if (!vtlApiInitialized)
  {
    printf("Error: The API has not been initialized.\n");
    return 1;
  }

  tube->initPiriformFossa(length_cm, volume_cm3);
  return 0;
}

int vtlTdsSetTubeAndRun(double *tubeLength_cm, double *tubeArea_cm2,
                        int *tubeArticulator, double incisorPos_cm,
                        double velumOpening_cm2, double tongueTipSideElevation,
                        bool filtering, int pressureSourceSection,
                        double pressureSourceAmp, double *secArea,
                        double *secLength, double *secR0, double *secR1,
                        double *secL, double *secC, double *secD, double *secE,
                        double *secAlpha, double *secBeta, double *secPressure,
                        double *bcMagnitude, double *mouthFlow,
                        double *nostrilFlow, double *skinFlow)
{
  if (!vtlApiInitialized)
  {
    printf("Error: The API has not been initialized.\n");
    return 1;
  }

  // Set tube geometry
  Tube::Articulator articulator[Tube::NUM_PHARYNX_MOUTH_SECTIONS];
  for (int i = 0; i < Tube::NUM_PHARYNX_MOUTH_SECTIONS; i++)
  {
    articulator[i] = (Tube::Articulator)tubeArticulator[i];
  }

  tube->setPharynxMouthGeometry(tubeLength_cm, tubeArea_cm2, articulator,
                                incisorPos_cm, tongueTipSideElevation);
  tube->setVelumOpening(velumOpening_cm2);

  // Load into TDS model
  tdsModel->setTube(tube, filtering);

  // Set pressure source
  if (pressureSourceSection >= 0)
  {
    tdsModel->setPressureSource(pressureSourceAmp, pressureSourceSection);
  }
  else
  {
    tdsModel->pressureSourceSection = -1;
    tdsModel->pressureSourceAmp = 0.0;
  }

  // Run one time step
  double mouth = 0.0, nostril = 0.0, skin = 0.0;
  tdsModel->proceedTimeStep(mouth, nostril, skin);

  // Copy out all internal state
  for (int i = 0; i < Tube::NUM_SECTIONS; i++)
  {
    secArea[i] = tdsModel->tubeSection[i].area;
    secLength[i] = tdsModel->tubeSection[i].length;
    secR0[i] = tdsModel->tubeSection[i].R[0];
    secR1[i] = tdsModel->tubeSection[i].R[1];
    secL[i] = tdsModel->tubeSection[i].L;
    secC[i] = tdsModel->tubeSection[i].C;
    secD[i] = tdsModel->tubeSection[i].D;
    secE[i] = tdsModel->tubeSection[i].E;
    secAlpha[i] = tdsModel->tubeSection[i].alpha;
    secBeta[i] = tdsModel->tubeSection[i].beta;
    secPressure[i] = tdsModel->tubeSection[i].pressure;
  }

  for (int i = 0; i < TdsModel::NUM_BRANCH_CURRENTS; i++)
  {
    bcMagnitude[i] = tdsModel->branchCurrent[i].magnitude;
  }

  *mouthFlow = mouth;
  *nostrilFlow = nostril;
  *skinFlow = skin;

  return 0;
}

// ****************************************************************************

// ****************************************************************************

int vtlGetTLIntermediateValues(
    double *tractParams, int numSpectrumSamples, TransferFunctionOptions *opts,
    int freqIndex, double *matrix_A_re, double *matrix_A_im,
    double *matrix_B_re, double *matrix_B_im, double *matrix_C_re,
    double *matrix_C_im, double *matrix_D_re, double *matrix_D_im,
    double *fossa_input_imp_re, double *fossa_input_imp_im,
    double *nose_rad_imp_re, double *nose_rad_imp_im, double *mouth_rad_imp_re,
    double *mouth_rad_imp_im)
{

  if (!vtlApiInitialized)
  {
    printf("Error: The API has not been initialized.\n");
    return 1;
  }

  TlModel *tlModel = new TlModel();

  // Set options
  if (opts != NULL)
  {
    TlModel::Options tlOpts;
    tlOpts.radiation = (TlModel::RadiationType)opts->radiationType;
    tlOpts.boundaryLayer = opts->boundaryLayer;
    tlOpts.heatConduction = opts->heatConduction;
    tlOpts.softWalls = opts->softWalls;
    tlOpts.hagenResistance = opts->hagenResistance;
    tlOpts.innerLengthCorrections = opts->innerLengthCorrections;
    tlOpts.lumpedElements = opts->lumpedElements;
    tlOpts.paranasalSinuses = opts->paranasalSinuses;
    tlOpts.piriformFossa = opts->piriformFossa;
    tlOpts.staticPressureDrops = opts->staticPressureDrops;
    tlModel->options = tlOpts;
  }

  // Calculate Tube
  for (int i = 0; i < VocalTract::NUM_PARAMS; i++)
  {
    vocalTract->param[i].x = tractParams[i];
  }
  vocalTract->calculateAll();
  vocalTract->getTube(&tlModel->tube);
  tlModel->tube.setGlottisArea(0.0);

  ComplexSignal dummySpectrum;
  tlModel->getSpectrum(TlModel::FLOW_SOURCE_TF, &dummySpectrum,
                       numSpectrumSamples, Tube::FIRST_PHARYNX_SECTION);

  if (freqIndex < 0 || freqIndex >= (numSpectrumSamples / 2))
  {
    delete tlModel;
    return 2;
  }

  for (int i = 0; i < Tube::NUM_SECTIONS; i++)
  {
    Matrix2x2 M = tlModel->getMatrixProduct(i, freqIndex);
    matrix_A_re[i] = M.A.real();
    matrix_A_im[i] = M.A.imag();
    matrix_B_re[i] = M.B.real();
    matrix_B_im[i] = M.B.imag();
    matrix_C_re[i] = M.C.real();
    matrix_C_im[i] = M.C.imag();
    matrix_D_re[i] = M.D.real();
    matrix_D_im[i] = M.D.imag();
  }

  ComplexValue fcc = tlModel->getFossaInputImpedance(freqIndex);
  *fossa_input_imp_re = fcc.real();
  *fossa_input_imp_im = fcc.imag();

  ComplexValue ncc = tlModel->getNoseRadiationImpedance(freqIndex);
  *nose_rad_imp_re = ncc.real();
  *nose_rad_imp_im = ncc.imag();

  ComplexValue mcc = tlModel->getMouthRadiationImpedance(freqIndex);
  *mouth_rad_imp_re = mcc.real();
  *mouth_rad_imp_im = mcc.imag();

  delete tlModel;
  return 0;
}

// ****************************************************************************
// Returns the 129 cross-section areas, positions, and articulators from the
// VocalTract model after calling calculateAll on the given tract parameters.
// ****************************************************************************

int vtlGetCrossSections(double *tractParams, double *crossSectionAreas,
                        double *crossSectionPositions,
                        int *crossSectionArticulators)
{
  if (!vtlApiInitialized)
  {
    printf("Error: The API has not been initialized.\n");
    return 1;
  }

  // Store the current control parameter values.
  vocalTract->storeControlParams();

  // Set the given vocal tract parameters.
  int i;
  for (i = 0; i < VocalTract::NUM_PARAMS; i++)
  {
    vocalTract->param[i].x = tractParams[i];
  }

  // Calculate the vocal tract shape.
  vocalTract->calculateAll();

  // Copy the cross-section data to the output arrays.
  for (i = 0; i < VocalTract::NUM_CENTERLINE_POINTS; i++)
  {
    crossSectionAreas[i] = vocalTract->crossSection[i].area;
    crossSectionPositions[i] = vocalTract->crossSection[i].pos;
    crossSectionArticulators[i] = (int)vocalTract->crossSection[i].articulator;
  }

  // Restore the previous control parameter values.
  vocalTract->restoreControlParams();

  return 0;
}

// ****************************************************************************
// Returns the upper and lower cross-sectional profiles at a specific
// centerline index for the given vocal tract parameters.
// ****************************************************************************

int vtlGetProfiles(double *tractParams, int centerlineIndex,
                   double *upperProfile, double *lowerProfile,
                   double *centerlineInfo)
{
  if (!vtlApiInitialized)
  {
    printf("Error: The API has not been initialized.\n");
    return 1;
  }

  if (centerlineIndex < 0 ||
      centerlineIndex >= VocalTract::NUM_CENTERLINE_POINTS)
  {
    printf("Error: centerlineIndex %d is out of range [0, %d).\n",
           centerlineIndex, VocalTract::NUM_CENTERLINE_POINTS);
    return 2;
  }

  // Store the current control parameter values.
  vocalTract->storeControlParams();

  // Set the given vocal tract parameters.
  int i;
  for (i = 0; i < VocalTract::NUM_PARAMS; i++)
  {
    vocalTract->param[i].x = tractParams[i];
  }

  // Calculate the vocal tract shape.
  vocalTract->calculateAll();

  // Get the cross profiles at the specified centerline index.
  double upProfile[VocalTract::NUM_PROFILE_SAMPLES];
  double loProfile[VocalTract::NUM_PROFILE_SAMPLES];
  Tube::Articulator articulator;

  vocalTract->getCrossProfiles(vocalTract->centerLine[centerlineIndex].point,
                               vocalTract->centerLine[centerlineIndex].normal,
                               upProfile, loProfile, true, articulator);

  // Copy profile data to output arrays.
  for (i = 0; i < VocalTract::NUM_PROFILE_SAMPLES; i++)
  {
    upperProfile[i] = upProfile[i];
    lowerProfile[i] = loProfile[i];
  }

  // Copy centerline info: point.x, point.y, normal.x, normal.y, area, pos
  centerlineInfo[0] = vocalTract->centerLine[centerlineIndex].point.x;
  centerlineInfo[1] = vocalTract->centerLine[centerlineIndex].point.y;
  centerlineInfo[2] = vocalTract->centerLine[centerlineIndex].normal.x;
  centerlineInfo[3] = vocalTract->centerLine[centerlineIndex].normal.y;
  centerlineInfo[4] = vocalTract->crossSection[centerlineIndex].area;
  centerlineInfo[5] = vocalTract->crossSection[centerlineIndex].pos;

  // Restore the previous control parameter values.
  vocalTract->restoreControlParams();

  return 0;
}

// ****************************************************************************
// Returns all 129 centerline points after computing vocal tract geometry.
// ****************************************************************************

int vtlGetCenterline(double *tractParams, double *centerlineData)
{
  if (!vtlApiInitialized)
  {
    printf("Error: The API has not been initialized.\n");
    return 1;
  }

  // Store the current control parameter values.
  vocalTract->storeControlParams();

  // Set the given vocal tract parameters.
  int i;
  for (i = 0; i < VocalTract::NUM_PARAMS; i++)
  {
    vocalTract->param[i].x = tractParams[i];
  }

  // Calculate the vocal tract shape.
  vocalTract->calculateAll();

  // Copy all 129 centerline points: x, y, normal_x, normal_y, pos
  for (i = 0; i < VocalTract::NUM_CENTERLINE_POINTS; i++)
  {
    centerlineData[i * 5 + 0] = vocalTract->centerLine[i].point.x;
    centerlineData[i * 5 + 1] = vocalTract->centerLine[i].point.y;
    centerlineData[i * 5 + 2] = vocalTract->centerLine[i].normal.x;
    centerlineData[i * 5 + 3] = vocalTract->centerLine[i].normal.y;
    centerlineData[i * 5 + 4] = vocalTract->centerLine[i].pos;
  }

  // Restore the previous control parameter values.
  vocalTract->restoreControlParams();

  return 0;
}

// ****************************************************************************
// Returns the 4 outlines used for centerline computation.
// ****************************************************************************

int vtlGetOutlines(double *tractParams, double *outlineData, int *outlineSizes)
{
  if (!vtlApiInitialized)
  {
    printf("Error: The API has not been initialized.\n");
    return 1;
  }

  vocalTract->storeControlParams();

  int i;
  for (i = 0; i < VocalTract::NUM_PARAMS; i++)
  {
    vocalTract->param[i].x = tractParams[i];
  }

  vocalTract->calculateAll();

  // Export 4 outlines: upper, lower, tongue, epiglottis
  // Each outline is stored as n points of (x, y) pairs
  int offset = 0;
  auto exportOutline = [&](LineStrip2D &outline, int outlineIdx)
  {
    int n = outline.getNumPoints();
    outlineSizes[outlineIdx] = n;
    for (int j = 0; j < n; j++)
    {
      outlineData[offset++] = outline.getControlPoint(j).x;
      outlineData[offset++] = outline.getControlPoint(j).y;
    }
  };

  exportOutline(vocalTract->upperOutline, 0);
  exportOutline(vocalTract->lowerOutline, 1);
  exportOutline(vocalTract->tongueOutline, 2);
  exportOutline(vocalTract->epiglottisOutline, 3);

  vocalTract->restoreControlParams();

  return 0;
}

// ****************************************************************************
// Returns tongue rib data (point, normal, leftSideHeight, rightSideHeight).
// ****************************************************************************

int vtlGetTongueRibData(double *tractParams, double *ribData, int *numRibs)
{
  if (!vtlApiInitialized)
  {
    printf("Error: The API has not been initialized.\n");
    return 1;
  }

  vocalTract->storeControlParams();

  int i;
  for (i = 0; i < VocalTract::NUM_PARAMS; i++)
  {
    vocalTract->param[i].x = tractParams[i];
  }

  vocalTract->calculateAll();

  *numRibs = VocalTract::NUM_TONGUE_RIBS;
  for (i = 0; i < VocalTract::NUM_TONGUE_RIBS; i++)
  {
    ribData[i * 6 + 0] = vocalTract->tongueRib[i].point.x;
    ribData[i * 6 + 1] = vocalTract->tongueRib[i].point.y;
    ribData[i * 6 + 2] = vocalTract->tongueRib[i].normal.x;
    ribData[i * 6 + 3] = vocalTract->tongueRib[i].normal.y;
    ribData[i * 6 + 4] = vocalTract->tongueRib[i].leftSideHeight;
    ribData[i * 6 + 5] = vocalTract->tongueRib[i].rightSideHeight;
  }

  vocalTract->restoreControlParams();

  return 0;
}

// ****************************************************************************
// Returns tongue rib width bounds (minX, maxX) after getCrossProfiles + lowpass.
// boundsData: flat array of size numRibs * 2 (minX, maxX per rib).
// ****************************************************************************

int vtlGetTongueWidthBounds(double *tractParams, double *boundsData,
                            int *numRibs)
{
  if (!vtlApiInitialized)
  {
    printf("Error: The API has not been initialized.\n");
    return 1;
  }

  vocalTract->storeControlParams();

  int i;
  for (i = 0; i < VocalTract::NUM_PARAMS; i++)
  {
    vocalTract->param[i].x = tractParams[i];
  }

  vocalTract->calculateAll();

  *numRibs = VocalTract::NUM_TONGUE_RIBS;
  for (i = 0; i < VocalTract::NUM_TONGUE_RIBS; i++)
  {
    boundsData[i * 2 + 0] = vocalTract->tongueRib[i].minX;
    boundsData[i * 2 + 1] = vocalTract->tongueRib[i].maxX;
  }

  vocalTract->restoreControlParams();

  return 0;
}

// ****************************************************************************
// Export all vertex positions for a specific surface.
// vertexData: flat array of size numRibs * numRibPoints * 3 (x, y, z per vertex)
// ****************************************************************************

int vtlGetSurfaceVertices(double *tractParams, int surfaceIndex,
                          double *vertexData, int *numRibs, int *numRibPoints)
{
  if (!vtlApiInitialized)
  {
    printf("Error: The API has not been initialized.\n");
    return 1;
  }
  if (surfaceIndex < 0 || surfaceIndex >= VocalTract::NUM_SURFACES)
  {
    printf("Error: Invalid surface index %d.\n", surfaceIndex);
    return 2;
  }

  vocalTract->storeControlParams();

  int i;
  for (i = 0; i < VocalTract::NUM_PARAMS; i++)
  {
    vocalTract->param[i].x = tractParams[i];
  }

  vocalTract->calculateAll();

  Surface *s = &vocalTract->surface[surfaceIndex];
  *numRibs = s->numRibs;
  *numRibPoints = s->numRibPoints;

  int idx = 0;
  for (int r = 0; r < s->numRibs; r++)
  {
    for (int p = 0; p < s->numRibPoints; p++)
    {
      Point3D v = s->getVertex(r, p);
      vertexData[idx++] = v.x;
      vertexData[idx++] = v.y;
      vertexData[idx++] = v.z;
    }
  }

  vocalTract->restoreControlParams();

  return 0;
}

// ****************************************************************************
// Returns raw triangle intersection cuts at a specific centerline index.
// cutData: flat array of size 2048 * 8 (P0.x, P0.y, P1.x, P1.y,
//          n.x, n.y, globalSurfaceIndex, localSurfaceIndex per cut).
// numCuts: output, number of cuts found.
// ****************************************************************************

int vtlGetCuts(double *tractParams, int centerlineIndex, double *cutData,
               int *numCuts)
{
  if (!vtlApiInitialized)
  {
    printf("Error: The API has not been initialized.\n");
    return 1;
  }

  if (centerlineIndex < 0 ||
      centerlineIndex >= VocalTract::NUM_CENTERLINE_POINTS)
  {
    printf("Error: centerlineIndex %d is out of range [0, %d).\n",
           centerlineIndex, VocalTract::NUM_CENTERLINE_POINTS);
    return 2;
  }

  vocalTract->storeControlParams();

  int i;
  for (i = 0; i < VocalTract::NUM_PARAMS; i++)
  {
    vocalTract->param[i].x = tractParams[i];
  }

  vocalTract->calculateAll();

  const int MAX_CUTS = 2048;
  *numCuts = vocalTract->getRawCuts(
      vocalTract->centerLine[centerlineIndex].point,
      vocalTract->centerLine[centerlineIndex].normal, cutData, MAX_CUTS);

  vocalTract->restoreControlParams();

  return 0;
}

// ****************************************************************************
// Apply anatomy parameters derived from age and gender to the loaded vocal
// tract.
// ****************************************************************************

int vtlSetAnatomyFromAge(int ageMonths, bool isMale)
{
  if (!vtlApiInitialized)
  {
    printf("Error: The API has not been initialized.\n");
    return 1;
  }

  AnatomyParams ap;
  ap.calcFromAge(ageMonths, isMale);
  ap.restrictParams();
  ap.setFor(vocalTract);

  return 0;
}

// ****************************************************************************
// Get the 13 anatomy parameters from the currently loaded vocal tract.
// ****************************************************************************

int vtlGetAnatomyParams(double *anatomyParams)
{
  if (!vtlApiInitialized)
  {
    printf("Error: The API has not been initialized.\n");
    return 1;
  }

  AnatomyParams ap;
  ap.getFrom(vocalTract);
  for (int i = 0; i < AnatomyParams::NUM_ANATOMY_PARAMS; i++)
  {
    anatomyParams[i] = ap.param[i].x;
  }

  return 0;
}

// ****************************************************************************
// Set the 13 anatomy parameters on the currently loaded vocal tract.
// ****************************************************************************

int vtlSetAnatomyParams(double *anatomyParams)
{
  if (!vtlApiInitialized)
  {
    printf("Error: The API has not been initialized.\n");
    return 1;
  }

  AnatomyParams ap;
  for (int i = 0; i < AnatomyParams::NUM_ANATOMY_PARAMS; i++)
  {
    ap.param[i].x = anatomyParams[i];
  }
  ap.restrictParams();
  ap.setFor(vocalTract);

  return 0;
}
