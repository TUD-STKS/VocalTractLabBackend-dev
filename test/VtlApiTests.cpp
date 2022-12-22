#include <gtest/gtest.h>

#include <iostream>
#include <string>
#include <vector>

#include <cstring>
#include <cstdlib>

#include "VocalTractLabApi.h"

#include <filesystem>

// Some files used in or produced by the various tests. Paths are relative to repo root (set executable working directory accordingly!)
const char* emaFile = "test/resources/example01.ema";
const char* exportedMeshPath = "test/resources/exportedMesh/";
const char* exportMeshBaseName = "mesh";
const char* gesFile = "test/resources/example01.ges";
const char* speakerFile = "resources/JD3.speaker";
const char* svgFile = "test/resources/ExportTractSvg.svg";
const char* tractSeqFile = "test/resources/TractSequence.seq";
const char* wavFile = "test/resources/example01.wav";


TEST(ApiTest, Version) 
{
    std::cout << std::filesystem::current_path() << std::endl;
	char version[64];
	vtlGetVersion(version);
	std::string s(version);
	EXPECT_EQ(s, std::string("API 2.4.2 ") + std::string(__DATE__));
}

TEST(ApiTest, Constants)
{

	std::cout << std::filesystem::current_path() << std::endl;
	vtlInitialize(speakerFile);

	int audioSamplingRate, numTubeSections, numVocalTractParams, numGlottisParams, numAudioSamplesPerTractState;
	double internalAudioSamplingRate;
	int ret = vtlGetConstants(&audioSamplingRate, &numTubeSections,
		&numVocalTractParams, &numGlottisParams, &numAudioSamplesPerTractState, &internalAudioSamplingRate);

	EXPECT_EQ(ret, 0);
	EXPECT_EQ(audioSamplingRate, 44100);
	EXPECT_EQ(numTubeSections, 40);
	EXPECT_EQ(numVocalTractParams, 19);
	EXPECT_EQ(numGlottisParams, 11);

	vtlClose();
}

TEST(ApiTest, TractParamInfo)
{
	vtlInitialize(speakerFile);

	int _, numVocalTractParams;
	double d_;
	vtlGetConstants(&_, &_, &numVocalTractParams, &_, &_, &d_);

	std::vector<char> names(numVocalTractParams * 10);
	std::vector<char> descriptions(numVocalTractParams * 50);
	std::vector<char> units(numVocalTractParams * 10);
	std::vector<double> paramMin(numVocalTractParams);
	std::vector<double> paramMax(numVocalTractParams);
	std::vector<double> paramStandard(numVocalTractParams);

	int ret = vtlGetTractParamInfo(&names[0], &descriptions[0], &units[0], &paramMin[0], &paramMax[0], &paramStandard[0]);
	EXPECT_EQ(ret, 0);
	
	vtlClose();
}

TEST(ApiTest, GlottisParamInfo)
{
	vtlInitialize(speakerFile);

	int _, numGlottisParams;
	double d_;
	vtlGetConstants(&_, &_,	&_, &numGlottisParams, &_, &d_);
	
	std::vector<char> names(numGlottisParams * 10);
	std::vector<char> descriptions(numGlottisParams * 50);
	std::vector<char> units(numGlottisParams * 10);
	std::vector<double> paramMin(numGlottisParams);
	std::vector<double> paramMax(numGlottisParams);
	std::vector<double> paramStandard(numGlottisParams);

	int ret = vtlGetGlottisParamInfo(&names[0], &descriptions[0], &units[0], &paramMin[0], &paramMax[0], &paramStandard[0]);
	EXPECT_EQ(ret, 0);

	vtlClose();
}

TEST(ApiTest, TractParams)
{
	vtlInitialize(speakerFile);

	int _, numVocalTractParams;
	double d_;
	vtlGetConstants(&_, &_, &numVocalTractParams, &_, &_, &d_);

	const char* shapeName = "a";
	std::vector<double> param(numVocalTractParams);
	int ret = vtlGetTractParams(shapeName, &param[0]);
	
	EXPECT_EQ(ret, 0);

	vtlClose();
}

TEST(ApiTest, ExportTractSvg)
{
	vtlInitialize(speakerFile);

	int _, numVocalTractParams;
	double d_;
	vtlGetConstants(&_, &_, &numVocalTractParams, &_, &_, &d_);
	const char* shapeName = "a";
	std::vector<double> param(numVocalTractParams);
	vtlGetTractParams(shapeName, &param[0]);

	int ret = vtlExportTractSvg(&param[0], svgFile);

	EXPECT_EQ(ret, 0);

	vtlClose();
}

TEST(ApiTest, TractToTube)
{
	vtlInitialize(speakerFile);

	int _, numTubeSections, numVocalTractParams;
	double d_;
	vtlGetConstants(&_, &numTubeSections, &numVocalTractParams, &_, &_, &d_);
	const char* shapeName = "a";
	std::vector<double> param(numVocalTractParams);
	vtlGetTractParams(shapeName, &param[0]);

	std::vector<double> tubeLength_cm(numTubeSections);
	std::vector<double> tubeArea_cm2(numTubeSections);
	std::vector<int> tubeArticulator(numTubeSections);
	std::vector<double> incisorPos_cm(numTubeSections);
	std::vector<double> tongueTipSideElevation(numTubeSections);
	std::vector<double> velumOpening_cm2(numTubeSections);

	int ret = vtlTractToTube(&param[0], &tubeLength_cm[0], &tubeArea_cm2[0], &tubeArticulator[0],
		&incisorPos_cm[0], &tongueTipSideElevation[0], &velumOpening_cm2[0]);

	EXPECT_EQ(ret, 0);

	vtlClose();
}

TEST(ApiTest, FastTractToTube)
{
	vtlInitialize(speakerFile);

	int _, numTubeSections, numVocalTractParams;
	double d_;
	vtlGetConstants(&_, &numTubeSections, &numVocalTractParams, &_, &_, &d_);
	const char* shapeName = "a";
	std::vector<double> param(numVocalTractParams);
	vtlGetTractParams(shapeName, &param[0]);

	std::vector<double> tubeLength_cm(numTubeSections);
	std::vector<double> tubeArea_cm2(numTubeSections);
	std::vector<int> tubeArticulator(numTubeSections);
	std::vector<double> incisorPos_cm(numTubeSections);
	std::vector<double> tongueTipSideElevation(numTubeSections);
	std::vector<double> velumOpening_cm2(numTubeSections);

	int ret = vtlFastTractToTube(&param[0], &tubeLength_cm[0], &tubeArea_cm2[0], &tubeArticulator[0],
		&incisorPos_cm[0], &tongueTipSideElevation[0], &velumOpening_cm2[0]);

	EXPECT_EQ(ret, 0);

	vtlClose();
}

TEST(ApiTest, DefaultTransferFunctionOptions)
{
	vtlInitialize(speakerFile);

	TransferFunctionOptions defaultOpts;
	int ret = vtlGetDefaultTransferFunctionOptions(&defaultOpts);

	EXPECT_EQ(ret, 0);
	EXPECT_EQ(defaultOpts.radiationType, PARALLEL_RADIATION);
	EXPECT_EQ(defaultOpts.boundaryLayer, true);
	EXPECT_EQ(defaultOpts.heatConduction, false);
	EXPECT_EQ(defaultOpts.softWalls, true);
	EXPECT_EQ(defaultOpts.hagenResistance, false);
	EXPECT_EQ(defaultOpts.lumpedElements, true);
	EXPECT_EQ(defaultOpts.innerLengthCorrections, false);
	EXPECT_EQ(defaultOpts.paranasalSinuses, true);
	EXPECT_EQ(defaultOpts.piriformFossa, true);
	EXPECT_EQ(defaultOpts.staticPressureDrops, true);
	EXPECT_EQ(defaultOpts.spectrumType, SPECTRUM_UU);
}

TEST(ApiTest, GetTransferFunction)
{
	vtlInitialize(speakerFile);

	TransferFunctionOptions opts;
	vtlGetDefaultTransferFunctionOptions(&opts);
	opts.spectrumType = SPECTRUM_PU;

	int _, numTubeSections, numVocalTractParams;
	double d_;
	vtlGetConstants(&_, &numTubeSections, &numVocalTractParams, &_, &_, &d_);
	
	std::vector<double> tractParams(numVocalTractParams);
	vtlGetTractParams("a", &tractParams[0]);
	
	const int numSamples = 4096;
	std::vector<double> magnitude(numSamples);
	std::vector<double> phase(numSamples);
	int ret = vtlGetTransferFunction(&tractParams[0], numSamples,
		&opts, &magnitude[0], &phase[0]);

	EXPECT_EQ(ret, 0);
}


TEST(ApiTest, GesToAudio_FileOut)
{

	std::cout << std::filesystem::current_path() << std::endl;

	vtlInitialize(speakerFile);
	
	int ret = vtlGesturalScoreToAudio(gesFile, wavFile, NULL, NULL, true);

	EXPECT_EQ(ret, 0);

	vtlClose();
}

TEST(ApiTest, GesToAudio_DataOut)
{
	vtlInitialize(speakerFile);

	std::vector<double> audio(44100 * 30);
	int numSamples;
	int ret = vtlGesturalScoreToAudio(gesFile, "", &audio[0], &numSamples, true);

	EXPECT_EQ(ret, 0);

	vtlClose();
}

TEST(ApiTest, GesToEma)
{
	vtlInitialize(speakerFile);

	int ret = vtlGesturalScoreToEma(gesFile, emaFile);

	EXPECT_EQ(ret, 0);

	vtlClose();
}

TEST(ApiTest, GesToEmaAndMesh)
{
	vtlInitialize(speakerFile);

	// Make sure the export path is exists but is empty
	std::filesystem::remove_all(exportedMeshPath);
	std::filesystem::create_directory(exportedMeshPath);

	int ret = vtlGesturalScoreToEmaAndMesh(gesFile, exportedMeshPath, exportMeshBaseName);

	EXPECT_EQ(ret, 0);

	vtlClose();
}

TEST(ApiTest, GesToTract)
{
	vtlInitialize(speakerFile);

	int ret = vtlGesturalScoreToTractSequence(gesFile, tractSeqFile);

	EXPECT_EQ(ret, 0);

	vtlClose();
}


TEST(ApiTest, TractToAudio_FileOut)
{
	vtlInitialize(speakerFile);

	int ret = vtlTractSequenceToAudio(tractSeqFile, wavFile, NULL, NULL);

	EXPECT_EQ(ret, 0);

	vtlClose();
}

TEST(ApiTest, TractToAudio_DataOut)
{
	vtlInitialize(speakerFile);

	std::vector<double> audio(44100 * 30);
	int numSamples;
	int ret = vtlTractSequenceToAudio(tractSeqFile, "", &audio[0], &numSamples);

	EXPECT_EQ(ret, 0);

	vtlClose();
}

TEST(ApiTest, TractToEmaAndMesh)
{
	vtlInitialize(speakerFile);

	// Make sure the export path is exists but is empty
	std::filesystem::remove_all(exportedMeshPath);
	std::filesystem::create_directory(exportedMeshPath);

	int numFrames = 10;
	int audioSamplingRate, numTubeSections, numVocalTractParams, numGlottisParams, numAudioSamplesPerTractState;
	double internalAudioSamplingRate;
	int ret = vtlGetConstants(&audioSamplingRate, &numTubeSections,
		&numVocalTractParams, &numGlottisParams, &numAudioSamplesPerTractState, &internalAudioSamplingRate);
	std::vector<double> tractParams(numVocalTractParams, 0.0);
	std::vector<double> tractSeq;
	std::vector<double> glottisParams(numGlottisParams, 0.0);
	std::vector<double> glottisSeq;

	const char* shapeName = "a";
	vtlGetTractParams(shapeName, &tractParams[0]);
	
	const char* glottisShapeName = "modal";
	vtlGetGlottisParams(glottisShapeName, &glottisParams[0]);

	for (int i = 0; i < numFrames; ++i)
	{
		tractSeq.insert(tractSeq.end(), tractParams.begin(), tractParams.end());
		glottisSeq.insert(glottisSeq.end(), glottisParams.begin(), glottisParams.end());
	}

	int numEmaPoints = 3;
	int surf[] = { 16,16,16 };   // Surface of selected EMA point
	int vert[] = { 115,225,335 };  // Vertex Index of selected EMA point

	ret = vtlTractSequenceToEmaAndMesh(tractSeq.data(), glottisSeq.data(),
		numVocalTractParams, numGlottisParams, numFrames, numEmaPoints,
		surf, vert,
		exportedMeshPath, exportMeshBaseName);

	EXPECT_EQ(ret, 0);

	vtlClose();
}

TEST(ApiTest, SaveSpeaker)
{	
	vtlInitialize(speakerFile);
	vtlSaveSpeaker("newSpeakerFile.speaker");
}
