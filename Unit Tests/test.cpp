#include "pch.h"

#include <iostream>
#include <string>
#include <vector>

#include "../VocalTractLabApi.h"

// Some files used in or produced by the various tests
const char* gesFile = "../../Unit Tests/example01.ges";
const char* speakerFile = "../../Unit Tests/JD2.speaker";
const char* svgFile = "../../Unit Tests/ExportTractSvg.svg";
const char* tractSeqFile = "../../Unit Tests/TractSequence.seq";
const char* wavFile = "../../Unit Tests/example01.wav";


TEST(ApiTest, Version) 
{
	char version[32];
	vtlGetVersion(version);
	std::string s(version);
	EXPECT_EQ(s, std::string(__DATE__));
}

TEST(ApiTest, Constants)
{
	vtlInitialize(speakerFile);

	int audioSamplingRate, numTubeSections, numVocalTractParams, numGlottisParams;
	int ret = vtlGetConstants(&audioSamplingRate, &numTubeSections,
		&numVocalTractParams, &numGlottisParams);

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
	vtlGetConstants(&_, &_, &numVocalTractParams, &_);

	std::vector<char> names(numVocalTractParams * 10);
	std::vector<double> paramMin(numVocalTractParams);
	std::vector<double> paramMax(numVocalTractParams);
	std::vector<double> paramNeutral(numVocalTractParams);

	int ret = vtlGetTractParamInfo(&names[0], &paramMin[0], &paramMax[0], &paramNeutral[0]);
	EXPECT_EQ(ret, 0);
	
	vtlClose();
}

TEST(ApiTest, GlottisParamInfo)
{
	vtlInitialize(speakerFile);

	int _, numGlottisParams;
	vtlGetConstants(&_, &_,	&_, &numGlottisParams);
	
	std::vector<char> names(numGlottisParams * 10);
	std::vector<double> paramMin(numGlottisParams);
	std::vector<double> paramMax(numGlottisParams);
	std::vector<double> paramNeutral(numGlottisParams);

	int ret = vtlGetGlottisParamInfo(&names[0], &paramMin[0], &paramMax[0], &paramNeutral[0]);
	EXPECT_EQ(ret, 0);

	vtlClose();
}

TEST(ApiTest, TractParams)
{
	vtlInitialize(speakerFile);

	int _, numVocalTractParams;
	vtlGetConstants(&_, &_, &numVocalTractParams, &_);

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
	vtlGetConstants(&_, &_, &numVocalTractParams, &_);
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
	vtlGetConstants(&_, &numTubeSections, &numVocalTractParams, &_);
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

TEST(ApiTest, DefaultTransferFunctionOptions)
{
	vtlInitialize(speakerFile);

	TransferFunctionOptions defaultOpts;
	int ret = vtlGetDefaultTransferFunctionOptions(&defaultOpts);

	EXPECT_EQ(ret, 0);
	EXPECT_EQ(defaultOpts.radiation, PARALLEL_RADIATION);
	EXPECT_EQ(defaultOpts.boundaryLayer, true);
	EXPECT_EQ(defaultOpts.heatConduction, false);
	EXPECT_EQ(defaultOpts.softWalls, true);
	EXPECT_EQ(defaultOpts.hagenResistance, false);
	EXPECT_EQ(defaultOpts.lumpedElements, true);
	EXPECT_EQ(defaultOpts.innerLengthCorrections, false);
	EXPECT_EQ(defaultOpts.paranasalSinuses, true);
	EXPECT_EQ(defaultOpts.piriformFossa, true);
	EXPECT_EQ(defaultOpts.staticPressureDrops, true);
	EXPECT_EQ(defaultOpts.type, SPECTRUM_UU);
}

TEST(ApiTest, GetTransferFunction)
{
	vtlInitialize(speakerFile);

	TransferFunctionOptions opts;
	vtlGetDefaultTransferFunctionOptions(&opts);
	opts.type = SPECTRUM_PU;

	int _, numVocalTractParams;
	vtlGetConstants(&_, &_, &numVocalTractParams, &_);
	
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