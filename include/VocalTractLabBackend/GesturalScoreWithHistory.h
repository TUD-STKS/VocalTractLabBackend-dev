// ****************************************************************************
// This file is part of VocalTractLab.
// Copyright (C) 2020, Peter Birkholz, Dresden, Germany
// www.vocaltractlab.de
// author: Simon Stone
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

#ifndef __GESTURAL_SCORE_WITH_HISTORY_H__
#define __GESTURAL_SCORE_WITH_HISTORY_H__

#include <vector>

#include "GesturalScore.h"

/// @brief This class is essentially a facade for the regular GesturalScore, but it additionally keeps track of all edits to the gestural score to allow undo and redo.
class GesturalScoreWithHistory
{
public:
	// Constructors
	GesturalScoreWithHistory() = default;
	GesturalScoreWithHistory(VocalTract* vocalTract, Glottis* glottis);
	// Copy and move constructors
	GesturalScoreWithHistory(const GesturalScoreWithHistory& other) = default;
	GesturalScoreWithHistory(GesturalScoreWithHistory&& other) = default;
	// Destructors
	~GesturalScoreWithHistory() = default;

	// Get the current gestural score
	std::vector<GesturalScore>::iterator get() const;

	bool CanUndo() const;
	bool CanRedo() const;

	void Undo();
	void Redo();

	/*
	 * Facade for GesturalScore interface
	 */
	void clear();
	void initTestScore();
	void createFromSegmentSequence(SegmentSequence* origSegmentSequence, bool enableConsoleOutput = true);

	void addClosingGesture(GesturalScore::GestureType gestureType, string gestureName,
		double closureBegin_s, double closureEnd_s, bool connectToPrevGesture);
	bool hasVocalTactClosure(GesturalScore::GestureType gestureType, string gestureName,
		double gestureBegin_s, double gestureEnd_s, double testTime_s) const;

	void addVelicOpeningGesture(double openingBegin_s, double openingEnd_s);
	bool hasVelicOpening(double gestureBegin_s, double gestureEnd_s, double testTime_s) const;

	bool changeGestureEnd(int gestureType, int gestureIndex, double newEnd_s, bool stretchNextGesture);
	bool changeTargetValue(int gestureType, int gestureIndex, double delta);
	bool deleteGesture(int gestureType, int gestureIndex);
	int insertGesture(int gestureType, double insertPos_s, int gestureIndex);

	bool setGestureValue(int gestureType, int gestureIndex, std::string newVal);
	bool setGestureValue(int gestureType, int gestureIndex, double newVal);
	bool setGestureDuration(int gestureType, int gestureIndex, double newDuration_s);
	bool setGestureNeutral(int gestureType, int gestureIndex, double isNeutral);

	bool loadGesturesXml(const string& fileName, bool& allValuesInRange);
	bool saveGesturesXml(const string& fileName) const;

	// MUST be called after any change to the score.
	void calcCurves() const;
	void getParams(double pos_s, double* vocalTractParams, double* glottisParams) const;
	double getScoreDuration_s() const;

	// Manipulation functions
	void changeF0Offset(double deltaF0_st);
	void changeF0Range(double factor);
	void changeF0TargetSlope(double deltaSlope_st_s);
	void substituteGlottalShapes(const string& oldShapeName, const string& newShapeName);
	void changeSubglottalPressure(double factor);
	void getF0Statistic(double& f0Mean_st, double& f0Sd_st, double& f0Mean_Hz, double& f0Sd_Hz) const;

	void changeDuration(double factor);
	void changeTimeConstants(double factor);

	GestureSequence* getGestures() const;
	Glottis* getGlottis() const;
	void setGlottis(Glottis* newGlottis);
	VocalTract* getVocalTract() const;
	void setVocalTract(VocalTract* newVocalTract);

	vector<Target>* getTractParamTargets();
	vector<Target>* getGlottisParamTargets();
	vector<double>* getTractParamCurve();
	vector<double>* getGlottisParamCurve();

	// **************************************************************************
	// Overwritten functions of the interface class.
	// **************************************************************************

	void getTube(Tube& tube) const;
	void getFlowSource(double& flow_cm3_s, int& section) const;
	void getPressureSource(double& pressure_dPa, int& section) const;

	void resetSequence() const;
	void incPos(const double pressure_dPa[]) const;
	int getDuration_pt() const;
	int getPos_pt() const;


private:
	[[nodiscard]] std::vector<GesturalScore>::iterator AddOperation();
	void ResetHistory();

private:
	// Memento design pattern: Keep a history of instances of the tracked object
	std::vector<GesturalScore> history_;
	std::vector<GesturalScore>::iterator current_;



public:
	// Operators
	GesturalScoreWithHistory& operator=(const GesturalScoreWithHistory& other) = default;
	GesturalScoreWithHistory& operator=(GesturalScoreWithHistory&& other) = default;	
};


#endif // __GESTURAL_SCORE_WITH_HISTORY_H__