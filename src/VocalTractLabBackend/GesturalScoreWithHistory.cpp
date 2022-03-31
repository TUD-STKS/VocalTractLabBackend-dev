#include "GesturalScoreWithHistory.h"

GesturalScoreWithHistory::GesturalScoreWithHistory(VocalTract* vocalTract, Glottis* glottis)
{
	history_.emplace_back(vocalTract, glottis);
	current_ = history_.begin();
}

std::vector<GesturalScore>::iterator GesturalScoreWithHistory::get() const
{
	return current_;
}

bool GesturalScoreWithHistory::CanUndo() const
{
	return current_ - history_.begin() > 0;
}

bool GesturalScoreWithHistory::CanRedo() const
{
	return current_ < history_.end() - 1;
}

void GesturalScoreWithHistory::Undo()
{
	if(CanUndo())
	{
		--current_;
	}
}

void GesturalScoreWithHistory::Redo()
{
	if (CanRedo())
	{
		++current_;
	}
}

void GesturalScoreWithHistory::clear()
{
	AddOperation()->clear();
}

void GesturalScoreWithHistory::initTestScore()
{
	AddOperation()->initTestScore();
}

void GesturalScoreWithHistory::createFromSegmentSequence(SegmentSequence* origSegmentSequence, bool enableConsoleOutput)
{
	AddOperation()->createFromSegmentSequence(origSegmentSequence, enableConsoleOutput);
}

void GesturalScoreWithHistory::addClosingGesture(GesturalScore::GestureType gestureType, string gestureName, double closureBegin_s, double closureEnd_s, bool connectToPrevGesture)
{
	AddOperation()->addClosingGesture(gestureType, gestureName, closureBegin_s, closureEnd_s, connectToPrevGesture);
}

bool GesturalScoreWithHistory::hasVocalTactClosure(GesturalScore::GestureType gestureType, string gestureName,
	double gestureBegin_s, double gestureEnd_s, double testTime_s) const
{
	return current_->hasVocalTactClosure(gestureType, gestureName, gestureBegin_s, gestureEnd_s, testTime_s);
}

void GesturalScoreWithHistory::addVelicOpeningGesture(double openingBegin_s, double openingEnd_s)
{
	AddOperation()->addVelicOpeningGesture(openingBegin_s, openingEnd_s);
}

bool GesturalScoreWithHistory::hasVelicOpening(double gestureBegin_s, double gestureEnd_s, double testTime_s) const
{
	return current_->hasVelicOpening(gestureBegin_s, gestureEnd_s, testTime_s);
}

bool GesturalScoreWithHistory::changeGestureEnd(int gestureType, int gestureIndex, double newEnd_s, bool stretchNextGesture)
{
	const bool success = AddOperation()->changeGestureEnd(gestureType, gestureIndex, newEnd_s, stretchNextGesture);
	if (!success)
	{
		// Nothing was changed, so don't keep a record
		history_.pop_back();
		current_ = history_.end() - 1;
	}
	return success;
}

bool GesturalScoreWithHistory::changeTargetValue(int gestureType, int gestureIndex, double delta)
{
	const bool success = AddOperation()->changeTargetValue(gestureType, gestureIndex, delta);
	if (!success)
	{
		// Nothing was changed, so don't keep a record
		history_.pop_back();
		current_ = history_.end() - 1;
	}
	return success;
}

bool GesturalScoreWithHistory::deleteGesture(int gestureType, int gestureIndex)
{
	const bool success = AddOperation()->deleteGesture(gestureType, gestureIndex);
	if (!success)
	{
		// Nothing was changed so don't keep a record
		history_.pop_back();
		current_ = history_.end() - 1;
	}
	return success;
}

int GesturalScoreWithHistory::insertGesture(int gestureType, double insertPos_s,
                                            int gestureIndex)
{
	return AddOperation()->insertGesture(gestureType, insertPos_s, gestureIndex);
}

bool GesturalScoreWithHistory::loadGesturesXml(const string& fileName, bool& allValuesInRange)
{
	return AddOperation()->loadGesturesXml(fileName, allValuesInRange);
}

bool GesturalScoreWithHistory::saveGesturesXml(const string& fileName) const
{
	return current_->saveGesturesXml(fileName);
}

void GesturalScoreWithHistory::calcCurves() const
{
	current_->calcCurves();
}

void GesturalScoreWithHistory::getParams(double pos_s, double* vocalTractParams, double* glottisParams) const
{
	return current_->getParams(pos_s, vocalTractParams, glottisParams);
}

double GesturalScoreWithHistory::getScoreDuration_s() const
{
	return current_->getScoreDuration_s();
}

void GesturalScoreWithHistory::changeF0Offset(double deltaF0_st)
{
	AddOperation()->changeF0Offset(deltaF0_st);
}

void GesturalScoreWithHistory::changeF0Range(double factor)
{
	AddOperation()->changeF0Range(factor);
}

void GesturalScoreWithHistory::changeF0TargetSlope(double deltaSlope_st_s)
{
	AddOperation()->changeF0TargetSlope(deltaSlope_st_s);
}

void GesturalScoreWithHistory::substituteGlottalShapes(const string& oldShapeName, const string& newShapeName)
{
	AddOperation()->substituteGlottalShapes(oldShapeName, newShapeName);
}

void GesturalScoreWithHistory::changeSubglottalPressure(double factor)
{
	AddOperation()->changeSubglottalPressure(factor);
}

void GesturalScoreWithHistory::getF0Statistic(double& f0Mean_st, double& f0Sd_st, double& f0Mean_Hz, double& f0Sd_Hz) const
{
	current_->getF0Statistic(f0Mean_st, f0Sd_st, f0Mean_Hz, f0Sd_Hz);
}

void GesturalScoreWithHistory::changeDuration(double factor)
{
	AddOperation()->changeDuration(factor);
}

void GesturalScoreWithHistory::changeTimeConstants(double factor)
{
	AddOperation()->changeTimeConstants(factor);
}

GestureSequence* GesturalScoreWithHistory::getGestures() const
{
	return current_->getGestures();
}

Glottis* GesturalScoreWithHistory::getGlottis() const
{
	return current_->getGlottis();
}

void GesturalScoreWithHistory::setGlottis(Glottis* newGlottis)
{
	AddOperation()->setGlottis(newGlottis);
}

VocalTract* GesturalScoreWithHistory::getVocalTract() const
{
	return current_->getVocalTract();
}

void GesturalScoreWithHistory::setVocalTract(VocalTract* newVocalTract)
{
	AddOperation()->setVocalTract(newVocalTract);
}

vector<Target>* GesturalScoreWithHistory::getTractParamTargets()
{
	return current_->getTractParamTargets();
}

vector<Target>* GesturalScoreWithHistory::getGlottisParamTargets()
{
	return current_->getGlottisParamTargets();
}

vector<double>* GesturalScoreWithHistory::getTractParamCurve()
{
	return current_->getTractParamCurve();
}

vector<double>* GesturalScoreWithHistory::getGlottisParamCurve()
{
	return current_->getGlottisParamCurve();
}

void GesturalScoreWithHistory::getTube(Tube& tube) const
{
	current_->getTube(tube);
}

void GesturalScoreWithHistory::getFlowSource(double& flow_cm3_s, int& section) const
{
	current_->getFlowSource(flow_cm3_s, section);
}

void GesturalScoreWithHistory::getPressureSource(double& pressure_dPa, int& section) const
{
	current_->getPressureSource(pressure_dPa, section);
}

void GesturalScoreWithHistory::resetSequence() const
{
	current_->resetSequence();
}

void GesturalScoreWithHistory::incPos(const double pressure_dPa[]) const
{
	current_->incPos(pressure_dPa);
}

int GesturalScoreWithHistory::getDuration_pt() const
{
	return current_->getDuration_pt();
}

int GesturalScoreWithHistory::getPos_pt() const
{
	return current_->getPos_pt();
}


std::vector<GesturalScore>::iterator GesturalScoreWithHistory::AddOperation()
{
	// First remove all steps from current position to end of vector
	history_.erase(current_ + 1, history_.end());
	// Then add a copy of the last state
	history_.emplace_back(*current_);
	current_ = history_.end() - 1;	
	return current_;
}

void GesturalScoreWithHistory::ResetHistory()
{
	history_.clear();
	current_ = history_.begin();
}
