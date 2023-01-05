// ****************************************************************************
// This file is part of VocalTractLab3D.
// Copyright (C) 2022, Peter Birkholz, Dresden, Germany
// www.vocaltractlab.de
// author: Peter Birkholz and Rémi Blandin
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

#ifndef _ACOUSTIC_3D_SIMULATION_
#define _ACOUSTIC_3D_SIMULATION_

#include "CrossSection2d.h"
#include "VocalTract.h"
#include "Signal.h"
#include <vector>

// for CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polygon_2<K>                            Polygon_2;
typedef CGAL::Point_3<K>                Point_3;

typedef Eigen::VectorXd                 Vec;

// **************************************************************************
// Identifiers for the glotis boundary condition
// **************************************************************************

enum openEndBoundaryCond {
  RADIATION,
  IFINITE_WAVGUIDE,
  HARD_WALL,
  ADMITTANCE_1,
  ZERO_PRESSURE
};

enum contourInterpolationMethod {
  AREA,
  BOUNDING_BOX,
  FROM_FILE
};

enum testType {
    MATRIX_E,
    DISCONTINUITY,
    ELEPHANT_TRUNK,
    SCALE_RAD_IMP
};

enum tfType {
  GLOTTAL,
  NOISE,
  INPUT_IMPED,
  RAD_IMPED_REAL_PART
};

class Acoustic3dSimulation
{
// **************************************************************************
/// Public functions
// **************************************************************************

public:

  ~Acoustic3dSimulation();
  static Acoustic3dSimulation *getInstance();

  // set simulation parameters
  void setBoundarySpecificAdmittance();
  void setSimulationParameters(double meshDensity, int secNoiseSource,
    struct simulationParameters simuParams,
    enum openEndBoundaryCond cond, enum contourInterpolationMethod scalingMethod);
  void setIdxSecNoiseSource(int idx) { m_idxSecNoiseSource = idx; }
  void setGeometryImported(bool isImported) { m_geometryImported = isImported; }
  void setGeometryFile(string fileName) { m_geometryFile = fileName; }
  void setContourInterpolationMethod(enum contourInterpolationMethod method);
  void requestReloadGeometry() { m_reloadGeometry = true; }
  void requestModesAndJunctionComputation() { m_simuParams.needToComputeModesAndJunctions = true; }
  void setNeedToComputeModesAndJunctions(bool val) { m_simuParams.needToComputeModesAndJunctions = val; }

  void generateLogFileHeader(bool cleanLog);

  // For geometry creation
  void extractContours(VocalTract* tract, 
    vector<vector<Polygon_2>>& contours, vector<vector<vector<int>>>& surfaceIdx,
    vector<Point2D>& centerLine, vector<Point2D>& normals);
  bool extractContoursFromCsvFile(
    vector<vector<Polygon_2>>& contours, vector<vector<vector<int>>>& surfaceIdx,
    vector<Point2D>& centerLine, vector<Point2D>& normals, 
    vector<pair<double, double>>& scalingFactors, bool simplifyContours);
  void addCrossSectionFEM(double areas, double spacing,
    Polygon_2 contours, vector<int> surfacesIdx,
    double length, Point2D ctrLinePt, Point2D normal,
    double scalingFactors[2]);
  void addCrossSectionRadiation(Point2D ctrLinePt, Point2D normal,
    double radius, double PMLThickness);
  bool createCrossSections(VocalTract* tract, bool createRadSection);
  void updateBoundingBox();
  void setBoundingBox(pair<Point2D, Point2D> &bbox);
  bool importGeometry(VocalTract* tract);

  // For transverse modes and junction matrices computation
  void computeMeshAndModes();
  void computeMeshAndModes(int segIdx);
  void computeJunctionMatrices(int segIdx);
  void computeJunctionMatrices(bool computeG);

  // For radiation impedance computation
  void preComputeRadiationMatrices(int nbRadFreqs, int idxRadSec);
  void initCoefInterpRadiationMatrices(int nbRadFreqs, int idxRadSec);
  void addRadMatToInterpolate(int nbRadFreqs, int idxRadSec, int idxRadFreq);
  void computeInterpCoefRadMat(int nbRadFreqs, int idxRadSec);

  void propagateImpedAdmitBranch(vector< Eigen::MatrixXcd> Q0, double freq,
    vector<int> startSections, vector<int> endSections, double direction);
  void propagateImpedAdmit(Eigen::MatrixXcd& startImped, Eigen::MatrixXcd& startAdmit, 
    double freq, int startSection, int endSection, std::chrono::duration<double> *time, int direction);
  void propagateImpedAdmit(Eigen::MatrixXcd& startImped, Eigen::MatrixXcd& startAdmit,
    double freq, int startSection, int endSection, std::chrono::duration<double> *time);
  void propagateVelocityPress(Eigen::MatrixXcd &startVelocity, Eigen::MatrixXcd &startPressure, 
    double freq, int startSection, int endSection, std::chrono::duration<double> *time, int direction);
  void propagateVelocityPress(Eigen::MatrixXcd& startVelocity, Eigen::MatrixXcd& startPressure,
    double freq, int startSection, int endSection, std::chrono::duration<double> *time);

  // For acoustic field and transfer function computation
  Point_3 movePointFromExitLandmarkToGeoLandmark(Point_3 pt);
  vector<Point_3> movePointFromExitLandmarkToGeoLandmark(vector<Point_3> pt);
  bool setTFPointsFromCsvFile(string fileName);
  void RayleighSommerfeldIntegral(vector<Point_3> points,
    Eigen::VectorXcd &radPress, double freq, int radSecIdx);
  void setAcousticFieldFreq(double freq) {m_simuParams.freqField = freq;}
  complex<double> acousticField(Point_3 queryPt);
  bool findSegmentContainingPoint(Point queryPt, int &idxSeg);
  Eigen::VectorXcd acousticField(vector<Point_3> queryPt);
  void prepareAcousticFieldComputation();
  void acousticFieldInLine(int idxLine);
  void acousticFieldInPlane();
  void precomputationsForTf();
  void solveWaveProblem(VocalTract* tract, double freq, bool precomputeRadImped,
    std::chrono::duration<double>& time, std::chrono::duration<double> *timeExp);
  void solveWaveProblem(VocalTract* tract, double freq,
    std::chrono::duration<double>& time, std::chrono::duration<double>* timeExp);
  void solveWaveProblemNoiseSrc(bool &needToExtractMatrixF, Matrix& F, double freq,
    std::chrono::duration<double>* time);
  void computeGlottalTf(int idxFreq, double freq);
  void computeNoiseSrcTf(int idxFreq);
  void generateSpectraForSynthesis(int tfIdx);
  void computeTransferFunction(VocalTract* tract);
  void computeAcousticField(VocalTract* tract);
  void coneConcatenationSimulation(string fileName);
  void runTest(enum testType tType, string fileName);
  void cleanAcousticField();
  double maxAmpField();
  double minAmpField();
  void setFieldImageComputation(bool val) { m_simuParams.computeFieldImage = val; }

  // for data exportation
  complex<double> interpolateTransferFunction(double freq, int idxPt, enum tfType type);
  void interpolateTransferFunction(vector<double> &freq, int idxPt, enum tfType type, 
    vector<complex<double>> &interpolatedValues);
  double interpolateAcousticField(Point querryPt);
  void interpolateAcousticField(Vec &coordX, Vec &coordY, Matrix &field);
  bool exportGeoInCsv(string fileName);
  bool exportTransferFucntions(string fileName, enum tfType type);
  bool exportAcousticField(string fileName);


// **************************************************************************
// accessors

  struct simulationParameters simuParams() const { return m_simuParams; }
  struct simulationParameters oldSimuParams() const { return m_oldSimuParams; }
  enum contourInterpolationMethod contInterpMeth() const {return m_contInterpMeth;}
  bool isReloadGeometryRequested() const { return m_reloadGeometry; }
  bool isGeometryImported() const { return m_geometryImported; }
  int sectionNumber() const;
  double soundSpeed() const;
  bool viscoThermalLosses() const { return m_simuParams.viscoThermalLosses; }
  CrossSection2d* crossSection(int csIdx) const;
  double meshDensity() const;
  double maxCutOnFreq() const;
  int numIntegrationStep() const;
  propagationMethod method() const { return m_simuParams.propMethod; }
  int spectrumLgthExponent() const { return m_simuParams.spectrumLgthExponent; }
  int oldSpectrumLgthExponent() const { return m_oldSimuParams.spectrumLgthExponent; }
  int idxSecNoiseSource() const { return m_idxSecNoiseSource; }
  pair<Point2D, Point2D> maxCSBoundingBox() const { return m_maxCSBoundingBox; }
  pair<Point2D, Point2D> bboxSagittalPlane() const 
  { 
    pair<Point2D, Point2D> bbox;
    bbox.first = Point2D(m_simuParams.bbox[0].x(),
      m_simuParams.bbox[0].y());
    bbox.second = Point2D(m_simuParams.bbox[1].x(),
      m_simuParams.bbox[1].y());
    return bbox;
  }
  int numberOfSegments() const { return m_crossSections.size(); }
  bool radImpedPrecomputed() const { return m_simuParams.radImpedPrecomputed; }
  openEndBoundaryCond mouthBoundaryCond() const {return m_mouthBoundaryCond;}
  double freqAcousticField() const {return m_simuParams.freqField;}
  enum physicalQuantity fieldPhysicalQuantity() const { return m_simuParams.fieldPhysicalQuantity; }
  bool showFieldAmplitude() const { return m_simuParams.showAmplitude; }
  bool fieldIndB() const { return m_simuParams.fieldIndB; }
  bool computeFieldImage() const { return m_simuParams.computeFieldImage; }
  int acousticFieldSize() const { return m_field.size(); }
  int numPtXField() const { return m_nPtx; }
  int numPtYField() const { return m_nPty; }
  int numFreqComputed() const { return m_numFreqComputed; }
  double freqSteps() const { return m_freqSteps; }
  double lastFreqComputed() const { return m_lastFreqComputed; }

// **************************************************************************
/// Public data
// **************************************************************************

public:
  ComplexSignal spectrum, spectrumNoise;

// **************************************************************************
/// Private data
// **************************************************************************

private:

  static Acoustic3dSimulation *instance;

  vector<unique_ptr<CrossSection2d>> m_crossSections;

  struct simulationParameters m_simuParams;
  struct simulationParameters m_oldSimuParams;

  // simulation parameters
  bool m_geometryImported;
  bool m_reloadGeometry;
  string m_geometryFile;
  contourInterpolationMethod m_contInterpMeth;
  double m_meshDensity;
  // the number of frequencies is 2 ^ (spectrumLgthExponent - 1)
  int m_numFreq;
  int m_numFreqPicture;
  double m_lastFreqComputed;
  int m_idxSecNoiseSource;
  openEndBoundaryCond m_glottisBoundaryCond;
  openEndBoundaryCond m_mouthBoundaryCond;
  vector<vector<vector<vector<double>>>> m_radiationMatrixInterp;
  vector<double> m_radiationFreqs;
  double m_freqSteps;
  int m_numFreqComputed;
  vector<double> m_tfFreqs;
  // transfer function points in the landmark of the geometry
  // (it is easier to give the coordinates of the transfer function points in the
  // the landmark of the exit segment, but one need their coordinate in the landmark 
  // of the geometry for the computation of the acoustic field)
  vector<Point_3> m_tfPoints;

  // parameters for field computation
  double m_lx;
  double m_ly;
  int m_nPtx;
  int m_nPty;

  // maximal bounding box of the cross-sections (for displaying mesh and modes)
  pair<Point2D, Point2D> m_maxCSBoundingBox;

  // simulation outputs
  Eigen::MatrixXcd m_glottalSourceTF;
  Eigen::MatrixXcd m_noiseSourceTF;
  Eigen::MatrixXcd m_planeModeInputImpedance;
  Eigen::MatrixXcd m_field;
  double m_maxAmpField;
  double m_minAmpField;
  double m_maxPhaseField;
  double m_minPhaseField;

// **************************************************************************
// Private functions.
// **************************************************************************

private:

  Acoustic3dSimulation();
  Point ctrLinePtOut(Point ctrLinePtIn, Vector normalIn, double circleArcAngle, 
    double curvatureRadius, double length);

  void createContour(double inputUpProf[VocalTract::NUM_PROFILE_SAMPLES],
    double inputLoProf[VocalTract::NUM_PROFILE_SAMPLES], 
    int upperProfileSurface[VocalTract::NUM_PROFILE_SAMPLES],
    int lowerProfileSurface[VocalTract::NUM_PROFILE_SAMPLES],
    vector<double> &areas, vector<double> &spacing,
    vector< Polygon_2> &contours, vector<vector<int>> &surfacesIdx);

  void createUniqueContour(double inputUpProf[VocalTract::NUM_PROFILE_SAMPLES],
    double inputLoProf[VocalTract::NUM_PROFILE_SAMPLES],
    int upperProfileSurface[VocalTract::NUM_PROFILE_SAMPLES],
    int lowerProfileSurface[VocalTract::NUM_PROFILE_SAMPLES],
    vector<double>& areas, vector<double>& spacing,
    vector< Polygon_2>& contours, vector<vector<int>>& surfacesIdx);

  void removeDuplicatedPoints(Polygon_2& contour);
  void mergeContours(vector<Polygon_2>& vecPoly, vector<vector<int>>& surfaceIdx);
  void makeContourConvexHull(Polygon_2& poly, vector<int>& surfaceIdx);

  void getCurvatureAngleShift(Point2D P1, Point2D P2,
    Point2D N1, Point2D N2, double& radius, double& angle, double& shift);
  
  // for radiation impedance 
  void interpolateRadiationImpedance(Eigen::MatrixXcd& imped, double freq, int idxRadSec); 
  void interpolateRadiationAdmittance(Eigen::MatrixXcd& admit, double freq, int idxRadSec);
  void radiationImpedance(Eigen::MatrixXcd& imped, double freq, double gridDensity, int idxRadSec);
  void getRadiationImpedanceAdmittance(Eigen::MatrixXcd& imped, Eigen::MatrixXcd& admit, 
    double freq, int idxRadSec);
};

#endif
