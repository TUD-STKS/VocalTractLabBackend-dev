#ifndef _ACOUSTIC_3D_SIMULATION_
#define _ACOUSTIC_3D_SIMULATION_

#include "CrossSection2d.h"
#include "VocalTract.h"
#include "Signal.h"
#include <vector>
#include <wx/wx.h>
#include <wx/progdlg.h>

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
  NOISE
};

class Acoustic3dSimulation
{
// **************************************************************************
/// Public functions
// **************************************************************************

public:

  ~Acoustic3dSimulation();
  static Acoustic3dSimulation *getInstance();

  void setBoundarySpecificAdmittance();
  void setSimulationParameters(double meshDensity,
    int secNoiseSource, int secConstriction,
    int expSpectrumLgth, struct simulationParameters simuParams,
    enum openEndBoundaryCond cond);
  void setIdxSecNoiseSource(int idx) { m_idxSecNoiseSource = idx; }
  void generateLogFileHeader(bool cleanLog);

  void setGeometryImported(bool isImported) {
    m_geometryImported = isImported;
  }
  void setGeometryFile(string fileName) {
    m_geometryFile = fileName;}
  void setContourInterpolationMethod(enum contourInterpolationMethod method);

  void addCrossSectionFEM(double areas, double spacing,
    Polygon_2 contours, vector<int> surfacesIdx,
    double length, Point2D ctrLinePt, Point2D normal, 
    double scalingFactors[2]);

  void addCrossSectionRadiation(Point2D ctrLinePt, Point2D normal,
    double radius, double PMLThickness);

  // For geometry creation
  void extractContours(VocalTract* tract, 
    vector<vector<Polygon_2>>& contours, vector<vector<vector<int>>>& surfaceIdx,
    vector<Point2D>& centerLine, vector<Point2D>& normals);
  bool extractContoursFromCsvFile(
    vector<vector<Polygon_2>>& contours, vector<vector<vector<int>>>& surfaceIdx,
    vector<Point2D>& centerLine, vector<Point2D>& normals, 
    vector<pair<double, double>>& scalingFactors, bool simplifyContours);
  bool createCrossSections(VocalTract* tract, bool createRadSection);
  void updateBoundingBox();
  void setBoundingBox(pair<Point2D, Point2D> &bbox);
  void importGeometry(VocalTract* tract);

  // For solving the wave problem 
  void computeMeshAndModes();
  void computeJunctionMatrices(bool computeG);
  void requestModesAndJunctionComputation() { m_simuParams.needToComputeModesAndJunctions = true; }
  void propagateImpedAdmitBranch(vector< Eigen::MatrixXcd> Q0, double freq,
    vector<int> startSections, vector<int> endSections, double direction);
  void propagateImpedAdmit(Eigen::MatrixXcd& startImped, Eigen::MatrixXcd& startAdmit, 
    double freq, int startSection, int endSection, std::chrono::duration<double> *time, int direction);
  void propagateImpedAdmit(Eigen::MatrixXcd& startImped, Eigen::MatrixXcd& startAdmit,
    double freq, int startSection, int endSection, std::chrono::duration<double> *time);
  void propagateAdmit(Eigen::MatrixXcd radImped, double freq);
  void propagateVelocityPress(Eigen::MatrixXcd &startVelocity, Eigen::MatrixXcd &startPressure, 
    double freq, int startSection, int endSection, std::chrono::duration<double> *time, int direction);
  void propagateVelocityPress(Eigen::MatrixXcd& startVelocity, Eigen::MatrixXcd& startPressure,
    double freq, int startSection, int endSection, std::chrono::duration<double> *time);
  void propagatePressure(Eigen::MatrixXcd startVelocity, double freq);
  //void propagateAcPressure(vector<Eigen::MatrixXcd> inputPressure, double freq);

  // For acoustic field and transfer function computation
  Point_3 movePointFromExitLandmarkToGeoLandmark(Point_3 pt);
  vector<Point_3> movePointFromExitLandmarkToGeoLandmark(vector<Point_3> pt);
  bool setTFPointsFromCsvFile(string fileName);
  void RayleighSommerfeldIntegral(vector<Point_3> points,
    Eigen::VectorXcd &radPress, double freq, int radSecIdx);
  complex<double> acousticField(Point_3 queryPt);
  bool findSegmentContainingPoint(Point queryPt, int &idxSeg);
  Eigen::VectorXcd acousticField(vector<Point_3> queryPt);
  void acousticFieldInPlane(Eigen::MatrixXcd& field);
  void solveWaveProblem(VocalTract* tract, double freq, bool precomputeRadImped,
    std::chrono::duration<double>& time, std::chrono::duration<double> *timeExp);
  void computeTransferFunction(VocalTract* tract);
  void computeAcousticField(VocalTract* tract);
  void coneConcatenationSimulation(string fileName);
  void runTest(enum testType tType, string fileName);
  void cleanAcousticField();

  // for data exportation
  complex<double> interpolateTransferFunction(double freq, int idxPt, enum tfType type);
  double interpolateAcousticField(Point querryPt);
  void interpolateAcousticField(Vec &coordX, Vec &coordY, Matrix &field);
  void exportGeoInCsv(string fileName);
  bool exportTransferFucntions(string fileName, enum tfType type);
  bool exportAcousticField(string fileName);


// **************************************************************************
// accessors

  struct simulationParameters simuParams() const { return m_simuParams; }
  bool isGeometryImported() const { return m_geometryImported; }
  int sectionNumber() const;
  double soundSpeed() const;
  bool freqDepLosses() const { return m_simuParams.freqDepLosses; }
  CrossSection2d* crossSection(int csIdx) const;
  double meshDensity() const;
  //int modeNumber() const;
  double maxCutOnFreq() const;
  int numIntegrationStep() const;
  propagationMethod method() const { return m_simuParams.propMethod; }
  int spectrumLgthExponent() const { return m_spectrumLgthExponent; }
  int idxSecNoiseSource() const { return m_idxSecNoiseSource; }
  int idxConstriction() const { return m_idxConstriction; }
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
  int numCrossSections(){ return m_crossSections.size(); }
  openEndBoundaryCond mouthBoundaryCond() const {return m_mouthBoundaryCond;}
  int acousticFieldSize() const { return m_field.size(); }
  double maxAmpField() const { return m_maxAmpField; }
  double minAmpField() const { return m_minAmpField; }
  double lastFreqComputed() const { return m_lastFreqComputed; }

// **************************************************************************
/// Public data
// **************************************************************************

public:
  ComplexSignal spectrum, spectrumNoise, spectrumConst;

// **************************************************************************
/// Private data
// **************************************************************************

private:

  static Acoustic3dSimulation *instance;

  wxGenericProgressDialog *progressDialog;
  vector<unique_ptr<CrossSection2d>> m_crossSections;

  struct simulationParameters m_simuParams;

  // simulation parameters
  bool m_geometryImported;
  string m_geometryFile;
  contourInterpolationMethod m_contInterpMeth;
  double m_meshDensity;
  // the number of frequencies is 2 ^ (m_spectrumLgthExponent - 1)
  int m_spectrumLgthExponent;
  int m_numFreq;
  int m_numFreqPicture;
  double m_lastFreqComputed;
  int m_idxSecNoiseSource;
  int m_idxConstriction;
  openEndBoundaryCond m_glottisBoundaryCond;
  openEndBoundaryCond m_mouthBoundaryCond;
  vector<vector<vector<vector<double>>>> m_radiationMatrixInterp;
  vector<double> m_radiationFreqs;
  vector<double> m_tfFreqs;

  // maximal bounding box of the cross-sections (for displaying mesh and modes)
  pair<Point2D, Point2D> m_maxCSBoundingBox;
  // for displaying segments and acoustic field
  //pair<Point2D, Point2D> m_bboxXZ;

  // simulation outputs
  Eigen::MatrixXcd m_glottalSourceTF;
  Eigen::MatrixXcd m_noiseSourceTF;
  Eigen::MatrixXcd m_field;
  double m_maxAmpField;
  double m_minAmpField;
  Eigen::VectorXcd m_planeModeInputImpedance;

// **************************************************************************
// Private functions.
// **************************************************************************

private:

  Acoustic3dSimulation();
  void createContour(double inputUpProf[VocalTract::NUM_PROFILE_SAMPLES],
    double inputLoProf[VocalTract::NUM_PROFILE_SAMPLES], 
    int upperProfileSurface[VocalTract::NUM_PROFILE_SAMPLES],
    int lowerProfileSurface[VocalTract::NUM_PROFILE_SAMPLES],
    vector<double> &areas, vector<double> &spacing,
    vector< Polygon_2> &contours, vector<vector<int>> &surfacesIdx);
  void getCurvatureAngleShift(Point2D P1, Point2D P2,
    Point2D N1, Point2D N2, double& radius, double& angle, double& shift);
  
  // for radiation impedance 
  void preComputeRadiationMatrices(int nbRadFreqs, int idxRadSec);
  void interpolateRadiationImpedance(Eigen::MatrixXcd& imped, double freq, int idxRadSec); 
  void interpolateRadiationAdmittance(Eigen::MatrixXcd& admit, double freq, int idxRadSec);
  void radiationImpedance(Eigen::MatrixXcd& imped, double freq, double gridDensity, int idxRadSec);
  void getRadiationImpedanceAdmittance(Eigen::MatrixXcd& imped, Eigen::MatrixXcd& admit, 
    double freq, int idxRadSec);
};

#endif
