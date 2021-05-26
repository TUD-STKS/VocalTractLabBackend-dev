#ifndef __CROSS_SECTION_2D_H__
#define __CROSS_SECTION_2D_H__

#include <string>
#include <chrono>		// to get the computation time
#include <ctime>  
#include <vector>
#include <boost/bimap.hpp>
#include "Geometry.h"

// for eigen
#include <Eigen/Dense>
#include <Eigen/Sparse>

// for CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include "Delaunay_mesh_vertex_base_with_info_2.h"
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Boolean_set_operations_2.h>

using namespace std;

// typedef for eigen
typedef Eigen::MatrixXd Matrix;
typedef Eigen::SparseMatrix<complex<double>> SparseMatC;

// typedef for CGAL
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_mesh_vertex_base_with_info_2<unsigned int, K>    Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Exact_intersections_tag                     Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds, Itag> CDT;
typedef CDT::Point									Point;

typedef CGAL::Polygon_2<K>                          Polygon_2;

enum propagationMethod {
	MAGNUS,
	STRAIGHT_TUBES
};

enum physicalQuantity {
	IMPEDANCE,
	ADMITTANCE,
	PRESSURE,
	VELOCITY
};

struct simulationParameters
{
	double temperature;
	double volumicMass;
	double sndSpeed;
	int numIntegrationStep;
	complex<double> viscousBndSpecAdm;
	complex<double> thermalBndSpecAdm;
	bool freqDepLosses;
	propagationMethod propMethod;
	double percentageLosses;
	bool wallLosses;
	bool curved;
	bool varyingArea;
	double maxComputedFreq;
};


/////////////////////////////////////////////////////////////////////////////
// classe Cross section 2d
/////////////////////////////////////////////////////////////////////////////

class CrossSection2d
{

public:

	CrossSection2d(double maxCutOnFreq, Point2D ctrLinePt, Point2D normal);
	~CrossSection2d();

	// cross section parameters
	//virtual void set(double area, double spacing,
	//	Polygon_2 contours, vector<int> surfacesIdx,
	//	double inMeshDensity, double maxCutOnFreq, double inLength,
	//	Point2D ctrLinePt, Point2D normal) {
	//	;
	//}
	virtual void setJunctionSection(bool junction) { ; }
	void setComputImpedance(bool imp) { m_computeImpedance = imp; }
	void setPreviousSection(int prevSec);
	void setPrevSects(vector<int> prevSects) { m_previousSections = prevSects; }
	void setNextSection(int nextSec);
	void setNextSects(vector<int> nextSects) { m_nextSections = nextSects; }
	void clearPrevSects() { m_previousSections.clear(); }
	void clearNextSects() { m_nextSections.clear(); }
	virtual void setCurvatureRadius(double& radius) { ; }
	virtual void setCurvatureAngle(double& angle) { ; }

	// cross section mesh and modes
	virtual void buildMesh() { ; }
	virtual void computeModes(struct simulationParameters simuParams) { ; }
	virtual Matrix interpolateModes(vector<Point> pts) { return Matrix(); }
	Matrix interpolateModes(vector<Point> pts, double scaling);

	// scatering  matrices 
	virtual void setMatrixF(vector<Matrix> & F) { ; }
	virtual void setIntersectionsArea(vector<double> areaInt) { ; }
	virtual void setMatrixGstart(Matrix Gs) { ; }
	virtual void setMatrixGend(Matrix Ge) { ; }

	// impedance, admittance, acoustic pressure and axial velocity
	void setImpedance(vector<Eigen::MatrixXcd> inputImped) { m_impedance = inputImped; }
	void setZin(Eigen::MatrixXcd imped);
	void setZout(Eigen::MatrixXcd imped);
	void clearImpedance() { m_impedance.clear(); }
	void setAdmittance(vector<Eigen::MatrixXcd> inputAdmit) { m_admittance = inputAdmit; }
	void setYin(Eigen::MatrixXcd admit);
	void setYout(Eigen::MatrixXcd admit);
	void clearAdmittance() { m_admittance.clear(); }
	virtual void characteristicImpedance(
		Eigen::MatrixXcd & characImped, double freq, struct simulationParameters simuParams) {;}
	virtual void characteristicAdmittance(
		Eigen::MatrixXcd& admit, double freq, struct simulationParameters simuParams) {;}
	virtual complex<double> getWallAdmittance( 
		struct simulationParameters simuParams, double freq) { return complex<double>(); }
	virtual void getSpecificBndAdm(struct simulationParameters simuParams, double freq, 
		Eigen::VectorXcd& bndSpecAdm) {;}
	void setAxialVelocity(vector<Eigen::MatrixXcd> inputVelocity) { m_axialVelocity = inputVelocity; }
	void clearAxialVelocity() { m_axialVelocity.clear(); }
	void setAcPressure(vector<Eigen::MatrixXcd> inputPressure) { m_acPressure = inputPressure; }
	void clearAcPressure() { m_acPressure.clear(); }

	// propagation 
	virtual void propagateMagnus(Eigen::MatrixXcd Q0, struct simulationParameters simuParams,
		double freq, double direction, enum physicalQuantity quant) {;}
	virtual void propagateImpedAdmiteRiccati(Eigen::MatrixXcd Z0,
		Eigen::MatrixXcd Y0, struct simulationParameters simuParams,
		double freq, double direction, std::chrono::duration<double>& time) {;}
	virtual void propagateImpedRiccati(Eigen::MatrixXcd Z0, double nextArea, double freq) {;}
	virtual void propagateAdmitRiccati(Eigen::MatrixXcd Y0, 
		struct simulationParameters simuParams, double nextArea, double freq, double direction) {;}
	virtual void propagatePressureVelocityRiccati(Eigen::MatrixXcd V0, Eigen::MatrixXcd P0, 
		struct simulationParameters simuParams, double nextArea, double freq, double direction) {;}
	virtual void propagatePressureRiccati(Eigen::MatrixXcd P0,
		struct simulationParameters simuParams, double nextArea, double freq) {;}
	virtual void propagateImpedAdmitStraight(Eigen::MatrixXcd Z0,
		Eigen::MatrixXcd Y0, double freq, struct simulationParameters simuParams,
		double prevArea, double nextArea) {;}
	virtual void propagatePressureVelocityStraight(Eigen::MatrixXcd V0,
		Eigen::MatrixXcd P0, double freq, struct simulationParameters simuParams, double nextArea) {;}
	// to get the amplitude of the pressure modes at a given distance from the exit
	virtual void radiatePressure(double distance, double freq,
		struct simulationParameters simuParams, Eigen::MatrixXcd& pressAmp) { ; }

	// **************************************************************************
	// accessors

	int numPrevSec() const;
	int numNextSec() const;
	int prevSec(int idx) const;
	vector<int> prevSections() const { return m_previousSections; }
	int nextSec(int idx) const;
	vector<int> nextSections() const { return m_nextSections; }
	bool computeImpedance() const { return m_computeImpedance; }
	Point2D ctrLinePt() const;
	Point2D normal() const;
	double area() const;
	int numberOfModes() const { return m_modesNumber; }

	virtual double scaleIn() const { return 1.;  }
	virtual double scaleOut() const { return 1.; }
	virtual double length() const { return 0.; }
	virtual vector<double> intersectionsArea() const { return vector<double>(0); }
	virtual double curvature() { return 0.; }
	virtual double circleArcAngle() const { return 0.; }
	virtual double spacing() const { return 0.; }
	virtual int numberOfVertices() const { return 0; }
	virtual int numberOfFaces() const { return 0; }
	virtual CDT triangulation() const { return CDT(); }
	virtual Polygon_2 contour() const { return Polygon_2(); }
	virtual bool isJunction() const { return bool(); }
	virtual vector<int> surfaceIdx() const { return vector<int>(); }
	virtual double eigenFrequency(int idxMode) const { return 0.; }
	virtual vector<array<double, 2>> getPoints() const { return vector<array<double, 2>>(); }
	virtual vector<array<int, 3>> getTriangles() const { return vector<array<int, 3>>(); }
	virtual Matrix getModes() const { return Matrix(); }
	virtual double getMaxAmplitude(int idxMode) const { return 0.; }
	virtual double getMinAmplitude(int idxMode) const { return 0.; }
	virtual vector<Matrix> getMatrixF() const { return vector<Matrix>(); }
	virtual Matrix getMatrixGStart() const { return Matrix(); }
	virtual Matrix getMatrixGEnd() const { return Matrix(); }
	virtual double curvRadius() const { return double(); }
	virtual double radius() const { return double(); }
	virtual double PMLThickness() const { return double(); }

	vector<Eigen::MatrixXcd> Z() const { return m_impedance; }
	Eigen::MatrixXcd Zin() const { return m_impedance.back(); }
	Eigen::MatrixXcd Zout() const { return m_impedance[0]; }
	vector<Eigen::MatrixXcd> Y() const { return m_admittance; }
	Eigen::MatrixXcd Yin() const { return m_admittance.back(); }
	Eigen::MatrixXcd Yout() const { return m_admittance[0]; }
	vector<Eigen::MatrixXcd> Q() const { return m_axialVelocity; }
	Eigen::MatrixXcd Qin() const { return m_axialVelocity[0]; }
	Eigen::MatrixXcd Qout() const { return m_axialVelocity.back(); }
	vector<Eigen::MatrixXcd> P() const { return m_acPressure; }
	Eigen::MatrixXcd Pin() const { return m_acPressure[0]; }
	Eigen::MatrixXcd Pout() const { return m_acPressure.back(); }

protected:

	vector<int> m_previousSections;
	vector<int> m_nextSections;
	Point2D m_ctrLinePt;
	Point2D m_normal;
	double m_area;
	int m_modesNumber;
	double m_maxCutOnFreq;
	vector<Eigen::MatrixXcd> m_impedance;
	vector<Eigen::MatrixXcd> m_admittance;
	vector<Eigen::MatrixXcd> m_axialVelocity;
	vector<Eigen::MatrixXcd> m_acPressure;
	bool m_computeImpedance;

};

/////////////////////////////////////////////////////////////////////////////
// classe Cross section 2d FEM
/////////////////////////////////////////////////////////////////////////////

class CrossSection2dFEM : public CrossSection2d
{
  // **************************************************************************
  // Public functions.
  // **************************************************************************

public:

	CrossSection2dFEM(double maxCutOnFreq, Point2D ctrLinePt, Point2D normal,
		double area, double spacing,
		Polygon_2 contour, vector<int> surfacesIdx,
		double inMeshDensity, double inLength, double scalingFactors[2]);
  ~CrossSection2dFEM();

  void setJunctionSection(bool junction);
  void setCurvatureRadius(double& radius);
  void setCurvatureAngle(double& angle);

  void buildMesh();
  void computeModes(struct simulationParameters simuParams);
  Matrix interpolateModes(vector<Point> pts);
  void setMatrixF(vector<Matrix> & F) { m_F = F; }
  // Set the area of the intersection with the following contour
  void setIntersectionsArea(vector<double> areaInt) { m_intersectionsArea = areaInt; }
  void setMatrixGstart(Matrix & Gs){ m_Gstart = Gs; }
  void setMatrixGend(Matrix Ge) {m_Gend = Ge;}
  void characteristicImpedance(Eigen::MatrixXcd& characImped, 
	  double freq, struct simulationParameters simuParams);
  void characteristicAdmittance(Eigen::MatrixXcd& admit, 
	  double freq, struct simulationParameters simuParams);
  complex<double> getWallAdmittance(struct simulationParameters simuParams, double freq);
  void getSpecificBndAdm(struct simulationParameters simuParams, double freq, 
	  Eigen::VectorXcd& bndSpecAdm);

  // propagation
  double curvature(bool curved);
  double scaling(double tau);
  void propagateMagnus(Eigen::MatrixXcd Q0, struct simulationParameters simuParams,
	  double freq, double direction, enum physicalQuantity quant);
  void propagateImpedAdmiteRiccati(Eigen::MatrixXcd Z0, Eigen::MatrixXcd Y0, 
	  struct simulationParameters simuParams,
	  double freq, double direction, std::chrono::duration<double>& time);
  void propagateAdmitRiccati(Eigen::MatrixXcd Y0, struct simulationParameters simuParams, double nextArea,
	  double freq, double direction);
  void propagatePressureVelocityRiccati(Eigen::MatrixXcd V0, Eigen::MatrixXcd P0, 
	  struct simulationParameters simuParams, double nextArea, double freq, double direction);
  void propagatePressureRiccati(Eigen::MatrixXcd P0,
	  struct simulationParameters simuParams, double nextArea, double freq);
  void propagateImpedAdmitStraight(Eigen::MatrixXcd Z0,
	  Eigen::MatrixXcd Y0, double freq, struct simulationParameters simuParams,
	  double prevArea, double nextArea);
  void propagatePressureVelocityStraight(Eigen::MatrixXcd V0,
	  Eigen::MatrixXcd P0, double freq, struct simulationParameters simuParams,
	  double nextArea);

  // **************************************************************************
  // accessors

  double scaleIn() const { return m_scalingFactors[0];  }
  double scaleOut() const { return m_scalingFactors[1]; }
  double length() const;
  double curvRadius() const { return m_curvatureRadius; }
  vector<double> intersectionsArea() const;
  double circleArcAngle() const { return m_circleArcAngle; }
  double spacing() const;
  int numberOfVertices() const;
  int numberOfFaces() const;
  CDT triangulation() const;
  Polygon_2 contour() const;
  bool isJunction() const;
  vector<int> surfaceIdx() const;
  double eigenFrequency(int idxMode) const;
  vector<array<double, 2>> getPoints() const;
  vector<array<int, 3>> getTriangles() const;
  Matrix getModes() const;
  double getMaxAmplitude(int idxMode) const;
  double getMinAmplitude(int idxMode) const;
  vector<Matrix> getMatrixF() const;
  Matrix getMatrixGStart() const;
  Matrix getMatrixGEnd() const;


  // **************************************************************************
  // Private data.
  // **************************************************************************

private:

	double m_scalingFactors[2];
	double m_curvatureRadius;
	double m_circleArcAngle;
	CDT m_mesh;
	vector<array<double, 2>> m_points;
	vector<array<int, 3>> m_triangles;
	vector<array<int, 2>> m_meshContourSeg;
	Polygon_2 m_contour;
	double m_perimeter;
	bool m_junctionSection;
	vector<int> m_surfaceIdx;				// surface indexes of cont pts
	vector<int> m_surfIdxList;				// list of different surf idx
	double m_length;
	double m_meshDensity;
	vector<double> m_intersectionsArea;
	double m_spacing;
	vector<double> m_eigenFreqs;
	Matrix m_modes;
	vector<double> m_maxAmplitude;
	vector<double> m_minAmplitude;
	vector<Matrix> m_F;
	Matrix m_Gstart;
	Matrix m_Gend;
	Matrix m_C;
	Matrix m_DN;
	//vector<Matrix> m_DR;
	Matrix m_E;
	vector<Matrix> m_KR2;


  // **************************************************************************
  // Private functions.
  // **************************************************************************

private:

	//void initializeCrossSection(double inputUpProf[VocalTract::NUM_PROFILE_SAMPLES],
	//	double inputLoProf[VocalTract::NUM_PROFILE_SAMPLES], int upperProfileSurface[VocalTract::NUM_PROFILE_SAMPLES],
	//	int lowerProfileSurface[VocalTract::NUM_PROFILE_SAMPLES]);

};

/////////////////////////////////////////////////////////////////////////////
// classe Cross section 2d radiation
/////////////////////////////////////////////////////////////////////////////

class CrossSection2dRadiation : public CrossSection2d
{
// **************************************************************************
// Public functions.
// **************************************************************************

public:

	CrossSection2dRadiation(double maxCutOnFreq, Point2D ctrLinePt, Point2D normal,
		double radius, double PMLThickness);
	~CrossSection2dRadiation() { ; }

	//void set(double maxCutOnFreq, Point2D ctrLinePt, Point2D normal,
	//	double radius, double PMLThickness);

	void computeModes(struct simulationParameters simuParams);
	Matrix interpolateModes(vector<Point> pts);

	void characteristicImpedance(
		Eigen::MatrixXcd& characImped, double freq, struct simulationParameters simuParams);
	void characteristicAdmittance(Eigen::MatrixXcd &admit, 
		double freq, struct simulationParameters simuParams);
	complex<double> getWallAdmittance(struct simulationParameters simuParams,
		double freq) {return complex<double>(); }
	void getSpecificBndAdm(struct simulationParameters simuParams, double freq,
		Eigen::VectorXcd& bndSpecAdm) {;}

	// propagation
	void propagateMagnus(Eigen::MatrixXcd Q0, struct simulationParameters simuParams,
		double freq, double direction, enum physicalQuantity quant);
	void propagateImpedAdmiteRiccati(Eigen::MatrixXcd Z0, Eigen::MatrixXcd Y0, 
		struct simulationParameters simuParams,
		double freq, double direction, std::chrono::duration<double>& time);
	void propagateImpedRiccati(Eigen::MatrixXcd Z0, double nextArea, double freq);
	void propagateAdmitRiccati(Eigen::MatrixXcd Y0, struct simulationParameters simuParams, double nextArea,
		double freq, double direction);
	void propagatePressureVelocityRiccati(Eigen::MatrixXcd V0, Eigen::MatrixXcd P0, 
		struct simulationParameters simuParams, double nextArea,
		double freq, double direction);
	void propagatePressureRiccati(Eigen::MatrixXcd P0,
		struct simulationParameters simuParams, double nextArea, double freq);
	void propagateImpedAdmitStraight(Eigen::MatrixXcd Z0,
		Eigen::MatrixXcd Y0, double freq, struct simulationParameters simuParams,
		double prevArea, double nextArea);
	void propagatePressureVelocityStraight(Eigen::MatrixXcd V0,
		Eigen::MatrixXcd P0, double freq, struct simulationParameters simuParams,
		double nextArea);
	// to get the amplitude of the pressure modes at a given distance from the exit
	void radiatePressure(double distance, double freq, 
		struct simulationParameters simuParams, Eigen::MatrixXcd& pressAmp);

	// **************************************************************************
	// Accessors

	double scaleIn() const { return 1.; }
	double scaleOut() const { return 1.; }
	double radius() const { return m_radius; }
	double PMLThickness() const { return m_PMLThickness; }
	double BesselZero(int m) const { return m_BesselZeros[m]; }
	int BesselOrder(int m) const { return m_BesselOrder[m]; }


// **************************************************************************
// Private data.
// **************************************************************************

private:

	double m_radius;
	double m_PMLThickness;
	vector<double> m_BesselZeros;
	vector<int> m_BesselOrder;
	vector<bool> m_degeneration;
	vector<double> m_normModes;
	SparseMatC m_CPML;
	SparseMatC m_DPML;
	Eigen::MatrixXcd m_eigVec, m_invEigVec;
	Eigen::VectorXcd m_eigVal;

// **************************************************************************
// Private functions.
// **************************************************************************

	void setBesselParam(double soundSpeed);
};



#endif