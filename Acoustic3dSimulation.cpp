#include "Acoustic3dSimulation.h"
#include "Constants.h"
#include "Dsp.h"
#include "TlModel.h"
#include "TdsModel.h"
#include <algorithm>
#include <chrono>		// to get the computation time
#include <ctime>  

// for Eigen
#include <Eigen/Dense>

// for CGAL
//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include "Delaunay_mesh_vertex_base_with_info_2.h"
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Aff_transformation_2.h>
#include <CGAL/Delaunay_mesher_no_edge_refinement_2.h>
#include <CGAL/Polyline_simplification_2/simplify.h>

// type for Eigen
typedef Eigen::MatrixXd Matrix;

// Types for CGAL
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_mesh_vertex_base_with_info_2<unsigned int, K>    Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Exact_intersections_tag                     Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds, Itag> CDT;
typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;
typedef CDT::Point										Point;
typedef CGAL::Point_3<K>								Point_3;
typedef CGAL::Vector_2<K>								Vector;
typedef CGAL::Polygon_2<K>                            Polygon_2;
typedef CGAL::Polygon_with_holes_2<K>                 Polygon_with_holes_2;
typedef std::list<Polygon_with_holes_2>               Pwh_list_2;
typedef CGAL::Aff_transformation_2<K>				Transformation;
typedef CGAL::Delaunay_mesher_no_edge_refinement_2<CDT, Criteria> MesherNoRefine;
typedef CGAL::Polyline_simplification_2::Stop_below_count_ratio_threshold Stop;
typedef CGAL::Polyline_simplification_2::Squared_distance_cost Cost;

// ****************************************************************************
// Independant functions
// ****************************************************************************

// ****************************************************************************
// Generate the coordinates of the points necessary for Gauss integration
// in each triangle of mesh

void gaussPointsFromMesh(vector<Point> &pts, vector<double> & areaFaces, const CDT &cdt)
{
	int numPts(cdt.number_of_vertices());
	double quadPtCoord[3][2]{ {1. / 6., 1. / 6.}, {2. / 3., 1. / 6.}, {1. / 6., 2. / 3.} };

	pts.clear();
	pts.reserve(numPts);
	areaFaces.clear();
	areaFaces.reserve(cdt.number_of_faces());

	for (auto itF = cdt.finite_faces_begin(); itF != cdt.finite_faces_end(); itF++)
	{
		// compute the area of the face
		areaFaces.push_back(abs(
			itF->vertex(0)->point().x() * (itF->vertex(1)->point().y() - itF->vertex(2)->point().y())
			+ itF->vertex(1)->point().x() * (itF->vertex(2)->point().y() - itF->vertex(0)->point().y())
			+ itF->vertex(2)->point().x() * (itF->vertex(0)->point().y() - itF->vertex(1)->point().y())) / 2);
	
		// create the Gauss integration points
		for (int g(0); g < 3; g++)
		{
			pts.push_back(Point(
			(1 - quadPtCoord[g][0] - quadPtCoord[g][1]) * itF->vertex(0)->point().x()
				+ quadPtCoord[g][0] * itF->vertex(1)->point().x()
				+ quadPtCoord[g][1] * itF->vertex(2)->point().x(),
				(1 - quadPtCoord[g][0] - quadPtCoord[g][1]) * itF->vertex(0)->point().y()
				+ quadPtCoord[g][0] * itF->vertex(1)->point().y()
				+ quadPtCoord[g][1] * itF->vertex(2)->point().y()));
		}
	}
}

// ****************************************************************************
// Constructor.
// ****************************************************************************

Acoustic3dSimulation::Acoustic3dSimulation()
// initialise the physical constants
	: m_geometryImported(false),
	m_meshDensity(10.),
	m_maxCutOnFreq(25000.),
	m_spectrumLgthExponent(10),
	m_idxSecNoiseSource(46), // for /sh/ 212, for vowels 46
	m_idxConstriction(40),
	m_glottisBoundaryCond(IFINITE_WAVGUIDE),
	m_contInterpMeth(BOUNDING_BOX)
{
	m_simuParams.temperature = 31.4266;
	m_simuParams.volumicMass = STATIC_PRESSURE_CGS * MOLECULAR_MASS / (GAS_CONSTANT *
		(m_simuParams.temperature + KELVIN_SHIFT));
	m_simuParams.numIntegrationStep = 3;
	m_simuParams.propMethod = MAGNUS;
	m_simuParams.freqDepLosses = true;
	m_simuParams.wallLosses = true;
	m_simuParams.sndSpeed = (sqrt(ADIABATIC_CONSTANT * STATIC_PRESSURE_CGS / m_simuParams.volumicMass));
	m_crossSections.reserve(2 * VocalTract::NUM_CENTERLINE_POINTS);
	m_simuParams.percentageLosses = 1.;
	m_simuParams.curved = true;
	m_simuParams.varyingArea = true;
	m_simuParams.maxComputedFreq = (double)SAMPLING_RATE/2.;
	m_numFreq = 1 << (m_spectrumLgthExponent - 1);
	spectrum.setNewLength(2 * m_numFreq);
	spectrumNoise.setNewLength(2 * m_numFreq);
	spectrumConst.setNewLength(2 * m_numFreq);

	setBoundarySpecificAdmittance();
}

// ****************************************************************************
// Set the boundary specific admittance depending if frequency dependant losses
// are taken into account or not

void Acoustic3dSimulation::setBoundarySpecificAdmittance()
{
	if (m_simuParams.freqDepLosses)
	{
		//********************************************
		// compute boundary specific admittances 
		// for visco-thermal losses
		//********************************************

		// characteristic viscous boundary layer length
		double lv(AIR_VISCOSITY_CGS / m_simuParams.volumicMass / m_simuParams.sndSpeed);

		// characteristic thermal boundary layer length
		double lt(HEAT_CONDUCTION_CGS * MOLECULAR_MASS / m_simuParams.volumicMass / m_simuParams.sndSpeed /
		SPECIFIC_HEAT_CGS);

		// viscous boundary specific admittance
		m_simuParams.viscousBndSpecAdm = complex<double>(1., 1.) * sqrt(M_PI * lv / m_simuParams.sndSpeed);

		// thermal boundary specific admittance
		m_simuParams.thermalBndSpecAdm = complex<double>(1., 1.) * sqrt(M_PI * lt / m_simuParams.sndSpeed)
		* (ADIABATIC_CONSTANT - 1.);
	}
	else
	{

		//********************************************
		// Simple boundary specific admittance
		//********************************************

		m_simuParams.viscousBndSpecAdm = complex<double>(0., 0.);
		m_simuParams.thermalBndSpecAdm = complex<double>(0.005, 0.);

	}
}

// ****************************************************************************
/// Destructor.
// ****************************************************************************

Acoustic3dSimulation::~Acoustic3dSimulation()
{
	//delete m_crossSections;
}

// ****************************************************************************
// Static data.
// ****************************************************************************

Acoustic3dSimulation *Acoustic3dSimulation::instance = NULL;


// ****************************************************************************
/// Returns the one instance of this class.
// ****************************************************************************

Acoustic3dSimulation *Acoustic3dSimulation::getInstance()
{
	if (instance == NULL)
	{
		instance = new Acoustic3dSimulation();
	}
	return instance;
}

// ****************************************************************************
// ****************************************************************************
// Set computation parameters

void Acoustic3dSimulation::setSimulationParameters(double meshDensity, double maxCutOnFreq,
	int secNoiseSource, int secConstriction, int expSpectrumLgth, struct simulationParameters simuParams)
{
	m_meshDensity = meshDensity;
	m_maxCutOnFreq = maxCutOnFreq;
	m_idxSecNoiseSource = secNoiseSource;
	m_idxConstriction = secConstriction;
	m_spectrumLgthExponent = expSpectrumLgth;
	m_simuParams = simuParams;

	m_numFreq = 1 << (m_spectrumLgthExponent - 1);
	spectrum.reset(2 * m_numFreq);
	spectrumNoise.reset(2 * m_numFreq);
	spectrumConst.reset(2 * m_numFreq);

	setBoundarySpecificAdmittance();
}

// ****************************************************************************
// Add a cross-section at the end

void Acoustic3dSimulation::addCrossSectionFEM(double area,
	double spacing,
	Polygon_2 contours, vector<int> surfacesIdx,
	double length, Point2D ctrLinePt, Point2D normal, double scalingFactors[2])
{
	m_crossSections.push_back(unique_ptr< CrossSection2d>(new CrossSection2dFEM(m_maxCutOnFreq, 
		ctrLinePt, normal, area, spacing, contours, surfacesIdx,
		m_meshDensity, length, scalingFactors)));

}

void Acoustic3dSimulation::addCrossSectionRadiation(Point2D ctrLinePt, Point2D normal,
	double radius, double PMLThickness)
{
	m_crossSections.push_back(unique_ptr< CrossSection2d>(new
		CrossSection2dRadiation(m_maxCutOnFreq, ctrLinePt, normal,
			radius, PMLThickness)));
}

// ****************************************************************************
// create the meshes and compute the propagation modes
void Acoustic3dSimulation::computeMeshAndModes()
{
	// Create the progress dialog
	progressDialog = new wxGenericProgressDialog("Modes computation progress",
		"Wait until the modes computation finished or press [Cancel]",
		m_crossSections.size(), NULL,
		wxPD_CAN_ABORT | wxPD_AUTO_HIDE | wxPD_ELAPSED_TIME);

	ofstream mesh;
	ofstream log;
	CDT cdt;
	log.open("log.txt", ofstream::app);
	//log << "Start computing modes" << endl;

	auto start = std::chrono::system_clock::now();
	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;

	for (int i(0); i < m_crossSections.size(); i++)
	{
		start = std::chrono::system_clock::now();
		m_crossSections[i]->buildMesh();
		end = std::chrono::system_clock::now();
		elapsed_seconds = end - start;
		log << "Sec " << i << " mesh, nb vertices: "
			<< m_crossSections[i]->numberOfVertices()
			<< " time: " << elapsed_seconds.count() << " s ";

		//// export a specific mesh as text file
		//if (i == 210)
		//{
		//	mesh.open("mesh.txt");
		//	cdt = m_crossSections[i]->triangulation();
		//	for (auto it = cdt.faces_begin(); it != cdt.faces_end(); it++)
		//	{
		//		for (int v(0); v < 3; v++)
		//		{
		//			mesh << it->vertex(v)->point().x() << "  "
		//				<< it->vertex(v)->point().y() << endl;
		//		}
		//		mesh << it->vertex(0)->point().x() << "  "
		//			<< it->vertex(0)->point().y() << endl;
		//		mesh << "nan  nan" << endl;
		//	}
		//	mesh.close();
		//}

		start = std::chrono::system_clock::now();
		m_crossSections[i]->computeModes(m_simuParams);
		end = std::chrono::system_clock::now();
		elapsed_seconds = end - start;
		log << m_crossSections[i]->numberOfModes() 
			<< " modes computed, time: "
			<< elapsed_seconds.count() << " s" << endl;

		// stop simulation if [Cancel] is pressed
		if (progressDialog->Update(i) == false)
		{
			progressDialog->Destroy();
			progressDialog = NULL;

			return;
		}
	}

	log.close();

	progressDialog->Update(VocalTract::NUM_CENTERLINE_POINTS);

	// destroy progress dialog
	progressDialog->Destroy();
	progressDialog = NULL;
}

// ****************************************************************************
// Compute the junction matrices between the different cross-sections
void Acoustic3dSimulation::computeJunctionMatrices(bool computeG)
{
	Polygon_2 contour, nextContour, intersecCont, prevContour;
	Pwh_list_2 intersections, differences;
	vector<Point> pts;
	vector<double> areaFaces, areaInt;
	CDT cdt;
	double spacing, scaling[2], areaDiff;
	Vector u, v;
	Matrix interpolation1, interpolation2;
	int nModes, nModesNext, nextSec, prevSec;
	double quadPtCoord[3][2]{ {1. / 6., 1. / 6.}, {2. / 3., 1. / 6.}, {1. / 6., 2. / 3.} };
	double quadPtWeight = 1. / 3.;
	vector<Matrix> matrixF;
	int tmpNextcontained, tmpPrevContained;
	Pwh_list_2 differencesLoc;
	vector<Matrix> matrixGStart;
	vector<Matrix> matrixGEnd;
	vector<int> nextContained;
	vector<int> prevContained;
	vector<Point> seeds;
	seeds.push_back(Point(0., 0.));

	// Create the progress dialog
	progressDialog = new wxGenericProgressDialog("Junction matrices computation progress",
		"Wait until the junction matrices computation finished or press [Cancel]",
		m_crossSections.size(), NULL,
		wxPD_CAN_ABORT | wxPD_AUTO_HIDE | wxPD_ELAPSED_TIME);

	ofstream log;
	log.open("log.txt", ofstream::app);
	log << "Start computing junction matrices" << endl;

	// loop over the cross-section
	for (int i(0); i < m_crossSections.size(); i++)
	{
		//log << "\nSection " << i << " num next Sec " 
		//	<< m_crossSections[i]->numNextSec() << endl;

		if (m_crossSections[i]->numNextSec() > 0)
		{
			nModes = m_crossSections[i]->numberOfModes();

			matrixF.clear();
			areaInt.clear();

			// get the contour of the current cross-section
			contour.clear();
			scaling[0] = m_crossSections[i]->scaleOut();
			Transformation scale(CGAL::SCALING, scaling[0]);
			contour = transform(scale, m_crossSections[i]->contour());
			//contour = m_crossSections[i]->contour();
			log << "Contour sec " << i << " scaling out " << scaling[0] << endl;

			if (computeG) {
				nextContained.clear();
				nextContained.reserve(1);
				prevContained.clear();
				prevContained.reserve(m_crossSections[i]->numNextSec());
				differences.clear();
				differencesLoc.clear();
				differencesLoc.push_back(Polygon_with_holes_2(contour));
				nextContained.push_back(-1);
			}

			// loop over the folowing connected cross-sections
			for (int ns(0); ns < m_crossSections[i]->numNextSec(); ns++)
			{
				nextSec = m_crossSections[i]->nextSec(ns);
				nModesNext = m_crossSections[nextSec]->numberOfModes();

				//log << "Sec " << i << "nModes " << nModes 
				//	<< " next sec " << nextSec << " nModesNext " << nModesNext << endl;

				Matrix F(Matrix::Zero(nModes, nModesNext));
				areaInt.push_back(0.);

				// get the next contour
				nextContour.clear();
				scaling[1] = m_crossSections[nextSec]->scaleIn();
				Transformation scale(CGAL::SCALING, scaling[1]);
				nextContour = transform(scale, m_crossSections[nextSec]->contour());
				//nextContour = m_crossSections[nextSec]->contour();
				//log << "Next contour, sec " << nextSec << " scaling in " << scaling[1] << endl;

				//////////////////////////////////////////////////////////////
				// Compute the intersections of the contours
				//////////////////////////////////////////////////////////////

				// compute the intersection between the contours of the current

				//log << "Next sec " << typeid(*m_crossSections[nextSec]).name()
				//	<< " is radiation "
				//	<< (typeid(*m_crossSections[nextSec]) == typeid(CrossSection2dRadiation))
				//	<< endl;

				intersections.clear();
				if ((typeid(*m_crossSections[nextSec]) == typeid(CrossSection2dRadiation))
					|| m_crossSections[i]->isJunction())
				{
					intersections.push_back(Polygon_with_holes_2(contour));
					//log << "Contour copied" << endl;
				}
				else if (m_crossSections[nextSec]->isJunction())
				{
					intersections.push_back(Polygon_with_holes_2(nextContour));
					//log << "Next contour copied" << endl;
				}
				else
				{
					CGAL::intersection(contour, nextContour, std::back_inserter(intersections));
					//log << "Intersection computed" << endl;
				}


				if (computeG) {
					//log.open("log.txt", ofstream::app);
					//log << "Before compute difference Gend" << endl;
					// compute the difference between this contour and the next one
					differences.clear();
					for (auto itD = differencesLoc.begin();
						itD != differencesLoc.end(); itD++)
					{
						// check if the next contour is contained inside without intersecting
						tmpNextcontained = nextSec;
						// loop over the points of the next contour
						for (auto itP = nextContour.begin(); itP != nextContour.end(); itP++)
						{
							if ((itD->outer_boundary()).has_on_unbounded_side(*itP))
							{
								tmpNextcontained = -1;
								break;
							}
						}
						if (tmpNextcontained != -1) { nextContained.back() = tmpNextcontained; }

						CGAL::difference(itD->outer_boundary(), nextContour, std::back_inserter(differences));

					}
					differencesLoc = differences;
					//log << "Difference Gend computed" << endl;
				}

				spacing = min(scaling[0] * m_crossSections[i]->spacing(), 
					m_crossSections[nextSec]->spacing());

				log << "Nb intersections " << intersections.size() << endl;

				if (intersections.size() > 0)
				{

					// loop over the intersection contours
					for (auto it = intersections.begin(); it != intersections.end(); ++it)
					{
						// add the area of the intersection contour to the total 
						// area of the intersection between both sections
						areaInt.back() += it->outer_boundary().area();

						//////////////////////////////////////////////////////////////
						// Mesh the intersection surfaces and generate integration points
						//////////////////////////////////////////////////////////////

						// mesh the intersection contours
						intersecCont.clear();
						pts.clear();
						areaFaces.clear();
						cdt.clear();
						intersecCont = it->outer_boundary();
						cdt.insert_constraint(intersecCont.begin(), intersecCont.end(), true);
						CGAL::refine_Delaunay_mesh_2(cdt, Criteria(0.125, spacing));

						// remove the faces which lies outside of the contour
						for (auto itF = cdt.finite_faces_begin();
							itF != cdt.finite_faces_end(); ++itF)
						{
							if (!itF->is_in_domain())
							{
								cdt.delete_face(itF);
							}
						}

						gaussPointsFromMesh(pts, areaFaces, cdt);

						//log << pts.size() << " integration points created" << endl;

						//////////////////////////////////////////////////////////////
						// Interpolate modes
						//////////////////////////////////////////////////////////////

						// interpolate the modes of the first cross-section
						interpolation1 = m_crossSections[i]->interpolateModes(pts
							,1. / scaling[0]
						);

						log << "interpolation of first section done " 
							<< interpolation1.rows() << "  " 
							<< interpolation1.cols() << endl;

						// interpolate the modes of the next cross-section
						interpolation2 = m_crossSections[nextSec]->interpolateModes(pts
							,1. / scaling[1]
						);

						log << "interpolation of second section done "
							<< interpolation2.rows() << "  "
							<< interpolation2.cols() << endl;

						//////////////////////////////////////////////////////////////
						// Compute scatering matrix F
						//////////////////////////////////////////////////////////////

						// loop over faces to integrate the modes product
						for (int f(0); f < areaFaces.size(); f++)
						{
							if (areaFaces[f] != 0)
							{
								// loop over modes of the first cross-section
								for (int m(0); m < nModes; m++)
								{
									// loop over modes of the next cross-section
									for (int n(0); n < nModesNext; n++)
									{
										// loop over the Gauss integration points
										for (int g(0); g < 3; g++)
										{
											F(m, n) += areaFaces[f] * interpolation1(f * 3 + g,m) *
												interpolation2(f * 3 + g,n) * quadPtWeight
												/ scaling[0] / scaling[1]
												;
										}
									}
								}
							}
						}
					}
				}

				log << "\nF\n" << F << endl << endl;

				matrixF.push_back(F);
			}
			m_crossSections[i]->setIntersectionsArea(areaInt);
			m_crossSections[i]->setMatrixF(matrixF);
		}

		//log << "matrix F computed" << endl;



		if (computeG) {
			//////////////////////////////////////////////////////////////
			// Compute matrices G.
			//////////////////////////////////////////////////////////////
			

			//////////////////////////////////////////////////////////////
			// MATRIX G CORRESPONDING TO THE END OF THE CURRENT SECTION
			//////////////////////////////////////////////////////////////

			if (m_crossSections[i]->numNextSec() > 0)
			{
				matrixGEnd.clear();
				Matrix Ge(nModes, nModes);
				Ge = Matrix::Zero(nModes, nModes);
				areaDiff = 0.;

				//////////////////////////////////////////////////////////////
				// Mesh differences between contours
				//////////////////////////////////////////////////////////////

				// loop over the difference contours
				for (auto it = differencesLoc.begin();
					it != differencesLoc.end(); it++)
				{
					// add the area of this part of the difference to the 
					// total area of the difference
					areaDiff += it->outer_boundary().area();

					cdt.clear();
					cdt.insert_constraint((it->outer_boundary()).begin(),
						(it->outer_boundary()).end(), true);

					//log.open("log.txt", ofstream::app);
					//log << it->outer_boundary().area() << endl;
					//log << "next contained " << nextContained[0] << endl;
					//log.close();
					if (nextContained[0] != -1)
					{
						// if the next contour is completely contained in this one
						nextContour.clear();
						nextContour = m_crossSections[nextContained[0]]->contour();
						cdt.insert_constraint(nextContour.begin(), nextContour.end(), true);
						MesherNoRefine mesher(cdt, Criteria(0., spacing));

						// generate seed to avoid meshing inside the next contour
						seeds.back() = Point(nextContour.begin()->x() + VocalTract::PROFILE_SAMPLE_LENGTH / 2,
							nextContour.begin()->y());
						mesher.set_seeds(seeds.begin(), seeds.end(), false);
						mesher.refine_mesh();
					}
					else
					{
						MesherNoRefine mesher(cdt, Criteria(0., spacing));
						mesher.refine_mesh();
					}

					// remove the faces which lies outside of the contour
					for (auto itF = cdt.finite_faces_begin();
						itF != cdt.finite_faces_end(); ++itF)
					{
						if (!itF->is_in_domain())
						{
							cdt.delete_face(itF);
						}
					}

					gaussPointsFromMesh(pts, areaFaces, cdt);

					//log << "Before interpolation " << endl;

					interpolation1 = m_crossSections[i]->interpolateModes(pts, scaling[0]);

					//log << "after interpolation " << endl;

					//////////////////////////////////////////////////////////////
					// Compute matrix G
					//////////////////////////////////////////////////////////////

					// loop over faces to integrate the modes product
					for (int f(0); f < areaFaces.size(); f++)
					{
						if (areaFaces[f] != 0)
						{
							// loop over modes
							for (int m(0); m < nModes; m++)
							{
								// loop over modes
								for (int n(m); n < nModes; n++)
								{
									// loop over the Gauss integration points
									for (int g(0); g < 3; g++)
									{
										Ge(m, n) += areaFaces[f] * interpolation1(f * 3 + g, m) *
											interpolation1(f * 3 + g, n) * quadPtWeight;
									}
									if (m != n) { Ge(n, m) = Ge(m, n); }
								}
							}
						}
					}
				}

				Ge = (Matrix::Identity(nModes, nModes) - Ge).fullPivLu().inverse();

				//log.open("Ge.txt", ofstream::app);
				//log << "Section " << i << endl;
				//log << Ge << endl << endl;
				//log.close();

				m_crossSections[i]->setMatrixGend(Ge);
			}

			if (m_crossSections[i]->numPrevSec() > 0)
			{
				//////////////////////////////////////////////////////////////
				// COMPUTE THE DIFFERENCE WITH THE PREVIOUS SECTION
				//////////////////////////////////////////////////////////////

				differences.clear();
				differencesLoc.clear();
				differencesLoc.push_back(Polygon_with_holes_2(contour));
				prevContained.clear();
				prevContained.push_back(-1);

				//log << "Before compute difference Gstart" << endl;

				// loop over the preceding connected cross-sections
				for (int ps(0); ps < m_crossSections[i]->numPrevSec(); ps++)
				{
					prevSec = m_crossSections[i]->prevSec(ps);
					prevContour = m_crossSections[prevSec]->contour();

					// compute the difference between the contour and the previous one
					differences.clear();
					for (auto it = differencesLoc.begin();
						it != differencesLoc.end(); it++)
					{
						// check if the previous contour is contained inside without 
						// intersecting
						tmpPrevContained = prevSec;
						// loop over the points of the previous contour
						for (auto itP = prevContour.begin();
							itP != prevContour.end(); itP++)
						{
							if ((it->outer_boundary()).has_on_unbounded_side(*itP))
							{
								tmpPrevContained = -1;
								break;
							}
						}
						if (tmpPrevContained != -1) { prevContained.back() = tmpPrevContained; }

						CGAL::difference(it->outer_boundary(), prevContour,
							std::back_inserter(differences));
					}
					differencesLoc = differences;

					//log << "Difference Gstart computed" << endl;
				}

				//////////////////////////////////////////////////////////////
				// MATRIX G CORRESPONDING TO THE BEGINING OF THE SECTION
				//////////////////////////////////////////////////////////////

				matrixGStart.clear();
				Matrix Gs(nModes, nModes);
				Gs = Matrix::Zero(nModes, nModes);
				areaDiff = 0.;

				// loop over the difference contours
				for (auto it = differencesLoc.begin();
					it != differencesLoc.end(); it++)
				{
					// add the area of this part of the difference to the 
					// total area of the difference
					areaDiff += it->outer_boundary().area();

					// mesh the surface
					cdt.clear();
					cdt.insert_constraint((it->outer_boundary()).begin(),
						(it->outer_boundary()).end(), true);

					//log << "Prev contained: " << prevContained[cn] << endl;

					if (prevContained[0] != -1)
					{
						// if the previous contour is completely contained inside 
						// the current one
						prevContour = m_crossSections[prevContained[0]]->contour();
						cdt.insert_constraint(prevContour.begin(), prevContour.end(), true);
						MesherNoRefine mesher(cdt, Criteria(0., spacing));

						// generate seed to avoid meshing inside the previous contour
						seeds.back() = Point(prevContour.begin()->x() + VocalTract::PROFILE_SAMPLE_LENGTH / 2,
							prevContour.begin()->y());
						mesher.set_seeds(seeds.begin(), seeds.end(), false);
						mesher.refine_mesh();
					}
					else
					{
						MesherNoRefine mesher(cdt, Criteria(0., spacing));
						mesher.refine_mesh();
					}

					// remove the faces which are outside of the domain
					for (auto itF = cdt.finite_faces_begin();
						itF != cdt.finite_faces_end(); itF++)
					{
						if (!itF->is_in_domain())
						{
							cdt.delete_face(itF);
						}
					}

					// generate the integration points
					gaussPointsFromMesh(pts, areaFaces, cdt);

					// interpolate the propagation modes
					interpolation2 = m_crossSections[i]->interpolateModes(pts);

					// loop over the faces to integrate the modes product
					for (int f(0); f < areaFaces.size(); f++)
					{
						if (areaFaces[f] != 0.)
						{
							// loop over modes
							for (int m(0); m < nModes; m++)
							{
								// loop over modes
								for (int n(m); n < nModes; n++)
								{
									// loop over Gauss integration points
									for (int g(0); g < 3; g++)
									{
										Gs(m, n) += areaFaces[f] * interpolation2(f * 3 + g, m)
											* interpolation2(f * 3 + g,n) * quadPtWeight;
									}
									if (m != n) { Gs(n, m) = Gs(m, n); }
								}
							}
						}
					}
					//log.open("Gs1.txt", ofstream::app);
					//log << "Section " << i << endl;
					//log << Gs << endl << endl;
					//log.close();
				}
				Gs = (Matrix::Identity(nModes, nModes) - Gs).fullPivLu().inverse();

				//log.open("Gs.txt", ofstream::app);
				//log << "Section " << i << endl;
				//log << Gs << endl << endl;
				//log.close();
				m_crossSections[i]->setMatrixGstart(Gs);
			}
		}

		if (progressDialog->Update(i) == false)
		{
			progressDialog->Destroy();
			progressDialog = NULL;

			return;
		}

	}
	log.close();

	progressDialog->Update(VocalTract::NUM_CENTERLINE_POINTS);

	// destroy progress dialog
	progressDialog->Destroy();
	progressDialog = NULL;
}

// **************************************************************************
// Propagate the impedance and admittance up to the other end of the geometry
// taking into account branches

void Acoustic3dSimulation::propagateImpedAdmitBranch(vector< Eigen::MatrixXcd> Q0, double freq,
	vector<int> startSections, vector<int> endSections, double direction)
{
	bool addSegsToList, isNotEndSeg, isNotInList;
	int m, n, idx, mn, ns;
	vector<int> prevSegs, nextSegs;
	vector<vector<int>> segToProp;
	vector<Matrix> Ftmp;
	Eigen::MatrixXcd Qout, Qini;
	Matrix F;
	complex<double> wallInterfaceAdmit(1i * 2. * M_PI * freq *
		m_simuParams.thermalBndSpecAdm / m_simuParams.sndSpeed);

	ofstream log("log.txt", ofstream::app);
	log << "Start branches" << endl;

	// initialise the list of segment lists to propagate
	for (auto it : startSections)
	{
		segToProp.push_back(vector<int>());
		segToProp.back().push_back(it);
	}

	// while there are non propagated segments in the list
	ns = 0;
	while (ns < segToProp.size())
	{
		log << "ns = " << ns << endl;
		log << "Segs ";
		for (auto it : segToProp[ns])
		{
			log << it << "  ";
		}
		log << endl;
		//for (auto it : segToProp)
		//{
		//	log << "{";
		//	for (auto seg : it)
		//	{
		//		log << seg << " ";
		//	}
		//	log << "} ";
		//}
		//log << endl;

		// if the segment is an initial segment
		if (ns < startSections.size())
		{
			if (m_crossSections[segToProp[ns][0]]->computeImpedance())
			{
				//log << "Propagate impedance" << endl;
				m_crossSections[segToProp[ns][0]]->propagateMagnus(Q0[ns], m_simuParams, freq, direction, IMPEDANCE);
			}
			else
			{
				m_crossSections[segToProp[ns][0]]->propagateMagnus(Q0[ns], m_simuParams, freq, direction, ADMITTANCE);
			}
			//log << "Impedance propagated" << endl;
		}
		else
		{
			// get the list of previous sections
			if (direction > 0)
			{
				prevSegs = m_crossSections[segToProp[ns][0]]->prevSections();
			}
			else
			{
				prevSegs = m_crossSections[segToProp[ns][0]]->nextSections();
			}
			log << "Prevsegs: ";
			for (auto it : prevSegs)
			{
				log << it << "  ";
			}
			log << endl;

			//***********************************
			// if the previous segment is larger
			//***********************************

			if (m_crossSections[segToProp[ns][0]]->area() < m_crossSections[prevSegs[0]]->area())
			{
				// in this case there can be only one segment connected to the current segment

				//log << "previous segment is larger" << endl;

				// Get the mode matching matrices
				if (direction > 0)
				{
					Ftmp = m_crossSections[prevSegs[0]]->getMatrixF();
				}
				else
				{
					// for each segment of the segment group
					Ftmp.clear();
					for (auto it : segToProp[ns])
					{
						Ftmp.push_back(m_crossSections[it]->getMatrixF()[0].transpose());
					}
				}
				// determine the dimension m,n of the concatenated mode matching matrix
				m = Ftmp[0].rows();
				n = 0;
				for (auto it : Ftmp) { n += it.cols(); }
				F.resize(m, n);
				// concatenate the mode matching matrices
				for (auto it : Ftmp) { F << it; }

				// get the output impedance of the previous segment
				// if the admittance have been computed in this segment
				if (!m_crossSections[prevSegs[0]]->computeImpedance())
				{
					// compute the corresponding impedance
					Qout = m_crossSections[prevSegs[0]]->Yin().fullPivLu().inverse();
				}
				else
				{
					Qout = m_crossSections[prevSegs[0]]->Zin();
				}

				log << "F\n" << F << endl << endl;

				// Compute the input impedance
				Qini = F.transpose() * Qout * F;

				log << "Qini\n" << Qini.cwiseAbs() << endl << endl;

				// propagate the impedance in each of the connected tube
				idx = 0;
				for (auto it : segToProp[ns])
				{
					// get the number of modes
					mn = m_crossSections[it]->numberOfModes();
					// propagate the impedance, the initial impedance of each connected
					// tube is a submatrix of Qini
					m_crossSections[it]->propagateMagnus(Qini.block(idx, idx, mn, mn),
						m_simuParams, freq, direction, IMPEDANCE);
					m_crossSections[it]->setComputImpedance(true);
					idx += mn;
				}
			}

			//***************************************
			// if the previous segment(s) is smaller
			//***************************************

			else 
			{
				//log << "previous segment(s) is smaller" << endl;

				//get the mode matching matrices
				Ftmp.clear();
				if (direction > 0)
				{
					for (auto it : prevSegs)
					{
						Ftmp.push_back(m_crossSections[it]->getMatrixF()[0].transpose());
					}
				}
				else
				{
					Ftmp = m_crossSections[segToProp[ns][0]]->getMatrixF();
				}
				// determine the dimension m,n of the concatenated mode matching matrix
				m = Ftmp[0].rows();
				n = 0;
				for (auto it : Ftmp) { n += it.cols(); }
				F.resize(m, n);
				// concatenate the mode matching matrices
				for (auto it : Ftmp) { F << it; }

				log << "F\n" << F << endl << endl;

				// build the output admittance matrix of all the previous segments
				Qout.setZero(n, n);
				idx = 0;
				for (auto it : prevSegs)
				{
					// get the mode number of the previous segment
					mn = m_crossSections[it]->numberOfModes();
					if (m_crossSections[it]->computeImpedance())
					{
						Qout.block(idx, idx, mn, mn) =
							m_crossSections[it]->Zin().fullPivLu().inverse();
					}
					else
					{
						Qout.block(idx, idx, mn, mn) = m_crossSections[it]->Yin();
					}
					idx += mn;
				}

				log << "Qout\n" << Qout.cwiseAbs() << endl << endl;

				// compute the input admittance matrix
				Qini = F * Qout * F.transpose();
				m_crossSections[segToProp[ns][0]]->propagateMagnus(Qini, m_simuParams, freq,
					direction, ADMITTANCE);
				m_crossSections[segToProp[ns][0]]->setComputImpedance(false);
			}
		}

		//***********************************************************************
		// add the following connected segments to list of segments to propagate
		//***********************************************************************

		for (auto it : segToProp[ns])
		{
			//log << "Check connection segment " << it << endl;

			// check if the segment is an end segment
			isNotEndSeg = true;
			for (auto endSec : endSections)
			{
				if (it == endSec) { isNotEndSeg = false; break; }
			}

			// get the indexes of the next connected segments
			if (direction > 0)
			{
				nextSegs = m_crossSections[it]->nextSections();
			}
			else
			{
				nextSegs = m_crossSections[it]->prevSections();
			}

			// check if there are connected segments 
			if (nextSegs.size() > 0) {

				// Check if the next segments are already in the list of segments 
				// to propagate
				isNotInList = true;
				for (int i(ns); i < segToProp.size(); i++)
				{
					if (segToProp[i][0] == nextSegs[0]) { isNotInList = false; break; }
				}

				if (isNotEndSeg && isNotInList)
				{
					//log << "Next segs ";
					//for (auto seg : nextSegs)
					//{
					//	log << seg;
					//}
					//log << endl;
					// check if the connected segments can be added to the list
					//
					// it is necessary to check only for the first segment of the 
					// list since either there is only one segment connected to several
					// others, or there is several segments, but connected to the same

					if (nextSegs.size() > 1)
					{
						addSegsToList = true;
					}
					else
					{
						if (direction > 0)
						{
							prevSegs = m_crossSections[nextSegs[0]]->prevSections();
						}
						else
						{
							prevSegs = m_crossSections[nextSegs[0]]->nextSections();
						}

						for (auto prevSeg : prevSegs)
						{
							//log << "Prev seg " << prevSeg << endl;

							addSegsToList = false;
							// check if this index is among all the segment indexes which 
							// have already been propagated
							for (int i(ns); i > -1; i--)
							{

								//log << "i = " << i << endl;
								for (auto propSeg : segToProp[i])
								{
									//log << "propSeg " << propSeg << endl;
									if (propSeg == prevSeg) { addSegsToList = true; break; }
								}
								if (addSegsToList) { break; }
							}
							if (!addSegsToList) { break; }
						}
						//log << "Propagated " << propagated << endl;

					}
					// if all the connected segments have been propagated
					// and the segments are not already in the list
					if (addSegsToList)
					{
						segToProp.push_back(nextSegs);
					}
					//log << "Next segmensts added" << endl;
				}
			}
		}
		if (ns < segToProp.size()) { ns++; }
	}
	log.close();
}

// **************************************************************************
// Propagate the impedance and admittance up to the other end of the geometry

void Acoustic3dSimulation::propagateImpedAdmit(Eigen::MatrixXcd & startImped,
	Eigen::MatrixXcd & startAdmit, double freq, int startSection, int endSection)
{
	Eigen::MatrixXcd prevImped;
	Eigen::MatrixXcd prevAdmit;
	vector<Matrix> F;
	//vector<double> areaInt;
	Matrix G;
	int numSec(m_crossSections.size()), nI, nPs;
	int direction, prevSec;
	double areaRatio;
	complex<double> wallInterfaceAdmit(1i*2.*M_PI*freq* 
		m_simuParams.thermalBndSpecAdm/m_simuParams.sndSpeed);
	vector<Eigen::MatrixXcd> inputImped;

	std::chrono::duration<double> time;

	// determine the direction of the propagation
	if (startSection > endSection)
	{
		direction = -1;
	}
	else
	{
		direction = 1;
	}
	

	//ofstream log;
	//log.open("log.txt", ofstream::app);
	//log << "start propagating admittance" << endl;
	//log << "Direction " << direction << endl;

	// set the initial impedance and admittance matrices
	m_crossSections[startSection]->clearImpedance();
	m_crossSections[startSection]->clearAdmittance();

	//log << "Admit imped cleared" << endl;
	
	switch(m_simuParams.propMethod)
	{
	case MAGNUS:
		//m_crossSections[startSection]->propagateImpedAdmiteRiccati(startImped, startAdmit, 
		//	m_simuParams, freq, (double)direction, time);
		//m_crossSections[startSection]->setComputImpedance(true);
		m_crossSections[startSection]->propagateMagnus(startAdmit, m_simuParams,
			freq, (double)direction, ADMITTANCE);
		inputImped.clear();
		for (auto it : m_crossSections[startSection]->Y())
		{
			inputImped.push_back(it.fullPivLu().inverse());
		}
		m_crossSections[startSection]->setImpedance(inputImped);
		break;
	case STRAIGHT_TUBES:
		m_crossSections[startSection]->propagateImpedAdmitStraight(startImped, startAdmit,
			freq, m_simuParams, 100.,
			m_crossSections[startSection+direction]->area());
		break;
	}
	
	//log << "First admittance propagated" << endl;

	// loop over sections
	for (int i(startSection + direction); i != (endSection + direction); i += direction)
	{
		//log << "sec " << i << " " << typeid(*m_crossSections[i]).name() << endl;
		//log << "sec " << i << endl;

		prevSec = i - direction;
		m_crossSections[i]->clearImpedance();
		m_crossSections[i]->clearAdmittance();

		//log << "Prev sec " << prevSec << endl;

		nI = m_crossSections[i]->numberOfModes();
		nPs = m_crossSections[prevSec]->numberOfModes();

		// Extract the scaterring matrix and its complementary
		F.clear();
		if (direction == -1)
		{
			F = m_crossSections[i]->getMatrixF();
			if (m_crossSections[i]->area() > m_crossSections[prevSec]->area())
			{
				G = Matrix::Identity(nI, nI) - F[0] * F[0].transpose();
			}
			else
			{
				G = Matrix::Identity(nPs, nPs) - F[0].transpose() * F[0];
			}
		}
		else
		{
			F = m_crossSections[prevSec]->getMatrixF();
			if (m_crossSections[i]->area() > m_crossSections[prevSec]->area())
			{
				G = Matrix::Identity(nI, nI) - F[0].transpose() * F[0];
			}
			else
			{
				G = Matrix::Identity(nPs, nPs) - F[0] * F[0].transpose();
			}
		}

		prevImped = Eigen::MatrixXcd::Zero(m_crossSections[i]->numberOfModes(),
			m_crossSections[i]->numberOfModes());
		prevAdmit = Eigen::MatrixXcd::Zero(m_crossSections[i]->numberOfModes(),
			m_crossSections[i]->numberOfModes());
		
		switch (m_simuParams.propMethod)
		{
		case MAGNUS:
			// case of a contraction: area(i) > area(ps)
			if (m_crossSections[i]->area() > m_crossSections[prevSec]->area())
			{
				if (direction == -1) // ps = i + 1, n_X > 0
				{ 
					//log << "area(i) > area(ps), dir = -1" << endl;
					prevAdmit += F[0] *
						m_crossSections[prevSec]->Yin()
						* (F[0].transpose())
						//- wallInterfaceAdmit*
						//////m_crossSections[i]->curvature() *
						//G
						;
					//m_crossSections[i]->setComputImpedance(false);
				}
				else // ps = i - 1, n_X < 0
				{
					//log << "area(i) > area(ps), dir = 1" << endl;
					prevAdmit += (F[0].transpose()) *
						m_crossSections[prevSec]->Yin() * F[0]
						//+ wallInterfaceAdmit * 
						//////m_crossSections[i]->curvature() * 
						//G
						;
					//m_crossSections[i]->setComputImpedance(false);
				}
				prevImped += prevAdmit.fullPivLu().inverse();
			}
			// case of an expansion: area(i) < area(ps)
			else
			{
				if (direction == -1) // ps = i + 1; n_X < 0
				{
					prevImped += F[0] * m_crossSections[prevSec]->Zin()
						//*(Matrix::Identity(nPs, nPs) - 
						//	wallInterfaceAdmit * 
						//	////m_crossSections[prevSec]->curvature() * 
						//	G*m_crossSections[prevSec]->Zin()).inverse()
						* (F[0].transpose());
					//m_crossSections[i]->setComputImpedance(false);
				}
				else // ps = i - 1; n_X > 0
				{
					//log << "area(i) < area(ps), dir = 1" << endl;
					prevImped += (F[0].transpose()) * m_crossSections[prevSec]->Zin()
						//*(Matrix::Identity(nPs, nPs) + 
						//	wallInterfaceAdmit *
						//	////m_crossSections[prevSec]->curvature() * 
						//	G*m_crossSections[prevSec]->Zin()).inverse()
						* F[0];
					//m_crossSections[i]->setComputImpedance(false);
				}
				prevAdmit += prevImped.fullPivLu().inverse();
			}

			break;
		case STRAIGHT_TUBES:

			areaRatio = max(m_crossSections[prevSec]->area(),
				m_crossSections[i]->area()) /
				min(m_crossSections[prevSec]->area(),
					m_crossSections[i]->area());
			// case of a contraction
			if (m_crossSections[i]->area() > m_crossSections[prevSec]->area())
			{
				prevAdmit += areaRatio * F[0] *
					m_crossSections[prevSec]->Yin()
					* (F[0].transpose());
				prevImped += prevAdmit.fullPivLu().inverse();
			}
			// case of an expansion
			else
			{
				prevImped += areaRatio * F[0] *
					m_crossSections[prevSec]->Zin()
					* (F[0].transpose());
				prevAdmit += prevImped.fullPivLu().inverse();
			}
			break;
		}

		//log << "prevImped computed" << endl;

		// propagate admittance in the section
		switch (m_simuParams.propMethod) {
		case MAGNUS:
			//m_crossSections[i]->propagateImpedAdmiteRiccati(prevImped, prevAdmit, m_simuParams,
			//	freq, (double)direction, time);
			//m_crossSections[i]->setComputImpedance(false);
			m_crossSections[i]->propagateMagnus(prevAdmit, m_simuParams,
				freq, (double)direction, ADMITTANCE);
			inputImped.clear();
			for (auto it : m_crossSections[i]->Y())
			{
				inputImped.push_back(it.fullPivLu().inverse());
			}
			m_crossSections[i]->setImpedance(inputImped);
			break;
		case STRAIGHT_TUBES:
			m_crossSections[i]->propagateImpedAdmitStraight(prevImped, prevAdmit,
				freq, m_simuParams, m_crossSections[prevSec]->area(),
				m_crossSections[i + direction]->area());
			break;
		}

		//log << "section " << i << " propagated" << endl;
	}
	//log.close();
}

// **************************************************************************
// Propagate the admittance up to the other end of the geometry

void Acoustic3dSimulation::propagateAdmit(Eigen::MatrixXcd radAdmit, double freq)
{
	Eigen::MatrixXcd prevAdmit;
	vector<Matrix> F;
	Matrix invF;
	int numSec(m_crossSections.size());
	double areaRatio;

	//ofstream log;
	//log.open("log.txt", ofstream::app);
	//log << "Start propagate admit" << endl;

	// compute admittance in the last section
	m_crossSections[numSec - 1]->clearAdmittance();
	m_crossSections[numSec - 1]->propagateAdmitRiccati(radAdmit, m_simuParams,
		m_crossSections[numSec - 1]->area(), freq, -1.);

	//log << "First admit propagated" << endl;

	// loop over sections
	for (int i(numSec - 2); i > 0; i--)
	{
		//log << "Section " << i << endl;

		m_crossSections[i]->clearAdmittance();

		// get scatering matrix
		F.clear();
		F = m_crossSections[i]->getMatrixF();

		// case of a contraction
		if (m_crossSections[i]->area() > m_crossSections[i + 1]->area())
		{
			//log << "Contraction " << endl;
			prevAdmit = F[0] * m_crossSections[i + 1]->Yin()
				* (F[0].transpose());
		}
		// case of an expansion
		else
		{
			//log << "Expension " << endl;
			invF = F[0].inverse();
			//log << "Size iF " << invF.rows() << " " << invF.cols() << endl;
			//log << "Size iFt " << invF.transpose().rows() << " "
			//	<< invF.transpose().cols() << endl;
			//log << "Size Y " << m_crossSections[i + 1]->startAdmittance().rows()
			//	<< " " << m_crossSections[i + 1]->startAdmittance().cols() << endl;
			prevAdmit = invF * 
				m_crossSections[i + 1]->Yin() * invF.transpose();
			//log << "Size prevadmit " << prevAdmit.rows() << " " << prevAdmit.cols() << endl;
		}
		//log << "Scatering computed" << endl;
		m_crossSections[i]->propagateAdmitRiccati(prevAdmit, m_simuParams,
			m_crossSections[i + 1]->area(), freq, -1.);
		//log << "Admittance propagated" << endl;
	}
	//log.close();
}

// **************************************************************************
// Propagate the axial velocity up to the other end of the geometry

void Acoustic3dSimulation::propagateVelocityPress(Eigen::MatrixXcd &startVelocity,
	Eigen::MatrixXcd& startPressure, double freq, int startSection, 
	int endSection)
{
	Eigen::MatrixXcd prevVelo(startVelocity), prevPress(startPressure);
	vector<Eigen::MatrixXcd> tmpQ, P, Y;
	vector<Matrix> F;
	Matrix G;
	int numSec(m_crossSections.size());
	int numX(m_simuParams.numIntegrationStep);
	int direction, nextSec, nI, nNs;
	Eigen::MatrixXcd pressure;
	double areaRatio;
	complex<double> wallInterfaceAdmit(1i * 2. * M_PI * freq * 
		m_simuParams.thermalBndSpecAdm / m_simuParams.sndSpeed);

	//ofstream ar;
	//ofstream log;
	//log.open("log.txt", ofstream::app);
	//log << "\n\nStart velocity propagation" << endl;

	// determine the direction of the propagation
	if (startSection > endSection)
	{
		direction = -1;
	}
	else
	{
		direction = 1;
	}

	//log << "Direction " << direction << endl;
	
	// loop over sections
	for (int i(startSection); i != (endSection); i += direction)
	{

		nextSec = i + direction;

		//log << "section " << i << " next sec " << nextSec << endl;

		m_crossSections[i]->clearAxialVelocity();
		m_crossSections[i]->clearAcPressure();
		nI = m_crossSections[i]->numberOfModes();
		nNs = m_crossSections[nextSec]->numberOfModes();

		// propagate axial velocity and acoustic pressure in the section
		switch (m_simuParams.propMethod) {
		case MAGNUS:
			//m_crossSections[i]->propagatePressureVelocityRiccati(prevVelo, prevPress,m_simuParams, 
			//	m_crossSections[nextSec]->area(), freq, (double)direction);
			m_crossSections[i]->propagateMagnus(prevPress, m_simuParams,
				freq, (double)direction, PRESSURE);
			tmpQ.clear(); P.clear(); Y.clear();
			P = m_crossSections[i]->P();
			Y = m_crossSections[i]->Y();
			numX = Y.size();
			for (int pt(0); pt < numX; pt++ )
			{
				tmpQ.push_back(Y[numX - 1 - pt] * P[pt]);
			}
			m_crossSections[i]->setAxialVelocity(tmpQ);
			break;
		case STRAIGHT_TUBES:
			m_crossSections[i]->propagatePressureVelocityStraight(prevVelo,
				prevPress, freq, m_simuParams,
				m_crossSections[nextSec]->area());
			break;
		}

		// get the scattering matrix 
		F.clear();
		if (direction == 1)
		{
			F = m_crossSections[i]->getMatrixF();
			if (m_crossSections[i]->area() > m_crossSections[nextSec]->area())
			{
				G = Matrix::Identity(nI, nI) - F[0] * F[0].transpose();
			}
			else
			{
				G = Matrix::Identity(nNs, nNs) - F[0].transpose() * F[0];
			}
		}
		else
		{
			F = m_crossSections[nextSec]->getMatrixF();
			if (m_crossSections[i]->area() > m_crossSections[nextSec]->area())
			{
				G = Matrix::Identity(nI, nI) - F[0].transpose() * F[0];
			}
			else
			{
				G = Matrix::Identity(nNs, nNs) - F[0] * F[0].transpose();
			}
		}

		prevVelo = Eigen::MatrixXcd::Zero(m_crossSections[nextSec]->numberOfModes(), 1);
		prevPress = Eigen::MatrixXcd::Zero(m_crossSections[nextSec]->numberOfModes(), 1);
		switch (m_simuParams.propMethod)
		{
		case MAGNUS:
			// if the section contracts: area(i) > area(ns)
			if (m_crossSections[i]->area() >
				m_crossSections[nextSec]->area())
			{
				if (direction == -1)
				{
					prevPress += F[0] *
						m_crossSections[i]->Pout();
				}
				else
				{
					prevPress += (F[0].transpose()) *
						m_crossSections[i]->Pout();
				}
				prevVelo +=
					m_crossSections[nextSec]->Yin() * prevPress;
			}
			// if the section expends: area(i) < area(ns)
			else
			{
				if (direction == -1)
				{
					prevVelo +=
						//(Matrix::Identity(nNs, nNs)
						//	+ wallInterfaceAdmit *
						//	////m_crossSections[nextSec]->curvature()*
						//	G * m_crossSections[nextSec]->Zin()).inverse() *
						F[0] * m_crossSections[i]->Qout();
				}
				else
				{
					prevVelo +=
						//(Matrix::Identity(nNs, nNs)
						//	- wallInterfaceAdmit *
						//	////m_crossSections[nextSec]->curvature() * 
						//	G * m_crossSections[nextSec]->Zin()).inverse() *
							(F[0].transpose()) * m_crossSections[i]->Qout();
				}
				prevPress +=
					m_crossSections[nextSec]->Zin() * prevVelo;
			}

			break;
		case STRAIGHT_TUBES:
			//areaRatio = sqrt(m_crossSections->at(i + 1).area() /
			//	m_crossSections->at(i).area());

			areaRatio = sqrt(max(m_crossSections[nextSec]->area(),
				m_crossSections[i]->area()) /
				min(m_crossSections[nextSec]->area(),
					m_crossSections[i]->area()));
			// if the section expends
			if (m_crossSections[nextSec]->area() >
				m_crossSections[i]->area())
			{
				prevVelo += areaRatio * (F[0].transpose()) *
					m_crossSections[i]->Qout();
				prevPress +=
					m_crossSections[nextSec]->Zin() * prevVelo;
			}
			// if the section contracts
			else
			{
				prevPress += areaRatio * (F[0].transpose()) *
					m_crossSections[i]->Pout();
				prevVelo +=
					m_crossSections[nextSec]->Yin() * prevPress;
			}
			break;
		}
	}

	// propagate in the last section
	m_crossSections[endSection]->clearAxialVelocity();
	m_crossSections[endSection]->clearAcPressure();
	switch (m_simuParams.propMethod) {
	case MAGNUS:
		//m_crossSections[endSection]->propagatePressureVelocityRiccati(prevVelo, prevPress,
		//	m_simuParams, 2.*m_crossSections[endSection]->area(), freq,
		//	(double)direction);
		m_crossSections[endSection]->propagateMagnus(prevPress, m_simuParams,
			freq, (double)direction, PRESSURE);
		tmpQ.clear(); P.clear(); Y.clear();
		Y = m_crossSections[endSection]->Y();
		P = m_crossSections[endSection]->P();
		numX = Y.size();
		for (int pt(0); pt < numX; pt++)
		{
			tmpQ.push_back(Y[numX - 1 - pt] * P[pt]);
		}
		m_crossSections[endSection]->setAxialVelocity(tmpQ);
		break;
	case STRAIGHT_TUBES:
		m_crossSections[endSection]->propagatePressureVelocityStraight(prevVelo,
			prevPress, freq, m_simuParams, 100.);
		break;
	}

	//log.close();
}

void Acoustic3dSimulation::propagatePressure(Eigen::MatrixXcd startPress, double freq)
{
	Eigen::MatrixXcd prevPress(startPress);
	vector<Matrix> F;
	int numSec(m_crossSections.size());

	ofstream log;
	log.open("log.txt", ofstream::app);

	// loop over sections
	for (int i(1); i < numSec - 1; i++)
	{
		// propagate pressure in the section
		m_crossSections[i]->clearAxialVelocity();
		m_crossSections[i]->propagatePressureRiccati(prevPress, m_simuParams,
			m_crossSections[i + 1]->area(), freq);

		// get scatering matrix
		F.clear();
		F = m_crossSections[i]->getMatrixF();

		// if the section expends
		if (m_crossSections[i + 1]->area() >
			m_crossSections[i]->area())
		{
			//log << "Size F " << F[0].rows() << " " << F[0].cols() << endl;
			//log << "Size iF " << F[0].inverse().rows() << " "
			//	<< F[0].inverse().cols() << endl;
			//log << "Size Pend " << m_crossSections[i]->endAcPressure().rows()
			//	<< " " << m_crossSections[i]->endAcPressure().cols() << endl;
			//Eigen::MatrixXcd tF = (Eigen::MatrixXcd)(F[0].transpose());
			//Eigen::MatrixXcd Pa = m_crossSections[i]->endAcPressure();
			//prevPress = tF.fullPivLu().solve(Pa);
			//prevPress = F[0].transpose().fullPivLu().inverse() * m_crossSections[i]->endAcPressure();
			prevPress = m_crossSections[i+1]->Yin().fullPivLu().inverse() * 
				F[0].transpose() * m_crossSections[i]->Yout() * 
				m_crossSections[i]->Pout();

			//log << "Size prevpress " << prevPress.rows() << " "
			//	<< prevPress.cols() << endl;
		}
		// if the section contracts
		else
		{
			prevPress = (F[0].transpose()) *
				m_crossSections[i]->Pout();
		}
	}
	// propagate in the last section
	m_crossSections[numSec - 1]->clearAcPressure();
	m_crossSections[numSec - 1]->propagatePressureRiccati(prevPress,
		m_simuParams, m_crossSections[numSec - 1]->area(), freq);
}

// **************************************************************************
// Run a static simulation

void Acoustic3dSimulation::staticSimulation(VocalTract* tract)
{
	int numSec; 
	double freq, freqSteps((double)SAMPLING_RATE/2./(double)m_numFreq);
	int numFreqComputed((int)ceil(m_simuParams.maxComputedFreq / freqSteps));
	complex<double> pressure, velocity;
	Eigen::MatrixXcd endPressAmpl, endVelAmpl, startPressure,
		radAdmit, radImped, upStreamImpAdm, totalImped, prevVelo, prevPress
		,pp ,pm;
	Matrix F;


	ofstream prop;
	ofstream log("log.txt", ofstream::out | ofstream::trunc);
	log.close();
	log.open("log.txt", ofstream::app);
	log << "Start simulation\n" << endl;
	log << "Air volumic mass: " << m_simuParams.volumicMass << " g/cm^3" << endl;
	log << "Sound speed: " << m_simuParams.sndSpeed << " cm/s" << endl;
	log << "viscous boundary specific admittance " << m_simuParams.viscousBndSpecAdm 
		<< " g.cm^-2 .s^-1"<< endl;
	log << "thermal boundary specific admittance " << m_simuParams.thermalBndSpecAdm
		<< " g.cm^-2 .s^-1" << endl;
	log << "Percentage losses " << m_simuParams.percentageLosses * 100. << " %" << endl;
	log << "Frequency dependant losses: ";
	if (m_simuParams.freqDepLosses)
	{
		log << "yes" << endl;
	}
	else
	{
		log << "no" << endl;
	}
	log << "Wall losses: ";
	if (m_simuParams.wallLosses)
	{
		log << "yes" << endl;
	}
	else
	{
		log << "no" << endl;
	}
	log << "glottis boundary condition: ";
	switch (m_glottisBoundaryCond)
	{
	case HARD_WALL:
		log << "HARD_WALL" << endl;
		break;
	case IFINITE_WAVGUIDE:
		log << "IFINITE_WAVGUIDE" << endl;
		break;
	}
	log << "Mesh density: " << m_meshDensity << endl;
	log << "Max cut-on frequency: " << m_maxCutOnFreq << " Hz" << endl;
	log << "Number of integration steps: " << m_simuParams.numIntegrationStep << endl;
	log << "Index of noise source section: " << m_idxSecNoiseSource << endl;
	log << "Index of constriction section: " << m_idxConstriction << endl;
	log << "Maximal computed frequency: " << m_simuParams.maxComputedFreq 
		<< " Hz" << endl;
	log << "Number of simulated frequencies: " << numFreqComputed << endl << endl;
	log << "Varying cross-sectional area: ";
	if (m_simuParams.varyingArea)
	{
		log << "yes" << endl;
	}
	else
	{
		log << "no" << endl;
	}
	log << "Propagation mmethod: ";
	switch (m_simuParams.propMethod)
	{
	case MAGNUS:
		log << "MAGNUS" << endl;
		break;
	case STRAIGHT_TUBES:
		log << "STRAIGHT_TUBES" << endl;
		break;
	}

	auto startTot = std::chrono::system_clock::now();

	//******************************
	// create the cross-sections
	//******************************

	if (m_geometryImported)
	{
		log << "Geometry imported from csv file " << m_geometryFile << endl;
	}
	else
	{
		log << "Geometry is from vocal tract lab" << endl;
	}
	if (createCrossSections(tract, true))
	{
		log << "Geometry successfully imported" << endl;
	}
	else
	{
		log << "Importation failed" << endl;
	}

	//if (m_geometryImported)
	//{
	//	if (createCSWithIntermediateCSFromCSVFile(false))
	//	{
	//		log << "Geometry successfully imported" << endl;
	//	}
	//	else
	//	{
	//		log << "Importation failed" << endl;
	//	}
	//}
	//else
	//{
	//	log << "Geometry is from vocal tract lab" << endl;
	//	createCSWithIntermediateCSGeneral(tract, true);
	//}
	
	numSec = m_crossSections.size();

	log << "Number of sections: " << numSec << endl;

	//// export cross-sections parameters
	//prop.open("sec.txt");
	//// check cross-section properties
	//for (int i(0); i < numSec - 2; i++)
	//{
	//	prop << m_crossSections[i]->ctrLinePt().x
	//		<< "  " << m_crossSections[i]->ctrLinePt().y
	//		<< "  " << m_crossSections[i]->normal().x
	//		<< "  " << m_crossSections[i]->normal().y
	//		<< "  " << m_crossSections[i]->area()
	//		<< "  " << m_crossSections[i]->length()
	//		<< "  " << m_crossSections[i]->curvature()
	//		<< "  " << m_crossSections[i]->circleArcAngle()
	//		<< "  " << m_crossSections[i]->isJunction()
	//		<< endl;

	//}
	//prop.close();

	//// check connexions of sections
	//for (int i(0); i < m_crossSections.size(); i++)
	//{
	//	log << "sec " << i << " length " << m_crossSections[i]->length() 
	//		<< endl;
	//}
	//for (int i(1); i < m_crossSections.size(); i++)
	//{
	//	log << "sec " << i
	//		<< " prev sec " << m_crossSections[i]->prevSec(0) << endl;
	//}
	//for (int i(0); i < m_crossSections.size() - 1; i++)
	//{
	//	log << "sec " << i
	//		<< " next sec " << m_crossSections[i]->nextSec(0) << endl;
	//}

	//// exctract some contours
	//ostringstream os;
	//for (int i(97); i < 103; i++)
	//{
	//	os.str("");
	//	os.clear();
	//	os << "c" << i << ".txt";
	//	prop.open(os.str());
	//	for (auto it : m_crossSections[i]->contour())
	//	{
	//		prop << it.x() << "  " << it.y() << endl;
	//	}
	//	prop.close();
	//}

	// FIXME: update the values of m_idxSecNoiseSource and m_idxConstriction
	// in the 3D simu properties dialog when they are changed
	// check if the noise source index is within the indexes range
	if (m_idxSecNoiseSource >= numSec)
	{
		m_idxSecNoiseSource = numSec - 2;
	}
	// check if the constriction location is within the indexes range
	if (m_idxConstriction >= numSec)
	{
		m_idxConstriction = numSec - 2;
	}

	//log << "before compute mode" << endl;

	// create the mesh and compute modes
	computeMeshAndModes();
	log << "Modes computed" << endl;

	// compute junction matrices
	auto start = std::chrono::system_clock::now();
	computeJunctionMatrices(false);
	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - startTot;
	log << "Junction matrices computed " << 
		elapsed_seconds.count() << " s" << endl;

	// extract the junction matrix of the noise source section
	F = m_crossSections[m_idxSecNoiseSource]->getMatrixF()[0];

	log << "noise source junction matrix extracted" << endl;

	// generate source matrix
	Eigen::MatrixXcd inputVelocity(Eigen::MatrixXcd::Zero(
	m_crossSections[0]->numberOfModes(), 1));
	//log << "inputVelocity initialized" << endl;
	Eigen::MatrixXcd inputVelocityNoise(Eigen::MatrixXcd::Zero(
		m_crossSections[0]->numberOfModes(), 1));
	//log << "inputVelocityNoise initialized" << endl;
	Eigen::MatrixXcd inputPressure(Eigen::MatrixXcd::Zero(
		m_crossSections[0]->numberOfModes(), 1));
	//log << "inputPressure initialzed" << endl;
	Eigen::MatrixXcd inputPressureNoise(Eigen::MatrixXcd::Zero(
		m_crossSections[m_idxSecNoiseSource]->numberOfModes(), 1));
	//log << "inputPressureNoise initialized" << endl;
	complex<double> V0(1.0, 0.0);
	// for a constant input velocity q = -j * w * rho * v 
	inputVelocity(0, 0) = -1i * 2. * M_PI * m_simuParams.volumicMass * V0;
	inputPressureNoise(0, 0) = V0;

	log << "source generated" << endl;

	//// compute the modal amplitude at the center of the last cross-section
	//vector<Point> pts;
	//pts.push_back(Point(0., 0.));

	//Matrix modeAmplitudes(
	//m_crossSections[numSec-2]->interpolateModes(pts));

	//*****************************************************************************
	//  Compute points for acoustic field computation
	//*****************************************************************************

	vector<Point_3> radPts;

	//// generate points distributed on a plane
	//int nPtsX(50), nPtsY(100);
	//double dX(5./(double)(nPtsX));
	//for (int i(0); i < nPtsX; i++)
	//{
	//	for (int j(0); j < nPtsY; j++)
	//	{
	//		radPts.push_back(Point_3((double)(i + 1) * dX, 
	//			((double)(j) - (double)(nPtsY)/2.) * dX, 0.));
	//	}
	//}

	//// generate points on a cicle arc
	//double radius(50.);
	//int nPts(91);
	//for (int p(0); p < nPts; p++)
	//{
	//	radPts.push_back(Point_3(radius * sin((double)p * M_PI / (double)(nPts-1)),
	//		radius * cos((double)p * M_PI / (double)(nPts-1)), 0.));
	//}

	// generate one point in front
	radPts.push_back(Point_3(25., 0., 0.));

	// print radiation points coordinates in the log file
	log << "Radiation point(s) coordinates (m):" << endl;
	for (auto const& it : radPts)
	{
		log << 0.01 * it.x() << "  "
			<< 0.01 * it.y() << "  "
			<< 0.01 * it.z() << endl;
	}

	// create vector containing the radiated field
	Eigen::VectorXcd radPress;

	// Precompute the radiation impedance and admittance at a few frequencies
	preComputeRadiationMatrices(16);

	end = std::chrono::system_clock::now();
	elapsed_seconds = end - startTot;
	log << "Time precomputations " << elapsed_seconds.count() << " s" << endl;

	//*****************************************************************************

	// Create the progress dialog
	progressDialog = new wxGenericProgressDialog("Transfer function computation",
		"Wait until the transfer function computation finished or press [Cancel]",
		numFreqComputed, NULL,
		wxPD_CAN_ABORT | wxPD_AUTO_HIDE | wxPD_ELAPSED_TIME);

	//// save the frequencies of the transfer function
	//prop.open("fSimu.txt");
	//for (int i(0); i < m_numFreq; i++)
	//{
	//	prop << max(10., (double)i* freqSteps) << endl;
	//}
	//prop.close();

	start = std::chrono::system_clock::now();

	//int mn(m_crossSections[numSec - 2]->numberOfModes());

	// loop over frequencies
	for (int i(0); i < numFreqComputed; i++)
	{
		freq = max(0.1, (double)i * freqSteps);
		log << "################################" << endl;
		log << "frequency " << i+1 << "/" << numFreqComputed << " f = " << freq
			<< " Hz" << endl;

		// generate radiation impedance and admittance
		start = std::chrono::system_clock::now();
		// m_crossSections[numSec - 1]->characteristicImpedance(radImped,
		//	 freq, m_soundSpeed, m_volumicMass, 0);
		//m_crossSections[numSec - 1]->characteristicAdmittance(radAdmit,
		//	freq, m_soundSpeed, m_volumicMass, 0);

		interpolateRadiationImpedance(radImped, freq);
		interpolateRadiationAdmittance(radAdmit, freq);

		//radImped.setConstant(mn , mn,complex<double>(0., 0.));
		//radAdmit.setConstant(mn, mn, complex<double>(1.e+10, 0.));

		end = std::chrono::system_clock::now();
		elapsed_seconds = end - start;
		log << "Time charac imped/admit " << elapsed_seconds.count() << " s" << endl;

		// propagate impedance and admittance
		start = std::chrono::system_clock::now();
		//propagateAdmit(radAdmit, freq, m_method);
		propagateImpedAdmit(radImped, radAdmit, freq, numSec-2, 0);
		
		end = std::chrono::system_clock::now();
		elapsed_seconds = end - start;
		log << "Time impedance " << elapsed_seconds.count() << " s" << endl;

		//// extract admittance
		//prop.open("adm.txt", ofstream::app);
		//for (int c(0); c < numSec - 1; c++)
		//{
		//	prop << abs(m_crossSections[c]->Yin()(0, 0)) << "  ";
		//}
		//prop << endl;
		//prop.close();

		//// print radiation impedance
		//prop.open("radR.txt", ofstream::app);
		//prop << m_crossSections[numSec-2]->endImpedance().real() << endl;
		//prop.close();
		//prop.open("radI.txt", ofstream::app);
		//prop << m_crossSections[numSec - 2]->endImpedance().imag() << endl;
		//prop.close();
		//prop.open("admR.txt", ofstream::app);
		//prop << m_crossSections[numSec - 2]->endAdmittance().real() << endl;
		//prop.close();
		//prop.open("admI.txt", ofstream::app);
		//prop << m_crossSections[numSec - 2]->endAdmittance().imag() << endl;
		//prop.close();

		// propagate axial velocity and pressure
		start = std::chrono::system_clock::now();
		//startPressure = m_crossSections[1]->startAdmittance().fullPivLu().inverse()
		//	* inputVelocity;
		//propagatePressure(startPressure, freq, m_method);
		log << "Compute impedance in first section " <<
			m_crossSections[0]->computeImpedance() << endl;
		inputVelocity(0, 0) = -1i * 2. * M_PI * freq * m_simuParams.volumicMass * 
			m_crossSections[0]->area();
		inputPressure = m_crossSections[0]->Zin() * inputVelocity;
		propagateVelocityPress(inputVelocity, inputPressure, freq, 0, numSec -2);
		end = std::chrono::system_clock::now();
		elapsed_seconds = end - start;
		log << "Time velocity pressure " << elapsed_seconds.count() << " s" << endl;

		// extract pressure
		prop.open("pamp.txt", ofstream::app);
		for (int c(0); c < numSec - 1; c++)
		{
			prop << abs(m_crossSections[c]->Pin()(0, 0)) << "  ";
		}
		prop << endl;
		prop.close();

		// extract admittance
		prop.open("adm.txt", ofstream::app);
		for (int c(0); c < numSec - 1; c++)
		{
			prop << abs(m_crossSections[c]->Yout()(0, 0)) << "  ";
		}
		prop << endl;
		prop.close();

		// Extract end volume velocity
		prop.open("uend.txt", ofstream::app);
		prop << abs(m_crossSections[numSec - 2]->Qout()(0, 0)) << endl;
		prop.close();
		
		start = std::chrono::system_clock::now();

		//*****************************************************************************
		//  Compute radiated field
		//*****************************************************************************

		RayleighSommerfeldIntegral(radPts, radPress, freq);

		//// save radiated field
		//prop.open("radR.txt", ofstream::app);
		//prop << radPress.transpose().real() << endl;
		//prop.close();
		//prop.open("radI.txt", ofstream::app);
		//prop << radPress.transpose().imag() << endl;
		//prop.close();

		//// compute the acoustic pressure at the center of the exit
		//pressure = (0., 0.);
		//velocity = (0., 0.);
		//endPressAmpl = m_crossSections[numSec - 2]->endAcPressure();
		//endVelAmpl = m_crossSections[numSec - 2]->endAxialVelocity();
		//for (int m(0); m < m_crossSections[numSec-2]->numberOfModes(); m++)
		//{
		//	pressure += endPressAmpl(m, 0) * modeAmplitudes(0, m);
		//	velocity += endVelAmpl(m, 0) * modeAmplitudes(0, m);
		//}

		/*prop << freq << "  " << abs(radPress(0, 0)) << "  " << arg(radPress(0, 0)) << endl;*/

		spectrum.setValue(i, radPress(0,0));

		log << "Radiated pressure computed" << endl;

		//*****************************************************************************
		//  Compute transfer function of the noise source
		//*****************************************************************************

		//// get the transfer function from the glottis to the constriction
		//// location for vowels

		//spectrumConst.setValue(i, m_crossSections[m_idxConstriction]->Pout()(0, 0));

		//// save the input impedance of the upstream part
		//// if the section expends
		////log << "Start save impedance" << endl;
		////log << "area next sec" << m_crossSections[m_idxSecNoiseSource + 1]->area() << endl;
		////log << "area " << m_crossSections[m_idxSecNoiseSource]->area() << endl;
		//if (m_crossSections[m_idxSecNoiseSource + 1]->area() >
		//	m_crossSections[m_idxSecNoiseSource]->area())
		//{
		//	upStreamImpAdm = m_crossSections[m_idxSecNoiseSource]->Zout();
		//}
		//// if the section contracts
		//else
		//{
		//	upStreamImpAdm = m_crossSections[m_idxSecNoiseSource]->Yout();
		//}
		////log << "upstream input impedance saved" << endl;

		//// set glottis boundary condition
		//
		//switch (m_glottisBoundaryCond)
		//{
		//case HARD_WALL:
		//	radImped.setZero();
		//	radImped.diagonal().setConstant(100000.);
		//	radAdmit.setZero();
		//	radAdmit.diagonal().setConstant(1. / 100000.);
		//	break;
		//case IFINITE_WAVGUIDE:
		//	m_crossSections[0]->characteristicImpedance(radImped, freq, m_simuParams);
		//	m_crossSections[0]->characteristicAdmittance(radAdmit, freq, m_simuParams);
		//	break;
		//}

		//// propagate impedance and admittance from the glottis to the location
		//// of the second sound source
		//propagateImpedAdmit(radImped, radAdmit, freq, 0, m_idxSecNoiseSource);

		//log << "Imped admit propagated" << endl;

		//// compute the pressure and the velocity at the entrance of the next section
		//// if the section expends
		//if (m_crossSections[m_idxSecNoiseSource +1]->area() >
		//	m_crossSections[m_idxSecNoiseSource]->area())
		//{
		//	prevVelo = (F.transpose()) * ((freq * upStreamImpAdm - freq * 
		//		m_crossSections[m_idxSecNoiseSource]->Zin()).householderQr()
		//		.solve(inputPressureNoise));
		//	prevPress = freq *
		//		m_crossSections[m_idxSecNoiseSource +1]->Zin() * prevVelo;
		//}
		//// if the section contracts
		//else
		//{
		//	prevPress = (F.transpose()) * ((upStreamImpAdm -
		//		m_crossSections[m_idxSecNoiseSource]->Yin()).householderQr()
		//		.solve( -m_crossSections[m_idxSecNoiseSource]->Yin() *
		//			inputPressureNoise));

		//	prevVelo =
		//		m_crossSections[m_idxSecNoiseSource +1]->Yin()* prevPress;
		//}

		//log << "Previous pressure/velocity computed" << endl;


		////for (int j(0); j < (numSec - 1); j++)
		////{
		////	log << m_crossSections[j]->startAcPressure()(0, 0) << "  "
		////		<< m_crossSections[j]->endAcPressure()(0, 0) << "  ";
		////}
		////log << endl;

		////log << "start sec: " << min(m_idxSecNoiseSource + 1, numSec - 2)
		////	<< " end sec: " << numSec - 2 << endl;
		//
		//// propagate the pressure and the velocity in the upstream part
		//propagateVelocityPress(prevVelo, prevPress, freq, 
		//	min(m_idxSecNoiseSource +1, numSec - 2), numSec - 2);

		////for (int j(0); j < (numSec - 1); j++)
		////{
		////	log << m_crossSections[j]->startAcPressure()(0, 0) << "  "
		////		<< m_crossSections[j]->endAcPressure()(0, 0) << "  ";
		////}
		////log << endl;

		//log << "Velocity and pressure propagated" << endl;

		//RayleighSommerfeldIntegral(radPts, radPress, freq);

		////// save radiated field
		////prop.open("radR.txt", ofstream::app);
		////prop << radPress.transpose().real() << endl;
		////prop.close();
		////prop.open("radI.txt", ofstream::app);
		////prop << radPress.transpose().imag() << endl;
		////prop.close();

		////// compute the acoustic pressure at the center of the exit
		////pressure = (0., 0.);
		////velocity = (0., 0.);
		////endPressAmpl = m_crossSections[numSec - 2]->endAcPressure();
		////endVelAmpl = m_crossSections[numSec - 2]->endAxialVelocity();
		////for (int m(0); m < m_crossSections[numSec - 2]->numberOfModes(); m++)
		////{
		////	pressure += endPressAmpl(m, 0) * modeAmplitudes(0, m);
		////	velocity += endVelAmpl(m, 0) * modeAmplitudes(0, m);
		////}


		end = std::chrono::system_clock::now();

		spectrumNoise.setValue(i, radPress(0,0));
		elapsed_seconds = end - start;
		log << "Remaining time " << elapsed_seconds.count() << " s" << endl;

		//*****************************************************************************

		if (progressDialog->Update(i) == false)
		{
			progressDialog->Destroy();
			progressDialog = NULL;

			return;
		}
	}

	// generate spectra values for negative frequencies
	log << "\n Generate spectrum negative part" << endl;
	for (int i(m_numFreq); i < 2*m_numFreq; i++)
	{
		spectrum.re[i] = spectrum.re[2 * m_numFreq - i - 1];
		spectrum.im[i] = -spectrum.im[2 * m_numFreq - i - 1];
		//spectrumNoise.re[i] = spectrumNoise.re[2 * m_numFreq - i - 1];
		//spectrumNoise.im[i] = -spectrumNoise.im[2 * m_numFreq - i - 1];
		//spectrumConst.re[i] = spectrumConst.re[2 * m_numFreq - i - 1];
		//spectrumConst.im[i] = -spectrumConst.im[2 * m_numFreq - i - 1];
	}

	// create the file containing the transfer function
	prop.open("press.txt");
	//prop << "num_points: " << spectrum.N << endl;
	//prop << "frequency_Hz  magnitude  phase_rad" << endl;
	for (int i(0); i < spectrum.N; i++)
	{
		freq = max(0.1, (double)i * freqSteps);
		prop << freq << "  " << spectrum.getMagnitude(i) << "  " 
			<< spectrum.getPhase(i) << endl;
	}
	prop.close();

	//prop.open("spec.txt");
	//for (int i(0); i < spectrum.N; i++)
	//{
	//	prop << spectrum.getRealPart(i) << "  "
	//		<< spectrum.getImaginaryPart(i) << endl;
	//}
	//prop.close();
	//
	//prop.open("specN.txt");
	//for (int i(0); i < spectrumNoise.N; i++)
	//{
	//	prop << spectrumNoise.getRealPart(i) << "  "
	//		<< spectrumNoise.getImaginaryPart(i) << endl;
	//}
	//prop.close();

	//prop.open("specC.txt");
	//for (int i(0); i < spectrumConst.N; i++)
	//{
	//	prop << spectrumConst.getRealPart(i) << "  "
	//		<< spectrumConst.getImaginaryPart(i) << endl;
	//}
	//prop.close();

	end = std::chrono::system_clock::now();
	elapsed_seconds = end - startTot;
	int hours(floor(elapsed_seconds.count() / 3600.));
	int minutes(floor((elapsed_seconds.count() - hours * 3600.) / 60.));
	double seconds(elapsed_seconds.count() - (double)hours * 3600. -
		(double)minutes * 60.);
	log << "\nTransfer function time "
		<< hours << " h " << minutes << " m " << seconds << " s" << endl;

	log.close();

	progressDialog->Update(VocalTract::NUM_CENTERLINE_POINTS);

	// destroy progress dialog
	progressDialog->Destroy();
	progressDialog = NULL;
}

// ****************************************************************************
// Compute the Rayleigh-Sommerfeld integral to compute radiated pressure

void Acoustic3dSimulation::RayleighSommerfeldIntegral(vector<Point_3> points,
	Eigen::VectorXcd& radPress, double freq)
{
	int nbPts(points.size());
	double quadPtCoord[3][2]{ {1. / 6., 1. / 6.}, {2. / 3., 1. / 6.}, {1. / 6., 2. / 3.} };
	double quadPtWeight = 1. / 3.;
	vector<Point> gaussPts;
	vector<double> areaFaces;
	double r, k(2.*M_PI*freq/m_simuParams.sndSpeed);
	complex<double> J(0., 1.);
	// the radiating cross-section is the section before the last one
	int radSecIdx(m_crossSections.size() - 2);

	//ofstream log;
	//log.open("log.txt", ofstream::app);

	// get velocity mode amplitude (v_x = j * q / w / rho)
	auto Vm = m_crossSections[radSecIdx]->Qout();

	// get quadrature points
	gaussPointsFromMesh(gaussPts, areaFaces, 
		m_crossSections[radSecIdx]->triangulation());

	// interpolate modes
	Matrix interpolatedModes = m_crossSections[radSecIdx]->
		interpolateModes(gaussPts);

	// compute Guauss integral
	radPress = Eigen::VectorXcd::Zero(nbPts);
	// loop over the faces of the mesh
	for (int f(0); f < areaFaces.size(); f++)
	{
		// loop over the modes
		for (int m(0); m < m_crossSections[radSecIdx]->numberOfModes(); m++)
		{
			// loop over points
			for (int p(0); p < nbPts; p++)
			{
				// loop over Gauss integration points
				for (int g(0); g < 3; g++)
				{
					// compute distance r between gauss point and the point
					r = sqrt(CGAL::squared_distance(points[p], 
						Point_3(0.,gaussPts[f * 3 + g].x(), gaussPts[f * 3 + g].y())));

					// compute integral
					radPress(p) -= areaFaces[f] * quadPtWeight *
						Vm(m) * interpolatedModes(f * 3 + g, m) * exp(-J * k * r) / r;
					//radPress(p) += J * areaFaces[f] * quadPtWeight * 
					//	(Vm(m).real() - J*Vm(m).imag()) *
					//	interpolatedModes(f * 3 + g, m) * exp(-J * k * r) / r;
				}
			}
		}
	}
	radPress /= 2 * M_PI;
	//log.close();
}

// ****************************************************************************
// **************************************************************************
// accessors

int Acoustic3dSimulation::sectionNumber() const { return((int)m_crossSections.size()); }
double Acoustic3dSimulation::soundSpeed() const { return(m_simuParams.sndSpeed); }
CrossSection2d * Acoustic3dSimulation::crossSection(int csIdx) const
{
	return m_crossSections[csIdx].get();
}
double Acoustic3dSimulation::meshDensity() const { return m_meshDensity; }
//int Acoustic3dSimulation::modeNumber() const { return m_modeNumber; }
double Acoustic3dSimulation::maxCutOnFreq() const { return m_maxCutOnFreq; }
int Acoustic3dSimulation::numIntegrationStep() const { return m_simuParams.numIntegrationStep; }

// **************************************************************************
// Private functions.
// **************************************************************************
//
// Create contours polygons and surface indexes lists from the upper and lower 
// profiles, and the corresponding upper and lower surface indexes
// generated by the articulatory models

void Acoustic3dSimulation::createContour(double inputUpProf[VocalTract::NUM_PROFILE_SAMPLES],
	double inputLoProf[VocalTract::NUM_PROFILE_SAMPLES], 
	int upperProfileSurface[VocalTract::NUM_PROFILE_SAMPLES],
	int lowerProfileSurface[VocalTract::NUM_PROFILE_SAMPLES],
	vector<double> &areas, vector<double> &spacing, 
	vector< Polygon_2> &contours, vector<vector<int>> &surfacesIdx)
{
	int idxContour(0);
	double temporaryUpProf[VocalTract::NUM_PROFILE_SAMPLES];
	std::fill_n(temporaryUpProf, VocalTract::NUM_PROFILE_SAMPLES, VocalTract::INVALID_PROFILE_SAMPLE);
	double temporaryLoProf[VocalTract::NUM_PROFILE_SAMPLES];
	std::fill_n(temporaryLoProf, VocalTract::NUM_PROFILE_SAMPLES, VocalTract::INVALID_PROFILE_SAMPLE);
	double tempArea(0.);
	Polygon_2 tempPoly;
	vector<int> tempSufIdx;
	tempSufIdx.reserve(3 * VocalTract::NUM_PROFILE_SAMPLES);
	double dist;
	int nIntermPts;
	Vector vecNextPt;
	Vector vecInsertPt;
	double alpha;
	bool toNewSurf(false);
	bool toNewSurfTeeth(true);

	// identify the samples between two contours as invalid samples
	for (int i(1); i < VocalTract::NUM_PROFILE_SAMPLES - 1; i++)
	{
		if ((inputUpProf[i - 1] == inputLoProf[i - 1]) && (inputUpProf[i + 1] == inputLoProf[i + 1]))
		{
			inputUpProf[i] = VocalTract::INVALID_PROFILE_SAMPLE;
			inputLoProf[i] = VocalTract::INVALID_PROFILE_SAMPLE;
		}
	}

	// initialize the first elements of the temporary profiles
	temporaryUpProf[0] = inputUpProf[0];
	temporaryLoProf[0] = inputLoProf[0];

	// create the meshes corresponding to the potentially multiple closed contours
	for (int i(1); i < VocalTract::NUM_PROFILE_SAMPLES; i++)
	{
		// store samples in the temporary profiles
		temporaryUpProf[i] = inputUpProf[i];
		temporaryLoProf[i] = inputLoProf[i];

		// compute area
		tempArea += 0.5 * (inputUpProf[i - 1] + inputUpProf[i] - inputLoProf[i - 1] - inputLoProf[i]) *
			VocalTract::PROFILE_SAMPLE_LENGTH;

		// if the contour ends, create the corresponding mesh and compute the
		// corresponding modes
		if ((inputUpProf[i - 1] != inputLoProf[i - 1])
			&& (inputUpProf[i] == inputLoProf[i])
			&& (inputUpProf[min(i + 1, VocalTract::NUM_PROFILE_SAMPLES - 1)]
				== inputLoProf[min(i + 1, VocalTract::NUM_PROFILE_SAMPLES - 1)]))
		{
			areas.push_back(tempArea);
			spacing.push_back(sqrt(tempArea) / m_meshDensity);

			// *********************************************************************
			// create the polygon corresponding to the upper part of the contour
			// *********************************************************************

			int idxBig;
			for (int p(0); p < (i + 1); p++)
			{
				if (temporaryUpProf[p] != VocalTract::INVALID_PROFILE_SAMPLE)
				{
					idxBig = p;
					for (int pt(idxBig); pt < (i + 1); pt++)
					{
						// insert new point in the polygon
						tempPoly.push_back(Point(pt * VocalTract::PROFILE_SAMPLE_LENGTH -
							VocalTract::PROFILE_LENGTH / 2, temporaryUpProf[pt]));
						tempSufIdx.push_back(upperProfileSurface[pt]);

						// check if it is necessary to add an intermediate point
						if ((pt != (VocalTract::NUM_PROFILE_SAMPLES - 1)) &&
							(temporaryUpProf[pt + 1] != VocalTract::INVALID_PROFILE_SAMPLE))
						{
							// compute the distance with the next point if it exists
							dist = sqrt(pow((temporaryUpProf[pt] - temporaryUpProf[pt + 1]), 2) +
								pow(VocalTract::PROFILE_SAMPLE_LENGTH, 2));

							// compute the number of necessary subdivision of contour segment
							nIntermPts = (int)floor(dist / VocalTract::PROFILE_SAMPLE_LENGTH / 2.) + 1;

							// if the next surface is different from the previous one
							if (upperProfileSurface[pt] !=
								upperProfileSurface[min(pt + 1, VocalTract::NUM_PROFILE_SAMPLES)])
							{
								toNewSurf = !toNewSurf;

								// if the exited surface is a  tooth
								if ((upperProfileSurface[pt] == 0) || (upperProfileSurface[pt] == 1))
								{
									toNewSurf = toNewSurfTeeth;
								}
								// if the next surface is a tooth
								else if ((upperProfileSurface[min(pt + 1, VocalTract::NUM_PROFILE_SAMPLES)] == 0)
									|| (upperProfileSurface[min(pt + 1, VocalTract::NUM_PROFILE_SAMPLES)] == 1))
								{
									// flip this value
									toNewSurfTeeth = !toNewSurfTeeth;
									toNewSurf = toNewSurfTeeth;
								}

							}

							// if necessary, add intermediate points
							if (nIntermPts > 1)
							{
								// next point on the upper profile
								vecNextPt = Vector(Point(0., 0.), Point((pt + 1) * VocalTract::PROFILE_SAMPLE_LENGTH -
									VocalTract::PROFILE_LENGTH / 2, temporaryUpProf[pt + 1]));

								for (int n(1); n < nIntermPts; n++)
								{
									alpha = 1. / (double)(nIntermPts - n + 1);
									vecInsertPt = alpha * vecNextPt +
										(1. - alpha) * Vector(Point(0., 0.), *(tempPoly.vertices_end() - 1));

									// add point to the upper profile
									tempPoly.push_back(Point(vecInsertPt.x(), vecInsertPt.y()));
									if (toNewSurf)
									{
										tempSufIdx.push_back(upperProfileSurface[min(pt + 1, VocalTract::NUM_PROFILE_SAMPLES)]);
									}
									else
									{
										tempSufIdx.push_back(tempSufIdx.back());
									}
								}
							}
						}
					}
					break;
				}
			}


			toNewSurf = true;

			// *********************************************************************
			// create the polygon corresponding to the lower part of the contour
			// *********************************************************************

			for (int p(i - 1); p > idxBig; p--)
			{
				vecNextPt = Vector(Point(0., 0.), Point(p * VocalTract::PROFILE_SAMPLE_LENGTH -
					VocalTract::PROFILE_LENGTH / 2, temporaryLoProf[p]));

				// compute the distance with the next point if it exists
				dist = sqrt(pow(vecNextPt.y() - (tempPoly.vertices_end() - 1)->y(), 2) +
					pow(VocalTract::PROFILE_SAMPLE_LENGTH, 2));

				// compute the number of necessary subdivision of contour segment
				nIntermPts = (int)floor(dist / VocalTract::PROFILE_SAMPLE_LENGTH / 2.) + 1;

				// if the next surface is different from the previous one
				if (tempSufIdx.back() != lowerProfileSurface[p])
				{
					toNewSurf = !toNewSurf;

					// if the exited surface is a  tooth
					if ((tempSufIdx.back() == 0) || (tempSufIdx.back() == 1))
					{
						toNewSurf = toNewSurfTeeth;
					}
					// if the next surface is a tooth
					else if ((lowerProfileSurface[p] == 0) || (lowerProfileSurface[p] == 1))
					{
						// flip this value
						toNewSurfTeeth = !toNewSurfTeeth;
						toNewSurf = toNewSurfTeeth;
					}

				}

				// if necessary, add intermediate points
				if (nIntermPts > 1)
				{

					for (int n(1); n < nIntermPts; n++)
					{
						alpha = 1. / (double)(nIntermPts - n + 1);
						vecInsertPt = alpha * vecNextPt +
							(1. - alpha) * Vector(Point(0., 0.), *(tempPoly.vertices_end() - 1));

						// add point to the upper profile
						tempPoly.push_back(Point(vecInsertPt.x(), vecInsertPt.y()));
						if (toNewSurf)
						{
							tempSufIdx.push_back(lowerProfileSurface[p]);
						}
						else
						{
							tempSufIdx.push_back(tempSufIdx.back());
						}
					}
				}

				// insert the next point
				tempPoly.push_back(Point(vecNextPt.x(), vecNextPt.y()));
				tempSufIdx.push_back(lowerProfileSurface[p]);
			}

			// check if it is necessary to add intermediary points in the last intervall
			vecNextPt = Vector(Point(0., 0.), Point(idxBig * VocalTract::PROFILE_SAMPLE_LENGTH -
				VocalTract::PROFILE_LENGTH / 2, temporaryLoProf[idxBig]));

			// compute the distance with the next point if it exists
			dist = sqrt(pow(vecNextPt.y() - (tempPoly.vertices_end() - 1)->y(), 2) +
				pow(VocalTract::PROFILE_SAMPLE_LENGTH, 2));

			// compute the number of necessary subdivision of contour segment
			nIntermPts = (int)floor(dist / VocalTract::PROFILE_SAMPLE_LENGTH / 2.) + 1;

			// if the next surface is different from the previous one
			if (tempSufIdx.back() != lowerProfileSurface[idxBig])
			{
				toNewSurf = !toNewSurf;

				// if the exited surface is a  tooth
				if ((tempSufIdx.back() == 0) || (tempSufIdx.back() == 1))
				{
					toNewSurf = toNewSurfTeeth;
				}
				// if the next surface is a tooth
				else if ((lowerProfileSurface[idxBig] == 0) || (lowerProfileSurface[idxBig] == 1))
				{
					// flip this value
					toNewSurfTeeth = !toNewSurfTeeth;
					toNewSurf = toNewSurfTeeth;
				}

			}

			// if necessary, add intermediate points
			if (nIntermPts > 1)
			{
				for (int n(1); n < nIntermPts; n++)
				{
					alpha = 1. / (double)(nIntermPts - n + 1);
					vecInsertPt = alpha * vecNextPt +
						(1. - alpha) * Vector(Point(0., 0.), *(tempPoly.vertices_end() - 1));

					// add point to the upper profile
					tempPoly.push_back(Point(vecInsertPt.x(), vecInsertPt.y()));
					if (toNewSurf)
					{
						tempSufIdx.push_back(lowerProfileSurface[idxBig]);
					}
					else
					{
						tempSufIdx.push_back(tempSufIdx.back());
					}
				}
			}

			// Add the created contour and the corresponding surface indexes
			// to the contour list and the surfaces index list
			contours.push_back(tempPoly);
			surfacesIdx.push_back(tempSufIdx);

			// reinitialize the temporary profiles and area
			tempArea = 0.;
			std::fill_n(temporaryUpProf, VocalTract::NUM_PROFILE_SAMPLES, VocalTract::INVALID_PROFILE_SAMPLE);
			std::fill_n(temporaryLoProf, VocalTract::NUM_PROFILE_SAMPLES, VocalTract::INVALID_PROFILE_SAMPLE);
			tempPoly.clear();
			tempSufIdx.clear();

			idxContour++;
		}
	}
}

//*****************************************************************************
// Determine the radius of curvatur, the angles of start and end of the 
// section and the necessary shift of the contour
//
// It is necessary to shift the contour because the radius obtained at the
// start and the end of the section are not exactely the same, so the 
// contour needs to be shifted to be at the average radius

void Acoustic3dSimulation::getCurvatureAngleShift(Point2D P1, Point2D P2, 
	Point2D N1, Point2D N2, double& radius, double& angle, double& shift)
{
	// inputs:
	//
	//	P1			centerline point of section 1
	//	P2			centerline point of section 2
	//	N1			centerline normal of section 1
	//	N2			centerline normal of section 2
	//
	// outputs:
	//
	//	radius		radius of circle arc joining sec 1 & 2
	//	angle		angle between sections 1 & 2
	//	shift		distance necessary to shift P2 along N2
	//				so that it is on the same circle arc as P1

	double radii[2], angles[2];

	// compute the radius of the circle arc whose center
	// is the intersection point of the normals and which
	// passes by P1, the point belonging to section i
	radii[0] = ((P2.y - P1.y) * N2.x - (P2.x - P1.x) * N2.y) /
		(N2.x * N1.y - N2.y * N1.x);

	// compute the radius of the circle arc whose center
	// is the intersection point of the normals and which
	// passes by P2, the point belonging to the next section i+1
	radii[1] = ((P2.y - P1.y) * N1.x - (P2.x - P1.x) * N1.y) /
		(N2.x * N1.y - N2.y * N1.x);

	// compute the angle corresponding to the position of P1
	angles[0] = atan2(-N1.y, -N1.x);
	angles[1] = atan2(-N2.y, -N2.x);

	// compute the average radius
	radius = (radii[0] + radii[1]) / 2.;

	// compute the vertical shift to apply to the contour
	shift = radii[0] - radius;

	// compute the angle between both angular positions
	angle = max(angles[0], angles[1]) - min(angles[0], angles[1]);
}

//*****************************************************************************
// Extract the contours, the surface indexes, the centerline and the normals
// from the VocalTractLab geometry

void Acoustic3dSimulation::extractContours(VocalTract* tract, vector<vector<Polygon_2>>& contours,
	vector<vector<vector<int>>>& surfaceIdx, vector<Point2D>& centerLine, vector<Point2D>& normals)
{
	// for cross profile data extraction
	double upperProfile[VocalTract::NUM_PROFILE_SAMPLES];
	double lowerProfile[VocalTract::NUM_PROFILE_SAMPLES];
	int upperProfileSurface[VocalTract::NUM_PROFILE_SAMPLES];
	int lowerProfileSurface[VocalTract::NUM_PROFILE_SAMPLES];
	Tube::Articulator articulator;

	// for contour data
	vector<double> areas, spacings;
	vector<Polygon_2> tmpContours;
	vector<vector<int>> tmpSurfacesIdx;


	for (int i(0); i < VocalTract::NUM_CENTERLINE_POINTS; i++)
	{
		// extract the data of the cross-section
		tract->getCrossProfiles(tract->centerLine[i].point, tract->centerLine[i].normal,
			upperProfile, lowerProfile, upperProfileSurface, lowerProfileSurface, true, articulator);

		// create the corresponding contours
		tmpContours.clear();
		tmpSurfacesIdx.clear();
		createContour(upperProfile, lowerProfile, upperProfileSurface,
			lowerProfileSurface, areas, spacings, tmpContours, tmpSurfacesIdx);

		contours.push_back(tmpContours);
		surfaceIdx.push_back(tmpSurfacesIdx);
		centerLine.push_back(tract->centerLine[i].point);
		normals.push_back(tract->centerLine[i].normal);
	}
}

//*****************************************************************************
// Extract the contours, the surface indexes, the centerline and the normals
// from a csv file 

bool Acoustic3dSimulation::extractContoursFromCsvFile(
	vector<vector<Polygon_2>>& contours, vector<vector<vector<int>>>& surfaceIdx,
	vector<Point2D>& centerLine, vector<Point2D>& normals,
	bool simplifyContours)
{
	//******************************************************************
	// Extract the centerline, centerline normals and the contours
	// from the csv file
	//******************************************************************

	vector<Polygon_2> tmpCont;
	vector<int> tmpIdx;
	vector<vector<int>> tmpVecIdx;
	string line, coordX, coordY;
	char separator(';');
	Cost cost;								// for contour simplification
	int idxCont(0);

	//ofstream log("log.txt", ofstream::app);
	//log << "Start geometry importation" << endl;

	ifstream geoFile(m_geometryFile);

	// initialize maximal bounding box
	m_maxCSBoundingBox.first = Point2D(0., 0.);
	m_maxCSBoundingBox.second = Point2D(0., 0.);

	// check if the file is opened
	if (!geoFile.is_open())
	{
		//log << "Cannot open " << m_geometryFile << endl;
		return false;
	}
	else
	{
		//log << "File successfully opened " << endl;

		while (getline(geoFile, line))
		{
			// initialize centerline point, centerline normal and contour
			tmpCont.clear();
			tmpCont.push_back(Polygon_2());

			// extract the line corresponding to the x components
			stringstream lineX(line);
			// extract the line corresponding to the y components
			getline(geoFile, line);
			stringstream lineY(line);

			// extract the centerline point
			getline(lineX, coordX, separator);
			getline(lineY, coordY, separator);
			centerLine.push_back(Point2D(stod(coordX), stod(coordY)));

			// extract x component of the normal
			getline(lineX, coordX, separator);
			getline(lineY, coordY, separator);
			normals.push_back(Point2D(stod(coordX), stod(coordY)));

			// extract contour
			while (getline(lineX, coordX, separator) && (coordX.length() > 0))
			{
				getline(lineY, coordY, separator);
				tmpCont.back().push_back(Point(stod(coordX), stod(coordY)));
			}

			// remove the last point if it is identical to the first point
			auto itFirst = tmpCont.back().vertices_begin();
			auto itLast = tmpCont.back().vertices_end() - 1;
			if (*itFirst == *itLast)
			{
				tmpCont.back().erase(itLast);
			}

			// compute the maximal bounding box of the created contour 
			// (this is used after to scale the image of the mesh and modes)
			m_maxCSBoundingBox.first.x = min(tmpCont.back().bbox().xmin(),
				m_maxCSBoundingBox.first.x);
			m_maxCSBoundingBox.first.y = min(tmpCont.back().bbox().ymin(),
				m_maxCSBoundingBox.first.y);
			m_maxCSBoundingBox.second.x = max(tmpCont.back().bbox().xmax(),
				m_maxCSBoundingBox.second.x);
			m_maxCSBoundingBox.second.y = max(tmpCont.back().bbox().ymax(),
				m_maxCSBoundingBox.second.y);

			// if requested, simplify the contour removing points which are close
			if (simplifyContours)
			{
				tmpCont.back() = CGAL::Polyline_simplification_2::simplify(
					tmpCont.back(), cost, Stop(0.5));
			}

			// Add the contour to the contour list
			contours.push_back(tmpCont);

			// generate the surface indexes (all zero since there is no clue 
			// about the surface type)
			tmpIdx.clear();
			tmpIdx.assign(tmpCont.back().size(), 0);
			tmpVecIdx.clear();
			tmpVecIdx.push_back(tmpIdx);
			surfaceIdx.push_back(tmpVecIdx);

			//log << "Contour " << idxCont << " extracted" << endl;
			idxCont++;
		}
		return true;
	}
}

//*************************************************************************
// Create cross-sections from VocalTractLab current geometry
// adding intermediate 0 length section where 
// one of the section is not exactely contained in the other

bool Acoustic3dSimulation::createCrossSections(VocalTract* tract,
	bool createRadSection)
{
	//*******************************************
	// Extract the contours and the centerline
	//*******************************************

	// variables for contour extraction
	vector<vector<Polygon_2>> contours;
	vector<vector<vector<int>>> surfaceIdx;
	vector<Point2D> centerLine;
	vector<Point2D> normals;
	vector<double> totAreas;
	vector<array<double, 4>> bboxes;
	array<double, 4> arrayZeros = { 0., 0., 0., 0. };

	ofstream log("log.txt", ofstream::app);

	if (m_geometryImported)
	{
		if ( !extractContoursFromCsvFile(contours, surfaceIdx, centerLine, 
			normals, false))
		{
			return false;
		}
	}
	else
	{
		extractContours(tract, contours, surfaceIdx, centerLine, normals);
	}

	// compute the total area and the bounding box of each contour groups
	for (auto conts : contours)
	{
		totAreas.push_back(0.);
		bboxes.push_back(arrayZeros);
		for (auto cont : conts){
			totAreas.back() += abs(cont.area());
			
			bboxes.back()[0] = min(bboxes.back()[0], cont.bbox().xmin());
			bboxes.back()[1] = max(bboxes.back()[1], cont.bbox().xmax());
			bboxes.back()[2] = min(bboxes.back()[2], cont.bbox().ymin());
			bboxes.back()[3] = max(bboxes.back()[3], cont.bbox().ymax());
		}
	}

	// export total area
	ofstream ar("area.txt");
	for (auto it : totAreas) { ar << it << endl; }
	ar.close();

	//*********************************************************
	// Create lambda expression to compute the scaling factors
	//*********************************************************

	auto getScalingFactor = [&](int idx1, int idx2)
	{
		log << "Start compute scaling factor idx1 " 
			<< idx1 << " idx2 " << idx2 << endl;
		double meanX((abs(bboxes[idx1][0]) + abs(bboxes[idx1][1])
			+ abs(bboxes[idx2][0]) + abs(bboxes[idx2][1])));
		double meanY((abs(bboxes[idx1][2]) + abs(bboxes[idx1][3])
			+ abs(bboxes[idx2][2]) + abs(bboxes[idx2][3])));

		log << "mean X " << meanX << " mean Y " << meanY << endl;
		if (meanX > meanY)
		{
			log << "X" << endl;
			return min(bboxes[idx2][0] / bboxes[idx1][0],
				bboxes[idx2][1] / bboxes[idx1][1]);
		}
		else
		{
			log << "Y" << endl;
			return min(bboxes[idx2][2] / bboxes[idx1][2],
				bboxes[idx2][3] / bboxes[idx1][3]);
		}
	};

	//*******************************************
	// Create the cross-sections
	//*******************************************

	// variables for cross-section creation
	double prevCurvRadius, curvRadius, prevAngle, angle, shift, area, length, radius;
	double prevScalingFactors[2] = { 1., 1. };
	double scalingFactors[2] = { 1., 1. };
	double array1[2] = { 1., 1. };
	const double MINIMAL_AREA(0.15);
	vector<int> tmpPrevSection, prevSecInt, listNextCont, tmpSurf;
	vector<vector<int>> prevSections, intSurfacesIdx;
	int secIdx(0), intSecIdx(0), nextSecIdx;
	vector<Polygon_2> intContours;
	Pwh_list_2 intersections;
	bool sidePrev, side;
	

	double propChg(1.);

	// clear cross-sections before 
	m_crossSections.clear();

	// compute the curvatur parameters of the first cross-section
	getCurvatureAngleShift(centerLine[0], centerLine[1], normals[0], normals[1],
		prevCurvRadius, prevAngle, shift);

	// shift the contours
	Transformation translate(CGAL::TRANSLATION, Vector(0., shift));
	for (auto cont : contours[0]){cont = transform(translate, cont);}

	// initialize the previous section index list
	for (int c(0); c < contours[0].size(); c++){prevSections.push_back(tmpPrevSection);}

	// initialise the scaling factors
	if (m_simuParams.varyingArea)
	{
		prevScalingFactors[0] = 1.;
		switch (m_contInterpMeth)
		{
		case AREA:
			prevScalingFactors[1] = (1. - propChg) +
				propChg * sqrt((totAreas[0] + totAreas[1]) / 2. / totAreas[0]);
			break;
		case BOUNDING_BOX:
			// determine the lagest dimension 
			prevScalingFactors[1] = (1. + getScalingFactor(0, 1)) / 2.;
			break;
		}
		
	}

	// create the cross-sections
	for (int i(1); i < contours.size(); i++)
	{
		log << "\ni= " << i << endl;
		//**********************************
		// Create previous cross-sections
		//**********************************

		// devide the angle and the length of the 2 first sections by 2 since they 
		// share they are defined between the same two points
		if (i < 3) { 
			prevAngle /= 2.; 
			length = centerLine[0].getDistanceFrom(centerLine[1]) / 2.;
		}
		else
		{
			length = centerLine[i - 1].getDistanceFrom(centerLine[i]);
		}

		// compute the scaling factors
		if (m_simuParams.varyingArea)
		{
			if (i == 1)
			{
				switch (m_contInterpMeth)
				{
				case AREA:
					scalingFactors[0] = (1. - propChg) +
						propChg * sqrt((totAreas[0] + totAreas[1]) / 2. / totAreas[1]);
					break;
				case BOUNDING_BOX:
					scalingFactors[0] = (1. + getScalingFactor(1, 0)) / 2.;
					break;
				}
			}
			else
			{
				switch (m_contInterpMeth)
				{
				case AREA:
					scalingFactors[0] = (1. - propChg) + propChg * sqrt(totAreas[i - 1] / totAreas[i]);
					break;
				case BOUNDING_BOX:
					log << "Before compute scaling factor i " << i 
						<< " i-1 " << i-1 << endl;
					scalingFactors[0] = getScalingFactor(i, i-1);
					break;
				}
			}
			scalingFactors[1] = 1.;
		}
		log << "scalingFactors[0] " << scalingFactors[0] << endl;
		log << "scalingFactors[1] " << scalingFactors[1] << endl;

		log << "i = " << i << " scaling factors " << prevScalingFactors[0]
			<< "  " << prevScalingFactors[1] << " area " 
			<< totAreas[i-1]<< " radius " << prevCurvRadius
			<< " angle " << prevAngle << endl;

		// loop over the created contours
		for (int c(0); c < contours[i-1].size(); c++)
		{
			area = abs(contours[i - 1][c].area());
			addCrossSectionFEM(area, sqrt(area)/m_meshDensity, contours[i - 1][c],
				surfaceIdx[i-1][c], length, centerLine[i - 1], normals[i - 1], 
				prevScalingFactors);

			// set the connections with the previous cross-sections if necessary
			if (prevSections[c].size() > 0)
			{
				for (int cn(0); cn < prevSections[c].size(); cn++)
				{
					m_crossSections[secIdx]->
						setPreviousSection(prevSections[c][cn]);
				}
			}

			// set the curvature radius
			m_crossSections[secIdx]->setCurvatureRadius(prevCurvRadius);

			// set the curvature angle
			m_crossSections[secIdx]->setCurvatureAngle(prevAngle);

			log << "Section " << secIdx << " created" << endl;
			secIdx++;
		}

		//******************************************
		// extract current cross-section parameters
		//******************************************

		// compute the curvatur parameters of the section
		getCurvatureAngleShift(centerLine[i - 1], centerLine[i], normals[i - 1], normals[i],
			curvRadius, angle, shift);

		// shift the contours
		Transformation translate(CGAL::TRANSLATION, Vector(0., shift));
		for (auto cont : contours[i]) {cont = transform(translate, cont);}

		log << "Contour shifted" << endl;

		// if the area is equal to the minimal area, the contours
		// are defined as the scaled previous contours
		//
		// Compute the sum of the areas of all the contours of the slice
		// FIXME: normally already computed
		area = 0.;
		for (auto cont : contours[i]) { area += abs(cont.area());
		}
		if (area <= MINIMAL_AREA)
		{
			// copy the previous contours scaling them so that they 
			// have the minimal area

			// compute the sum of the areas of all the contours 
			// of the previous slice
			area = 0.;
			for (auto cont : contours[i-1]) { area += abs(cont.area()); }
			scalingFactors[0] = MINIMAL_AREA / area;
			Transformation scale(CGAL::SCALING, scalingFactors[0]);
			scalingFactors[0] = 1.;
			scalingFactors[1] = 1.;

			// FIXME: check if the copy works as intended
			contours[i] = contours[i-1];
			surfaceIdx[i] = surfaceIdx[i-1];
			//scale the contours
			for (auto cont : contours[i]){cont = transform(scale, cont);
			}
		}

		log << "Minimal area checked" << endl;

		//**********************************
		// Create intermediate 0 length
		// sections if necessary
		//**********************************

		// FIXME A particular case of all points of one polygon beeing
		// Either outside or inside but an edge cutting a part of the other
		// polygon is not taken into account
		//
		// This can be solved by applying the check to the other polygon
		// if no intersection is found
		// 
		// However, this is not very likely to be found and it is more likely 
		// when the polygons have few points

		// clear intermediate contour list
		intContours.clear();
		intSurfacesIdx.clear();
		intersections.clear();
		prevSecInt.clear();
		listNextCont.clear();
		prevSections.clear();
		intSecIdx = 0;

		// loop over the contours of the current cross-section
		for (int c(0); c < contours[i].size(); c++)
		{
			// clean the temporary previous section list
			tmpPrevSection.clear();

			// extract the contour and scale it
			Transformation scale(CGAL::SCALING, scalingFactors[0]);
			Polygon_2 cont(transform(scale, contours[i][c]));

			log << "Contour extracted and scaled sc = " 
				<< scalingFactors[0] << endl;

			// loop over the contours of the previous cross-section
			for (int cp(0); cp < contours[i-1].size(); cp++)
			{
				// extract the previous contour and scale it
				Transformation scale(CGAL::SCALING, prevScalingFactors[1]);
				Polygon_2 prevCont(transform(scale, contours[i - 1][cp]));

				log << "Prev Contour extracted and scaled sc = " 
					<< prevScalingFactors[1] << endl;

				// Check if the current and previous contour intersect:
				// check if the first point of the previous contour
				// is on the bounded side of the current contour
				auto itP = prevCont.begin();
				sidePrev = cont.has_on_bounded_side(*itP);

				log << "Sideprev " << sidePrev << endl;

				// loop over the points of the previous contour
				for (; itP != prevCont.end(); itP++)
				{
					side = cont.has_on_bounded_side(*itP);
					log << "Side " << side << endl;
					// if the previous point and the next point are on
					// different sides
					if (side != sidePrev)
					{
						log << "side != sidePrev" << endl;

						ofstream os("cont.txt");
						for (auto pt : cont) { os << pt.x() << "  " << pt.y() << endl; }
						os.close();
						os.open("pcont.txt");
						for (auto pt : prevCont) { os << pt.x() << "  " << pt.y() << endl; }
						os.close();

						// compute the intersections of both contours
						intersections.clear();
						CGAL::intersection(prevCont, cont, back_inserter(intersections));

						log << intersections.size() << " intersections computed" << endl;

						// loop over the intersection polygons created
						for (auto pol = intersections.begin();
							pol != intersections.end(); pol++)
						{
							// add the corresponding previous section index
							// to the previous section index list for the 
							// intermediate section which will be created
							prevSecInt.push_back(secIdx - contours[i - 1].size() + cp);

							// add the corresponding index of the next contour
							listNextCont.push_back(c);

							// add the index of the intermediate section 
							// which will be created from this contour 
							// to the previous section index list
							tmpPrevSection.push_back(secIdx + intSecIdx);

							// add it the intermediate contour list
							intContours.push_back(pol->outer_boundary());

							// create a surface index vector (the value does not
							// mater since they are not used after)
							tmpSurf.clear();
							tmpSurf.assign(intContours.back().size(), 0);
							intSurfacesIdx.push_back(tmpSurf);

							intSecIdx++;
						}
						break;
					}
					else
					{
						sidePrev = side;
					}
				}
				// if no intersection has been found, check if one contour is 
				// completely contained inside the other
				if ((sidePrev == side) &&
					CGAL::do_intersect(contours[i][c], contours[i - 1][cp]))
				{
					tmpPrevSection.push_back(secIdx - contours[i - 1].size() + cp);
				}
			}
			// add the list of previous section to connect to the 
			// current section
			prevSections.push_back(tmpPrevSection);
		}

		log << "Intersection computed" << endl;

		// set the next section indexes to the previous section 
		nextSecIdx = secIdx + intSecIdx; // index of the first next section
		for (int c(0); c < prevSections.size(); c++)
		{
			if (prevSections[c].size() > 0)
			{
				// loop over the previous sections
				for (int cp(0); cp < prevSections[c].size(); cp++)
				{
					// if the section is not an intermediate one
					if (prevSections[c][cp] < secIdx)
					{
						m_crossSections[prevSections[c][cp]]->setNextSection(
							nextSecIdx + c
						);
					}
				}
			}
		}

		log << "Connections set" << endl;

		// if intermediate contours have been created, add corresponding
		// cross-sections
		if (intContours.size() != 0)
		{
			for (int c(0); c < intContours.size(); c++) {

				area = abs(intContours[c].area());
				addCrossSectionFEM(area, sqrt(area)/m_meshDensity, intContours[c],
					intSurfacesIdx[c], 0., centerLine[i], normals[i], array1);

				// define as a junction section
				m_crossSections[secIdx]->setJunctionSection(true);
				// set the previous section index
				m_crossSections[secIdx]->setPreviousSection(prevSecInt[c]);
				// set this section as the next section of its previous one
				m_crossSections[prevSecInt[c]]->setNextSection(secIdx);
				// set the next section
				m_crossSections[secIdx]->setNextSection(nextSecIdx +
					listNextCont[c]);
				//log << "Intermediate sec " << secIdx << " created" << endl;
				secIdx++;
			}
		}

		// set current cross-section as previous cross-section
		log << "Copy data" << endl;
		std::copy(begin(scalingFactors), end(scalingFactors), begin(prevScalingFactors));
		prevCurvRadius = curvRadius;
		prevAngle = angle;
		log << "Data copied" << endl;
	}

	//********************************
	// create last cross-sections
	//********************************

	log << "Create last cross-section" << endl;

	radius = 0.;	// initalise the radius of the radiation cross-section
	tmpPrevSection.clear(); // the list of section connected to the radiation section
	for (int c(0); c < contours.back().size(); c++)
	{
		// add the index of the created cross-section to the list of cross-sections
		// to connect to the radiation cross-section
		tmpPrevSection.push_back(secIdx);

		area = abs(contours.back()[c].area());
		length = centerLine.rbegin()[1].getDistanceFrom(centerLine.back());
		addCrossSectionFEM(area, sqrt(area)/m_meshDensity, contours.back()[c],
			surfaceIdx.back()[c], length, centerLine.back(), normals.back(), 
			prevScalingFactors);

		// set the connections with the previous cross-sections if necessary
		if (prevSections[c].size() > 0)
		{
			for (int cn(0); cn < prevSections[c].size(); cn++)
			{
				m_crossSections[secIdx]->
					setPreviousSection(prevSections[c][cn]);
			}
		}
		// set the curvature radius
		m_crossSections[secIdx]->setCurvatureRadius(prevCurvRadius);

		// set the curvature angle
		m_crossSections[secIdx]->setCurvatureAngle(prevAngle);

		// set the radius of the radiation cross-section so that the last 
		// cross-sections are contained inside
		radius = max(radius, max({ contours.back()[c].bbox().xmax(),
			contours.back()[c].bbox().ymax(),
			abs(contours.back()[c].bbox().xmin()),
			abs(contours.back()[c].bbox().ymin()) }));

		secIdx++;
	}

	//***************************************************************
	// Create radiation cross-section
	//***************************************************************

	if (createRadSection)
	{
		log << "Create radiation cross-section" << endl;

		double PMLThickness = radius;
		radius *= 2.1;

		addCrossSectionRadiation(centerLine.back(), normals.back(), 
			radius, PMLThickness);

		// Connect with the previous cross-sections
		for (int i(0); i < tmpPrevSection.size(); i++)
		{
			m_crossSections[secIdx]->setPreviousSection(tmpPrevSection[i]);
			m_crossSections[tmpPrevSection[i]]->setNextSection(secIdx);
		}
	}
	return true;
	log.close();
}

//*************************************************************************
// Export the geometry extracted as CSV file

void Acoustic3dSimulation::exportGeoInCsv(string fileName)
{
	ofstream of(fileName);
	stringstream strX, strY;
	string separator(";");
	Point Pt;
	Vector N;
	double theta, thetaN;

	ofstream log("log.txt", ofstream::app);

	for (int i(0); i < m_crossSections.size(); i++)
	{
		if (m_crossSections[i]->length() > 0.)
		{
			//**********************************************************
			// contour in

			if (i == 0)
			{
				Pt = Point(m_crossSections[i]->ctrLinePt().x,
					m_crossSections[i]->ctrLinePt().y);
				N = Vector(m_crossSections[i]->normal().x,
					m_crossSections[i]->normal().y);
			}
			else
			{
				Pt = Point(m_crossSections[i]->ctrLinePt().x,
					m_crossSections[i]->ctrLinePt().y);
				N = Vector(m_crossSections[i]->normal().x,
					m_crossSections[i]->normal().y);
				theta = m_crossSections[i]->circleArcAngle()/2.;
				Transformation rotate(CGAL::ROTATION, sin(M_PI / 2. - theta),
					cos(M_PI / 2. - theta));
				Transformation translate(CGAL::TRANSLATION,
					2. * abs(m_crossSections[i]->curvRadius()) * sin(theta) * rotate(N));
				log << "Radius " << abs(m_crossSections[i]->curvRadius()) << endl;
				log << "theta " << theta << endl;
				log << "Rotated N " << rotate(N) << endl;
				log << "Vector translate " <<
					2. * abs(m_crossSections[i]->curvRadius()) * sin(theta) * rotate(N) << endl;
				log << "Pt before " << Pt << endl;
				Pt = translate(Pt);
				if (signbit(m_crossSections[i]->curvRadius()))
				{
					thetaN = 2. * theta;
				}
				else
				{
					thetaN = -2. * theta;
				}
				Transformation rotateN(CGAL::ROTATION, sin(thetaN), cos(thetaN));
				log << "Pt after " << Pt << endl;
				N = rotateN(N);
			}

			// write centerline point coordinates
			strX << Pt.x() << separator;
			strY << Pt.y() << separator;

			// write normal coordinates
			strX << N.x() << separator;
			strY << N.y() << separator;

			// write contour
			for (auto pt : m_crossSections[i]->contour())
			{
				strX << m_crossSections[i]->scaleIn() * pt.x() << separator;
				strY << m_crossSections[i]->scaleIn() * pt.y() << separator;
			}

			of << strX.str() << endl << strY.str() << endl;

			strX.str("");
			strX.clear();
			strY.str("");
			strY.clear();

			//********************************************************
			// contour out

			if (i == 0)
			{
				Pt = Point(m_crossSections[i]->ctrLinePt().x,
					m_crossSections[i]->ctrLinePt().y);
				N = Vector(m_crossSections[i]->normal().x,
					m_crossSections[i]->normal().y);
				theta = m_crossSections[i]->circleArcAngle() / 2.;
				Transformation rotate(CGAL::ROTATION, sin(-M_PI / 2. + theta),
					cos(-M_PI / 2. + theta));
				Transformation translate(CGAL::TRANSLATION,
					2. * abs(m_crossSections[i]->curvRadius()) * sin(theta) * rotate(N));
				log << "Radius " << abs(m_crossSections[i]->curvRadius()) << endl;
				log << "theta " << theta << endl;
				log << "Rotated N " << rotate(N) << endl;
				log << "Vector translate " <<
					2. * abs(m_crossSections[i]->curvRadius()) * sin(theta) * rotate(N) << endl;
				log << "Pt before " << Pt << endl;
				Pt = translate(Pt);
				log << "Pt after " << Pt << endl;
				Transformation rotateN(CGAL::ROTATION, sin(2. * theta), cos(2. * theta));
				N = rotateN(N);
			}
			else
			{
				Pt = Point(m_crossSections[i]->ctrLinePt().x,
					m_crossSections[i]->ctrLinePt().y);
				N = Vector(m_crossSections[i]->normal().x,
					m_crossSections[i]->normal().y);
			}

			// write centerline point coordinates
			strX << Pt.x() << separator;
			strY << Pt.y() << separator;

			// write normal coordinates
			strX << N.x() << separator;
			strY << N.y() << separator;

			// write contour
			for (auto pt : m_crossSections[i]->contour())
			{
				strX << m_crossSections[i]->scaleOut() * pt.x() << separator;
				strY << m_crossSections[i]->scaleOut() * pt.y() << separator;
			}

			of << strX.str() << endl << strY.str() << endl;

			strX.str("");
			strX.clear();
			strY.str("");
			strY.clear();
		}
	}
	of.close();
	log.close();
}

//*************************************************************************
// Interpolate the radiation and admittance matrices with splines

void Acoustic3dSimulation::preComputeRadiationMatrices(int nbRadFreqs) {

	int numSec(m_crossSections.size());
	int mn(m_crossSections[numSec - 2]->numberOfModes());
	double freq;
	double radFreqSteps((double)SAMPLING_RATE / 2. / (double)(nbRadFreqs - 1));
	vector<double> stepRadFreqs;
	m_radiationFreqs.clear();
	m_radiationFreqs.reserve(nbRadFreqs);
	stepRadFreqs.reserve(nbRadFreqs - 1);
	Matrix A(Matrix::Zero(nbRadFreqs - 2, nbRadFreqs - 2)), imped(mn, mn);
	Eigen::VectorXd R(Eigen::VectorXd::Zero(nbRadFreqs - 2)), B;
	Eigen::MatrixXcd radImped(mn,mn), radAdmit(mn,mn);
	vector<double>* a, * b, * c, * d;

	//ofstream prop;

	// initialize coefficient structures
	//        a       |        b       |        c       |       d       
	// Zr  Zi  Ir  Ii | Zr  Zi  Ir  Ii | Zr  Zi  Ir  Ii | Zr  Zi  Ir  Ii  	
	// 0   1   2   3  | 4   5   6   7  | 8   9   10  11 | 12  13  14  15
	m_radiationMatrixInterp.clear();
	for (int m(0); m < 16; m++)
	{
		m_radiationMatrixInterp.push_back(vector<vector<vector<double>>>());
		for (int i(0); i < mn; i++) {

			m_radiationMatrixInterp.back().push_back(vector<vector<double>>());

			for (int j(0); j < mn; j++)
			{
				m_radiationMatrixInterp.back().back().push_back(vector<double>());
				m_radiationMatrixInterp.back().back().back().reserve(nbRadFreqs);
			}
		}
	}

	//prop.open("rad.txt");

	// loop over a few frequencies
	for (int i(0); i < nbRadFreqs; i++)
	{
		freq = max(500., (double)i * radFreqSteps);
		m_radiationFreqs.push_back(freq);

		//// Compute radiation impedance matrix at this frequency
		//m_crossSections[numSec - 1]->characteristicImpedance(radImped,
		//	freq, m_soundSpeed, m_volumicMass, 0);
		//m_crossSections[numSec - 1]->characteristicAdmittance(radAdmit,
		//	freq, m_soundSpeed, m_volumicMass, 0);
		//propagateImpedAdmit(radImped, radAdmit, freq, numSec - 1, numSec - 2, m_method);

		//radImped = m_crossSections[numSec - 2]->endImpedance();
		//radAdmit = m_crossSections[numSec - 2]->endAdmittance();

		radiationImpedance(radImped, freq, 15.);
		radAdmit = radImped.inverse();

		//prop << radAdmit.imag() << endl;

		// Create first spline coefficient
		for (int m(0); m < mn; m++)
		{
			for (int n(0); n < mn; n++)
			{
				m_radiationMatrixInterp[0][m][n].push_back(radImped(m, n).real());
				m_radiationMatrixInterp[1][m][n].push_back(radImped(m, n).imag());
				m_radiationMatrixInterp[2][m][n].push_back(radAdmit(m, n).real());
				m_radiationMatrixInterp[3][m][n].push_back(radAdmit(m, n).imag());
			}
		}
	}
	//prop.close();
	//prop.open("fRad.txt");
	// compute frequency steps
	for (int i(0); i < nbRadFreqs - 1; i++)
	{
		//prop << m_radiationFreqs[i] << endl;
		stepRadFreqs.push_back(m_radiationFreqs[i + 1] - m_radiationFreqs[i]);
	}
	//prop << m_radiationFreqs.back() << endl;
	//prop.close();

	// compute the spline coefficients
	for (int m(0); m < 4; m++)
	{
		for (int i(0); i < mn; i++)
		{
			for (int j(0); j < mn; j++)
			{
				A.setZero(nbRadFreqs - 2, nbRadFreqs - 2);
				R.setZero(nbRadFreqs - 2);

				a = &m_radiationMatrixInterp[m][i][j];
				b = &m_radiationMatrixInterp[m + 4][i][j];
				c = &m_radiationMatrixInterp[m + 8][i][j];
				d = &m_radiationMatrixInterp[m + 12][i][j];

				// build matrice A and vector R to solve the equation A * c = R
				// in order to find the coefficient c of the spline
				A(0, 0) = 2 * (stepRadFreqs[0] + stepRadFreqs[1]);
				A(0, 1) = stepRadFreqs[1];
				R(0) = 3. * ((*a)[2] - (*a)[1]) / stepRadFreqs[1]
					- 3 * ((*a)[1] - (*a)[0]) / stepRadFreqs[0];

				for (int f(1); f < nbRadFreqs - 3; f++)
				{
					A(f, f - 1) = stepRadFreqs[f];
					A(f, f) = 2. * (stepRadFreqs[f] + stepRadFreqs[f + 1]);
					A(f, f + 1) = stepRadFreqs[f + 1];
					R(f) = 3. * ((*a)[f + 2] - (*a)[f + 1]) / stepRadFreqs[f + 1]
						- 3. * ((*a)[f + 1] - (*a)[f]) / stepRadFreqs[f];
				}

				A(nbRadFreqs - 3, nbRadFreqs - 4) = stepRadFreqs[nbRadFreqs - 3];
				A(nbRadFreqs - 3, nbRadFreqs - 3) = 2. * (stepRadFreqs[nbRadFreqs - 3] + stepRadFreqs[nbRadFreqs - 2]);
				R(nbRadFreqs - 3) = 3. * ((*a)[nbRadFreqs - 1] - (*a)[nbRadFreqs - 2])
					/ stepRadFreqs[nbRadFreqs - 2]
					- 3. * ((*a)[nbRadFreqs - 2] - (*a)[nbRadFreqs - 3])
					/ stepRadFreqs[nbRadFreqs - 3];

				// compute c coefficient
				B = A.householderQr().solve(R);
				(*c).push_back(0.);
				for (int f(0); f < B.rows(); f++)
				{
					(*c).push_back(B(f));
				}
				(*c).push_back(0.);

				// compute b coefficient
				for (int f(0); f < nbRadFreqs - 1; f++)
				{
					(*b).push_back(((*a)[f + 1] - (*a)[f])
						/ stepRadFreqs[f] - stepRadFreqs[f] *
						((*c)[f + 1] + 2. * (*c)[f]) / 3.);
					(*d).push_back(((*c)[f + 1] - (*c)[f]) / 3. / stepRadFreqs[f]);
				}
			}
		}
	}
}

//*************************************************************************
// Interpolate the radiation impedance matrix at a given frequency

void Acoustic3dSimulation::interpolateRadiationImpedance(Eigen::MatrixXcd& imped, double freq)
{
	// find the index corresponding to the coefficient to use for this frequency
	int nbRadFreqs(m_radiationMatrixInterp[0][0][0].size());
	int idx(nbRadFreqs - 2);
	while (m_radiationFreqs[idx] > freq) { idx--; }
	idx = max(0, idx);

	// get number of modes of the radiating section
	int nbCS(m_crossSections.size());
	int mn(m_crossSections[nbCS - 2]->numberOfModes());

	imped.setZero(mn, mn);
	for (int m(0); m < mn; m++)
	{
		for (int n(0); n < mn; n++)
		{
			imped(m, n) = complex<double>(m_radiationMatrixInterp[0][m][n][idx] +
				m_radiationMatrixInterp[4][m][n][idx] * (freq - m_radiationFreqs[idx]) +
				m_radiationMatrixInterp[8][m][n][idx] * pow(freq - m_radiationFreqs[idx], 2) +
				m_radiationMatrixInterp[12][m][n][idx] * pow(freq - m_radiationFreqs[idx], 3), 
				(m_radiationMatrixInterp[1][m][n][idx] +
				m_radiationMatrixInterp[5][m][n][idx] * (freq - m_radiationFreqs[idx]) +
				m_radiationMatrixInterp[9][m][n][idx] * pow(freq - m_radiationFreqs[idx], 2) +
				m_radiationMatrixInterp[13][m][n][idx] * pow(freq - m_radiationFreqs[idx], 3)));
		}
	}
}

//*************************************************************************
// Interpolate the radiation admittance matrix at a given frequency

void Acoustic3dSimulation::interpolateRadiationAdmittance(Eigen::MatrixXcd& admit, double freq)
{
	// find the index corresponding to the coefficient to use for this frequency
	int nbRadFreqs(m_radiationMatrixInterp[0][0][0].size());
	int idx(nbRadFreqs - 2);
	while (m_radiationFreqs[idx] > freq) { idx--; }
	idx = max(0, idx);

	// get number of modes of the radiating section
	int nbCS(m_crossSections.size());
	int mn(m_crossSections[nbCS - 2]->numberOfModes());

	admit.setZero(mn, mn);
	for (int m(0); m < mn; m++)
	{
		for (int n(0); n < mn; n++)
		{
			admit(m, n) = complex<double>(m_radiationMatrixInterp[2][m][n][idx] +
				m_radiationMatrixInterp[6][m][n][idx] * (freq - m_radiationFreqs[idx]) +
				m_radiationMatrixInterp[10][m][n][idx] * pow(freq - m_radiationFreqs[idx], 2) +
				m_radiationMatrixInterp[14][m][n][idx] * pow(freq - m_radiationFreqs[idx], 3),
				(m_radiationMatrixInterp[3][m][n][idx] +
				m_radiationMatrixInterp[7][m][n][idx] * (freq - m_radiationFreqs[idx]) +
				m_radiationMatrixInterp[11][m][n][idx] * pow(freq - m_radiationFreqs[idx], 2) +
				m_radiationMatrixInterp[15][m][n][idx] * pow(freq - m_radiationFreqs[idx], 3)));
		}
	}
}

//*************************************************************************
// Compute the radiation impedance according to Blandin et al 2019
// Multimodal radiation impedance of a waveguide with arbitrary
// cross - sectional shape terminated in an infinite baffle

void Acoustic3dSimulation::radiationImpedance(Eigen::MatrixXcd& imped, double freq, double gridDensity)
{
	int numSec(m_crossSections.size());
	int mn(m_crossSections[numSec - 2]->numberOfModes());
	complex<double> J(0., 1.);
	double area2(pow(m_crossSections[numSec - 2]->area(),2));

	imped = Eigen::MatrixXcd::Zero(mn, mn);
	Eigen::MatrixXcd integral2(mn, mn);

	//ofstream log("log.txt", ofstream::app);
	//log << "\nStart computing radiation impedance" << endl;

	//******************************
	// generate cartesian grid
	//******************************

	// very good precision is obtained with gridDensity = 30
	//  good precision is obtained with gridDensity = 15;
	double spacing(sqrt(m_crossSections[numSec - 2]->area()) / gridDensity);

	//log << "Spacing: " << spacing << endl;

	Polygon_2 contour(m_crossSections[numSec - 2]->contour());
	vector<Point> cartGrid;
	Point pt;
	double xmin(contour.bbox().xmin());
	double ymin(contour.bbox().ymin());
	int nx(ceil((contour.bbox().xmax() - xmin) / spacing));
	int ny(ceil((contour.bbox().ymax() - ymin) / spacing));
	for (int i(0); i < nx; i++)
	{
		for (int j(0); j < ny; j++)
		{
			pt = Point(xmin + i * spacing, ymin + j * spacing);
			if (contour.has_on_bounded_side(pt))
			{
				cartGrid.push_back(pt);
			}
		}
	}

	//log << "Cartesian grid generated: " << cartGrid.size() << " points" << endl;

	// Interpolate the propagation modes on the cartesian grid
	Matrix intCartGrid(m_crossSections[numSec - 2]->interpolateModes(cartGrid));

	//// export grid
	//ofstream out;
	//out.open("grid.txt");
	//for (auto it : cartGrid)
	//{
	//	out << it.x() << "  " << it.y() << endl;
	//}
	//out.close();
	//// export contour
	//out.open("cont.txt");
	//for (auto it : contour)
	//{
	//	out << it.x() << "  " << it.y() << endl;
	//}
	//out.close();

	// loop over the points of the cartesian grid
	for (int c(0); c < cartGrid.size(); c++)
	{

		//******************************
		// generate polar grid
		//******************************

		//FIXME: don't work for lines intersecting several times the polygon

		double dist;
		int numDirections;
		double angleSpacing, direction;
		vector<Point> polGrid;
		vector<double> radius;
		int cnt, nbPts(0);
		double r, sumH;
		Point ptToAdd;

		// get center point from the cartesian grid
		pt = cartGrid[c];

		// estimate the ratio [number of direction] / [number of point]
		numDirections = 50;
		angleSpacing = 2. * M_PI / (double)numDirections;
		for (int i(0); i < numDirections; i++)
		{
			direction = (double)i * angleSpacing - M_PI;
			cnt = 0;
			r = (0.5 + (double)cnt) * spacing;
			ptToAdd = Point(r * cos(direction) + pt.x(), r * sin(direction) + pt.y());
			while (contour.has_on_bounded_side(ptToAdd))
			{
				nbPts++;
				cnt++;
				r = (0.5 + (double)cnt) * spacing;
				ptToAdd = Point(r * cos(direction) + pt.x(), r * sin(direction) + pt.y());
			}
		}

		//log << "Nb points estimate " << nbPts << endl;


		// Rough estimate of the number of needed directions
		numDirections = cartGrid.size() * numDirections / nbPts;
		//numDirections = cartGrid.size() * contour.size() / nbPts;
		//numDirections = cartGrid.size() / (int)(avDist / spacing);

		// generate angles of the polar grid
		angleSpacing = 2. * M_PI / (double)numDirections;
		for (int i(0); i < numDirections; i++)
		{
			direction = (double)i * angleSpacing - M_PI;

			// generate the points of the polar grid for each direction
			cnt = 0;
			r = (0.5 + (double)cnt) * spacing;
			ptToAdd = Point(r * cos(direction) + pt.x(), r * sin(direction) + pt.y());
			while (contour.has_on_bounded_side(ptToAdd))
			{
				polGrid.push_back(ptToAdd);
				radius.push_back(r);
				cnt++;
				r = (0.5 + (double)cnt) * spacing;
				ptToAdd = Point(r * cos(direction) + pt.x(), r * sin(direction) + pt.y());
			} 
		}

		//log << "Polar grid generated"
		//	<< "\nspacing:\t" << spacing
		//	<< "\nNum directions:\t\t" << numDirections
		//	<< "\nNum points:\t\t" << polGrid.size() << endl;

		// interpolate the polar grid
		Matrix intPolGrid(m_crossSections[numSec - 2]->interpolateModes(polGrid));

		//// export polar grid
		//if (c == 0) {
		//	out.open("pGrid.txt");
		//	for (auto it : polGrid)
		//	{
		//		out << it.x() << "  " << it.y() << endl;
		//	}
		//	out.close();
		//}

		//******************************
		// Compute first integral
		//******************************

		sumH = 0;
		integral2.setZero();

		// loop over the points of the polar grid
		for (int p(0); p < polGrid.size(); p++)
		{
			sumH += radius[p];

			// loop over the modes
			for (int m(0); m < mn; m++)
			{
				// loop over the modes
				for (int n(0); n < mn; n++)
				{
					integral2(m, n) += intPolGrid(p, m) * intCartGrid(c, n) *
						exp(-J * 2. * M_PI * freq * radius[p] / m_simuParams.sndSpeed);
				}
			}
		}

		imped += - integral2 / sumH / 2. / M_PI / cartGrid.size();
	}

	imped *= pow(m_crossSections[numSec - 2]->area(),2);

	//log.close();
}