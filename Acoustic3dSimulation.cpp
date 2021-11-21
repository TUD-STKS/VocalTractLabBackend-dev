#include "Acoustic3dSimulation.h"
#include "Constants.h"
#include "Dsp.h"
#include "TlModel.h"
#include "TdsModel.h"
#include <algorithm>
#include <chrono>    // to get the computation time
#include <ctime>  
#include <string>
#include <regex>

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
typedef CDT::Point                    Point;
typedef CGAL::Point_3<K>                Point_3;
typedef CGAL::Vector_2<K>                Vector;
typedef CGAL::Polygon_2<K>                            Polygon_2;
typedef CGAL::Polygon_with_holes_2<K>                 Polygon_with_holes_2;
typedef std::list<Polygon_with_holes_2>               Pwh_list_2;
typedef CGAL::Aff_transformation_2<K>        Transformation;
typedef CGAL::Aff_transformation_3<K>         Transformation3;
typedef CGAL::Delaunay_mesher_no_edge_refinement_2<CDT, Criteria> MesherNoRefine;
typedef CGAL::Polyline_simplification_2::Stop_below_count_ratio_threshold Stop;
typedef CGAL::Polyline_simplification_2::Squared_distance_cost Cost;

//const double MINIMAL_DISTANCE = 1e-14;

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
// Check if 2 contours are similar with a distance criterion

bool similarContours(Polygon_2& cont1, Polygon_2& cont2, double minDist)
{
  bool similar(false);
  if (cont1.size() == cont2.size())
  {
    int nPt(cont1.size());
    similar = true;
    for (int i(0); i < nPt; i++)
    {
      if ((abs(cont1[i].x() - cont2[i].x()) > minDist)
        || (abs(cont1[i].y() - cont2[i].y()) > minDist))
      { 
        similar = false;
        break;
      }
    }
  }
  return similar;
}

// ****************************************************************************
// Constructor.
// ****************************************************************************

Acoustic3dSimulation::Acoustic3dSimulation()
// initialise the physical constants
  : m_geometryImported(false),
  m_meshDensity(15.),
  m_maxCutOnFreq(20000.),
  m_spectrumLgthExponent(10),
  m_idxSecNoiseSource(46), // for /sh/ 212, for vowels 46
  m_idxConstriction(40),
  m_glottisBoundaryCond(IFINITE_WAVGUIDE),
  m_mouthBoundaryCond(ADMITTANCE_1),
  m_contInterpMeth(AREA)
{
  //m_simuParams.temperature = 31.4266; // for 350 m/s
  m_simuParams.temperature = 21.0735; // for 344 m/s 
  m_simuParams.volumicMass = STATIC_PRESSURE_CGS * MOLECULAR_MASS / (GAS_CONSTANT *
    (m_simuParams.temperature + KELVIN_SHIFT));
  m_simuParams.numIntegrationStep = 25;
  m_simuParams.orderMagnusScheme = 2;
  m_simuParams.propMethod = MAGNUS;
  m_simuParams.freqDepLosses = false;
  m_simuParams.wallLosses = false;
  m_simuParams.sndSpeed = (sqrt(ADIABATIC_CONSTANT * STATIC_PRESSURE_CGS / m_simuParams.volumicMass));
  m_crossSections.reserve(2 * VocalTract::NUM_CENTERLINE_POINTS);
  m_simuParams.percentageLosses = 1.;
  m_simuParams.curved = false;
  m_simuParams.varyingArea = true;
  m_simuParams.maxComputedFreq = 10000.; // (double)SAMPLING_RATE / 2.;
  m_simuParams.freqField = 5000.;
  m_simuParams.bboxField[0] = Point(-5., -10.);
  m_simuParams.bboxField[1] = Point(10., 5.);
  m_simuParams.fieldResolution = 30;
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
// Clean the log file and print the simulation parameters

void Acoustic3dSimulation::generateLogFileHeader(bool cleanLog) {

  double freq, freqSteps((double)SAMPLING_RATE / 2. / (double)m_numFreq);
  int numFreqComputed((int)ceil(m_simuParams.maxComputedFreq / freqSteps));

  ofstream log;
  if (cleanLog) {
    log.open("log.txt", ofstream::out | ofstream::trunc);
    log.close();
  }

  log.open("log.txt", ofstream::app);

  // print the date of the simulation
  time_t start_time = std::chrono::system_clock::to_time_t(chrono::system_clock::now());
  log << ctime(&start_time) << endl;
  log << "Start simulation\n" << endl;
  log << "Air volumic mass: " << m_simuParams.volumicMass << " g/cm^3" << endl;
  log << "Sound speed: " << m_simuParams.sndSpeed << " cm/s" << endl;
  log << "viscous boundary specific admittance " << m_simuParams.viscousBndSpecAdm
    << " g.cm^-2 .s^-1" << endl;
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
  log << "mouth boundary condition: ";
  switch (m_mouthBoundaryCond)
  {
  case RADIATION:
    log << "RADIATION" << endl;
    break;
  case IFINITE_WAVGUIDE:
    log << "IFINITE_WAVGUIDE" << endl;
    break;
  case HARD_WALL:
    log << "HARD_WALL" << endl;
    break;
  case ADMITTANCE_1:
    log << "ADMITTANCE_1" << endl;
    break;
  }
  log << "Mesh density: " << m_meshDensity << endl;
  log << "Max cut-on frequency: " << m_maxCutOnFreq << " Hz" << endl;
  log << "Number of integration steps: " << m_simuParams.numIntegrationStep << endl;
  log << "Index of noise source section: " << m_idxSecNoiseSource << endl;
  log << "Index of constriction section: " << m_idxConstriction << endl;
  log << "Maximal computed frequency: " << m_simuParams.maxComputedFreq
    << " Hz" << endl;
  log << "Number of simulated frequencies: " << numFreqComputed << endl;
  log << "Varying cross-sectional area: ";
  if (m_simuParams.varyingArea)
  {
    log << "yes\nscaling factor computation method: ";
    switch (m_contInterpMeth)
    {
      case AREA:
        log << "AREA" << endl;
        break;
      case BOUNDING_BOX:
        log << "BOUNDING_BOX" << endl;
        break;
      case FROM_FILE:
        log << "FROM_FILE" << endl;
        break;
    }
  }
  else
  {
    log << "no" << endl;
  }
  log << "Take into account curvature: ";
  if (m_simuParams.curved) { log << "yes" << endl; }
  else { log << "no" << endl; }
  log << "Propagation mmethod: ";
  switch (m_simuParams.propMethod)
  {
  case MAGNUS:
    log << "MAGNUS order " << m_simuParams.orderMagnusScheme << endl;
    break;
  case STRAIGHT_TUBES:
    log << "STRAIGHT_TUBES" << endl;
    break;
  }

  if (m_geometryImported)
  {
    log << "Geometry imported from csv file:\n  " << m_geometryFile << endl;
  }
  else
  {
    log << "Geometry is from vocal tract lab" << endl;
  }
  log << "Acoustic field computation at " << m_simuParams.freqField
  << " Hz with " << m_simuParams.fieldResolution << " points per cm" << endl;
  log << "Bounding box: min x " << m_simuParams.bboxField[0].x() 
    << " max x " << m_simuParams.bboxField[1].x() 
    << " min y " << m_simuParams.bboxField[0].y() 
    << " max y " << m_simuParams.bboxField[1].y() << endl;
  log.close();
}

// ****************************************************************************

void Acoustic3dSimulation::setContourInterpolationMethod(enum contourInterpolationMethod method)
{
  m_contInterpMeth = method;
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

  //ofstream mesh;
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
    //  mesh.open("mesh.txt");
    //  cdt = m_crossSections[i]->triangulation();
    //  for (auto it = cdt.faces_begin(); it != cdt.faces_end(); it++)
    //  {
    //    for (int v(0); v < 3; v++)
    //    {
    //      mesh << it->vertex(v)->point().x() << "  "
    //        << it->vertex(v)->point().y() << endl;
    //    }
    //    mesh << it->vertex(0)->point().x() << "  "
    //      << it->vertex(0)->point().y() << endl;
    //    mesh << "nan  nan" << endl;
    //  }
    //  mesh.close();
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
  int idxMinArea;
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

  //ofstream log;
  //log.open("log.txt", ofstream::app);
  //log << "Start computing junction matrices" << endl;

  // loop over the cross-section
  for (int i(0); i < m_crossSections.size(); i++)
  {
    //log << "Num next sec " << m_crossSections[i]->numNextSec() << endl;

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
      //log << "Contour sec " << i << " scaling out " << scaling[0] 
      //  << " area " << contour.area() << endl;
      //if (m_crossSections[i]->isJunction()) { log << "Junction section" << endl; }

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
        //  << " next sec " << nextSec << " nModesNext " << nModesNext << endl;

        Matrix F(Matrix::Zero(nModes, nModesNext));
        areaInt.push_back(0.);

        // get the next contour
        nextContour.clear();
        scaling[1] = m_crossSections[nextSec]->scaleIn();
        Transformation scale(CGAL::SCALING, scaling[1]);
        nextContour = transform(scale, m_crossSections[nextSec]->contour());
        if (contour.area() >= nextContour.area())
        {
          idxMinArea = 1;
        }
        else
        {
          idxMinArea = 0;
        }
        //nextContour = m_crossSections[nextSec]->contour();
        //log << "Next contour, sec " << nextSec << " scaling in " << scaling[1] 
        //  << " area " << nextContour.area() << endl;
        //if (m_crossSections[nextSec]->isJunction()) { log << "Junction section" << endl; }

        //////////////////////////////////////////////////////////////
        // Compute the intersections of the contours
        //////////////////////////////////////////////////////////////

        // compute the intersection between the contours of the current

        //log << "Next sec " << typeid(*m_crossSections[nextSec]).name()
        //  << " is radiation "
        //  << (typeid(*m_crossSections[nextSec]) == typeid(CrossSection2dRadiation))
        //  << endl;

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
          if (!similarContours(contour, nextContour, MINIMAL_DISTANCE_DIFF_POLYGONS))
          {
            CGAL::intersection(contour, nextContour, std::back_inserter(intersections));
            //log << "Intersection computed" << endl;
          }
          else
          {
            intersections.push_back(Polygon_with_holes_2(nextContour));
            //log << "Next contour copied" << endl;
          }
        }

        if (computeG) {
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

        //log << "Nb intersections " << intersections.size() << endl;

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

            //log << "interpolation of first section done " 
            //  << interpolation1.rows() << "  " 
            //  << interpolation1.cols() << endl;

            // interpolate the modes of the next cross-section
            interpolation2 = m_crossSections[nextSec]->interpolateModes(pts
              ,1. / scaling[1]
            );

            //log << "interpolation of second section done "
            //  << interpolation2.rows() << "  "
            //  << interpolation2.cols() << endl;

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
                         / pow(scaling[idxMinArea], 2)
                             /// scaling[0] / scaling[1]
                        ;
                    }
                  }
                }
              }
            }
          }
        }

        //log << "\nF\n" << F << endl << endl;

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
  //log.close();

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
    //  log << "{";
    //  for (auto seg : it)
    //  {
    //    log << seg << " ";
    //  }
    //  log << "} ";
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
          //  log << seg;
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

void Acoustic3dSimulation::propagateImpedAdmit(Eigen::MatrixXcd& startImped,
  Eigen::MatrixXcd& startAdmit, double freq, int startSection, int endSection)
{
  int direction;
  if (startSection > endSection)
  {
    direction = -1;
  }
  else
  {
    direction = 1;
  }
  propagateImpedAdmit(startImped, startAdmit, freq, startSection, endSection,
    direction);
}

void Acoustic3dSimulation::propagateImpedAdmit(Eigen::MatrixXcd & startImped,
  Eigen::MatrixXcd & startAdmit, double freq, int startSection, int endSection,
  int direction)
{
  Eigen::MatrixXcd prevImped;
  Eigen::MatrixXcd prevAdmit;
  vector<Matrix> F;
  //vector<double> areaInt;
  Matrix G;
  int numSec(m_crossSections.size()), nI, nPs;
  int prevSec;
  double areaRatio;
  complex<double> wallInterfaceAdmit(1i*2.*M_PI*freq* 
    m_simuParams.thermalBndSpecAdm/m_simuParams.sndSpeed);
  vector<Eigen::MatrixXcd> inputImped;

  std::chrono::duration<double> time;

 /* ofstream log;
  log.open("log.txt", ofstream::app);
  log << "start propagating admittance" << endl;*/
  //log << "Direction " << direction << endl;

  // set the initial impedance and admittance matrices
  m_crossSections[startSection]->clearImpedance();
  m_crossSections[startSection]->clearAdmittance();

  // set the propagation direction of the first section
  m_crossSections[startSection]->setZdir(direction);
  m_crossSections[startSection]->setYdir(direction);
  
  switch(m_simuParams.propMethod)
  {
  case MAGNUS:
    //m_crossSections[startSection]->propagateImpedAdmiteRiccati(startImped, startAdmit, 
    //  m_simuParams, freq, (double)direction, time);
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
      m_crossSections[max(0, min(numSec, startSection+direction))]->area());
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
    m_crossSections[i]->setZdir(direction);
    m_crossSections[i]->setYdir(direction);

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

    //log << "Size F " << F.size() << endl;
    //log << "Size F " << F[0].rows() << "  " << F[0].cols() << endl;

    prevImped = Eigen::MatrixXcd::Zero(m_crossSections[i]->numberOfModes(),
      m_crossSections[i]->numberOfModes());
    prevAdmit = Eigen::MatrixXcd::Zero(m_crossSections[i]->numberOfModes(),
      m_crossSections[i]->numberOfModes());
    
    switch (m_simuParams.propMethod)
    {
    case MAGNUS:
      if (direction == -1)
      {
      // case of a contraction: area(i) > area(ps)
        if ((m_crossSections[i]->area() * 
                    pow(m_crossSections[i]->scaleOut(), 2)) >
           (m_crossSections[prevSec]->area() * 
            pow(m_crossSections[prevSec]->scaleIn(), 2)))
        {
          prevAdmit += 
           //   (pow(m_crossSections[i]->scaleOut(), 2)/
           //pow(m_crossSections[prevSec]->scaleIn(), 2))*
            F[0] * m_crossSections[prevSec]->Yin()
            * (F[0].transpose())
            - wallInterfaceAdmit*
            //////m_crossSections[i]->curvature() *
            G
            ;
        }
      // case of an expansion: area(i) < area(ps)
        else
        {
          prevImped += 
            //(pow(m_crossSections[prevSec]->scaleIn(),2) /
            //pow(m_crossSections[i]->scaleOut(), 2))* 
            F[0] * m_crossSections[prevSec]->Zin()
            *(Matrix::Identity(nPs, nPs) -
              wallInterfaceAdmit * 
            //  ////m_crossSections[prevSec]->curvature() * 
              G*m_crossSections[prevSec]->Zin()).inverse()
            * (F[0].transpose())
            ;
        prevAdmit += prevImped.fullPivLu().inverse();
        }
      }
      else
      {
      // case of a contraction: area(i) > area(ps)
        if ((m_crossSections[i]->area() * 
                    pow(m_crossSections[i]->scaleIn(), 2)) >
           (m_crossSections[prevSec]->area() * 
            pow(m_crossSections[prevSec]->scaleOut(), 2)))
        {
          prevAdmit += 
              (pow(m_crossSections[i]->scaleIn(), 2)/
            pow(m_crossSections[prevSec]->scaleOut(), 2)) * 
            (F[0].transpose()) *
            m_crossSections[prevSec]->Yout() * F[0]
            //+ wallInterfaceAdmit * 
            //////m_crossSections[i]->curvature() * 
            //G
            ;
        }
      // case of an expansion: area(i) < area(ps)
        else
        {
          prevImped += 
              (pow(m_crossSections[prevSec]->scaleOut(),2) /
            pow(m_crossSections[i]->scaleIn(),2)) * 
            (F[0].transpose()) * m_crossSections[prevSec]->Zout()
            //*(Matrix::Identity(nPs, nPs) + 
            //  wallInterfaceAdmit *
            //  ////m_crossSections[prevSec]->curvature() * 
            //  G*m_crossSections[prevSec]->Zout()).inverse()
            * F[0]
            ;
        prevAdmit += prevImped.fullPivLu().inverse();
        }
      }

      break;
    case STRAIGHT_TUBES:

      areaRatio = max(m_crossSections[prevSec]->area(),
        m_crossSections[i]->area()) /
        min(m_crossSections[prevSec]->area(),
          m_crossSections[i]->area());
      if (direction == -1)
      {
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
      }
      else
      {
        // case of a contraction
        if (m_crossSections[i]->area() > m_crossSections[prevSec]->area())
        {
          prevAdmit += areaRatio * F[0] *
            m_crossSections[prevSec]->Yout()
            * (F[0].transpose());
          prevImped += prevAdmit.fullPivLu().inverse();
        }
        // case of an expansion
        else
        {
          prevImped += areaRatio * F[0] *
            m_crossSections[prevSec]->Zout()
            * (F[0].transpose());
          prevAdmit += prevImped.fullPivLu().inverse();
        }
      }
      break;
    }

    //log << "prevImped computed" << endl;

    // propagate admittance in the section
    switch (m_simuParams.propMethod) {
    case MAGNUS:
      //m_crossSections[i]->propagateImpedAdmiteRiccati(prevImped, prevAdmit, m_simuParams,
      //  freq, (double)direction, time);
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
        m_crossSections[max(0,min(numSec, i + direction))]->area());
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
      //  << invF.transpose().cols() << endl;
      //log << "Size Y " << m_crossSections[i + 1]->startAdmittance().rows()
      //  << " " << m_crossSections[i + 1]->startAdmittance().cols() << endl;
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

void Acoustic3dSimulation::propagateVelocityPress(Eigen::MatrixXcd& startVelocity,
  Eigen::MatrixXcd& startPressure, double freq, int startSection, int endSection)
{
  int direction;
  if (startSection > endSection)
  {
    direction = -1;
  }
  else
  {
    direction = 1;
  }
  propagateVelocityPress(startVelocity, startPressure, freq, startSection, endSection,
    direction);
}

void Acoustic3dSimulation::propagateVelocityPress(Eigen::MatrixXcd &startVelocity,
  Eigen::MatrixXcd& startPressure, double freq, int startSection, 
  int endSection, int direction)
{
  Eigen::MatrixXcd prevVelo(startVelocity), prevPress(startPressure);
  vector<Eigen::MatrixXcd> tmpQ, P, Y;
  vector<Matrix> F;
  Matrix G;
  int numSec(m_crossSections.size());
  int numX(m_simuParams.numIntegrationStep);
  int nextSec, nI, nNs;
  Eigen::MatrixXcd pressure;
  double areaRatio, tau, scaling;
  complex<double> wallInterfaceAdmit(1i * 2. * M_PI * freq * 
    m_simuParams.thermalBndSpecAdm / m_simuParams.sndSpeed);

  //ofstream log;
  //log.open("log.txt", ofstream::app);
  //log << "\n\nStart velocity propagation" << endl;
  //log << "Direction " << direction << endl;
  
  // loop over sections
  for (int i(startSection); i != (endSection); i += direction)
  {

    nextSec = i + direction;

    //log << "section " << i << " next sec " << nextSec << endl;

    m_crossSections[i]->clearAxialVelocity();
    m_crossSections[i]->clearAcPressure();
    m_crossSections[i]->setQdir(direction);
    m_crossSections[i]->setPdir(direction);
    nI = m_crossSections[i]->numberOfModes();
    nNs = m_crossSections[nextSec]->numberOfModes();

    // propagate axial velocity and acoustic pressure in the section
    switch (m_simuParams.propMethod) {
    case MAGNUS:
      //m_crossSections[i]->propagatePressureVelocityRiccati(prevVelo, prevPress,m_simuParams, 
      //  m_crossSections[nextSec]->area(), freq, (double)direction);
      m_crossSections[i]->propagateMagnus(prevPress, m_simuParams,
        freq, (double)direction, PRESSURE);
      tmpQ.clear(); P.clear(); Y.clear();
      P = m_crossSections[i]->P();
      Y = m_crossSections[i]->Y();
      numX = Y.size();
      for (int pt(0); pt < numX; pt++ )
      {
        if (numX > 1) {
          if (direction == 1) {
            tau = (double)pt / (double)(numX - 1);
          }
          else
          {
            tau = (double)(numX - 1 - pt) / (double)(numX - 1);
          }
        }
        else
        {
          tau = 1.;
        }
        scaling = m_crossSections[i]->scaling(tau);
        tmpQ.push_back(Y[numX - 1 - pt] * P[pt]);
        //tmpQ.push_back(Y[numX - 1 - pt] * P[pt] / pow(scaling, 2));
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

    //log << "i = " << i << " nextSec = " << nextSec << endl;

    prevVelo = Eigen::MatrixXcd::Zero(m_crossSections[nextSec]->numberOfModes(), 1);
    prevPress = Eigen::MatrixXcd::Zero(m_crossSections[nextSec]->numberOfModes(), 1);
    switch (m_simuParams.propMethod)
    {
    case MAGNUS:
      if (direction == -1)
      {
          // if the section contracts: area(i) > area(ns)
          if ((m_crossSections[i]->area()*
                     pow(m_crossSections[i]->scaleIn(), 2)) >
            (m_crossSections[nextSec]->area() *
             pow(m_crossSections[nextSec]->scaleOut(), 2)))
              {
            //log << "area(i) > area(ns) compute pressure" << endl;
                prevPress += F[0] *
                    m_crossSections[i]->Pin()
                  / m_crossSections[i]->scaleIn()
                  / m_crossSections[nextSec]->scaleOut();
                prevVelo +=
                  m_crossSections[nextSec]->Yout() * prevPress;
              }
          // if the section expends: area(i) < area(ns)
          else
          {
            //log << "area(i) < area(ns) compute velocity" << endl;
            prevVelo +=
            //(Matrix::Identity(nNs, nNs)
            //  + wallInterfaceAdmit *
            //  ////m_crossSections[nextSec]->curvature()*
            //  G * m_crossSections[nextSec]->Zin()).inverse() *
            F[0] * m_crossSections[i]->Qin()
              * m_crossSections[i]->scaleIn()
              * m_crossSections[nextSec]->scaleOut();
            prevPress +=
              m_crossSections[nextSec]->Zout() * prevVelo;
          }
      }
      else
      {
          // if the section contracts: area(i) > area(ns)
          if ((m_crossSections[i]->area()*
                     pow(m_crossSections[i]->scaleOut(), 2)) >
            (m_crossSections[nextSec]->area() *
             pow(m_crossSections[nextSec]->scaleIn(), 2)))
          {
            //log << "area(i) > area(ns), compute pressure" << endl;
            prevPress += 
                (F[0].transpose()) *
              m_crossSections[i]->Pout()
              //* m_crossSections[i]->scaleOut()
              /// m_crossSections[nextSec]->scaleIn()
              ;
            prevVelo +=
              m_crossSections[nextSec]->Yin() * prevPress;
          }
          // if the section expends: area(i) < area(ns)
          else
          {
            //log << "area(i) < area(ns), compute velocity" << endl;
          prevVelo +=
            (Matrix::Identity(nNs, nNs)
              - wallInterfaceAdmit *
            //  ////m_crossSections[nextSec]->curvature() * 
              G * m_crossSections[nextSec]->Zin()).inverse() *
              (F[0].transpose()) * m_crossSections[i]->Qout()
            //* m_crossSections[nextSec]->scaleIn()
            /// m_crossSections[i]->scaleOut()
            ;
          prevPress += m_crossSections[nextSec]->Zin() * prevVelo;
          }
      }

      break;
    case STRAIGHT_TUBES:
      //areaRatio = sqrt(m_crossSections->at(i + 1).area() /
      //  m_crossSections->at(i).area());

      areaRatio = sqrt(max(m_crossSections[nextSec]->area(),
        m_crossSections[i]->area()) /
        min(m_crossSections[nextSec]->area(),
          m_crossSections[i]->area()));
      if (direction == -1)
      {
        // if the section expends
        if (m_crossSections[nextSec]->area() >
          m_crossSections[i]->area())
        {
          prevVelo += areaRatio * (F[0].transpose()) *
            m_crossSections[i]->Qin();
          prevPress +=
            m_crossSections[nextSec]->Zout() * prevVelo;
        }
        // if the section contracts
        else
        {
          prevPress += areaRatio * (F[0].transpose()) *
            m_crossSections[i]->Pin();
          prevVelo +=
            m_crossSections[nextSec]->Yout() * prevPress;
        }
      }
      else
      {
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
      }
      break;
    }
  }

  // propagate in the last section
  m_crossSections[endSection]->clearAxialVelocity();
  m_crossSections[endSection]->clearAcPressure();
  m_crossSections[endSection]->setQdir(direction);
  m_crossSections[endSection]->setPdir(direction);
  switch (m_simuParams.propMethod) {
  case MAGNUS:
    //m_crossSections[endSection]->propagatePressureVelocityRiccati(prevVelo, prevPress,
    //  m_simuParams, 2.*m_crossSections[endSection]->area(), freq,
    //  (double)direction);
    m_crossSections[endSection]->propagateMagnus(prevPress, m_simuParams,
      freq, (double)direction, PRESSURE);
    tmpQ.clear(); P.clear(); Y.clear();
    Y = m_crossSections[endSection]->Y();
    P = m_crossSections[endSection]->P();
    numX = Y.size();
    for (int pt(0); pt < numX; pt++)
    {
      if (numX > 1) {
        if (direction == 1) {
          tau = (double)pt / (double)(numX - 1);
        }
        else
        {
          tau = (double)(numX - 1 - pt) / (double)(numX - 1);
        }
      }
      else
      {
        tau = 1.;
      }
      scaling = m_crossSections[endSection]->scaling(tau);
      tmpQ.push_back(Y[numX - 1 - pt] * P[pt]);
      //tmpQ.push_back(Y[numX - 1 - pt] * P[pt] / pow(scaling, 2));
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
      //  << F[0].inverse().cols() << endl;
      //log << "Size Pend " << m_crossSections[i]->endAcPressure().rows()
      //  << " " << m_crossSections[i]->endAcPressure().cols() << endl;
      //Eigen::MatrixXcd tF = (Eigen::MatrixXcd)(F[0].transpose());
      //Eigen::MatrixXcd Pa = m_crossSections[i]->endAcPressure();
      //prevPress = tF.fullPivLu().solve(Pa);
      //prevPress = F[0].transpose().fullPivLu().inverse() * m_crossSections[i]->endAcPressure();
      prevPress = m_crossSections[i+1]->Yin().fullPivLu().inverse() * 
        F[0].transpose() * m_crossSections[i]->Yout() * 
        m_crossSections[i]->Pout();

      //log << "Size prevpress " << prevPress.rows() << " "
      //  << prevPress.cols() << endl;
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
// Extract the internal acoustic field at a point 
// (given in cartesian coordinates)

complex<double> Acoustic3dSimulation::acousticFieldInside(Point_3 queryPt)
{
  bool ptFound(false);
  int numSec(m_crossSections.size());
  Point_3 outPt;
  complex<double> field;

  //ofstream log("log.txt", ofstream::app);
  //log << "Start acoustic field inside numSec " << numSec << endl;

  for (int s(0); s < numSec; s++)
  {
    if (m_crossSections[s]->getCoordinateFromCartesianPt(queryPt, outPt))
    {
      //log << "Pt found " << outPt << endl;
      ptFound = true;
      field = m_crossSections[s]->interiorField(outPt, m_simuParams, PRESSURE);
      //log << "field computed" << endl;
      break;
    }
  }
  if (!ptFound)
  {
    field = complex<double>(NAN, NAN);
  }
  //log.close();

  return(field);
}

// **************************************************************************
// Extract the internal acoustic field in a plane

void Acoustic3dSimulation::acousticFieldInPlane(Eigen::MatrixXcd& field)
{
  double lx(m_simuParams.bboxField[1].x() - m_simuParams.bboxField[0].x());
  double ly(m_simuParams.bboxField[1].y() - m_simuParams.bboxField[0].y());
  int nPtx(round(lx * (double)m_simuParams.fieldResolution));
  int nPty(round(ly * (double)m_simuParams.fieldResolution));
  int cnt(0);
  int numSec(m_crossSections.size());
  Point_3 queryPt, outPt;
  bool ptFound;
  ofstream log("log.txt", ofstream::app);

  field.resize(nPty, nPtx);

  for (int i(0); i < nPtx; i++)
  {
    log << "Point " << cnt << " out of " << nPtx * nPty << "  "
      << queryPt << endl;
    for (int j(0); j < nPty; j++)
    {
      cnt++;

      //log << "Point " << cnt << " out of " << nPtx * nPty << "  "
  //<< queryPt << endl;

      // generate cartesian coordinates point to search
      queryPt = Point_3(lx * (double)i / (double)(nPtx - 1) + m_simuParams.bboxField[0].x(), 0.,
        ly * (double)j / (double)(nPty - 1) + m_simuParams.bboxField[0].y());

      //log << "Query point coordinate computed " << queryPt << endl;

      field(j, i) = acousticFieldInside(queryPt);

      //log << "Field computed" << endl;
    }
  }
  log.close();
}

// **************************************************************************
// Run a static simulation

void Acoustic3dSimulation::staticSimulation(VocalTract* tract)
{
  int numSec, lastSec; 
  double freq, freqSteps((double)SAMPLING_RATE/2./(double)m_numFreq);
  int numFreqComputed((int)ceil(m_simuParams.maxComputedFreq / freqSteps));
  int mn; // mode number in the last section
  complex<double> pressure, velocity;
  Eigen::MatrixXcd endPressAmpl, endVelAmpl, startPressure,
    radAdmit, radImped, upStreamImpAdm, totalImped, prevVelo, prevPress
    ,pp ,pm;
  Matrix F;
  Point_3 internalFieldPt;
  Point internalFieldPt2D;
  ofstream prop;

  generateLogFileHeader(true);
  ofstream log("log.txt", ofstream::app);

  auto startTot = std::chrono::system_clock::now();

  //******************************
  // create the cross-sections
  //******************************

  if (createCrossSections(tract, false))
  {
    log << "Geometry successfully imported" << endl;
  }
  else
  {
    log << "Importation failed" << endl;
  }

  numSec = m_crossSections.size();
  lastSec = numSec - 1;

  log << "Number of sections: " << numSec << endl;

  // export cross-sections parameters
  for (int i(0); i < numSec; i++)
  {
    log << "Scetion " << i << endl;
    log << *m_crossSections[i] << endl;
  }

  //// exctract some contours
  //ostringstream os;
  //for (int i(97); i < 103; i++)
  //{
  //  os.str("");
  //  os.clear();
  //  os << "c" << i << ".txt";
  //  prop.open(os.str());
  //  for (auto it : m_crossSections[i]->contour())
  //  {
  //    prop << it.x() << "  " << it.y() << endl;
  //  }
  //  prop.close();
  //}

  // FIXME: update the values of m_idxSecNoiseSource and m_idxConstriction
  // in the 3D simu properties dialog when they are changed
  // check if the noise source index is within the indexes range
  if (m_idxSecNoiseSource >= numSec)
  {
    m_idxSecNoiseSource = numSec - 1;
  }
  // check if the constriction location is within the indexes range
  if (m_idxConstriction >= numSec)
  {
    m_idxConstriction = numSec - 1;
  }

  //log << "before compute mode" << endl;

  // create the mesh and compute modes
  computeMeshAndModes();
  log << "Modes computed" << endl;
  mn = m_crossSections[lastSec]->numberOfModes();

  // compute junction matrices
  auto start = std::chrono::system_clock::now();
  computeJunctionMatrices(false);
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - startTot;
  log << "Junction matrices computed " << 
    elapsed_seconds.count() << " s" << endl;

  // check junction matrix first term
  //prop.open("chkF.txt");
  //double F00, S, Sn, A, An, sc, scn;
  //Matrix modes;
  //for (int i(0); i < numSec - 2; i++)
  //{
    //sc = m_crossSections[i]->scaleOut();
    //scn = m_crossSections[i + 1]->scaleIn();
    //S = m_crossSections[i]->area();
    //Sn = m_crossSections[i + 1]->area();
    
    //F00 = min(S*pow(sc,2), Sn*pow(scn, 2)) / sqrt(S * Sn)/sc/scn;

    //modes = m_crossSections[i]->getModes();
    //A = modes(0, 0);
    //modes = m_crossSections[i + 1]->getModes();
    //An = modes(0, 0);

    //prop << m_crossSections[i]->getMatrixF()[0](0, 0) << "  "
      //<< F00 << "  " << S << "  " << Sn << "  "
      //<< A << "  " << An << "  " << sc << "  " << scn << endl;
  //}
  //prop.close();
  
  //// export scaling factors
  //prop.open("sc.txt");
  //for (int i(0); i < numSec - 2; i++)
  //{
  //  prop << m_crossSections[i]->scaleIn() << "  "
  //    << m_crossSections[i]->scaleOut() << endl;
  //}
  //prop.close();

  // generate source matrix
  Eigen::MatrixXcd inputVelocity(Eigen::MatrixXcd::Zero(
  m_crossSections[0]->numberOfModes(), 1));
  //log << "inputVelocity initialized" << endl;
  Eigen::MatrixXcd inputPressure(Eigen::MatrixXcd::Zero(
    m_crossSections[0]->numberOfModes(), 1));
  //log << "inputPressure initialzed" << endl;
  complex<double> V0(1.0, 0.0);
  // for a constant input velocity q = -j * w * rho * v 
  inputVelocity(0, 0) = -1i * 2. * M_PI * m_simuParams.volumicMass * V0;
  Eigen::MatrixXcd inputPressureNoise;

  if (m_idxSecNoiseSource < lastSec)
  {
    // extract the junction matrix of the noise source section
    F = m_crossSections[m_idxSecNoiseSource]->getMatrixF()[0];
    log << "noise source junction matrix extracted" << endl;
    inputPressureNoise.setZero(
        m_crossSections[m_idxSecNoiseSource]->numberOfModes(), 1);
    inputPressureNoise(0, 0) = V0;
  }

  log << "source generated" << endl;

  //*****************************************************************************
  //  Compute points for acoustic field computation
  //*****************************************************************************

  vector<Point_3> radPts;

  //// generate points distributed on a plane
  //int nPtsX(50), nPtsY(100);
  //double dX(5./(double)(nPtsX));
  //for (int i(0); i < nPtsX; i++)
  //{
  //  for (int j(0); j < nPtsY; j++)
  //  {
  //    radPts.push_back(Point_3((double)(i + 1) * dX, 
  //      ((double)(j) - (double)(nPtsY)/2.) * dX, 0.));
  //  }
  //}

  //// generate points on a cicle arc
  //double radius(50.);
  //int nPts(91);
  //for (int p(0); p < nPts; p++)
  //{
  //  radPts.push_back(Point_3(radius * sin((double)p * M_PI / (double)(nPts-1)),
  //    radius * cos((double)p * M_PI / (double)(nPts-1)), 0.));
  //}

  // generate one point in front
  radPts.push_back(Point_3(25., 0., 0.));

  // get the coordinates of the point where to compute the internal field
  Transformation rotate(CGAL::ROTATION, sin(M_PI / 2.),
    cos(M_PI / 2.));
  Transformation translate(CGAL::TRANSLATION, 0.3 * rotate(
    m_crossSections[lastSec]->normalOut()
  ));
  //internalFieldPt2D = translate(m_crossSections[lastSec]->ctrLinePtOut());
  internalFieldPt2D = m_crossSections[lastSec]->ctrLinePtOut();
  log << "Internal field point 2D " << internalFieldPt2D << endl;
  internalFieldPt = Point_3(internalFieldPt2D.x(), 0., internalFieldPt2D.y());
  log << "Internal field point " << internalFieldPt << endl;

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

  if (m_mouthBoundaryCond == RADIATION)
  {
    // Precompute the radiation impedance and admittance at a few frequencies
    preComputeRadiationMatrices(16, lastSec);
  }

  end = std::chrono::system_clock::now();
  elapsed_seconds = end - startTot;
  log << "Time precomputations " << elapsed_seconds.count() << " s" << endl;

  //*****************************************************************************

  // Create the progress dialog
  progressDialog = new wxGenericProgressDialog("Transfer function computation",
    "Wait until the transfer function computation finished or press [Cancel]",
    numFreqComputed, NULL,
    wxPD_CAN_ABORT | wxPD_AUTO_HIDE | wxPD_ELAPSED_TIME);

  start = std::chrono::system_clock::now();

  // loop over frequencies
  for (int i(0); i < numFreqComputed; i++)
  {
    freq = max(0.1, (double)i * freqSteps);
    log << "################################" << endl;
    log << "frequency " << i+1 << "/" << numFreqComputed << " f = " << freq
      << " Hz" << endl;

    // generate radiation impedance and admittance
    start = std::chrono::system_clock::now();
    // m_crossSections[lastSec]->characteristicImpedance(radImped,
    //   freq, m_soundSpeed, m_volumicMass, 0);
    //m_crossSections[lastSec]->characteristicAdmittance(radAdmit,
    //  freq, m_soundSpeed, m_volumicMass, 0);

    switch (m_mouthBoundaryCond)
    {
    case RADIATION:
      interpolateRadiationImpedance(radImped, freq, lastSec);
      interpolateRadiationAdmittance(radAdmit, freq, lastSec);
      break;
    case ADMITTANCE_1:
      radAdmit.setZero(mn, mn);
      radAdmit.diagonal() = Eigen::VectorXcd::Constant(mn, complex<double>(
        pow(m_crossSections[lastSec]->scaleOut(), 2), 0.));
      radImped = radAdmit.inverse();
      break;
    }

    end = std::chrono::system_clock::now();
    elapsed_seconds = end - start;
    log << "Time charac imped/admit " << elapsed_seconds.count() << " s" << endl;

    // propagate impedance and admittance
    start = std::chrono::system_clock::now();
    //propagateAdmit(radAdmit, freq, m_method);
    propagateImpedAdmit(radImped, radAdmit, freq, lastSec, 0);
    
    end = std::chrono::system_clock::now();
    elapsed_seconds = end - start;
    log << "Time impedance " << elapsed_seconds.count() << " s" << endl;

    //// extract admittance
    //prop.open("adm.txt", ofstream::app);
    //for (int c(0); c < lastSec; c++)
    //{
    //  prop << abs(m_crossSections[c]->Yin()(0, 0)) << "  ";
    //}
    //prop << endl;
    //prop.close();

    //// extract impedance
    //prop.open("imp.txt", ofstream::app);
    //for (int c(0); c < lastSec; c++)
    //{
    //  prop << abs(m_crossSections[c]->Zin()(0, 0)) << "  ";
    //}
    //prop << endl;
    //prop.close();

    //// print radiation impedance
    //prop.open("radR.txt", ofstream::app);
    //prop << m_crossSections[lastSec]->endImpedance().real() << endl;
    //prop.close();
    //prop.open("radI.txt", ofstream::app);
    //prop << m_crossSections[lastSec]->endImpedance().imag() << endl;
    //prop.close();
    //prop.open("admR.txt", ofstream::app);
    //prop << m_crossSections[lastSec]->endAdmittance().real() << endl;
    //prop.close();
    //prop.open("admI.txt", ofstream::app);
    //prop << m_crossSections[lastSec]->endAdmittance().imag() << endl;
    //prop.close();

    // propagate axial velocity and pressure
    start = std::chrono::system_clock::now();
    //startPressure = m_crossSections[1]->startAdmittance().fullPivLu().inverse()
    //  * inputVelocity;
    //propagatePressure(startPressure, freq, m_method);
    //log << "Compute impedance in first section " <<
    //  m_crossSections[0]->computeImpedance() << endl;
    inputVelocity(0, 0) = -1i * 2. * M_PI * freq * m_simuParams.volumicMass * 
      m_crossSections[0]->area();
    //log << "input velocity\n" << inputVelocity << endl;
    inputPressure = m_crossSections[0]->Zin() * inputVelocity;
    propagateVelocityPress(inputVelocity, inputPressure, freq, 0, lastSec);
    end = std::chrono::system_clock::now();
    elapsed_seconds = end - start;
    log << "Time velocity pressure " << elapsed_seconds.count() << " s" << endl;

    //// extract pressure
    //prop.open("pamp.txt", ofstream::app);
    //for (int c(0); c < lastSec; c++)
    //{
    //  prop << abs(m_crossSections[c]->Pin()(0,0)) << "  "
    //    << abs(m_crossSections[c]->Pout()(0, 0)) << "  ";
    //}
    //prop << endl;
    //prop.close();

    //// extract potential
    //prop.open("qamp.txt", ofstream::app);
    //for (int c(0); c < lastSec; c++)
    //{
    //  prop << abs(m_crossSections[c]->Qin()(0, 0)) << "  "
    //    << abs(m_crossSections[c]->Qout()(0, 0)) << "  ";
    //}
    //prop.close();

    start = std::chrono::system_clock::now();

    //*****************************************************************************
    //  Compute acoustic pressure inside and outside
    //*****************************************************************************

    RayleighSommerfeldIntegral(radPts, radPress, freq, lastSec);

    //// save radiated field
    //prop.open("radR.txt", ofstream::app);
    //prop << radPress.transpose().real() << endl;
    //prop.close();
    //prop.open("radI.txt", ofstream::app);
    //prop << radPress.transpose().imag() << endl;
    //prop.close();

    spectrum.setValue(i, radPress(0,0));

    //pressure = acousticFieldInside(internalFieldPt);
    pressure = m_crossSections[lastSec]->pout(Point(0., 0.));
    spectrumConst.setValue(i, pressure);

    log << "Radiated pressure computed" << endl;

    //*****************************************************************************
    //  Compute transfer function of the noise source
    //*****************************************************************************

    if (m_idxSecNoiseSource < lastSec)
    {
      //// get the transfer function from the glottis to the constriction
      //// location for vowels
      //spectrumConst.setValue(i, m_crossSections[m_idxConstriction]->Pout()(0, 0));

      // save the input impedance of the upstream part
      // if the section expends
      if ((pow(m_crossSections[m_idxSecNoiseSource + 1]->scaleIn(), 2) *
        m_crossSections[m_idxSecNoiseSource + 1]->area()) >
        (pow(m_crossSections[m_idxSecNoiseSource]->scaleOut(), 2) * 
        m_crossSections[m_idxSecNoiseSource]->area()))
      {
        upStreamImpAdm = m_crossSections[m_idxSecNoiseSource]->Zout();
      }
      // if the section contracts
      else
      {
        upStreamImpAdm = m_crossSections[m_idxSecNoiseSource]->Yout();
      }
      log << "upstream input impedance saved" << endl;

      // set glottis boundary condition
      
      switch (m_glottisBoundaryCond)
      {
      case HARD_WALL:
        radImped.setZero();
        radImped.diagonal().setConstant(100000.);
        radAdmit.setZero();
        radAdmit.diagonal().setConstant(1. / 100000.);
        break;
      case IFINITE_WAVGUIDE:
        m_crossSections[0]->characteristicImpedance(radImped, freq, m_simuParams);
        m_crossSections[0]->characteristicAdmittance(radAdmit, freq, m_simuParams);
        break;
      }

      // propagate impedance and admittance from the glottis to the location
      // of the second sound source
      propagateImpedAdmit(radImped, radAdmit, freq, 0, m_idxSecNoiseSource);

      log << "Imped admit propagated" << endl;

      // compute the pressure and the velocity at the entrance of the next section
      // if the section expends
      if ((pow(m_crossSections[m_idxSecNoiseSource + 1]->scaleIn(), 2) * 
        m_crossSections[m_idxSecNoiseSource +1]->area()) >
        (pow(m_crossSections[m_idxSecNoiseSource]->scaleOut(), 2) * 
        m_crossSections[m_idxSecNoiseSource]->area()))
      {
        prevVelo = (F.transpose()) * ((freq * upStreamImpAdm - freq * 
          m_crossSections[m_idxSecNoiseSource]->Zout()).householderQr()
          .solve(inputPressureNoise));
        prevPress = freq *
          m_crossSections[m_idxSecNoiseSource +1]->Zin() * prevVelo;
      }
      // if the section contracts
      else
      {
        prevPress = (F.transpose()) * ((upStreamImpAdm -
          m_crossSections[m_idxSecNoiseSource]->Yout()).householderQr()
          .solve( -m_crossSections[m_idxSecNoiseSource]->Yout() *
            inputPressureNoise));

        prevVelo =
          m_crossSections[m_idxSecNoiseSource +1]->Yin()* prevPress;
      }

      log << "Previous pressure/velocity computed" << endl;
      
      // propagate the pressure and the velocity in the upstream part
      propagateVelocityPress(prevVelo, prevPress, freq, 
        min(m_idxSecNoiseSource +1, lastSec), lastSec);

      log << "Velocity and pressure propagated" << endl;

      RayleighSommerfeldIntegral(radPts, radPress, freq, lastSec);

      end = std::chrono::system_clock::now();

      spectrumNoise.setValue(i, radPress(0,0));
    }
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
    spectrumNoise.re[i] = spectrumNoise.re[2 * m_numFreq - i - 1];
    spectrumNoise.im[i] = -spectrumNoise.im[2 * m_numFreq - i - 1];
    spectrumConst.re[i] = spectrumConst.re[2 * m_numFreq - i - 1];
    spectrumConst.im[i] = -spectrumConst.im[2 * m_numFreq - i - 1];
  }

  // create the file containing the transfer function
  prop.open("press.txt");
  //prop << "num_points: " << spectrum.N << endl;
  //prop << "frequency_Hz  magnitude  phase_rad" << endl;
  for (int i(0); i < spectrum.N; i++)
  {
    freq = max(0.1, (double)i * freqSteps);
    prop << freq << "  " 
      << spectrum.getMagnitude(i) << "  " 
      << spectrum.getPhase(i) << "  "
      << spectrumNoise.getMagnitude(i) << "  "
      << spectrumNoise.getPhase(i) << "  "
      << spectrumConst.getMagnitude(i) << "  "
      << spectrumConst.getPhase(i)
      << endl;
  }
  prop.close();

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
// Run a simulation at a specific frequency and compute the aoustic field 

void Acoustic3dSimulation::computeAcousticField(VocalTract* tract)
{
  int numSec, lastSec; 
  double freq(m_simuParams.freqField);

  generateLogFileHeader(true);
  ofstream log("log.txt", ofstream::app);

  auto start = std::chrono::system_clock::now();

  // create the cross-sections
  if (createCrossSections(tract, false))
  {
    log << "Geometry successfully imported" << endl;
  }
  else
  {
    log << "Importation failed" << endl;
  }
  numSec = m_crossSections.size();
  lastSec = numSec - 1;
  log << "Number of sections: " << numSec << endl;

  // Export cross-sections parameters
  for (int i(0); i < numSec; i++)
  {
    log << "\nSection " << i << endl;
    log << *m_crossSections[i] << endl;
  }

  // create the mesh and compute modes
  computeMeshAndModes();
  log << "Modes computed" << endl;

  // Compute junction matrices
  computeJunctionMatrices(false);
  log << "Junction matrices computed " << endl;

  // generate source matrix
  Eigen::MatrixXcd inputVelocity(Eigen::MatrixXcd::Zero(
  m_crossSections[0]->numberOfModes(), 1));
  Eigen::MatrixXcd inputPressure(Eigen::MatrixXcd::Zero(
    m_crossSections[0]->numberOfModes(), 1));
  complex<double> V0(1.0, 0.0);
  // for a constant input velocity q = -j * w * rho * v 
  inputVelocity(0, 0) = -1i * 2. * M_PI * m_simuParams.volumicMass * V0;
  log << "source generated" << endl;

  // Compute the radiation impedance and admittance 
  Eigen::MatrixXcd radImped, radAdmit;
  int mn(m_crossSections[lastSec]->numberOfModes());
  switch (m_mouthBoundaryCond)
  {
  case RADIATION:
    radiationImpedance(radImped, freq, 15., lastSec);
    radAdmit = radImped.inverse();
    break;
  case ADMITTANCE_1:
    radAdmit = Eigen::MatrixXcd::Zero(mn, mn);
    radAdmit.diagonal().setConstant(complex<double>(
          pow(m_crossSections[lastSec]->scaleOut(), 2), 0.));
    radImped = radAdmit.inverse();
  }
  log << "End impedance computed" << endl;

  // propagate impedance and admittance
  propagateImpedAdmit(radImped, radAdmit, freq, lastSec, 0);
  log << "Impedance and admittance computed" << endl;

  // propagate velocity and pressure
  inputVelocity(0, 0) = -1i * 2. * M_PI * freq * m_simuParams.volumicMass * 
    m_crossSections[0]->area();
  inputPressure = m_crossSections[0]->Zin() * inputVelocity;
  propagateVelocityPress(inputVelocity, inputPressure, freq, 0, lastSec);
  log << "Pressure and velocity computed" << endl;

  //******************************************************
  // Extract acoustic field
  //******************************************************

  ofstream ofs;
  ofs.open("sec.txt");
  // extract parameters of the sections
  for (int s(0); s < numSec; s++)
  {
    if (m_crossSections[s]->length() > 0.)
    {
      // entrance features
      ofs << m_crossSections[s]->ctrLinePt().x << "  "
        << m_crossSections[s]->ctrLinePt().y << "  "
        << m_crossSections[s]->normal().x << "  "
        << m_crossSections[s]->normal().y << "  "
        << m_crossSections[s]->curvRadius() << "  "
        << m_crossSections[s]->circleArcAngle() << endl;

      // exist features
      ofs << m_crossSections[s]->ctrLinePtOut().x() << "  "
        << m_crossSections[s]->ctrLinePtOut().y() << "  "
        << m_crossSections[s]->normalOut().x() << "  "
        << m_crossSections[s]->normalOut().y() << "  "
        << m_crossSections[s]->curvRadius() << "  "
        << m_crossSections[s]->circleArcAngle() << endl;
    }
  }
  ofs.close();

  Eigen::MatrixXcd field;
  acousticFieldInPlane(field);
  ofs.open("field.txt");
  stringstream txtField;
  txtField << field.cwiseAbs();
  ofs << regex_replace(txtField.str(), regex("-nan\\(ind\\)"), "nan");
  ofs.close();

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  int hours(floor(elapsed_seconds.count() / 3600.));
  int minutes(floor((elapsed_seconds.count() - hours * 3600.) / 60.));
  double seconds(elapsed_seconds.count() - (double)hours * 3600. -
    (double)minutes * 60.);
  log << "\nAcoustic field computation time "
    << hours << " h " << minutes << " m " << seconds << " s" << endl;

  log.close();
}

// ****************************************************************************
// Run a simulation for a concatenation of cylinders

void Acoustic3dSimulation::coneConcatenationSimulation(string fileName)
{
  // for data extraction
  double endAdmit, freqField;
  vector<int> vIdx;
  vector<double> rads, shifts, scaleIn, scaleOut, lengths, curvAngles;
  ifstream ifs;
  string line, str;
  char separator(';');
  stringstream strSt;

  // for cross-section creation
  Polygon_2 contour;
  int nbAngles(100), nbSec;
  double angle, area, length, inAngle, inRadius, maxRad;
  double scalingFactors[2];
  vector<int> surfaceIdx(nbAngles, 0);

  // for solving the wave problem 
  bool reverse(false);
  int mn;
  double freq, freqMax, x, y, totalLength;
  Eigen::MatrixXcd radImped, radAdmit, inputVelocity, inputPressure;
  complex<double> pout, vout, yin;
  string::size_type idxStr;
  ofstream ofs, ofs2, ofs3, ofs4, ofs5, ofs6;
  Point ptOut;
  Point_3 pointComputeField;

  m_geometryImported = true; // to have the good bounding box for modes plot
  //m_simuParams.sndSpeed = 34400;

  generateLogFileHeader(true);
  ofstream log("log.txt", ofstream::app);
  log << "\nStart cylinder concatenation simulation" << endl;
  if (reverse) { log << "Propagation direction reversed" << endl; }
  log << "Geometry from file " << fileName << endl;

  //***************************
  // load geometry parameters
  //***************************

  ifs.open(fileName);
  if (!ifs.is_open())
  {
    log << "failed to opened parameters file" << endl;
  }
  else
  {
    // extract number of modes
    vIdx.clear();
    getline(ifs, line); // to remove comment line
    getline(ifs, line);
    strSt.clear();
    strSt.str(line);
    while (getline(strSt, str, separator))
    {
      vIdx.push_back(stoi(str));
    }

    // extract radius
    rads.clear();
    getline(ifs, line); // to remove comment line
    getline(ifs, line);
    strSt.clear();
    strSt.str(line);
    while (getline(strSt, str, separator))
    {
      rads.push_back(stod(str));
    }

    // extract contour shift
    shifts.clear();
    getline(ifs, line); // to remove comment line
    getline(ifs, line);
    strSt.clear();
    strSt.str(line);
    while (getline(strSt, str, separator))
    {
      shifts.push_back(stod(str));
    }

    // extract scaling in
    scaleIn.clear();
    getline(ifs, line); // to remove comment line
    getline(ifs, line);
    strSt.clear();
    strSt.str(line);
    while (getline(strSt, str, separator))
    {
      scaleIn.push_back(stod(str));
    }

    // extract scaling out
    scaleOut.clear();
    getline(ifs, line); // to remove comment line
    getline(ifs, line);
    strSt.clear();
    strSt.str(line);
    while (getline(strSt, str, separator))
    {
      scaleOut.push_back(stod(str));
    }

    // extract length
    lengths.clear();
    getline(ifs, line); // to remove comment line
    getline(ifs, line);
    strSt.clear();
    strSt.str(line);
    while (getline(strSt, str, separator))
    {
      lengths.push_back(stod(str));
    }

    // extract curvature angle
    curvAngles.clear();
    getline(ifs, line); // to remove comment line
    getline(ifs, line);
    strSt.clear();
    strSt.str(line);
    while (getline(strSt, str, separator))
    {
      curvAngles.push_back(stod(str));
    }

    // extract end admittance
    getline(ifs, line); // to remove comment line
    getline(ifs, line);
    strSt.clear();
    strSt.str(line);
    getline(strSt, str, separator);
    endAdmit =  stod(str); // * pow(scaleOut[1], 2); // open end admittance
    getline(strSt, str, separator);
    // wall admittance
    m_simuParams.thermalBndSpecAdm = complex<double>(stod(str), 0.);

    // extract index frequency field extraction
    getline(ifs, line); // to remove comment line
    getline(ifs, line);
    freqField = stod(line);

  }

  log << "Geometry parameters extracted" << endl;

  //***********************
  // Create cross-sections
  //***********************

  maxRad = 0.;
  m_crossSections.clear();
  for (int s(0); s < vIdx.size(); s++)
  { 
    maxRad = max(maxRad, rads[s]);
    // Generate a circular contour
    contour.clear();
    for (int i(0); i < nbAngles; i++)
    {
      angle = 2. * M_PI * (double)(i) / (double)(nbAngles);
      contour.push_back(Point(rads[s] * cos(angle), rads[s] * (sin(angle) + shifts[s])));
    }

    area = pow(rads[s], 2) * M_PI;
    scalingFactors[0] = scaleIn[s];
    scalingFactors[1] = scaleOut[s];
    length = lengths[s];
    inAngle = curvAngles[s];
    inRadius = length / inAngle;
    addCrossSectionFEM(area, sqrt(area) / m_meshDensity, contour,
      surfaceIdx, length, Point2D(0., 0.), Point2D(0., 1.),
      scalingFactors);
    m_crossSections[s]->setCurvatureRadius(inRadius);
    m_crossSections[s]->setCurvatureAngle(inAngle);
    m_crossSections[s]->setModesNumber(vIdx[s]);
    if (s > 0){m_crossSections[s]->setPreviousSection(s - 1);}
    if (s < vIdx.size() -1) {m_crossSections[s]->setNextSection(s+1);}

    log << "Section " << s << " created" << endl;
  }
  nbSec = m_crossSections.size();
  // define bounding box for modes and mesh display
  m_maxCSBoundingBox.first = Point2D(-2. * maxRad, -2. * maxRad);
  m_maxCSBoundingBox.second = Point2D(2. * maxRad, 2. * maxRad);

  log << nbSec << " sections created" << endl;

  // Export cross-sections parameters
  for (int i(0); i < nbSec; i++)
  {
    log << "\nSection " << i << endl;
    log << *m_crossSections[i] << endl;
  }

  //*********************
  // solve wave problem
  //*********************

  computeMeshAndModes();
  log << "Modes computed" << endl;

  computeJunctionMatrices(false);
  log << "Junctions computed" << endl;

  // initialize input pressure and velocity vectors
  mn = m_crossSections[0]->numberOfModes();
  inputPressure.setZero(mn, 1);
  inputVelocity.setZero(mn, 1);

  freqMax = 10000.;
  m_numFreq = 501;
  idxStr = fileName.find_last_of("/\\");
  str = fileName.substr(0, idxStr + 1) + "tfMM.txt"; 
  ofs.open(str);
  // define the coordinate of the point where the acoustic field is computed
  // for transfer fucntion computation
  if (reverse) { ptOut = Point(0., shifts[0] * rads[0]); }
  else { ptOut = Point(0., shifts.back() * rads.back()); }
  log << "Point for transfer function computation " << ptOut << endl;
  for (int i(0); i < m_numFreq; i++)
  {
    freq = max(0.1, freqMax * (double)i / (double)(m_numFreq - 1));
    log << "f = " << freq << " Hz" << endl;

    if (reverse)
    {
      radAdmit.setZero(vIdx[0], vIdx[0]);
      radAdmit.diagonal() = Eigen::VectorXcd::Constant(vIdx[0], 
        complex<double>(pow(m_crossSections[0]->scaleIn(), 2) * endAdmit, 0.));
    }
    else
    {
      radAdmit.setZero(vIdx.back(), vIdx.back());
      radAdmit.diagonal() = Eigen::VectorXcd::Constant(vIdx.back(), 
        complex<double>(pow(m_crossSections[nbSec - 1]->scaleOut(), 2) * endAdmit, 0.));
    }
    radImped = radAdmit.inverse();

    if (m_simuParams.propMethod == STRAIGHT_TUBES)
    {
      radImped *= -1i * 2. * M_PI * freq * m_simuParams.volumicMass
        / m_crossSections[nbSec - 1]->area();
      radAdmit /= -1i * 2. * M_PI * freq * m_simuParams.volumicMass
        / m_crossSections[nbSec - 1]->area();
    }

    if (reverse) { propagateImpedAdmit(radImped, radAdmit, freq, 0, nbSec - 1, 1); }
    else { propagateImpedAdmit(radImped, radAdmit, freq, nbSec - 1, 0, -1); }

    //log << "Impedance propagated" << endl;

    if (reverse)
    {
      inputVelocity(0, 0) = -1i * 2. * M_PI * freq * m_simuParams.volumicMass
        * pow(m_crossSections[nbSec - 1]->scaleOut(), 2)
        * sqrt(m_crossSections[nbSec - 1]->area());
      inputPressure = m_crossSections[nbSec - 1]->Zout() * inputVelocity;
      propagateVelocityPress(inputVelocity, inputPressure, freq, nbSec - 1, 0, -1);
      //log << "Velocity and pressure propagated" << endl;
      // compute the acoustic pressure and the particle velocity at the 
      // center of the exit surface
      pout = m_crossSections[0]->pin(ptOut); // reversed order
      vout = -m_crossSections[0]->qin(ptOut)
        / 1i / 2. / M_PI / freq / m_simuParams.volumicMass;
      yin = m_crossSections[0]->interiorField(
        Point_3(m_crossSections[0]->length(), ptOut.x(), ptOut.y()), 
        m_simuParams, ADMITTANCE);
    }
    else
    {
      inputVelocity(0, 0) = -1i * 2. * M_PI * freq * m_simuParams.volumicMass
        * pow(m_crossSections[0]->scaleIn(), 2)
        * sqrt(m_crossSections[0]->area());
      inputPressure = m_crossSections[0]->Zin() * inputVelocity;
      propagateVelocityPress(inputVelocity, inputPressure, freq, 0, nbSec - 1, 1);
      //log << "Velocity and pressure propagated" << endl;
      // compute the acoustic pressure and the particle velocity at the 
      // center of the exit surface
      pout = m_crossSections[nbSec - 1]->pout(ptOut);
      vout = -m_crossSections[nbSec - 1]->qout(ptOut)
        / 1i / 2. / M_PI / freq / m_simuParams.volumicMass;
      yin = m_crossSections[nbSec - 1]->interiorField(
        Point_3(0., ptOut.x(), ptOut.y()), m_simuParams, ADMITTANCE);

    }

    // write result in a text file
    ofs << freq << "  "
      << "  " << abs(vout) // modulus of particle velocity
      << "  " << arg(vout) // phase of particle velocity
      << "  " << abs(pout) // modulus of acoustic pressure
      << "  " << arg(pout) // phase of acoustic pressure
      << "  " << abs(yin)  // modulus of the input impedance
      << "  " << arg(yin)  // phase of the input impedance
      << endl;
  }
  ofs.close();

  //******************************************************
  // Extract the acoustic field at a specified frequency
  //******************************************************

  if (freqField > 0.)
  {
    log << "Compute acoustic field at the frequency " << freqField << endl;

    freq = freqField;

    if (reverse)
    {
      radAdmit.setZero(vIdx[0], vIdx[0]);
      radAdmit.diagonal() = Eigen::VectorXcd::Constant(vIdx[0], complex<double>(
        pow(m_crossSections[0]->scaleIn(), 2) * endAdmit, 0.));
    }
    else
    {
      radAdmit.setZero(vIdx.back(), vIdx.back());
      radAdmit.diagonal() = Eigen::VectorXcd::Constant(vIdx.back(), complex<double>(
        pow(m_crossSections[nbSec - 1]->scaleOut(), 2) * endAdmit, 0.));
    }
    radImped = radAdmit.inverse();

    if (m_simuParams.propMethod == STRAIGHT_TUBES)
    {
      radImped *= -1i * 2. * M_PI * freq * m_simuParams.volumicMass
        / m_crossSections[nbSec - 1]->area();
      radAdmit /= -1i * 2. * M_PI * freq * m_simuParams.volumicMass
        / m_crossSections[nbSec - 1]->area();
    }

    if (reverse) { propagateImpedAdmit(radImped, radAdmit, freq, 0, nbSec - 1, 1); }
    else { propagateImpedAdmit(radImped, radAdmit, freq, nbSec - 1, 0, -1); }

    log << "Impedance propagated" << endl;

    if (reverse)
    {
      inputVelocity(0, 0) = -1i * 2. * M_PI * freq * m_simuParams.volumicMass
        * pow(m_crossSections[nbSec - 1]->scaleOut(), 2)
        * sqrt(m_crossSections[nbSec - 1]->area());
      inputPressure = m_crossSections[nbSec - 1]->Zout() * inputVelocity;
      propagateVelocityPress(inputVelocity, inputPressure, freq, nbSec - 1, 0, -1);
    }
    else
    {
      inputVelocity(0, 0) = -1i * 2. * M_PI * freq * m_simuParams.volumicMass
        * pow(m_crossSections[0]->scaleIn(), 2)
        * sqrt(m_crossSections[0]->area());
      inputPressure = m_crossSections[0]->Zin() * inputVelocity;
      propagateVelocityPress(inputVelocity, inputPressure, freq, 0, nbSec - 1, 1);
    }
    log << "Velocity and pressure propagated" << endl;

    // export the field in text files
    str = fileName.substr(0, idxStr + 1) + "q.txt"; 
    ofs2.open(str);
    str = fileName.substr(0, idxStr + 1) + "Y.txt"; 
    ofs3.open(str);
    str = fileName.substr(0, idxStr + 1) + "p.txt"; 
    ofs4.open(str);
    str = fileName.substr(0, idxStr + 1) + "cx.txt";
    ofs5.open(str);
    str = fileName.substr(0, idxStr + 1) + "cy.txt";
    ofs6.open(str);
    int cnt(0);
    totalLength = 0.;
    for (int s(0); s < nbSec; s++)
    {
      if (m_crossSections[s]->length() > 0.)
      {
        cnt = 0;
        for (int px(0); px < 100; px++)
        {
          for (int pz(0); pz < 51; pz++)
          {
            x = lengths[s] * (double)(px) / 99.;
            ofs5 << x + totalLength << "  ";
            y = rads[s] * (2. * (double)(pz) / 50. - 1.);
            ofs6 << y * m_crossSections[s]->scaling((double)(px) / 99.) << "  ";
            pointComputeField = Point_3(x, 0., y);
            ofs3 << abs(m_crossSections[s]->interiorField(pointComputeField, 
                  m_simuParams, ADMITTANCE)) << "  ";
            ofs2 << abs(m_crossSections[s]->q(pointComputeField, m_simuParams)) << "  ";
            ofs4 << abs(m_crossSections[s]->p(pointComputeField, m_simuParams)) << "  ";
            cnt++;
            log << "Point " << cnt << " over " << 100 * 51 << " " << pointComputeField << endl;
          }
          ofs2 << endl;
          ofs3 << endl;
          ofs4 << endl;
          ofs5 << endl;
          ofs6 << endl;
        }
        totalLength += m_crossSections[s]->length();
      }
    }
    ofs4.close();
    ofs3.close();
    ofs2.close();
  }

  log.close();
}

// ****************************************************************************
// Run specific simulation tests

void Acoustic3dSimulation::runTest(enum testType tType, string fileName)
{
  ofstream ofs, ofs2, ofs3, ofs4;
  ifstream ifs;
  string line, str;
  char separator(';');
  stringstream strs;
  string::size_type idxStr;
  ofstream log("log.txt", ofstream::app);
  log << "\nStart test" << endl;

  Polygon_2 contour;
  Point_3 pointComputeField;
  double radius, angle, area, length, inRadius, inAngle, sc, r;
  int nbAngles(100);
  vector<int> surfaceIdx(nbAngles, 0), slectedModesIdx;
  double scalingFactors[2] = { 1., 1. };
  double a(5.5), b(3.);
  double freq, freqMax, freqField, endAdmit;
  complex<double> result;
  int nbFreqs, mn, idxField, cnt;
  Eigen::MatrixXcd characImped, radImped, radAdmit, inputVelocity, inputPressure;
  vector<Point_3> radPts;
  vector<array<double, 2>> pts;
  vector<array<int, 3>> triangles;
  Eigen::VectorXcd radPress;
  vector<Eigen::MatrixXcd> Q;
  vector<Matrix> KR2;
  // to determine the rectangle mode indexes
  vector<array<int, 2>> modeIdxs;
  vector<int> vIdx;
  vector<double> k2, rads, shifts, scaleIn, scaleOut, lengths, curvAngles;
  int nCombinations, ie, je, me, ne;
  Matrix analyticE, E, modes;
  double E1y, E1z, E2y, E2z;
  bool radiationCondition;

  auto start = std::chrono::system_clock::now();
  auto stop = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsedTime;
  int hours;
  int minutes;
  double seconds;

  switch(tType){
  case MATRIX_E:
    /////////////////////////////////////////////////////////////////////////////////
    //*********************************************
    // Matrix E 
    //*********************************************

    log << "Start matrix E test" << endl;

    // create contour
    contour.push_back(Point(0., 0.));
    contour.push_back(Point(a, 0.));
    contour.push_back(Point(a, b));
    contour.push_back(Point(0., b));
    m_maxCSBoundingBox.first = Point2D(-1.2*a, -1.2*a);
    m_maxCSBoundingBox.second = Point2D(1.2*a, 1.2*a);
    m_geometryImported = true; // to have the good bounding box for modes plot

    m_crossSections.clear();
    area = a*b;
    length = 20.;
    surfaceIdx.resize(4);
    addCrossSectionFEM(area, sqrt(area) / m_meshDensity, contour,
      surfaceIdx, length, Point2D(0., 0.), Point2D(0., 1.),
      scalingFactors);

    m_crossSections[0]->buildMesh();

    //*******************
    // export mesh
    //*******************

    // export points
    strs << "points_rec_" << m_meshDensity << ".txt";
    ofs.open(strs.str());
    pts = m_crossSections[0]->getPoints();
    for (auto it : pts)
    {
      ofs << it[0] << "  " << it[1] << endl;
    }
    ofs.close();

    // export triangles
    strs.str("");
    strs << "triangles_rec_" << m_meshDensity << ".txt";
    ofs.open(strs.str());
    triangles = m_crossSections[0]->getTriangles();
    for (auto it : triangles)
    {
      for (int i(0); i < 3; i++) { ofs << it[i] + 1 << "  "; }
      ofs << endl;
    }
    ofs.close();

    m_crossSections[0]->computeModes(m_simuParams);

    //***************************************
    // Export modes and propagation matrices
    //***************************************

    modes = m_crossSections[0]->getModes();
    strs.str("");
    strs << "modes_rec_" << m_meshDensity << ".txt";
    ofs.open(strs.str());
    ofs << modes << endl;
    ofs.close();

    strs.str("");
    strs << "eigen_freqs_rec_" << m_meshDensity << ".txt";
    ofs.open(strs.str());
    for (int i(0); i < 100; i++)
    {
      ofs << m_crossSections[0]->eigenFrequency(i) << endl;
    }
    ofs.close();

    modes = m_crossSections[0]->getMatrixC();
    strs.str("");
    strs << "C_rec_" << m_meshDensity << ".txt";
    ofs.open(strs.str());
    ofs << modes << endl;
    ofs.close();

    modes = m_crossSections[0]->getMatrixD();
    strs.str("");
    strs << "D_rec_" << m_meshDensity << ".txt";
    ofs.open(strs.str());
    ofs << modes << endl;
    ofs.close();

    modes = m_crossSections[0]->getMatrixE();
    strs.str("");
    strs << "E_rec_" << m_meshDensity << ".txt";
    ofs.open(strs.str());
    ofs << modes << endl;
    ofs.close();

    modes = m_crossSections[0]->getMatrixKR2(0);
    strs.str("");
    strs << "KR2_rec_" << m_meshDensity << ".txt";
    ofs.open(strs.str());
    ofs << modes << endl;
    ofs.close();

    //**********************************************
    // Compute matrix E from analytical expression
    //**********************************************

    // generate modes indexes
    mn = m_crossSections[0]->numberOfModes();
    nCombinations = 10000;
    modeIdxs.reserve(nCombinations);
    k2.reserve(nCombinations);
    vIdx.resize(nCombinations);
    iota(vIdx.begin(), vIdx.end(), 0);
    for (int m(0); m < 100; m++)
    {
      for (int n(0); n < 100; n++)
      {
        k2.push_back(pow((double)m * b / a, 2) + pow((double)n, 2));
        modeIdxs.push_back({ m,n });
      }
    }
    sort(vIdx.begin(), vIdx.end(), [&](int i, int j){return k2[i]<k2[j];});

/*    for (int i(0); i < 20; i ++)
    {
      log << k2[vIdx[i]] << "  " << modeIdxs[vIdx[i]][0] << "  " 
          << modeIdxs[vIdx[i]][1] << endl;
    }*/

    analyticE.setZero(mn, mn);
    for (int m(0); m < mn; m++)
    {
      for (int n(0); n < mn; n++)
      {
        // extract the indexes of the modes
        ie = modeIdxs[vIdx[m]][0];
        je = modeIdxs[vIdx[m]][1];
        me = modeIdxs[vIdx[n]][0];
        ne = modeIdxs[vIdx[n]][1];

        // compute E1y
        if (me == 0) {E1y = 0.;}
        if ((ie == 0) && (me != 0)) {E1y = sqrt(2.)*cos((double)me*M_PI);}
        if ((ie == me) && (me != 0)) {E1y = 0.5;}
        if ((ie != me) && (ie != 0) && (me != 0)){E1y = (double)me*(
          cos((double)(ie + me)*M_PI)/(double)(ie + me) -
          cos((double)(ie - me)*M_PI)/(double)(ie - me));}

        // compute E1z
        if (je == ne) {E1z = 1.;} else {E1z = 0.;}

        // compute E2y
        if (ie == me){E2y = 1.;} else {E2y = 0.;}

        // compute E2z
        if (ne == 0) {E2z = 0.;}
        if ((je == 0) && (ne != 0)) {E2z = sqrt(2.)*cos((double)(ne)*M_PI);}
        if ((je == ne) && (ne != 0)) {E2z = 0.5;}
        if ((je != ne) && (je != 0) && (ne != 0)) {E2z = (double)ne*(
          cos((double)(je + ne)*M_PI)/(double)(je + ne) -
          cos((double)(je - ne)*M_PI)/(double)(je - ne));} 
        
        // Compute E(m,n)
        analyticE(m,n) = E1y*E1z + E2y*E2z;
      }
    }

    ofs.open("anE.txt");
    ofs << analyticE << endl;
    ofs.close();
    ofs.open("nuE.txt");
    ofs << m_crossSections[0]->getMatrixE() << endl;
    ofs.close();

    break;
    /////////////////////////////////////////////////////////////////////////////////
  case DISCONTINUITY:

    //*********************************************
    // Discontinuity
    //*********************************************
  
    // Set the proper simulation parameters
    m_simuParams.numIntegrationStep = 165;
    m_meshDensity = 20.;
    m_simuParams.percentageLosses = 0.;
    m_simuParams.wallLosses = false;
    m_simuParams.curved = false;
    m_maxCutOnFreq = 30000.;
    m_geometryImported = true; // to have the good bounding box for modes plot
    m_simuParams.sndSpeed = 34400;

    generateLogFileHeader(true);
  
    // Generate a circular contour
    radius = 4.;
    m_maxCSBoundingBox.first = Point2D(-radius, -radius);
    m_maxCSBoundingBox.second = Point2D(radius, radius);
    for (int i(0); i < nbAngles; i++)
    {
      angle = 2.*M_PI*(double)(i) / (double)(nbAngles);
      contour.push_back(Point(radius * cos(angle), radius * sin(angle)));
    }
  
    //// Check contour
    //ofs.open("cont.txt");
    //for (auto it : contour) { ofs << it.x() << "  " << it.y() << endl; }
    //ofs.close();
  
    m_crossSections.clear();
    area = M_PI*pow(radius, 2);
    length = 30.;
    addCrossSectionFEM(area, sqrt(area) / m_meshDensity, contour,
      surfaceIdx, length, Point2D(0., 0.), Point2D(0., 1.),
      scalingFactors);
    m_crossSections[0]->setAreaVariationProfileType(GAUSSIAN);
  
    log << "Cross-section created" << endl;
  
    // Check the scaling factor
    ofs.open("sc.txt");
    for (int i(0); i < nbAngles; i++)
    {
      ofs << m_crossSections[0]->scaling((double)(i) / (double)(nbAngles - 1))
        << endl;
    }
    ofs.close();
  
    // Check the scaling factor derivative
    ofs.open("dsc.txt");
    for (int i(0); i < nbAngles; i++)
    {
      ofs << m_crossSections[0]->scalingDerivative((double)(i) / (double)(nbAngles - 1))
        << endl;
    }
    ofs.close();
  
    log << "Parameters set" << endl;
  
    // compute propagation modes
    m_crossSections[0]->buildMesh();
  
    log << "Mesh generated" << endl;
    //// extract mesh
    //CDT tr(m_crossSections[0]->triangulation());
    //ofs.open("mesh.txt");
    //for (auto it = tr.all_faces_begin(); it != tr.all_faces_end(); it++)
    //{
    //  for (int i(0); i < 3; i++)
    //  {
    //    ofs << it->vertex(i)->point().x() << "  "
    //  {
    //    ofs << it->vertex(i)->point().x() << "  "
    //      << it->vertex(i)->point().y() << endl;
    //  }
    //  ofs << "nan  nan" << endl;
    //}
    //ofs.close();
  
    m_crossSections[0]->computeModes(m_simuParams);
    log << m_crossSections[0]->numberOfModes() << " modes computed" << endl;

    slectedModesIdx.push_back(0);
    slectedModesIdx.push_back(5);
    slectedModesIdx.push_back(16);
    slectedModesIdx.push_back(31);
    slectedModesIdx.push_back(52);
    slectedModesIdx.push_back(106);

    m_crossSections[0]->selectModes(slectedModesIdx);

    // Extract value of the 1st mode
    E = m_crossSections[0]->getModes();
    log << "1st mode: " << E(0, 0) << endl;

    // modify matrix E
    /*E = m_crossSections[0]->getMatrixE();
    E(0, 5) = -E(0, 5);
    m_crossSections[0]->setMatrixE(E);*/

    //// Export matrix E
    //ofs.open("mE.txt");
    //ofs << m_crossSections[0]->getMatrixE() << endl;
    //ofs.close();

    preComputeRadiationMatrices(16, 0);
  
    freqMax = 2500.;
    nbFreqs = 1500;
    ofs.open("imp.txt");
    ofs2.open("freqs.txt");
    ofs3.open("rad.txt");
    for (int i(0); i < nbFreqs; i++)
    {
      // get the output impedance
      freq = max(0.1, freqMax*(double)i/(double)(nbFreqs - 1));
      log << "f = " << freq << " Hz" << endl;
      
      //m_crossSections[0]->characteristicImpedance(characImped, freq, m_simuParams);
      interpolateRadiationImpedance(radImped, freq, 0);
      ofs3 << radImped.imag() << endl;
  
      // Propagate impedance
      m_crossSections[0]->propagateMagnus(radImped, m_simuParams, freq, -1., IMPEDANCE);
  
      ofs2 << freq << endl;
      ofs << m_crossSections[0]->Zin().cwiseAbs() << endl;
      //  //-1i* 2. * M_PI * freq * m_simuParams.volumicMass*
      //  m_crossSections[0]->Zin()(0, 0)
      //  /// characImped(0,0)
      //) << "  "
      //  << arg(
      //    //-1i* 2. * M_PI * freq * m_simuParams.volumicMass*
      //    m_crossSections[0]->Zin()(0, 0)
      //    /// characImped(0, 0)
      //  ) << endl;

    }
    ofs3.close();
    ofs2.close();
    ofs.close();
    
    break;
    /////////////////////////////////////////////////////////////////////////////////
    case ELEPHANT_TRUNK:
    
    //*********************************************
    // Elephant trunk 
    //*********************************************

      start = std::chrono::system_clock::now();

      radiationCondition = true;
    
    // Set the proper simulation parameters
    m_simuParams.freqDepLosses = false;
//    m_simuParams.wallLosses = false;
      //m_simuParams.orderMagnusScheme = 4;
    //m_meshDensity = 15.;
    //m_maxCutOnFreq = 20000;
    m_simuParams.numIntegrationStep = 50;
    m_simuParams.maxComputedFreq = 10000;
    m_simuParams.curved = true;
    m_geometryImported = true; // to have the good bounding box for modes plot

    generateLogFileHeader(true);

    // Generate a circular contour
    radius = 3.;
    m_maxCSBoundingBox.first = Point2D(-2.*radius, -2.*radius);
    m_maxCSBoundingBox.second = Point2D(2.*radius, 2.*radius);
    for (int i(0); i < nbAngles; i++)
    {
      angle = 2.*M_PI*(double)(i) / (double)(nbAngles);
      contour.push_back(Point(radius * cos(angle), radius * sin(angle)));
    }

    // export contour
    ofs.open("cont.txt");
    for (auto pt : contour)
    {
      ofs << pt.x() << "  " << pt.y() << endl;
    }
    ofs.close();

    // create the cross-section
    m_crossSections.clear();
    area = pow(radius, 2) * M_PI;
    scalingFactors[0] = 0.25;
    scalingFactors[1] = 1.;
    inRadius = 1.25*4.*1.5;
    log << "inRadius " << inRadius << endl;
    inAngle = 2.26;
    length = inAngle * inRadius;
    log << "length: " << length << endl;
    addCrossSectionFEM(area, sqrt(area) / m_meshDensity, contour,
      surfaceIdx, length, Point2D(0., 0.), Point2D(-1., 0.),
      scalingFactors);
    m_crossSections[0]->setAreaVariationProfileType(ELEPHANT);
    m_crossSections[0]->setCurvatureRadius(-inRadius);
    m_crossSections[0]->setCurvatureAngle(inAngle);
    mn = 6;
    m_crossSections[0]->setModesNumber(mn);
  
    log << "Cross-section created" << endl;

    // test get coordinates
    //ofs.open("pt3.txt");
    //for (int i(0); i < 100; i++)
    //{
      //for (int j(0); j < 100; j++)
      //{
        //angle = inAngle * (double)i/99.;
        //r = inRadius + m_crossSections[0]->scaling(inRadius * angle / length)
          //* 2. * radius * ((double)j/99. - 0.5);
        //log << "angle = " << angle << endl;
        //log << "r = " << r << " inRadius " << inRadius << " scaling " 
          //<< m_crossSections[0]->scaling(inRadius * angle / length)  
          //<< " radius " << radius << " val " << ((double)j/99. - 0.5) 
          //<< " z " << m_crossSections[0]->scaling(inRadius * angle / length)
          //* ((double)j/99. - 0.5) << endl;
        //pointComputeField = m_crossSections[0]->getCoordinateFromCartesianPt(Point_3(
              //-r * cos(angle) + inRadius, 0., r * sin(angle)));
        //ofs << pointComputeField << "  " << -r * cos(angle) + inRadius
          //<< "  " << r * sin(angle) << endl;

         //m_crossSections[0]->getCoordinateFromCartesianPt(Point_3(
              //28.*(double)i/99. - 14., 0., 28.*(double)j/99. - 14.), pointComputeField);
        //ofs << pointComputeField << "  " << 28.*(double)i/99. - 14.
          //<< "  " << 28.*(double)j/99. - 14. << endl;
      //}
    //}
    //ofs.close();

    // Check the scaling factor
    ofs.open("sc.txt");
    for (int i(0); i < nbAngles; i++)
    {
      ofs << m_crossSections[0]->scaling((double)(i) / (double)(nbAngles - 1))
        << endl;
    }
    ofs.close();
  
    // Check the scaling factor derivative
    ofs.open("dsc.txt");
    for (int i(0); i < nbAngles; i++)
    {
      ofs << m_crossSections[0]->scalingDerivative((double)(i) / (double)(nbAngles - 1))
        << endl;
    }
    ofs.close();

    // compute propagation modes
    m_crossSections[0]->buildMesh();
    m_crossSections[0]->computeModes(m_simuParams);
    mn = m_crossSections[0]->numberOfModes();
    log << mn << " modes computed" << endl;

    // initialize input pressure and velocity vectors
    inputPressure.setZero(mn, 1);
    inputVelocity.setZero(mn, 1);

    // Open end boundary consition
    if (radiationCondition)
    { 
      preComputeRadiationMatrices(16, 0); 
    }
    else
    {
      radAdmit.setZero(mn, mn);
      radAdmit.diagonal() = Eigen::VectorXcd::Constant(mn, complex<double>(1000000000., 0.));
      radImped = radAdmit.inverse();
    }

    // define the position of the radiation point
    radPts.push_back(Point_3(3., 0., 0.));

    freqMax = m_simuParams.maxComputedFreq;
    ofs.open("elephant_ac_press_MM.txt");
    m_numFreq = 101;
    for (int i(0); i < m_numFreq; i++)
    {
      freq = max(0.1, freqMax*(double)i/(double)(m_numFreq - 1));
      log << "f = " << freq << " Hz" << endl;

      if (radiationCondition) { interpolateRadiationImpedance(radImped, freq, 0); }

      ofs2 << characImped.real() << endl;
      ofs3 << characImped.imag() << endl;

      //log << "Radiation impedance interpolated" << endl;
  
      if (radiationCondition) {
        // Propagate impedance
        m_crossSections[0]->propagateMagnus(radImped, m_simuParams, freq, -1., IMPEDANCE);
      }
      else
      {
        // Propagate admittance
        m_crossSections[0]->propagateMagnus(radAdmit, m_simuParams, freq, -1., ADMITTANCE);
        //log << "Admittance propagated" << endl;
      }

      // propagate velocity or pressure
      inputVelocity(0, 0) = -1i * 2. * M_PI * freq * m_simuParams.volumicMass
       //* m_crossSections[0]->area() * pow(m_crossSections[0]->scaling(0.), 2)
        ;

      // export result
      ofs << freq << "  ";
      if (radiationCondition)
      {
        m_crossSections[0]->propagateMagnus(inputVelocity, m_simuParams, freq, 1., VELOCITY);
        // compute radiated pressure
        RayleighSommerfeldIntegral(radPts, radPress, freq, 0);
        spectrum.setValue(i, radPress(0,0));
        ofs << abs(1e5*radPress(0)/2./M_PI) << "  "
          << arg(1e5*radPress(0)/2./M_PI) << "  " << endl;
      }
      else
      {
        inputPressure = m_crossSections[0]->Yin().inverse() * inputVelocity;
        m_crossSections[0]->propagateMagnus(inputPressure, m_simuParams, freq, 1., PRESSURE);
        //log << "Pressure propagated" << endl;
        result = m_crossSections[0]->area() * pow(m_crossSections[0]->scaleIn(), 2) * 1e5 *
            m_crossSections[0]->interiorField(
          Point_3(m_crossSections[0]->length() - 0.03, 0., 0.), m_simuParams, PRESSURE);
        ofs << abs(result) << "  "
          << arg(result) << "  " << endl;
      }
      
      // extract acoustic field
      cnt = 0;
      if (freq == 2400.)
      {
        ofs2.open("p_field.txt");
        Eigen::MatrixXcd field;
        acousticFieldInPlane(field);
        stringstream txtField;
        txtField << field.cwiseAbs();
        ofs2 << regex_replace(txtField.str(), regex("-nan\\(ind\\)"), "nan");

        //for (int i(0); i < 500; i++)
        //{
        //  for (int j(0); j < 500; j++)
        //  {
        //    if (m_crossSections[0]->getCoordinateFromCartesianPt(Point_3(
        //         15.*(double)i/499. - 1., 0., 15.*(double)j/499. - 1), pointComputeField))
        //    {
        //      result = m_crossSections[0]->interiorField(pointComputeField, m_simuParams, PRESSURE);
        //      cnt ++;
        //      log << "Pressure computed " << cnt << " over " << 250000 << endl;
        //      ofs2 << abs(result) << "  ";
        //    }
        //    else
        //    {
        //      ofs2 << "nan  ";
        //    }
        //  }
        //  ofs2 << endl;
        //}
        ofs2.close();
      }
    }
    ofs.close();

    //// extract the amplitude of the acoustic pressure along the guide
    //ofs.open("amp.txt");
    //Q = m_crossSections[0]->Z();
    //for (int i(0); i < Q.size(); i++)
    //{
    //  sc = m_crossSections[0]->scaling((double)i / (double)(Q.size() - 1));
    //  ofs << abs(Q[i](0, 0)) << "  " << sc << endl;
    //}

    // get total time of the simulation
    stop = std::chrono::system_clock::now();
    elapsedTime = stop - start;
    hours = floor(elapsedTime.count() / 3600.);
    minutes = floor((elapsedTime.count() - hours * 3600.) / 60.);
    seconds = elapsedTime.count() - (double)hours * 3600. -
      (double)minutes * 60.;
    log << "\nTotal time "
      << hours << " h " << minutes << " m " << seconds << " s" << endl;
    ofs.close();

    //// generate spectra values for negative frequencies
    //for (int i(m_numFreq); i < 2 * m_numFreq; i++)
    //{
    //  spectrum.re[i] = spectrum.re[2 * m_numFreq - i - 1];
    //  spectrum.im[i] = -spectrum.im[2 * m_numFreq - i - 1];
    //}

    break;
    /////////////////////////////////////////////////////////////////////////////////
    case SCALE_RAD_IMP:

    //*****************************************************
    // Test radiation impedance computation with scaling
    //*****************************************************

    // Generate reference radiation impedance
    //***************************************


      // Set the proper simulation parameters
      m_geometryImported = true; // to have the good bounding box for modes plot
      //m_maxCutOnFreq = 1.;

      generateLogFileHeader(true);
      log << "Start test scale rad imped" << endl;

      // Generate a circular contour
      radius = 3.;
      m_maxCSBoundingBox.first = Point2D(-2. * radius, -2. * radius);
      m_maxCSBoundingBox.second = Point2D(2. * radius, 2. * radius);
      for (int i(0); i < nbAngles; i++)
      {
        angle = 2. * M_PI * (double)(i) / (double)(nbAngles);
        contour.push_back(Point(radius * cos(angle), radius * sin(angle)));
      }

      // creeate cross section
      m_crossSections.clear();
      area = pow(radius, 2) * M_PI;
      scalingFactors[0] = 1.;
      scalingFactors[1] = 1.;
      length = 1.;
      addCrossSectionFEM(area, sqrt(area) / m_meshDensity, contour,
        surfaceIdx, length, Point2D(0., 0.), Point2D(0., 1.),
        scalingFactors);
      m_crossSections[0]->setCurvatureRadius(1.);
      m_crossSections[0]->setCurvatureAngle(1.);

      log << "Cross-section created" << endl;

      // compute propagation modes
      m_crossSections[0]->buildMesh();
      m_crossSections[0]->computeModes(m_simuParams);
      mn = m_crossSections[0]->numberOfModes();
      log << mn << " modes computed" << endl;

      // compute radiation impedance
      preComputeRadiationMatrices(16, 0);

      log << "Precomputation of rad imped done" << endl;
      ofs.open("imp.txt");
      for (int i(0); i < m_numFreq; i++)
      {
        freq = max(0.1, m_simuParams.maxComputedFreq * 
          (double)i / (double)(m_numFreq - 1));
        log << "f = " << freq << " Hz" << endl;

        interpolateRadiationImpedance(radImped, freq, 0);

        ofs << freq << "  " << radImped(0, 0).real() << "  "
          << radImped(0, 0).imag() << "  "
          << radImped(1, 1).real() << "  "
          << radImped(1, 1).imag() << endl;
      }
      ofs.close();

      // Generate radiation impedance computed with scaling
      //***************************************************

      // Generate a circular contour
      contour.clear();
      radius = 1.5;
      m_maxCSBoundingBox.first = Point2D(-2. * radius, -2. * radius);
      m_maxCSBoundingBox.second = Point2D(2. * radius, 2. * radius);
      for (int i(0); i < nbAngles; i++)
      {
        angle = 2. * M_PI * (double)(i) / (double)(nbAngles);
        contour.push_back(Point(radius * cos(angle), radius * sin(angle)));
      }

      // creeate cross section
      m_crossSections.clear();
      area = pow(radius, 2) * M_PI;
      scalingFactors[0] = 1.;
      scalingFactors[1] = 2.;
      length = 1.;
      addCrossSectionFEM(area, sqrt(area) / m_meshDensity, contour,
        surfaceIdx, length, Point2D(0., 0.), Point2D(0., 1.),
        scalingFactors);
      m_crossSections[0]->setCurvatureRadius(1.);
      m_crossSections[0]->setCurvatureAngle(1.);

      log << "Cross-section created" << endl;

      // compute propagation modes
      m_crossSections[0]->buildMesh();
      m_crossSections[0]->computeModes(m_simuParams);
      mn = m_crossSections[0]->numberOfModes();
      log << mn << " modes computed" << endl;

      // compute radiation impedance
      preComputeRadiationMatrices(16, 0);

      log << "Precomputation of rad imped done" << endl;
      ofs.open("impS.txt");
      for (int i(0); i < m_numFreq; i++)
      {
        freq = max(0.1, m_simuParams.maxComputedFreq *
          (double)i / (double)(m_numFreq - 1));
        log << "f = " << freq << " Hz" << endl;

        interpolateRadiationImpedance(radImped, freq, 0);

        ofs << freq << "  " << radImped(0, 0).real() << "  "
          << radImped(0, 0).imag() << "  " 
          << radImped(1, 1).real() << "  "
          << radImped(1, 1).imag() << endl;
      }
      ofs.close();

    break;
  }

  log.close();
}

// ****************************************************************************
// Compute the Rayleigh-Sommerfeld integral to compute radiated pressure

void Acoustic3dSimulation::RayleighSommerfeldIntegral(vector<Point_3> points,
  Eigen::VectorXcd& radPress, double freq, int radSecIdx)
{
  int nbPts(points.size());
  double quadPtCoord[3][2]{ {1. / 6., 1. / 6.}, {2. / 3., 1. / 6.}, {1. / 6., 2. / 3.} };
  double quadPtWeight = 1. / 3.;
  vector<Point> gaussPts;
  vector<double> areaFaces;
  double r, k(2.*M_PI*freq/m_simuParams.sndSpeed), scaling;

  //ofstream log;
  //log.open("log.txt", ofstream::app);

  // get scaling
  scaling = m_crossSections[radSecIdx]->scaleOut();

  // get velocity mode amplitude (v_x = j * q / w / rho)
  auto Vm = m_crossSections[radSecIdx]->Qout();

  // get quadrature points
  gaussPointsFromMesh(gaussPts, areaFaces, 
    m_crossSections[radSecIdx]->triangulation());

  // interpolate modes
  Matrix interpolatedModes = m_crossSections[radSecIdx]->
    interpolateModes(gaussPts);

  // scale the position of the integration points
  Transformation3 scale(CGAL::SCALING, 1./scaling);
  for (int i(0); i < points.size(); i++)
  {
    points[i] = scale(points[i]);
  }

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
            Vm(m) * interpolatedModes(f * 3 + g, m) * exp(-1i * k * scaling * r) / r;
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
  //  P1      centerline point of section 1
  //  P2      centerline point of section 2
  //  N1      centerline normal of section 1
  //  N2      centerline normal of section 2
  //
  // outputs:
  //
  //  radius    radius of circle arc joining sec 1 & 2
  //  angle    angle between sections 1 & 2
  //  shift    distance necessary to shift P2 along N2
  //        so that it is on the same circle arc as P1

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
  //radius = (radii[0] + radii[1]) / 2.;
  radius = radii[0];

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
  vector<pair<double, double>>& scalingFactors, bool simplifyContours)
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
  Cost cost;                // for contour simplification
  int idxCont(0);

  ofstream log("log.txt", ofstream::app);
  //log << "Start geometry importation" << endl;

  ifstream geoFile(m_geometryFile);

  // initialize maximal bounding box
  m_maxCSBoundingBox.first = Point2D(0., 0.);
  m_maxCSBoundingBox.second = Point2D(0., 0.);

  // check if the file is opened
  if (!geoFile.is_open())
  {
    log << "Cannot open " << m_geometryFile << endl;
    return false;
  }
  else
  {

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

      // extract the normal to the centerline
      getline(lineX, coordX, separator);
      getline(lineY, coordY, separator);
      normals.push_back(Point2D(stod(coordX), stod(coordY)));

      // extract the scaling factors
      getline(lineX, coordX, separator);
      getline(lineY, coordY, separator);
      scalingFactors.push_back(pair<double, double>(stod(coordX), stod(coordY)));

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
  log.close();
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
  vector<pair<double, double>> vecScalingFactors;
  vector<double> totAreas;
  vector<array<double, 4>> bboxes;
  array<double, 4> arrayZeros = { 0., 0., 0., 0. };

  //ofstream log("log.txt", ofstream::app);
  //log << "Start cross-section creation" << endl;

  if (m_geometryImported)
  {
    if ( !extractContoursFromCsvFile(contours, surfaceIdx, centerLine, 
      normals, vecScalingFactors, false))
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
  //ofstream ar("area.txt");
  //for (auto it : totAreas) { ar << it << endl; }
  //ar.close();

  //*********************************************************
  // Create lambda expression to compute the scaling factors
  //*********************************************************

  auto getScalingFactor = [bboxes](int idx1, int idx2)
  {
    double meanX((abs(bboxes[idx1][0]) + abs(bboxes[idx1][1])
      + abs(bboxes[idx2][0]) + abs(bboxes[idx2][1])));
    double meanY((abs(bboxes[idx1][2]) + abs(bboxes[idx1][3])
      + abs(bboxes[idx2][2]) + abs(bboxes[idx2][3])));

    // a factor 0.999 is added to avoid that points of successive
    // contours end being exactly the same which cause 
    // CGAL::intersection to seg fault
    if (meanX > meanY)
    {
      return 0.999*min(bboxes[idx2][0] / bboxes[idx1][0],
        bboxes[idx2][1] / bboxes[idx1][1]);
    }
    else
    {
      return 0.999*min(bboxes[idx2][2] / bboxes[idx1][2],
        bboxes[idx2][3] / bboxes[idx1][3]);
    }
  };

  // variables for cross-section creation
  double prevCurvRadius, curvRadius, prevAngle, angle, shift, area, length, radius;
  double prevScalingFactors[2] = { 1., 1. };
  double scalingFactors[2] = { 1., 1. };
  double array1[2] = { 1., 1. };
  //const double minDist(1.e-4);
  const double MINIMAL_AREA(0.15);
  vector<int> tmpPrevSection, prevSecInt, listNextCont, tmpSurf;
  vector<vector<int>> prevSections, intSurfacesIdx;
  int secIdx(0), intSecIdx(0), nextSecIdx, nbCont(contours.size());
  vector<Polygon_2> intContours;
  Pwh_list_2 intersections;
  bool sidePrev, side;

  //**********************************************************************
  // Add an intermediate centerline point and normal before the last ones
  //**********************************************************************

  centerLine.push_back(centerLine.back());
  normals.push_back(normals.back());
  int lastCtl(centerLine.size() - 1);
  //log << "Before last ctl " << centerLine[lastCtl - 2].x << "  "
  //  << centerLine[lastCtl - 2].y << " normal "
  //  << normals[lastCtl - 2].x << "  "
  //  << normals[lastCtl - 2].y << endl;
  //log << "Last ctl " << centerLine[lastCtl].x << "  "
  //  << centerLine[lastCtl].y << " normal "
  //  << normals[lastCtl].x << "  "
  //  << normals[lastCtl].y << endl;

  getCurvatureAngleShift(centerLine[lastCtl - 2], centerLine[lastCtl], 
    normals[lastCtl - 2], normals[lastCtl], curvRadius, angle, shift);

  //log << "Curv radius " << curvRadius << endl;

  Point pt(Point(centerLine.back().x, centerLine.back().y));
  Vector N(Vector(normals.back().x, normals.back().y));

  if (angle > MINIMAL_DISTANCE)
  {
    angle /= -4.;
    Transformation rotate(CGAL::ROTATION, sin(angle - M_PI / 2.),
      cos(angle - M_PI / 2.));
    Transformation translateInit(CGAL::TRANSLATION,
      2. * abs(curvRadius) * sin(angle) * rotate(N));

    pt = translateInit(pt);
    centerLine[lastCtl - 1].x = pt.x();
    centerLine[lastCtl - 1].y = pt.y();

    if (curvRadius > 0.)
    {
      angle *= 2.;
    }
    else
    {
      angle *= -2.;
    }
    Transformation rotateN(CGAL::ROTATION, sin(angle), cos(angle));
    N = rotateN(N);
    normals[lastCtl - 1].x = N.x();
    normals[lastCtl - 1].y = N.y();
  }
  else
  {
    Vector tr((centerLine[lastCtl - 2].x - centerLine[lastCtl].x)/2.,
      (centerLine[lastCtl - 2].y - centerLine[lastCtl].y)/2.);
    //log << "Vec tr " << tr << endl;
    Transformation translateInit(CGAL::TRANSLATION, tr);

    pt = translateInit(pt);
    centerLine[lastCtl - 1].x = pt.x();
    centerLine[lastCtl - 1].y = pt.y();
  }

  //log << "Trans ctl " << pt.x() << "  " << pt.y() << " normal "
  //  << N.x() << "  " << N.y() << endl;

  //*******************************************
  // Create the cross-sections
  //*******************************************

  // clear cross-sections before 
  m_crossSections.clear();

  // compute the curvatur parameters of the first cross-section
  getCurvatureAngleShift(centerLine[0], centerLine[1], normals[0], normals[1],
    prevCurvRadius, prevAngle, shift);

  //// shift the contours
  //Transformation translate(CGAL::TRANSLATION, Vector(0., shift));
  //for (auto cont : contours[0]){cont = transform(translate, cont);}
  // FEXME: one should shift also the ccenterline point probably

  // initialize the previous section index list
  for (int c(0); c < contours[0].size(); c++){prevSections.push_back(tmpPrevSection);}

  // initialise the scaling factors
  if (m_simuParams.varyingArea)
  {
    switch (m_contInterpMeth)
    {
    case AREA:
      prevScalingFactors[0] = 1.;
      prevScalingFactors[1] = sqrt(totAreas[1] / totAreas[0]);
      break;
    case BOUNDING_BOX:
      prevScalingFactors[0] = 1.;
      prevScalingFactors[1] = getScalingFactor(0, 1);
      break;
    case FROM_FILE:
      prevScalingFactors[0] = vecScalingFactors[0].first;
      prevScalingFactors[1] = vecScalingFactors[0].second;
      break;
    }
  }

  // create the cross-sections
  for (int i(1); i < nbCont; i++)
  {
    //log << "\ni= " << i << endl;

    //**********************************
    // Create previous cross-sections
    //**********************************

      length = centerLine[i-1].getDistanceFrom(centerLine[i]);
      //log << "length " << length << endl;

    // compute the scaling factors
    if (m_simuParams.varyingArea)
    {
      if (i < (nbCont - 2))
      {
        switch (m_contInterpMeth)
        {
        case AREA:
          scalingFactors[0] = 1.;
          scalingFactors[1] =  sqrt(totAreas[i + 1] / totAreas[i]);
          break;
        case BOUNDING_BOX:
          scalingFactors[0] = 1.;
          scalingFactors[1] = getScalingFactor(i, i + 1);
          break;
        case FROM_FILE:
          scalingFactors[0] = vecScalingFactors[i].first;
          scalingFactors[1] = vecScalingFactors[i].second;
          break;
        }
      }
      else if (i == nbCont - 2)
      {
        switch (m_contInterpMeth)
        {
        case AREA: 
          scalingFactors[0] = 1.;
          scalingFactors[1] = sqrt((totAreas[i] + totAreas[i + 1]) / 2. / totAreas[i]);
          break;
        case BOUNDING_BOX:
          scalingFactors[0] = 1.;
          scalingFactors[1] = getScalingFactor(i, i + 1) / 2.;
          break;
        case FROM_FILE:
          scalingFactors[0] = vecScalingFactors[i].first;
          scalingFactors[1] = vecScalingFactors[i].second;
          break;
        }
      }
      else if (i == nbCont - 1)
      {
        switch (m_contInterpMeth)
        {
        case AREA:
          scalingFactors[0] = sqrt((totAreas[i-1] + totAreas[i]) / 2. / totAreas[i]);
          scalingFactors[1] = 1.;
          break;
        case BOUNDING_BOX:
          scalingFactors[1] = getScalingFactor(i-1, i) / 2.;
          scalingFactors[1] = 1.;
          break;
        case FROM_FILE:
          scalingFactors[0] = vecScalingFactors[i].first;
          scalingFactors[1] = vecScalingFactors[i].second;
        }
      }
    }
    //log << "SCaling factors " << scalingFactors[0] << "  " << scalingFactors[1] << endl;

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
      //log << "Curv radius " << prevCurvRadius << endl;

      // set the curvature angle
      m_crossSections[secIdx]->setCurvatureAngle(prevAngle);
      //log << "angle " << prevAngle << endl;

      //log << "Section " << secIdx << " created" << endl;
      secIdx++;
    }

    //******************************************
    // extract current cross-section parameters
    //******************************************

    // compute the curvatur parameters of the section
    getCurvatureAngleShift(centerLine[i], centerLine[i+1], normals[i], normals[i+1],
      curvRadius, angle, shift);

    //// shift the contours
    //Transformation translate(CGAL::TRANSLATION, Vector(0., shift));
    //for (auto cont : contours[i]) {cont = transform(translate, cont);}

    //log << "Contour shifted" << endl;

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

    //log << "Minimal area checked" << endl;

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

      //log << "Contour extracted and scaled sc = " 
      //  << scalingFactors[0] << endl;

      // loop over the contours of the previous cross-section
      for (int cp(0); cp < contours[i-1].size(); cp++)
      {
        // extract the previous contour and scale it
        Transformation scale(CGAL::SCALING, prevScalingFactors[1]);
        Polygon_2 prevCont(transform(scale, contours[i - 1][cp]));

        //log << "Prev Contour extracted and scaled sc = " 
        //  << prevScalingFactors[1] << endl;

        if (!similarContours(cont, prevCont, MINIMAL_DISTANCE_DIFF_POLYGONS))
        {

          // Check if the current and previous contour intersect:
          // check if the first point of the previous contour
          // is on the bounded side of the current contour
          auto itP = prevCont.begin();
          sidePrev = cont.has_on_bounded_side(*itP);

          //log << "Sideprev " << sidePrev << endl;

          // loop over the points of the previous contour
          for (; itP != prevCont.end(); itP++)
          {
            side = cont.has_on_bounded_side(*itP);
            // log << "Side " << side << endl;
            // if the previous point and the next point are on
            // different sides
            if (side != sidePrev)
            {
              //log << "side != sidePrev" << endl;

              //ofstream os("cont.txt");
              //for (auto pt : cont) { os << pt.x() << "  " << pt.y() << endl; }
              //os.close();
              //os.open("pcont.txt");
              //for (auto pt : prevCont) { os << pt.x() << "  " << pt.y() << endl; }
              //os.close();

              //log << "Before intersection" << endl;

              // compute the intersections of both contours
              intersections.clear();
              CGAL::intersection(prevCont, cont, back_inserter(intersections));

              //log << intersections.size() << " intersections computed" << endl;

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
        else
        {
          tmpPrevSection.push_back(secIdx - contours[i - 1].size() + cp);
        }
      }
      // add the list of previous section to connect to the 
      // current section
      prevSections.push_back(tmpPrevSection);
    }

    //log << "Intersection computed" << endl;

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

    //log << "Connections set" << endl;

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
    std::copy(begin(scalingFactors), end(scalingFactors), begin(prevScalingFactors));
    prevCurvRadius = curvRadius;
    prevAngle = angle;
  }

  //********************************
  // create last cross-sections
  //********************************

  //log << "\nCreate last cross-section" << endl;

  radius = 0.;  // initalise the radius of the radiation cross-section
  tmpPrevSection.clear(); // the list of section connected to the radiation section
  for (int c(0); c < contours.back().size(); c++)
  {
    // add the index of the created cross-section to the list of cross-sections
    // to connect to the radiation cross-section
    tmpPrevSection.push_back(secIdx);

    area = abs(contours.back()[c].area());
    length = centerLine.rbegin()[1].getDistanceFrom(centerLine.back());

    //log << "length " << length << endl;
    //log << "Scaling factors " << prevScalingFactors[0] << "  "
    //  << prevScalingFactors[1] << endl;

    addCrossSectionFEM(area, sqrt(area)/m_meshDensity, contours.back()[c],
      surfaceIdx.back()[c], length, centerLine.rbegin()[1], normals.rbegin()[1],
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
    //log << "Curv radius " << prevCurvRadius << endl;

    // set the curvature angle
    m_crossSections[secIdx]->setCurvatureAngle(prevAngle);
    //log << "angle " << prevAngle << endl;

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
    //log << "Create radiation cross-section" << endl;

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
  //log.close();
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
    else
    {
      // export junction cross-section

      // write centerline point coordinates
      strX << m_crossSections[i]->ctrLinePt().x << separator;
      strY << m_crossSections[i]->ctrLinePt().y << separator;

      // write normal coordinates
      strX << m_crossSections[i]->normal().x << separator;
      strY << m_crossSections[i]->normal().y << separator;

      // write contour
      for (auto pt : m_crossSections[i]->contour())
      {
        strX << pt.x() << separator;
        strY << pt.y() << separator;
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

void Acoustic3dSimulation::preComputeRadiationMatrices(int nbRadFreqs, int idxRadSec) {

  int numSec(m_crossSections.size());
  int mn(m_crossSections[idxRadSec]->numberOfModes());
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

  ofstream log("log.txt", ofstream::app);
  log << "Start precompute radiation impedance" << endl;
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
    //  freq, m_soundSpeed, m_volumicMass, 0);
    //m_crossSections[numSec - 1]->characteristicAdmittance(radAdmit,
    //  freq, m_soundSpeed, m_volumicMass, 0);
    //propagateImpedAdmit(radImped, radAdmit, freq, numSec - 1, numSec - 2, m_method);

    //radImped = m_crossSections[numSec - 2]->endImpedance();
    //radAdmit = m_crossSections[numSec - 2]->endAdmittance();

    radiationImpedance(radImped, freq, 15., idxRadSec);
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

    log << "Freq " << freq << " Hz " << i << " over " << nbRadFreqs << endl;
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
  log.close();
}

//*************************************************************************
// Interpolate the radiation impedance matrix at a given frequency

void Acoustic3dSimulation::interpolateRadiationImpedance(Eigen::MatrixXcd& imped, 
   double freq, int idxRadSec)
{
  //ofstream log("log.txt", ofstream::app);
  //log << "Start interpolate rad imped " << endl;

  // find the index corresponding to the coefficient to use for this frequency
  int nbRadFreqs(m_radiationMatrixInterp[0][0][0].size());
  int idx(nbRadFreqs - 2);
  while (m_radiationFreqs[idx] > freq) { idx--; }
  idx = max(0, idx);

  //log << "Idx coef " << idx << endl;

  // get number of modes of the radiating section
  int mn(m_crossSections[idxRadSec]->numberOfModes());

  //log << "Number of modes " << mn << endl;
  //log << "Size m_radiationMatrixInterp " << m_radiationMatrixInterp.size() << endl;

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
  //log.close();
}

//*************************************************************************
// Interpolate the radiation admittance matrix at a given frequency

void Acoustic3dSimulation::interpolateRadiationAdmittance(Eigen::MatrixXcd& admit, 
  double freq, int idxRadSec)
{
  // find the index corresponding to the coefficient to use for this frequency
  int nbRadFreqs(m_radiationMatrixInterp[0][0][0].size());
  int idx(nbRadFreqs - 2);
  while (m_radiationFreqs[idx] > freq) { idx--; }
  idx = max(0, idx);

  // get number of modes of the radiating section
  int mn(m_crossSections[idxRadSec]->numberOfModes());

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

void Acoustic3dSimulation::radiationImpedance(Eigen::MatrixXcd& imped, double freq, 
  double gridDensity, int idxRadSec)
{
  int mn(m_crossSections[idxRadSec]->numberOfModes());

  imped = Eigen::MatrixXcd::Zero(mn, mn);
  Eigen::MatrixXcd integral2(mn, mn);

  //ofstream log("log.txt", ofstream::app);
  //log << "\nStart computing radiation impedance" << endl;

  //******************************
  // generate cartesian grid
  //******************************

  // very good precision is obtained with gridDensity = 30
  //  good precision is obtained with gridDensity = 15;
  double scaling(m_crossSections[idxRadSec]->scaleOut());
  double spacing(sqrt(m_crossSections[idxRadSec]->area()) 
    / gridDensity);

  //log << "scaling: " << scaling << endl;
  //log << "Spacing: " << spacing << endl;

  //Transformation scale(CGAL::SCALING, scaling);
  //Polygon_2 contour(transform(scale, m_crossSections[idxRadSec]->contour()));
  Polygon_2 contour(m_crossSections[idxRadSec]->contour());
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

  //// export grid
  //ofstream out;
  //out.open("grid.txt");
  //for (auto it : cartGrid)
  //{
  //  out << it.x() << "  " << it.y() << endl;
  //}
  //out.close();
  //// export contour
  //out.open("cont.txt");
  //for (auto it : contour)
  //{
  //  out << it.x() << "  " << it.y() << endl;
  //}
  //out.close();

  // Interpolate the propagation modes on the cartesian grid
  Matrix intCartGrid(m_crossSections[idxRadSec]->interpolateModes(cartGrid));

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
      //<< "\nspacing:\t" << spacing
      //<< "\nNum directions:\t\t" << numDirections
      //<< "\nNum points:\t\t" << polGrid.size() << endl;

    // interpolate the polar grid
    Matrix intPolGrid(m_crossSections[idxRadSec]->interpolateModes(polGrid));

    //// export polar grid
    //if (c == 0) {
    //  out.open("pGrid.txt");
    //  for (auto it : polGrid)
    //  {
    //    out << it.x() << "  " << it.y() << endl;
    //  }
    //  out.close();
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
            exp(-1i * 2. * M_PI * freq * scaling * radius[p] / m_simuParams.sndSpeed);
        }
      }
    }

    imped += - integral2 / sumH / 2. / M_PI / cartGrid.size() / scaling;
  }

  imped *= pow(m_crossSections[idxRadSec]->area(),2);

  //log.close();
}
